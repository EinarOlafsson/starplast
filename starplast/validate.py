"""Does a cluster-based annotation actually work? Measured by hiding labels we already have.

The workflow this exists to check: hold out a label, build a map, find a cluster that is mostly one
category, and propose the unlabeled members of that cluster as candidates. That produces a list of
genes and no error rate, and a candidate list with no error rate is a list of guesses.

The check is to make the same inference where the answer is already known. Hide a fraction of the
genes that DO carry a category, rebuild nothing -- the embedding is fixed and was built without the
label -- cluster it, and ask: of the hidden genes, how many land in the cluster the visible ones
picked out? That is precision and recall for exactly the inference being made, on genes whose true
label we withheld from ourselves.

Two things this is careful about.

The embedding must not have seen the label. If `compartment` fed the map, a cluster matching
compartment is circular and this measures nothing -- so the caller passes the columns the embedding
used and the run is refused when the target column is among them. The check is on the COLUMN: it
used to be on the category value, comparing "dense granules" against a list of column names, which
is never true -- and the application passed block names, so the guard never fired at all.

The hidden genes are hidden from the SCORING, not from the embedding. Re-embedding per fold is the
stricter test and costs a full embedding per fold; `refit` does it, and requires a `rebuild`
callable to do it with. Asking for it without one is an error rather than a flag, because `refit`
used to change only the sentence the result prints -- so a run that re-fit nothing described itself
as "re-embedded per fold". Whether it re-fit is recorded on every row either way, since the
difference is invisible in the numbers.

A candidate list is the output all of this exists to qualify, so `candidates` attaches the
composition of the cluster it came from, and `orthogonal_support` counts how many genes already
carrying the category share an orthogroup or a domain with each candidate -- agreement from
evidence the map never saw. A column that fed the embedding is refused there too: its agreement
would be the map agreeing with itself.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

# What counts as "no label". Shared with the held-out search so absence means one thing everywhere.
from .search import ABSENCE_LABELS


@dataclass
class FoldResult:
    """One fold: what was hidden, and how much of it came back."""
    category: str
    fold: int
    n_hidden: int
    n_visible: int
    cluster: int
    precision: float          # of the cluster's unlabelled-at-scoring-time members, how many were the category
    recall: float             # of the hidden genes, how many landed in that cluster
    f1: float


@dataclass
class Validation:
    """The answer, with the caveats that have to travel with it."""
    category: str
    folds: list = field(default_factory=list)
    refit: bool = False
    note: str = ""

    @property
    def precision(self) -> float:
        return float(np.mean([f.precision for f in self.folds])) if self.folds else float("nan")

    @property
    def recall(self) -> float:
        return float(np.mean([f.recall for f in self.folds])) if self.folds else float("nan")

    @property
    def f1(self) -> float:
        return float(np.mean([f.f1 for f in self.folds])) if self.folds else float("nan")

    def summary(self) -> str:
        """One line: the numbers, and whether the embedding was refit per fold."""
        if not self.folds:
            return f"{self.category}: not enough labeled genes to hide any"
        how = "re-embedded per fold" if self.refit else "one fixed embedding, labels hidden only from scoring"
        return (f"{self.category}: precision {self.precision:.2f}, recall {self.recall:.2f}, "
                f"F1 {self.f1:.2f} over {len(self.folds)} folds ({how})")


def circularity_error(truth: pd.Series, used_columns, target_column: str | None = None) -> str:
    """Why this label cannot be scored against this embedding, or "" if it can.

    Returned as a sentence rather than raised from three call sites, because the interface has to
    SHOW this: a tab that pops an exception dialog teaches people that validation is broken, when
    what happened is that validation correctly refused to produce a meaningless number.

    The test is on the label COLUMN, not on the category value. It used to be the value -- `category
    in used_columns` -- which compares "dense granules" against a list of column names and is
    therefore never true. In the application it was worse than never true: the panel passed the
    embedding's BLOCK names, so the guard was comparing a compartment against "expression_summary"
    and the tab would happily score a target the map was built from.
    """
    col = target_column or getattr(truth, "name", None)
    used = set(used_columns)
    if col is None:
        return ""
    if col in used:
        return (f"{col!r} is one of the columns this embedding was built from, so a cluster matching "
                f"it is circular by construction and the number would mean nothing. Rebuild the map "
                f"without it -- on the Data tab -- and validate that.")

    # The same test the search applies, by name alone. Naming the label column is not enough: a map
    # built on `lopit_prob_map` is a map built on hyperLOPIT's own output, and scoring hyperLOPIT's
    # compartment against it is the same circularity by a longer route. The search guard was taught
    # this and this one was not -- which would have put a validated-looking number on precisely the
    # maps that cannot carry one, in the tab whose whole job is to say how much to believe a cluster.
    from . import datasets
    from .search import SAME_QUANTITY
    same = datasets.provenance(col)
    experiment = {c for c in (same.columns if same else ()) if c != col} & used
    quantity = {q: {c for c in family if c != col} & used
                for q, family in SAME_QUANTITY.items() if col in family}
    quantity = {q: cols for q, cols in quantity.items() if cols}
    reasons = []
    if experiment:
        reasons.append(f"{', '.join(sorted(experiment))} -- produced by the same experiment "
                       f"({same.name})")
    for q, cols in quantity.items():
        rest = sorted(cols - experiment)
        if rest:
            reasons.append(f"{', '.join(rest)} -- another estimate of the same thing ({q})")
    if reasons:
        return (f"this embedding was built from {'; and '.join(reasons)}. A cluster matching "
                f"{col!r} is then not evidence about the map. Rebuild without those columns -- on "
                f"the Data tab -- and validate that.")
    return ""


def masked_recovery(labels: np.ndarray, truth: pd.Series, category: str, *,
                    folds: int = 5, hold_frac: float = 0.2, seed: int = 0,
                    used_columns=(), target_column: str | None = None,
                    refit: bool = False, rebuild=None) -> Validation:
    """Hide some of a category's genes and see whether the clustering puts them back.

    `labels` is a clustering of a FIXED embedding. `truth` is the held-out label column. For each
    fold, `hold_frac` of the genes carrying `category` are hidden; the cluster is chosen by which one
    holds most of the still-visible members -- exactly as a user would pick it by eye -- and then
    scored against the hidden ones, who had no say in the choice.

    Precision is over the cluster's members that were not visibly of this category, which is the set
    a user would actually annotate from. Recall is over the hidden genes.

    `refit` requires `rebuild`, a callable taking the fold number and returning that fold's labels
    -- a fresh embedding and clustering, so the estimate is not conditional on one particular map.
    Asking for it without supplying one is an error rather than a flag: `refit=True` used to change
    only the sentence the result prints, so a run that had re-fit nothing reported itself as
    "re-embedded per fold". A number that overstates how it was obtained is worse than a weaker
    number that describes itself accurately.
    """
    problem = circularity_error(truth, used_columns, target_column)
    if problem:
        raise ValueError(problem)
    if refit and rebuild is None:
        raise ValueError("refit=True needs `rebuild`, a callable returning each fold's labels from "
                         "a fresh embedding; without one nothing would be re-embedded and the "
                         "result would claim otherwise")

    labels = np.asarray(labels)
    v = truth.astype("object").where(truth.notna(), "").astype(str).to_numpy()
    is_cat = v == category
    idx = np.flatnonzero(is_cat)
    out = Validation(category=category, refit=refit)
    n_hold = int(round(len(idx) * hold_frac))
    if len(idx) < 4 or n_hold < 1:
        out.note = f"only {len(idx)} genes carry {category!r}; too few to hide any"
        return out

    rng = np.random.default_rng(seed)
    for k in range(folds):
        hidden = rng.choice(idx, size=n_hold, replace=False)
        hidden_mask = np.zeros(len(labels), bool)
        hidden_mask[hidden] = True
        visible = is_cat & ~hidden_mask
        if not visible.any():
            continue
        # A fresh map per fold when one was asked for. The label never fed the embedding either way,
        # so what re-fitting buys is that the answer is not conditional on one map and one
        # clustering -- and on this proteome that is the difference between an estimate and an
        # anecdote about a particular seed.
        fold_labels = np.asarray(rebuild(k)) if refit else labels
        if len(fold_labels) != len(labels):
            raise ValueError(f"fold {k}: rebuild returned {len(fold_labels)} labels for "
                             f"{len(labels)} genes")
        # The cluster a user would pick: the one holding most of the visible members. Noise (-1) is
        # never a candidate -- "it is in the noise" is not an annotation.
        real = fold_labels >= 0
        counts = {c: int(((fold_labels == c) & visible).sum()) for c in np.unique(fold_labels[real])}
        if not counts or max(counts.values()) == 0:
            continue
        best = max(counts, key=counts.get)
        in_cluster = (fold_labels == best)

        # Precision over exactly the set a user would annotate from: cluster members that do not
        # visibly carry the category. Some of those are the hidden ones, and those are the hits.
        candidates = in_cluster & ~visible
        n_cand = int(candidates.sum())
        hits = int((candidates & hidden_mask).sum())
        precision = hits / n_cand if n_cand else 0.0
        recall = hits / n_hold
        f1 = 0.0 if precision + recall == 0 else 2 * precision * recall / (precision + recall)
        out.folds.append(FoldResult(category=category, fold=k, n_hidden=n_hold,
                                    n_visible=int(visible.sum()), cluster=int(best),
                                    precision=precision, recall=recall, f1=f1))
    return out


def validate_all(labels: np.ndarray, truth: pd.Series, *, min_size: int = 15, log=print,
                 **kw) -> pd.DataFrame:
    """Run `masked_recovery` for every category with enough genes, and rank them.

    Per-category rather than one global number, deliberately. A method that recovers hidden dense
    granules but not hidden rhoptries is not "60% accurate"; it works for one compartment and not the
    other, and only the per-category table says so.

    `log` is called once per category, which is what makes the run reportable and, in the
    application, stoppable: the cancellation is raised from the reporting callable, so a run with
    `refit` on -- a full embedding per fold per category -- can be stopped between categories
    instead of at the end.
    """
    # Checked once, here, as well as per category. A table where every category was too small to
    # score would otherwise come back empty and clean from a target the embedding was built on,
    # which reads as "nothing to report" rather than as "this question cannot be asked of this map".
    problem = circularity_error(truth, kw.get("used_columns", ()), kw.get("target_column"))
    if problem:
        raise ValueError(problem)
    v = truth.astype("object").where(truth.notna(), "").astype(str)
    counts = v.value_counts()
    cols = ["category", "n_labelled", "n_folds", "precision", "recall", "f1", "refit", "note"]
    rows = []
    scorable = [(c, n) for c, n in counts.items()
                if str(c).lower() not in ABSENCE_LABELS and n >= min_size]
    log(f"validating {len(scorable)} categories of {getattr(truth, 'name', 'the target')!r} "
        f"({int(counts.sum()):,} genes, {len(counts) - len(scorable)} categories too small or absent)")
    for i, (cat, n) in enumerate(scorable, start=1):
        log(f"  {i}/{len(scorable)} {cat} ({int(n):,} labeled)")
        r = masked_recovery(labels, truth, str(cat), **kw)
        rows.append({"category": cat, "n_labelled": int(n), "n_folds": len(r.folds),
                     "precision": r.precision, "recall": r.recall, "f1": r.f1,
                     # Carried on every row: a score that did not re-embed per fold is a weaker
                     # claim than one that did, and the difference is invisible in the numbers.
                     "refit": bool(r.refit), "note": r.note})
    if not rows:
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(rows, columns=cols).sort_values("f1", ascending=False).reset_index(drop=True)


def best_cluster(labels: np.ndarray, truth: pd.Series, category: str) -> int | None:
    """The cluster holding most of a category's labeled genes, or None if none does.

    The same choice `masked_recovery` makes inside each fold, and the same one a user makes by eye,
    so the cluster a candidate list comes from is the cluster the error rate was measured on. Noise
    is never chosen: "it is in the noise" is not an annotation.
    """
    labels = np.asarray(labels)
    v = truth.astype("object").where(truth.notna(), "").astype(str).to_numpy()
    is_cat = v == category
    counts = {int(c): int(((labels == c) & is_cat).sum()) for c in np.unique(labels[labels >= 0])}
    if not counts or max(counts.values()) == 0:
        return None
    return max(counts, key=counts.get)


#: Grouping columns whose agreement is evidence the map did not use: shared orthogroup, shared Pfam,
#: shared InterPro. Multi-valued columns are `;`-separated in this table.
SUPPORT_COLUMNS = ("orthogroup", "pfam_id", "interpro_id")


def orthogonal_support(nodes: pd.DataFrame, candidate_ids, truth: pd.Series, category: str, *,
                       columns=SUPPORT_COLUMNS, used_columns=(), log=print) -> pd.DataFrame:
    """For each candidate, how many genes already labeled `category` it shares a group with.

    A candidate's whole claim is that it sits near genes of one category in a map. That is one piece
    of evidence, and the honest question about it is whether anything the map never saw agrees.
    Sharing an orthogroup or a domain with the labeled members is such a thing: not proof -- paralogs
    of a dense-granule protein need not be dense-granule proteins -- but independent, and countable.

    A column that fed the embedding is refused rather than counted, and the refusal is logged. Its
    agreement would be the map agreeing with itself, which is the failure this whole tab exists to
    put a number on.

    Zero support is the common answer and is reported as zero, never as a blank: most of this
    proteome has no orthogroup neighbour with a measured label, and a candidate with no independent
    support is exactly the one to be careful about.
    """
    ids = np.asarray(candidate_ids, dtype=object)
    out = pd.DataFrame({"gene_id": ids})
    v = truth.astype("object").where(truth.notna(), "").astype(str).to_numpy()
    labelled = np.flatnonzero(v == category)
    gene_id = nodes["gene_id"].to_numpy()
    pos = {g: i for i, g in enumerate(gene_id)}
    used = set(used_columns)
    for col in columns:
        if col not in nodes.columns:
            continue
        if col in used:
            log(f"orthogonal support: {col!r} fed the embedding, so its agreement is not "
                f"independent evidence and it is not counted")
            continue
        tokens = _tokens(nodes[col])
        # token -> the labelled genes carrying it, so a labelled gene sharing two domains with a
        # candidate is counted once rather than twice.
        holders = {}
        for i in labelled:
            for t in tokens[i]:
                holders.setdefault(t, set()).add(int(i))
        support = []
        for g in ids:
            i = pos.get(g)
            hits = set()
            if i is not None:
                for t in tokens[i]:
                    hits |= holders.get(t, set())
            support.append(len(hits))
        out[f"shares_{col}"] = support
    return out


def _tokens(values: pd.Series) -> list:
    """Each row's group memberships. `;`-separated where a gene has several domains."""
    v = values.astype("object").where(values.notna(), "").astype(str).to_numpy()
    return [[t for t in str(x).split(";") if t and t.lower() not in ABSENCE_LABELS] for x in v]


def candidates(labels: np.ndarray, truth: pd.Series, category: str, cluster: int,
               gene_ids) -> pd.DataFrame:
    """The genes a cluster would have you annotate: its members with no label of their own.

    Returned with the cluster's composition attached, so a candidate never travels without the
    number that says how much to believe it.
    """
    v = truth.astype("object").where(truth.notna(), "").astype(str).to_numpy()
    gene_ids = np.asarray(gene_ids)
    labels = np.asarray(labels)
    inc = labels == cluster
    absent = np.isin(np.char.lower(v.astype(str)), list(ABSENCE_LABELS))
    n_in = int(inc.sum())
    n_cat = int((inc & (v == category)).sum())
    n_other = int((inc & ~absent & (v != category)).sum())
    picked = inc & absent
    # How common the category is across the whole map, which is the number the cluster's own
    # fraction has to be read against. A cluster that is 10% dense granules sounds like something
    # until you notice that 9% of the map is, and then it is a cluster like any other.
    #
    # Over the same denominator as `cluster_frac_category` -- every gene, not only the labelled ones
    # -- because the two are meant to be divided by each other. Computed over labelled genes alone
    # it read 1.3x higher than the cluster fraction it is compared against, on this table, purely
    # from the 68% of the proteome that carries no label.
    prevalence = (int((v == category).sum()) / len(v)) if len(v) else 0.0
    frac = (n_cat / n_in) if n_in else 0.0
    return pd.DataFrame({
        "gene_id": gene_ids[picked],
        "proposed": category,
        "cluster": cluster,
        "cluster_size": n_in,
        # How much of the cluster already carries the category, and how much carries something else.
        # A cluster that is 90% the category with 2% contradictions is a different proposition from
        # one that is 30% the category with 40% contradictions, and the numbers must travel with it.
        "cluster_frac_category": frac,
        "cluster_frac_contradicting": (n_other / n_in) if n_in else 0.0,
        "overall_frac_category": prevalence,
        "enrichment": (frac / prevalence) if prevalence else float("nan"),
    })
