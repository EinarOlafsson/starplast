"""Does a cluster-based annotation actually work? Measured by hiding labels we already have.

The workflow this exists to check: hold out a label, build a map, find a cluster that is mostly one
category, and propose the unlabelled members of that cluster as candidates. That produces a list of
genes and no error rate, and a candidate list with no error rate is a list of guesses.

The check is to make the same inference where the answer is already known. Hide a fraction of the
genes that DO carry a category, rebuild nothing -- the embedding is fixed and was built without the
label -- cluster it, and ask: of the hidden genes, how many land in the cluster the visible ones
picked out? That is precision and recall for exactly the inference being made, on genes whose true
label we withheld from ourselves.

Two things this is careful about.

The embedding must not have seen the label. If `compartment` fed the map, a cluster matching
compartment is circular and this measures nothing -- so the caller passes the columns the embedding
used and `masked_recovery` refuses when the target is among them.

The hidden genes are hidden from the SCORING, not from the embedding. Re-embedding per fold would
be the stricter test and costs an hour per fold on this proteome; it is offered as `refit` for when
that is affordable, and its absence is recorded in the result rather than glossed over, because a
number that did not re-fit is a weaker claim than one that did.
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
        if not self.folds:
            return f"{self.category}: not enough labelled genes to hide any"
        how = "re-embedded per fold" if self.refit else "one fixed embedding, labels hidden only from scoring"
        return (f"{self.category}: precision {self.precision:.2f}, recall {self.recall:.2f}, "
                f"F1 {self.f1:.2f} over {len(self.folds)} folds ({how})")


def _labelled(truth: pd.Series):
    """Indices of genes carrying a real label, absence excluded."""
    v = truth.astype("object").where(truth.notna(), "").astype(str)
    return np.flatnonzero(~np.isin(np.char.lower(v.to_numpy().astype(str)),
                                   list(ABSENCE_LABELS)))


def masked_recovery(labels: np.ndarray, truth: pd.Series, category: str, *,
                    folds: int = 5, hold_frac: float = 0.2, seed: int = 0,
                    used_columns=(), refit: bool = False) -> Validation:
    """Hide some of a category's genes and see whether the clustering puts them back.

    `labels` is a clustering of a FIXED embedding. `truth` is the held-out label column. For each
    fold, `hold_frac` of the genes carrying `category` are hidden; the cluster is chosen by which one
    holds most of the still-visible members -- exactly as a user would pick it by eye -- and then
    scored against the hidden ones, who had no say in the choice.

    Precision is over the cluster's members that were not visibly of this category, which is the set
    a user would actually annotate from. Recall is over the hidden genes.
    """
    if category in set(used_columns):
        raise ValueError(
            f"{category!r} was among the columns the embedding was built from, so a cluster matching "
            f"it is circular by construction and this measures nothing")

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
        # The cluster a user would pick: the one holding most of the visible members. Noise (-1) is
        # never a candidate -- "it is in the noise" is not an annotation.
        real = labels >= 0
        counts = {c: int(((labels == c) & visible).sum()) for c in np.unique(labels[real])}
        if not counts or max(counts.values()) == 0:
            continue
        best = max(counts, key=counts.get)
        in_cluster = (labels == best)

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


def validate_all(labels: np.ndarray, truth: pd.Series, *, min_size: int = 15,
                 **kw) -> pd.DataFrame:
    """Run `masked_recovery` for every category with enough genes, and rank them.

    Per-category rather than one global number, deliberately. A method that recovers hidden dense
    granules but not hidden rhoptries is not "60% accurate"; it works for one compartment and not the
    other, and only the per-category table says so.
    """
    v = truth.astype("object").where(truth.notna(), "").astype(str)
    counts = v.value_counts()
    rows = []
    for cat, n in counts.items():
        if str(cat).lower() in ABSENCE_LABELS or n < min_size:
            continue
        r = masked_recovery(labels, truth, str(cat), **kw)
        rows.append({"category": cat, "n_labelled": int(n), "n_folds": len(r.folds),
                     "precision": r.precision, "recall": r.recall, "f1": r.f1,
                     "note": r.note})
    if not rows:
        return pd.DataFrame(columns=["category", "n_labelled", "n_folds",
                                     "precision", "recall", "f1", "note"])
    return pd.DataFrame(rows).sort_values("f1", ascending=False).reset_index(drop=True)


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
    return pd.DataFrame({
        "gene_id": gene_ids[picked],
        "proposed": category,
        "cluster": cluster,
        "cluster_size": n_in,
        # How much of the cluster already carries the category, and how much carries something else.
        # A cluster that is 90% the category with 2% contradictions is a different proposition from
        # one that is 30% the category with 40% contradictions, and the numbers must travel with it.
        "cluster_frac_category": (n_cat / n_in) if n_in else 0.0,
        "cluster_frac_contradicting": (n_other / n_in) if n_in else 0.0,
    })
