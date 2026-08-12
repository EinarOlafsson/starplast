"""What "good structure" means, made explicit and selectable.

F1 is one answer to a question that has several, and the right one depends on what the structure is
for. Hunting a single compartment is not the same job as organising the whole proteome, and a search
that only ever maximises one blended number cannot tell you it found the other thing.

The four objectives a user actually asks for, in their own words, and what each one is:

    "most clusters are one label"                  PRECISION, averaged over clusters
    "most labels are in one cluster"               RECALL, averaged over labels
    "one or more clusters is mostly one label"     MAX PRECISION over (cluster, label) pairs
    "one or more labels is in mostly its own"      MAX F1 over pairs

Precision is `|c and l| / |c|` -- of this cluster's members, what fraction share a label. Recall is
`|c and l| / |l|` -- of this label's genes, what fraction are in this cluster. Everything below is a
different way of aggregating those two numbers.

EVERY OBJECTIVE HERE HAS A DEGENERATE MAXIMISER, and the guards matter more than the metrics:

  - mean precision over clusters is maximised by many TINY clusters; a singleton is perfectly pure.
  - mean recall over labels is maximised by ONE GIANT cluster; recall is 1.0 for every label. This
    is not hypothetical. On the full proteome the positive control scored 0.675 -- the highest of
    four targets -- from a two-cluster solution whose per-label recalls were 1.000, 1.000, 1.000.
  - max precision is maximised by one tiny pure cluster: three co-located genes, precision 1.0,
    and nothing to annotate from.

So `min_cluster_size` is enforced here rather than left to the caller, and every objective reports
the coverage and cluster count beside its score. A number without those is not interpretable.

WEIGHTING IS A CHOICE, not a detail. `search.score_recovery` weights by label size, so on this
proteome nucleus-chromatin at 769 genes dominates and dense granules at 167 barely registers -- and
someone hunting dense granules is optimising against themselves. `macro` weights every label
equally and is the right default for finding a rare class.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .search import ABSENCE_LABELS

NOISE = -1

#: The explainer, shown in the application under Help. Kept here beside the implementation so the
#: words and the code cannot drift apart, and so the API reference carries it too.
EXPLANATION = """\
PRECISION and RECALL answer different questions, and which one you want depends on what you are
looking for.

For one cluster c and one label l:

    precision = |c and l| / |c|     of THIS CLUSTER's members, what fraction share the label
    recall    = |c and l| / |l|     of THIS LABEL's genes, what fraction are in the cluster
    F1        = harmonic mean       both at once; high only when neither is low

So:

    "most clusters are one label"                 precision, averaged over CLUSTERS
    "most labels are in one cluster"              recall, averaged over LABELS
    "one or more clusters is mostly one label"    the single best precision
    "one or more labels has mostly its own"       the single best F1 -- "its own cluster"
                                                  claims both directions at once

EVERY ONE OF THESE HAS A DEGENERATE SOLUTION THAT WINS IT OUTRIGHT. Measured on a synthetic case
with three labels of sixty genes:

                        perfect     one big cluster     all singletons
    mean precision         1.00                0.33               1.00
    mean recall            1.00                1.00               0.02
    mean F1                1.00                0.50               0.03
    best precision         1.00                0.33               1.00
    V-measure              1.00                0.00               0.35
    ARI / AMI              1.00                0.00               0.00

Any precision objective is won by shattering the map into singletons, because a cluster of one is
perfectly pure. Any recall objective is won by one giant cluster, because everything is then in its
best cluster. Neither is hypothetical: on the full proteome the positive control scored 0.675 -- the
highest of four targets -- from a two-cluster solution whose per-label recalls were 1.000, 1.000 and
1.000. A hyperparameter walk is an optimiser, and an optimiser finds exactly what you reward.

Three things follow, and they are enforced rather than advised:

  - a minimum cluster size, which is what stops the singleton exploit;
  - coverage and cluster count reported beside every score, because the attention control's winning
    configuration was 52% noise;
  - chance-corrected agreement (ARI, AMI) available alongside, since it is the one family immune to
    BOTH degenerate cases. A high objective next to an ARI near zero means the objective was gamed.

WEIGHTING IS ALSO A CHOICE. Averaging over labels can weight every label equally (macro) or every
gene equally (size). On this proteome nucleus-chromatin has 769 genes and dense granules 167, so
size weighting means a search for dense granules is dominated by the nucleus. Macro is the default
for that reason.

FINALLY, precision_at_recall is usually the objective an annotation actually wants: the purest
cluster that still holds enough of the label to be worth annotating from. Best precision alone
selects three co-located genes at precision 1.00 and gives you nothing to work with."""

#: name -> one-line description, for a UI to offer and for a result to record.
OBJECTIVES = {
    "mean_precision": "most clusters are one label (purity, averaged over clusters)",
    "mean_recall": "most labels are in one cluster (completeness, averaged over labels)",
    # Recommended because it is the only one of the four common phrasings that resists BOTH
    # degenerate solutions on its own: singletons win any precision objective, one giant cluster
    # wins any recall objective, and needing both at once is what rules each of them out.
    "mean_f1": "labels and clusters correspond, on average  (recommended)",
    "best_precision": "at least one cluster is mostly one label",
    "best_f1": "at least one label has mostly its own cluster",
    "n_recovered": "how many labels clear a threshold (immune to class prevalence)",
    "v_measure": "the whole clustering agrees with the whole labelling",
    "precision_at_recall": "the purest cluster that still holds enough of a label to annotate from "
                           "(recommended for annotation)",
}

#: How many shuffles the permutation null uses by default. Twenty is enough to separate a real
#: structure from chance without doubling the cost of a walk.
NULL_PERMUTATIONS = 20


def pairs(labels: np.ndarray, truth: pd.Series, *, min_label: int = 15,
          min_cluster: int = 10, exclude_labels=ABSENCE_LABELS) -> pd.DataFrame:
    """Precision, recall and F1 for every (cluster, label) pair worth considering.

    One table that every objective below is computed from, so they cannot disagree about the
    underlying numbers. Noise is dropped: "it is in the noise" is not a structure.
    """
    labels = np.asarray(labels)
    v = truth.astype("object").where(truth.notna(), "").astype(str).to_numpy()
    keep = (labels != NOISE) & (v != "")
    lab, tru = labels[keep], v[keep]
    if len(lab) == 0:
        return pd.DataFrame(columns=["cluster", "label", "n_cluster", "n_label", "n_both",
                                     "precision", "recall", "f1"])
    absent = {str(x).lower() for x in exclude_labels}
    sizes_c = pd.Series(lab).value_counts()
    sizes_l = pd.Series(tru).value_counts()
    rows = []
    for c in sizes_c.index:
        if sizes_c[c] < min_cluster:
            continue
        for t in sizes_l.index:
            if str(t).lower() in absent or sizes_l[t] < min_label:
                continue
            both = int(((lab == c) & (tru == t)).sum())
            if not both:
                continue
            p = both / int(sizes_c[c])
            r = both / int(sizes_l[t])
            rows.append({"cluster": int(c), "label": str(t), "n_cluster": int(sizes_c[c]),
                         "n_label": int(sizes_l[t]), "n_both": both, "precision": p, "recall": r,
                         "f1": 2 * p * r / (p + r) if (p + r) else 0.0})
    return pd.DataFrame(rows)


def score(labels: np.ndarray, truth: pd.Series, *, objective: str = "mean_f1",
          weighting: str = "macro", category=None,
          min_recall: float = 0.25, threshold: float = 0.5, **kw) -> dict:
    """Score a clustering under one named objective.

    `category` restricts scoring to one label or to several -- a string, or any sequence of them.
    "Precision for dense granules" and "recall for dense granules and rhoptries" are both ordinary
    questions, and they are different questions, so the objective and the label set are independent
    choices rather than one combined mode.

    With several, the aggregation follows the objective: the mean_* objectives average over the
    chosen labels and the best_* objectives take the best among them. Nothing else changes, so a
    score over two labels is comparable with a score over all of them.

    `weighting` is macro (every label equal) or size (every gene equal); macro is the default
    because the alternative buries exactly the rare classes worth hunting.
    """
    if objective not in OBJECTIVES:
        raise ValueError(f"unknown objective {objective!r}; choose from {sorted(OBJECTIVES)}")
    P = pairs(labels, truth, **kw)
    labels = np.asarray(labels)
    n_clusters = int(len(set(labels[labels != NOISE])))
    coverage = float((labels != NOISE).mean()) if len(labels) else 0.0
    out = {"objective": objective, "weighting": weighting, "n_clusters": n_clusters,
           "coverage": coverage, "score": 0.0, "detail": ""}
    if P.empty:
        out["detail"] = "no (cluster, label) pair passed the size floors"
        return out
    if category is not None:
        wanted = [category] if isinstance(category, str) else list(category)
        out["categories"] = wanted
        P = P[P.label.isin(wanted)]
        if P.empty:
            shown = ", ".join(map(str, wanted))
            out["detail"] = f"no cluster holds enough of: {shown}"
            return out

    def _avg(frame, col, by):
        best = frame.loc[frame.groupby(by)[col].idxmax()]
        if weighting == "size" and by == "label":
            w = best.n_label / best.n_label.sum()
            return float((best[col] * w).sum())
        return float(best[col].mean())

    if objective == "mean_precision":
        # Averaged over CLUSTERS: asks of each cluster whether its members share a label.
        out["score"] = _avg(P, "precision", "cluster")
    elif objective == "mean_recall":
        out["score"] = _avg(P, "recall", "label")
    elif objective == "mean_f1":
        out["score"] = _avg(P, "f1", "label")
    elif objective == "best_precision":
        row = P.loc[P.precision.idxmax()]
        out["score"] = float(row.precision)
        out["detail"] = (f"cluster {int(row.cluster)} is {row.precision:.0%} {row.label} "
                         f"and holds {row.recall:.0%} of it")
    elif objective == "best_f1":
        row = P.loc[P.f1.idxmax()]
        out["score"] = float(row.f1)
        out["detail"] = (f"{row.label}: cluster {int(row.cluster)} at precision "
                         f"{row.precision:.2f}, recall {row.recall:.2f}")
    elif objective == "n_recovered":
        best = P.loc[P.groupby("label").f1.idxmax()]
        out["score"] = float((best.f1 >= threshold).sum())
        out["detail"] = f"of {len(best)} labels scored, at F1 >= {threshold}"
    elif objective == "v_measure":
        out["score"] = v_measure(labels, truth)
    elif objective == "precision_at_recall":
        # The annotation objective: the purest cluster that still holds enough of the label to be
        # worth annotating from. Precision alone picks three co-located genes.
        ok = P[P.recall >= min_recall]
        if ok.empty:
            out["detail"] = f"no cluster holds at least {min_recall:.0%} of any label"
            return out
        row = ok.loc[ok.precision.idxmax()]
        out["score"] = float(row.precision)
        out["detail"] = (f"{row.label}: cluster {int(row.cluster)} is {row.precision:.0%} pure "
                         f"while holding {row.recall:.0%} of the label")
    return out


def v_measure(labels: np.ndarray, truth: pd.Series) -> float:
    """Homogeneity and completeness together, as their harmonic mean.

    The whole-clustering version of the user's first two objectives: homogeneity IS "most clusters
    are one label" and completeness IS "most labels are in one cluster", so V-measure is what asking
    for both at once means. Unlike either alone it cannot be won by one giant cluster (completeness
    1, homogeneity 0) or by singletons (homogeneity 1, completeness 0).
    """
    try:
        from sklearn.metrics import v_measure_score
    except ImportError:                                    # pragma: no cover - optional dependency
        return float("nan")
    labels = np.asarray(labels)
    v = truth.astype("object").where(truth.notna(), "").astype(str).to_numpy()
    absent = {str(x).lower() for x in ABSENCE_LABELS}
    keep = (labels != NOISE) & (v != "") & ~np.isin(np.char.lower(v.astype(str)), list(absent))
    if keep.sum() < 2 or len(set(labels[keep])) < 1:
        return float("nan")
    return float(v_measure_score(v[keep], labels[keep]))


def agreement(labels: np.ndarray, truth: pd.Series) -> dict:
    """Chance-corrected agreement between the clustering and the labelling.

    Reported alongside any objective because it is the one family immune to BOTH degenerate cases:
    one giant cluster and all-singletons each score about zero, since neither beats chance. A high
    objective with an ARI near zero means the objective was gamed.
    """
    out = {"adjusted_rand": float("nan"), "adjusted_mutual_info": float("nan")}
    try:
        from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score
    except ImportError:                                    # pragma: no cover - optional dependency
        return out
    labels = np.asarray(labels)
    v = truth.astype("object").where(truth.notna(), "").astype(str).to_numpy()
    absent = {str(x).lower() for x in ABSENCE_LABELS}
    keep = (labels != NOISE) & (v != "") & ~np.isin(np.char.lower(v.astype(str)), list(absent))
    if keep.sum() < 2:
        return out
    out["adjusted_rand"] = float(adjusted_rand_score(v[keep], labels[keep]))
    out["adjusted_mutual_info"] = float(adjusted_mutual_info_score(v[keep], labels[keep]))
    return out

def null_score(labels: np.ndarray, truth: pd.Series, *, n_permutations: int = 20,
               seed: int = 0, **kw) -> float:
    """What this objective scores on the SAME clustering with the labels shuffled.

    The baseline moves with the number of classes and with how uneven they are, and none of the raw
    precision/recall/F1 objectives know that. Mean F1 of 0.30 across three balanced classes is worse
    than guessing; the same 0.30 across twenty-four compartments of wildly different sizes is
    remarkable. A search that compares configurations across targets is comparing against different
    nulls without saying so.

    Permutation rather than a formula: shuffling the labels preserves exactly how many classes there
    are, how many genes are in each, and the sizes of the clusters they are being matched against,
    so the null answers "what would this objective give me for nothing, on this data". That is the
    same idea that makes ARI and AMI trustworthy, applied to whichever objective was chosen.
    """
    rng = np.random.default_rng(seed)
    v = truth.astype("object").where(truth.notna(), "").astype(str).to_numpy()
    out = []
    for _ in range(max(1, n_permutations)):
        out.append(score(labels, pd.Series(rng.permutation(v)), **kw)["score"])
    return float(np.mean(out))


def adjusted(labels: np.ndarray, truth: pd.Series, *, n_permutations: int = 20,
             seed: int = 0, **kw) -> dict:
    """A score, its permutation null, and the score corrected for it.

    `adjusted = (score - null) / (1 - null)`, so 0 is "no better than shuffled labels" and 1 is
    perfect, whatever the objective and however many classes there are. Negative means the structure
    is worse than chance, which is a real result and should not be clipped away.
    """
    r = score(labels, truth, **kw)
    r["null"] = null_score(labels, truth, n_permutations=n_permutations, seed=seed, **kw)
    denom = 1.0 - r["null"]
    r["adjusted"] = (r["score"] - r["null"]) / denom if abs(denom) > 1e-9 else float("nan")
    r["n_labels"] = int(pd.Series(
        truth.astype("object").where(truth.notna(), "").astype(str)).nunique())
    return r
