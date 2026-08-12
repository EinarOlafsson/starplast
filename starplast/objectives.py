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

#: name -> one-line description, for a UI to offer and for a result to record.
OBJECTIVES = {
    "mean_precision": "most clusters are one label (purity, averaged over clusters)",
    "mean_recall": "most labels are in one cluster (completeness, averaged over labels)",
    "mean_f1": "labels and clusters correspond, on average",
    "best_precision": "at least one cluster is mostly one label",
    "best_f1": "at least one label has mostly its own cluster",
    "n_recovered": "how many labels clear a threshold (immune to class prevalence)",
    "v_measure": "the whole clustering agrees with the whole labelling",
    "precision_at_recall": "the purest cluster that still holds enough of a label to annotate from",
}


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
          weighting: str = "macro", category: str | None = None,
          min_recall: float = 0.25, threshold: float = 0.5, **kw) -> dict:
    """Score a clustering under one named objective.

    `category` restricts to a single label, which is what "find me a map where the GRAs cluster"
    means. `weighting` is macro (every label equal) or size (every gene equal); macro is the default
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
        P = P[P.label == category]
        if P.empty:
            out["detail"] = f"no cluster holds enough {category!r}"
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
