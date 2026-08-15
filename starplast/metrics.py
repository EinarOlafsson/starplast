#!/usr/bin/env python3
"""Scoring a map by how well it would let you ANNOTATE, not by how tidy it looks.

The scores this project already had are agreement scores: F1 between clusters and labels, V-measure,
adjusted Rand. They answer "does this clustering correspond to that labelling", which is a question
about a partition and inherits every weakness of one -- change `min_cluster_size` and the number
moves, because a different partition of the same map is a different answer.

The ones here answer a different question: **if I ranked every gene by how much its neighbours
suggest it is a rhoptry protein, how good would that ranking be?** That is what guilt by association
actually does, and a ranking has its own well-worn statistics.

## Why AUPRC is the one to read, and AUROC the one to report

Both summarise a ranking; they disagree about what matters, and on this data the disagreement is
large.

**AUROC** asks: take a random rhoptry gene and a random non-rhoptry gene -- how often is the rhoptry
one ranked higher? It is prevalence-independent, which sounds like a virtue and here is a trap. With
120 rhoptry genes among 8,140, a ranking can put four thousand irrelevant genes above the first true
one and still score 0.9, because the metric spends most of its mass on comparisons nobody will ever
make. Reported because it is what everyone else reports and comparisons need it.

**AUPRC** (average precision) asks about the top of the list only: of the genes I would actually go
and test, what share are right. It moves with prevalence, so a raw value means nothing on its own --
which is why `lift` is computed alongside it, dividing by the prevalence a random ranking achieves.
An AUPRC of 0.30 for a class that is 1.5% of the proteome is a twentyfold enrichment over chance and
a genuinely useful shortlist; an AUPRC of 0.30 for a class that is 28% is worse than guessing.

## Leave-one-out, which is the whole reason these numbers are honest

A gene's score for class c is the share of its neighbours carrying c -- **excluding itself**. Without
that exclusion every labelled gene votes for its own label, every score is inflated by exactly one
vote, and small clusters (where one vote is a large share) are inflated most. The result is a map
that scores beautifully and predicts nothing, which is the same failure as reading an enrichment off
a layer that built the map.

## What leave-one-out costs, stated because it changes how the numbers read

Removing a gene's own vote biases the chance level DOWNWARD, and not by a rounding error. In a
cluster of n genes of which m carry the class, a positive gene scores (m-1)/(n-1) and a negative one
scores m/(n-1) -- so inside any cluster that carries no real information, every positive sits just
below every negative. Measured on four balanced classes in clusters of sixty, a randomly shuffled
clustering scores AUROC 0.36 rather than 0.50.

Three consequences worth carrying:

- **Do not read 0.5 as the floor.** `chance_level` measures the floor for a given clustering by
  permuting the labels, the same permutation null `objectives.null_score` uses, and that is what a
  score should be compared against.
- **The bias shrinks as clusters grow**, which is a mild pressure toward larger clusters among
  clusterings that are equally uninformative. It does not reach clusterings that carry real signal
  -- a pure cluster scores near 1 either way -- but it is a reason not to optimise on AUROC alone.
- **AUPRC lift is the more robust target**, being normalised by prevalence, which is why it is what
  `optimize` climbs.

## The neighbourhood metrics, and why they matter more than they look

`knn_purity` and `trustworthiness` score the EMBEDDING, with no clustering involved at all. That
makes them the right tools for the problem this project actually has -- a clustering that returns two
clusters over eight thousand genes tells you nothing about whether the map underneath it is good.
Separating "is the map good" from "did the clustering find what is in it" is the difference between
tuning two things at once and tuning them one at a time.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .clustering import NOISE
from .discovery import ABSENT

#: Neighbours consulted when scoring the embedding directly. 15 is the usual working value: small
#: enough that a neighbourhood is one biological group, large enough that a single mislabelled gene
#: does not swing it.
K = 15


def _classes(truth: pd.Series, min_count: int = 10):
    """The categories worth scoring: real labels, present often enough to have a curve at all."""
    s = truth.astype("object")
    text = s.map(lambda v: str(v).strip().lower() if v is not None and v == v else "")
    clean = s.where(~text.isin(ABSENT) & (text != ""), other=np.nan)
    counts = clean.value_counts()
    return clean, [c for c, n in counts.items() if n >= min_count]


def cluster_scores(labels: np.ndarray, truth: pd.Series, category) -> np.ndarray:
    """Each gene's score for one category: the share of its CLUSTER carrying it, minus itself.

    Noise points score zero rather than being dropped. A map that calls half the proteome noise has
    genuinely failed to rank those genes, and letting them out of the denominator would hide exactly
    that failure behind a better-looking curve.
    """
    labels = np.asarray(labels)
    is_c = (truth.to_numpy() == category).astype(float)
    out = np.zeros(len(labels), dtype=float)
    for k in {int(x) for x in labels}:
        if k == NOISE:
            continue
        here = labels == k
        n = int(here.sum())
        if n <= 1:
            continue
        # Leave-one-out: a gene's own label is removed from the vote that scores it.
        out[here] = (is_c[here].sum() - is_c[here]) / (n - 1)
    return out


def neighbour_scores(coords: np.ndarray, truth: pd.Series, category, k: int = K) -> np.ndarray:
    """The same score from the k nearest genes in the map, with no clustering in the way."""
    from sklearn.neighbors import NearestNeighbors
    coords = np.asarray(coords, dtype=float)
    k = int(min(max(k, 1), max(len(coords) - 1, 1)))
    nn = NearestNeighbors(n_neighbors=k + 1).fit(coords)
    idx = nn.kneighbors(coords, return_distance=False)[:, 1:]      # column 0 is the gene itself
    return (truth.to_numpy()[idx] == category).mean(axis=1)


def curve_scores(score: np.ndarray, hit: np.ndarray) -> dict:
    """AUROC, AUPRC and the lift of AUPRC over a random ranking, for one class.

    Undefined rather than zero where a class has no positives or no negatives: an AUROC over a
    column of ones is not 0.5, it is a question that was not asked, and averaging a fabricated 0.5
    into a summary quietly drags every real number toward the middle.
    """
    from sklearn.metrics import average_precision_score, roc_auc_score
    hit = np.asarray(hit).astype(bool)
    if hit.all() or not hit.any() or len(hit) < 3:
        return {"auroc": np.nan, "auprc": np.nan, "lift": np.nan, "prevalence": float(hit.mean())}
    prevalence = float(hit.mean())
    auprc = float(average_precision_score(hit, score))
    return {"auroc": float(roc_auc_score(hit, score)), "auprc": auprc,
            "lift": auprc / prevalence if prevalence > 0 else np.nan, "prevalence": prevalence}


def chance_level(labels: np.ndarray, truth: pd.Series, category, n_permutations: int = 20,
                 seed: int = 42) -> dict:
    """What this clustering scores for this category when the labels are shuffled.

    The floor is not 0.5 and not the prevalence -- leave-one-out moves both, by an amount that
    depends on how big the clusters are. Measuring it costs twenty cheap re-scorings and is the
    difference between "AUROC 0.61, promising" and "AUROC 0.61 against a floor of 0.58".
    """
    rng = np.random.default_rng(seed)
    clean, _ = _classes(truth, 1)
    known = clean.notna().to_numpy()
    out = []
    for _ in range(max(int(n_permutations), 1)):
        shuffled = pd.Series(rng.permutation(clean.to_numpy()), index=clean.index)
        # Scored against the SHUFFLED membership, not the real one. Keeping the real `hit` here
        # looks like the more obvious null and measures the wrong thing: the leave-one-out deduction
        # is tied to a gene's own label, so scoring shuffled labels against true membership breaks
        # exactly the coupling that causes the bias and reports a floor of 0.49 for a metric whose
        # floor is really 0.36.
        out.append(curve_scores(cluster_scores(labels, shuffled, category)[known],
                                (shuffled.to_numpy() == category)[known]))
    return {k: float(np.nanmean([o[k] for o in out])) for k in ("auroc", "auprc", "lift")}


def ranking(labels: np.ndarray, truth: pd.Series, coords: np.ndarray = None, k: int = K,
            min_count: int = 10) -> pd.DataFrame:
    """One row per category: how well this map ranks the genes that belong to it.

    Both rankings where coordinates are given -- by cluster and by neighbourhood -- because the gap
    between them is diagnostic. A map with good neighbourhood numbers and poor cluster numbers has
    the structure and lost it in the clustering, which is a different problem with a different fix.
    """
    if pd.api.types.is_numeric_dtype(truth):
        # A continuous measurement has no classes to rank. Treating each repeated float as a class
        # produced tiny but plausible-looking AUPRC values in saved searches; those numbers answer
        # no biological question. Keep the table shape stable and mark why it is empty so
        # `summarise` can distinguish this from a categorical layer with too little coverage.
        out = pd.DataFrame(columns=["category", "n", "auroc", "auprc", "lift", "prevalence"])
        out.attrs["continuous_truth"] = True
        return out
    clean, categories = _classes(truth, min_count)
    known = clean.notna().to_numpy()
    rows = []
    for category in categories:
        hit = (clean.to_numpy() == category)[known]
        row = {"category": str(category), "n": int(hit.sum())}
        row.update(curve_scores(cluster_scores(labels, clean, category)[known], hit))
        if coords is not None and len(coords) == len(clean):
            by_nn = curve_scores(neighbour_scores(coords, clean, category, k)[known], hit)
            row.update({f"nn_{a}": b for a, b in by_nn.items() if a != "prevalence"})
        rows.append(row)
    columns = ["category", "n", "auroc", "auprc", "lift", "prevalence"]
    if not rows:
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows).sort_values("auprc", ascending=False).reset_index(drop=True)


def summarise(table: pd.DataFrame) -> dict:
    """A ranking table as a handful of numbers, weighted and unweighted.

    Both, because they answer different questions and disagree loudly on this data. The plain mean
    is a statement about CATEGORIES -- a rare compartment counts as much as a common one, which is
    what a reader looking for somewhere to work wants. The weighted mean is a statement about GENES,
    and is dominated by whichever two categories are largest.
    """
    if table is not None and table.empty and table.attrs.get("continuous_truth"):
        return {"mean_auprc": np.nan, "mean_auroc": np.nan, "mean_lift": np.nan,
                "weighted_auprc": np.nan, "n_categories": 0}
    if table is None or table.empty:
        return {"mean_auprc": 0.0, "mean_auroc": 0.0, "mean_lift": 0.0,
                "weighted_auprc": 0.0, "n_categories": 0}
    w = table["n"].to_numpy(dtype=float)
    good = np.isfinite(table["auprc"].to_numpy(dtype=float))
    out = {"n_categories": int(good.sum())}
    for column in ("auprc", "auroc", "lift"):
        v = pd.to_numeric(table.get(column), errors="coerce").to_numpy(dtype=float)
        out[f"mean_{column}"] = float(np.nanmean(v)) if good.any() else 0.0
    out["weighted_auprc"] = (float(np.nansum(table.auprc.to_numpy(dtype=float) * w)
                                   / max(np.nansum(w[good]), 1e-9)) if good.any() else 0.0)
    for column in ("nn_auprc", "nn_auroc"):
        if column in table.columns:
            out[f"mean_{column}"] = float(np.nanmean(pd.to_numeric(table[column],
                                                                   errors="coerce")))
    return out


# --------------------------------------------------------------------------- the map, on its own
def knn_purity(coords: np.ndarray, truth: pd.Series, k: int = K) -> float:
    """What share of a labelled gene's k nearest neighbours share its label.

    No clustering anywhere in it, which is the point: this is the number to tune an embedding on
    before choosing how to cut it up. Chance level is the sum of squared class frequencies, not
    zero, so read it against `expected_purity` below.
    """
    clean, _ = _classes(truth, 1)
    known = clean.notna().to_numpy()
    if known.sum() < 3:
        return 0.0
    from sklearn.neighbors import NearestNeighbors
    coords = np.asarray(coords, dtype=float)[known]
    values = clean[known].to_numpy()
    k = int(min(max(k, 1), max(len(coords) - 1, 1)))
    idx = NearestNeighbors(n_neighbors=k + 1).fit(coords).kneighbors(
        coords, return_distance=False)[:, 1:]
    return float((values[idx] == values[:, None]).mean())


def expected_purity(truth: pd.Series) -> float:
    """What `knn_purity` would be if the map were shuffled: the sum of squared class frequencies."""
    clean, _ = _classes(truth, 1)
    freq = clean.value_counts(normalize=True).to_numpy(dtype=float)
    return float((freq ** 2).sum()) if len(freq) else 0.0


def trustworthiness(X: np.ndarray, coords: np.ndarray, k: int = K, sample: int = 2000,
                    seed: int = 42) -> float:
    """How much of the original neighbourhood structure survived the projection, in [0, 1].

    The one metric here that needs no labels at all, and the only defence against a map that looks
    beautifully clustered because the projection invented the clusters. Subsampled, because it is
    quadratic in the number of genes and the estimate is stable well before 8,140 of them.
    """
    from sklearn.manifold import trustworthiness as _t
    X, coords = np.asarray(X, dtype=float), np.asarray(coords, dtype=float)
    if len(X) != len(coords) or len(X) < k + 2:
        return float("nan")
    if len(X) > sample:
        take = np.random.default_rng(seed).choice(len(X), sample, replace=False)
        X, coords = X[take], coords[take]
    return float(_t(np.nan_to_num(X), coords, n_neighbors=int(min(k, len(X) - 1))))


def partition(labels: np.ndarray, truth: pd.Series, coords: np.ndarray = None) -> dict:
    """The agreement and shape statistics for a clustering, in one call.

    Silhouette is on the clustering's own geometry and says nothing about biology -- a map cut into
    tidy round groups that correspond to nothing scores well. It is here to be read ALONGSIDE the
    label metrics: high silhouette with low AUPRC is the signature of a clustering that found the
    shape of the projection rather than the shape of the data.
    """
    from sklearn.metrics import (adjusted_mutual_info_score, adjusted_rand_score,
                                 silhouette_score)
    labels = np.asarray(labels)
    clean, _ = _classes(truth, 1)
    known = clean.notna().to_numpy() & (labels != NOISE)
    sizes = pd.Series(labels[labels != NOISE]).value_counts().to_numpy(dtype=float)
    out = {
        "n_clusters": int(len(sizes)),
        "noise": float(np.mean(labels == NOISE)),
        # How lopsided the clustering is: 0 when every cluster is the same size, near 1 when one
        # cluster holds everything. The two-clusters-over-eight-thousand failure, as a number.
        "largest_share": float(sizes.max() / max(sizes.sum(), 1)) if len(sizes) else 1.0,
        "ari": float("nan"), "ami": float("nan"), "silhouette": float("nan"),
    }
    if known.sum() > 2 and len(set(labels[known])) > 1:
        out["ari"] = float(adjusted_rand_score(clean[known].astype(str), labels[known]))
        out["ami"] = float(adjusted_mutual_info_score(clean[known].astype(str), labels[known]))
    if coords is not None and len(set(labels[labels != NOISE])) > 1:
        keep = labels != NOISE
        out["silhouette"] = float(silhouette_score(np.asarray(coords)[keep], labels[keep]))
    return out


def report(labels: np.ndarray, truth: pd.Series, coords: np.ndarray = None, X: np.ndarray = None,
           k: int = K) -> dict:
    """Every number above, for one configuration, as a flat dict a results table can hold."""
    table = ranking(labels, truth, coords, k=k)
    out = {**summarise(table), **partition(labels, truth, coords)}
    if coords is not None:
        out["knn_purity"] = knn_purity(coords, truth, k)
        out["knn_purity_chance"] = expected_purity(truth)
        out["knn_lift"] = (out["knn_purity"] / out["knn_purity_chance"]
                           if out["knn_purity_chance"] > 0 else float("nan"))
    if X is not None and coords is not None:
        out["trustworthiness"] = trustworthiness(X, coords, k)
    return out
