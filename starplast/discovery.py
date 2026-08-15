#!/usr/bin/env python3
"""Three questions worth optimising a map for, and the instances of each that it finds.

Everything else in this project scores a structure by how well it RECOVERS something already known:
hold out the localisation labels, cluster, see how many come back. That is the right way to decide
whether a map is trustworthy and the wrong way to decide whether it is useful, because a map that
recovers what is known perfectly and says nothing about anything else has produced no work.

These are the two shapes of claim a map like this can actually make.

**Guilt by association.** A cluster whose labelled members are overwhelmingly one thing, and which
also holds genes nobody has labelled. The labelled part is the evidence, the unlabelled part is the
prediction, and the strength of the claim is the enrichment: 41 of 58 labelled genes in one cluster
being IMC is a statement about the 12 unlabelled ones sitting with them. The same shape works on a
measured quantity rather than a category -- a cluster whose screened members are uniformly essential
says something about the members of that cluster nobody screened.

**Layer disagreement.** A cluster that is homogeneous in one layer and split in another. Same
compartment, opposite fitness. Same cell-cycle phase, opposite stage expression. Agreement between
layers mostly re-derives what is already known; disagreement is where a category is hiding a
distinction, and it is the more interesting of the two because nothing about it is predictable from
either layer alone.

## What makes these honest rather than a slot machine

Three things, all of them enforced here rather than left to the reader.

**Multiplicity.** A run tests every cluster against every category. At 40 clusters and 12 categories
that is 480 tests, and at p < 0.05 two dozen come back by chance. Every p-value here is corrected
across the whole family of tests the run performed, and the corrected value is what ranks a finding.

**Circularity.** A cluster enriched for localisation, in a map built from localisation, is a
statement about arithmetic. Findings carry the columns that built the map, and one whose layer was
an input is marked rather than dropped -- it is still a fact about the map, it is just not evidence
about biology. The same guard as `search.excluded_for`, applied at the point of interpretation.

**Prevalence.** A cluster that is 90% "unassigned" in a proteome that is 60% unassigned has found
nothing. Enrichment is always against the background rate in the labelled population, never against
uniformity.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .clustering import NOISE, _bh

#: Kinds of layer, and what a claim about each one means. A category is something a gene IS; a
#: quantity is something it MEASURES AT. The two need different statistics and produce differently
#: worded claims, and keeping them apart here is what lets one optimiser drive both.
DISCRETE, CONTINUOUS = "discrete", "continuous"

#: Labels that mean "nobody has said", which is not a category a cluster can be enriched for. Shared
#: in spirit with `search.ABSENCE_LABELS`; kept separate because that set guards a scoring loop and
#: this one guards a claim.
ABSENT = frozenset({"unassigned", "unknown", "none", "nan", "", "unlabelled", "unlabeled",
                    "no data", "not assigned", "unclear"})

#: The smallest cluster, and the smallest count of a category within one, that may carry a claim.
#: Below this the statistics are formally fine and the biology is not: three genes agreeing is a
#: coincidence with a p-value attached.
MIN_CLUSTER = 10
MIN_CATEGORY = 5

#: The smallest number of genes a finding must be ABOUT to be worth reporting. A cluster that is
#: beautifully enriched and contains no unlabelled genes has confirmed something, not found it.
MIN_PREDICTIONS = 3

#: The three guards that decide whether an enrichment MEANS anything, as opposed to merely being
#: significant. Every one of them was added after watching this module report nonsense with a
#: p-value of 1e-32 attached, and each rules out a different way of getting one.
#:
#: `MIN_LIFT` -- how many times the background rate. A p-value is a statement about sample size as
#: much as effect: 395 of 1621 genes being cytosol where 15.9% of the proteome is cytosol is a lift
#: of 1.5 and a p of 1e-32, and it predicts nothing about anything. Below 2x, a cluster is a slice
#: of the proteome rather than a group of related genes.
#: `MIN_PURITY` -- what share of the labelled members agree. A cluster that is 10% ER is not an ER
#: cluster however unlikely that 10% was.
#: `MAX_SHARE` -- how much of the map one cluster may cover. A clustering that returns two clusters
#: over eight thousand genes will report both as enriched for something, and "these 2,085 unlabelled
#: genes are cytosol" is not a prediction, it is a clustering that failed stated as a result.
MIN_LIFT = 2.0
MIN_PURITY = 0.25
MAX_SHARE = 0.25

#: How far apart the two halves of a split must sit, in standard deviations of the layer ACROSS THE
#: WHOLE MAP. One is deliberately demanding: a subdivision that moves a gene by less than the
#: ordinary spread of the measurement is not a subdivision anybody can act on.
MIN_SPLIT = 1.0

#: The smallest standardised difference a quantity may show and still be called a finding. Half a
#: standard deviation is the conventional "small but real"; below it, a cluster large enough will
#: always differ from the rest of the map by something.
MIN_EFFECT = 0.5


def _clean(values: pd.Series) -> pd.Series:
    """A categorical layer with its absences removed."""
    s = values.astype("object")
    text = s.map(lambda v: str(v).strip().lower() if v is not None and v == v else "")
    return s.where(~text.isin(ABSENT) & (text != ""), other=np.nan)


def _hypergeom(hits: int, draws: int, successes: int, total: int) -> float:
    """P(at least `hits` of this category in a cluster this size), exactly.

    Fisher's exact test in the one-sided direction that matters. Not a chi-square: cluster counts are
    small and skewed, and a chi-square on an expected count of two is a number with no meaning.
    """
    from scipy.stats import hypergeom
    if min(draws, successes) <= 0 or total <= 0:
        return 1.0
    return float(hypergeom.sf(hits - 1, total, successes, draws))


def _effect(inside: np.ndarray, outside: np.ndarray) -> tuple:
    """(standardised difference, rank-sum p) for a quantity inside a cluster against the rest.

    Cohen's d on a pooled standard deviation, and Mann-Whitney rather than a t-test: fitness scores
    and expression are not normal, and the rank test asks the question actually being made -- are
    the genes in this cluster systematically higher or lower than the rest.
    """
    from scipy.stats import mannwhitneyu
    if len(inside) < 3 or len(outside) < 3:
        return 0.0, 1.0
    sd = np.sqrt((np.var(inside, ddof=1) * (len(inside) - 1)
                  + np.var(outside, ddof=1) * (len(outside) - 1))
                 / max(len(inside) + len(outside) - 2, 1))
    if sd <= 1e-12:                          # every value identical: no difference to test
        return 0.0, 1.0
    d = float((np.mean(inside) - np.mean(outside)) / sd)
    return d, float(mannwhitneyu(inside, outside, alternative="two-sided").pvalue)


def _split(values: np.ndarray, scale: float) -> tuple:
    """How strongly a set of measurements falls into two groups, as (separation, low, high, n_low).

    One-dimensional two-means, which on sorted data is exact rather than iterative: try every cut,
    take the one minimising within-group variance.

    The gap is measured against `scale` -- the spread of the layer across the WHOLE map -- and not
    against the spread within the cluster, which is what it did first and which is wrong in a way
    that manufactures findings. Standardised by its own spread, a cluster whose fitness scores are
    all -3.0 give or take 0.05 splits into "-2.97" and "-3.03" three standard deviations apart, and
    a set of genes that agree about everything is reported as a subdivision. Against the layer's own
    spread, that same cluster scores essentially zero, which is the right answer.

    A dip test would be the textbook answer and is not available without another dependency; on
    fifteen fitness scores it would also be underpowered to the point of decoration.
    """
    v = np.sort(np.asarray(values, dtype=float))
    scale = float(scale) if np.isfinite(scale) and scale > 1e-12 else 1.0
    if len(v) < 6 or np.allclose(v, v[0]):
        return 0.0, float(v[0]) if len(v) else 0.0, float(v[-1]) if len(v) else 0.0, 0
    cuts = np.arange(3, len(v) - 2)
    within = np.array([np.var(v[:c]) * c + np.var(v[c:]) * (len(v) - c) for c in cuts])
    c = int(cuts[int(np.argmin(within))])
    low, high = v[:c], v[c:]
    gap = float((np.mean(high) - np.mean(low)) / scale)
    return gap, float(np.mean(low)), float(np.mean(high)), c


# --------------------------------------------------------------------------- guilt by association
def guilt(nodes: pd.DataFrame, labels: np.ndarray, layer: str, kind: str = None,
          min_cluster: int = MIN_CLUSTER, min_category: int = MIN_CATEGORY,
          min_predictions: int = MIN_PREDICTIONS, min_lift: float = MIN_LIFT,
          max_share: float = MAX_SHARE, used_columns=()) -> pd.DataFrame:
    """Every cluster that says something about the genes in it nobody has measured.

    Returns one row per claim: which cluster, which layer, how strong the evidence, how many genes
    it is about, and which genes those are. Empty where the structure supports no claim at all,
    which is a real answer and the commonest one for a badly chosen embedding.
    """
    labels = np.asarray(labels)
    if layer not in nodes.columns or len(labels) != len(nodes):
        return pd.DataFrame()
    kind = kind or (CONTINUOUS if pd.api.types.is_numeric_dtype(nodes[layer]) else DISCRETE)
    ids = nodes.index.to_numpy()
    genes = (nodes["gene_id"] if "gene_id" in nodes.columns else pd.Series(ids, index=nodes.index))
    rows = []

    if kind == DISCRETE:
        truth = _clean(nodes[layer])
        known = truth.notna().to_numpy()
        total = int(known.sum())
        for k in sorted({int(x) for x in labels if int(x) != NOISE}):
            here = labels == k
            if here.sum() < min_cluster or not (here & known).any():
                continue
            if here.sum() > max_share * len(labels):
                continue
            unlabelled = int((here & ~known).sum())
            counts = truth[here & known].value_counts()
            for category, hits in counts.items():
                if hits < min_category:
                    continue
                background = int((truth == category).sum())
                purity = float(hits / max((here & known).sum(), 1))
                lift = purity / max(background / max(total, 1), 1e-9)
                if lift < min_lift or purity < MIN_PURITY:
                    continue
                rows.append({
                    "kind": "guilt", "layer": layer, "layer_kind": DISCRETE, "cluster": k,
                    "category": str(category), "n_cluster": int(here.sum()),
                    "n_known": int((here & known).sum()), "n_hits": int(hits),
                    "purity": purity,
                    "background": float(background / max(total, 1)),
                    "lift": float(lift),
                    "p": _hypergeom(int(hits), int((here & known).sum()), background, total),
                    "n_predicted": unlabelled,
                    "genes": list(genes[here & ~known].astype(str)[:60]),
                })
    else:
        values = pd.to_numeric(nodes[layer], errors="coerce")
        known = values.notna().to_numpy()
        for k in sorted({int(x) for x in labels if int(x) != NOISE}):
            here = labels == k
            if here.sum() < min_cluster or (here & known).sum() < min_category:
                continue
            if here.sum() > max_share * len(labels):
                continue
            inside = values[here & known].to_numpy(dtype=float)
            outside = values[~here & known].to_numpy(dtype=float)
            d, p = _effect(inside, outside)
            if abs(d) < MIN_EFFECT:
                continue
            rows.append({
                "kind": "guilt", "layer": layer, "layer_kind": CONTINUOUS, "cluster": k,
                "category": "", "n_cluster": int(here.sum()), "n_known": int((here & known).sum()),
                "n_hits": int((here & known).sum()), "purity": np.nan,
                "mean_in": float(np.mean(inside)), "mean_out": float(np.mean(outside)),
                "effect": d, "lift": abs(d), "p": p,
                "n_predicted": int((here & ~known).sum()),
                "genes": list(genes[here & ~known].astype(str)[:60]),
            })

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["q"] = _bh(out["p"].to_numpy())
    out["circular"] = _touches(layer, used_columns)
    return out[out.n_predicted >= min_predictions].sort_values("q").reset_index(drop=True)


# --------------------------------------------------------------------------- layer disagreement
def disagreement(nodes: pd.DataFrame, labels: np.ndarray, layer_a: str, layer_b: str,
                 min_cluster: int = MIN_CLUSTER, min_category: int = MIN_CATEGORY,
                 purity: float = 0.6, used_columns=()) -> pd.DataFrame:
    """Every cluster that agrees about one layer and disagrees about another.

    `layer_a` is the one that has to be homogeneous -- the thing these genes share -- and `layer_b`
    is where they part company. The claim is a subdivision: this compartment contains two kinds of
    gene, and here is the measurement that separates them and the members of each side.

    Not symmetric, and deliberately called twice with the layers swapped when both directions are
    wanted: "same compartment, split fitness" and "same fitness, split compartment" are different
    claims about different biology.
    """
    labels = np.asarray(labels)
    if layer_a not in nodes.columns or layer_b not in nodes.columns or len(labels) != len(nodes):
        return pd.DataFrame()
    genes = (nodes["gene_id"] if "gene_id" in nodes.columns
             else pd.Series(nodes.index.to_numpy(), index=nodes.index))
    a = _clean(nodes[layer_a])
    b_numeric = pd.api.types.is_numeric_dtype(nodes[layer_b])
    b = pd.to_numeric(nodes[layer_b], errors="coerce") if b_numeric else _clean(nodes[layer_b])
    rows = []
    for k in sorted({int(x) for x in labels if int(x) != NOISE}):
        here = (labels == k) & a.notna().to_numpy()
        if here.sum() < min_cluster or (labels == k).sum() > MAX_SHARE * len(labels):
            continue
        counts = a[here].value_counts()
        share = float(counts.iloc[0] / counts.sum())
        if share < purity:                    # not homogeneous in A: nothing to disagree WITH
            continue
        shared = str(counts.index[0])
        agree = here & (a == shared).to_numpy()
        vals = b[agree]
        if b_numeric:
            v = vals.dropna()
            if len(v) < min_category * 2:
                continue
            gap, low, high, n_low = _split(v.to_numpy(dtype=float),
                                           float(np.nanstd(b.to_numpy(dtype=float))))
            if gap < MIN_SPLIT:
                continue
            cut = (low + high) / 2.0
            side = v <= cut
            rows.append({
                "kind": "disagreement", "layer": layer_a, "layer_kind": DISCRETE,
                "other": layer_b, "other_kind": CONTINUOUS, "cluster": k, "category": shared,
                "n_cluster": int(here.sum()), "n_known": int(len(v)), "purity": share,
                "gap": gap, "lift": gap, "mean_low": low, "mean_high": high,
                "n_low": int(side.sum()), "n_high": int((~side).sum()),
                "p": _effect(v[side].to_numpy(dtype=float),
                             v[~side].to_numpy(dtype=float))[1],
                "n_predicted": int(min(side.sum(), (~side).sum())),
                "genes": list(genes[v[side].index].astype(str)[:40]),
                "genes_high": list(genes[v[~side].index].astype(str)[:40]),
            })
        else:
            v = vals.dropna()
            split = v.value_counts()
            split = split[split >= min_category]
            if len(split) < 2:                # everyone agrees in B as well: no disagreement
                continue
            rows.append({
                "kind": "disagreement", "layer": layer_a, "layer_kind": DISCRETE,
                "other": layer_b, "other_kind": DISCRETE, "cluster": k, "category": shared,
                "n_cluster": int(here.sum()), "n_known": int(len(v)), "purity": share,
                # Evenness of the split, in the units of the continuous case: a cluster that is half
                # one thing and half another is a stronger subdivision than one with a lone outlier.
                "gap": float(split.iloc[1] / split.iloc[0]) * float(len(split)),
                "lift": float(split.iloc[1] / split.iloc[0]) * float(len(split)),
                "groups": [str(x) for x in split.index[:6]],
                "counts": [int(x) for x in split.values[:6]],
                "p": _hypergeom(int(split.iloc[1]), int(len(v)),
                                int((b == split.index[1]).sum()), int(b.notna().sum())),
                "n_predicted": int(split.iloc[1]),
                "genes": list(genes[v[v == split.index[1]].index].astype(str)[:40]),
            })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["q"] = _bh(out["p"].to_numpy())
    out["circular"] = _touches(layer_a, used_columns) or _touches(layer_b, used_columns)
    return out.sort_values("q").reset_index(drop=True)


# --------------------------------------------------------------------------- crossed factors
def conjunction(nodes: pd.DataFrame, labels: np.ndarray, layer_a: str, layer_b: str,
                min_cluster: int = MIN_CLUSTER, min_category: int = MIN_CATEGORY,
                min_interaction: float = 1.25, max_share: float = MAX_SHARE,
                used_columns=()) -> pd.DataFrame:
    """Clusters enriched for a combination of two categorical layers.

    ``joint_lift`` says how enriched the pair is over its proteome-wide prevalence. The qualifying
    ``interaction_ratio`` divides that by the stronger single-layer lift: values above one mean the
    pair says more than either margin alone. ``product_ratio`` also records the stricter literal
    joint/(A-lift x B-lift) diagnostic; in a balanced, perfectly resolved A x B grid it is exactly
    one by construction, so using it as the discovery threshold would paradoxically reject the
    canonical crossed-factor fixture.
    """
    labels = np.asarray(labels)
    if (layer_a not in nodes.columns or layer_b not in nodes.columns or len(labels) != len(nodes)
            or pd.api.types.is_numeric_dtype(nodes[layer_a])
            or pd.api.types.is_numeric_dtype(nodes[layer_b])):
        return pd.DataFrame()
    a, b = _clean(nodes[layer_a]), _clean(nodes[layer_b])
    known = a.notna().to_numpy() & b.notna().to_numpy()
    total = int(known.sum())
    if total < min_cluster:
        return pd.DataFrame()
    genes = (nodes["gene_id"] if "gene_id" in nodes.columns
             else pd.Series(nodes.index.to_numpy(), index=nodes.index))
    a_counts, b_counts = a[known].value_counts(), b[known].value_counts()
    joint_counts = pd.DataFrame({"a": a[known], "b": b[known]}).value_counts()
    rows = []
    for cluster in sorted({int(x) for x in labels if int(x) != NOISE}):
        here = labels == cluster
        known_here = here & known
        draws = int(known_here.sum())
        if here.sum() < min_cluster or here.sum() > max_share * len(labels) or draws < min_category:
            continue
        pairs = pd.DataFrame({"a": a[known_here], "b": b[known_here]}).value_counts()
        for (category_a, category_b), hits in pairs.items():
            if hits < min_category:
                continue
            success = int(joint_counts.get((category_a, category_b), 0))
            p_joint_here = hits / draws
            p_a_here = float((a[known_here] == category_a).mean())
            p_b_here = float((b[known_here] == category_b).mean())
            p_a = float(a_counts.get(category_a, 0) / total)
            p_b = float(b_counts.get(category_b, 0) / total)
            p_joint = success / total
            lift_a = p_a_here / max(p_a, 1e-12)
            lift_b = p_b_here / max(p_b, 1e-12)
            joint_lift = p_joint_here / max(p_joint, 1e-12)
            interaction = joint_lift / max(lift_a, lift_b, 1e-12)
            product_ratio = joint_lift / max(lift_a * lift_b, 1e-12)
            unmeasured = here & ~known
            rows.append({
                "kind": "conjunction", "layer": layer_a, "layer_kind": DISCRETE,
                "other": layer_b, "other_kind": DISCRETE, "cluster": cluster,
                "category": str(category_a), "other_category": str(category_b),
                "n_cluster": int(here.sum()), "n_known": draws, "n_hits": int(hits),
                "purity": p_joint_here, "background": p_joint, "lift_a": lift_a,
                "lift_b": lift_b, "joint_lift": joint_lift,
                "interaction_ratio": interaction, "product_ratio": product_ratio,
                "lift": interaction,
                "p": _hypergeom(int(hits), draws, success, total),
                "n_predicted": int(unmeasured.sum()),
                "genes": list(genes[unmeasured].astype(str)[:60]),
            })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    # Correct before filtering: every cluster x A x B cell inspected belongs to the family.
    out["q"] = _bh(out["p"].to_numpy())
    out["circular"] = (_touches(layer_a, used_columns)
                       or _touches(layer_b, used_columns))
    out = out[(out.interaction_ratio >= float(min_interaction))
              & (out.joint_lift >= MIN_LIFT) & (out.purity >= MIN_PURITY)]
    return out.sort_values("q").reset_index(drop=True)


def _mean_pairwise_distance(matrix: np.ndarray) -> float:
    """Mean total-variation distance between rows, or zero with fewer than two rows."""
    matrix = np.asarray(matrix, dtype=float)
    if len(matrix) < 2:
        return 0.0
    return float(np.mean([0.5 * np.abs(matrix[i] - matrix[j]).sum()
                          for i in range(len(matrix)) for j in range(i + 1, len(matrix))]))


def explains_fragmentation(nodes: pd.DataFrame, labels: np.ndarray, layer_a: str, layer_b: str,
                           permutations: int = 200, seed: int = 42,
                           min_cluster: int = MIN_CLUSTER,
                           min_category: int = MIN_CATEGORY) -> pd.DataFrame:
    """Do sibling clusters enriched for the same A category separate on categorical layer B?

    Returns one row per A category with at least two enriched sibling clusters and a final
    ``__run__`` row. Distances are total-variation distances between B-composition vectors. The
    null shuffles B within the A category, preserving both margins while breaking its assignment to
    sibling clusters.
    """
    labels = np.asarray(labels)
    if (layer_a not in nodes.columns or layer_b not in nodes.columns or len(labels) != len(nodes)
            or pd.api.types.is_numeric_dtype(nodes[layer_a])
            or pd.api.types.is_numeric_dtype(nodes[layer_b])):
        return pd.DataFrame()
    a, b = _clean(nodes[layer_a]), _clean(nodes[layer_b])
    a_known = a.notna().to_numpy()
    total = int(a_known.sum())
    rng, rows = np.random.default_rng(seed), []
    for category, background_count in a[a_known].value_counts().items():
        background = float(background_count / max(total, 1))
        siblings = []
        for cluster in sorted({int(x) for x in labels if int(x) != NOISE}):
            here = labels == cluster
            measured = here & a_known
            hits = int((measured & (a == category).to_numpy()).sum())
            if (here.sum() >= min_cluster and here.sum() <= MAX_SHARE * len(labels)
                    and hits >= min_category):
                purity = hits / max(int(measured.sum()), 1)
                if purity >= MIN_PURITY and purity / max(background, 1e-12) >= MIN_LIFT:
                    siblings.append(cluster)
        if len(siblings) < 2:
            continue
        eligible = (a == category).to_numpy() & b.notna().to_numpy()
        categories_b = [value for value, count in b[eligible].value_counts().items()
                        if count >= min_category]
        if len(categories_b) < 2:
            continue

        def composition(values):
            return np.array([[np.mean(values[eligible & (labels == cluster)] == value)
                              if np.any(eligible & (labels == cluster)) else 0.0
                              for value in categories_b] for cluster in siblings])

        observed = _mean_pairwise_distance(composition(b.to_numpy(dtype=object)))
        null = []
        base = b.to_numpy(dtype=object).copy()
        where = np.flatnonzero(eligible)
        for _ in range(max(int(permutations), 1)):
            shuffled = base.copy()
            shuffled[where] = rng.permutation(shuffled[where])
            null.append(_mean_pairwise_distance(composition(shuffled)))
        p = (1.0 + sum(value >= observed for value in null)) / (len(null) + 1.0)
        rows.append({"category": str(category), "n_sibling_clusters": len(siblings),
                     "sibling_clusters": siblings, "observed_distance": observed,
                     "null_mean": float(np.mean(null)), "excess": observed - float(np.mean(null)),
                     "p": p})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["q"] = _bh(out.p.to_numpy())
    explained = out.q <= 0.05
    n_siblings = int(out.n_sibling_clusters.sum())
    n_explained = int(out.loc[explained, "n_sibling_clusters"].sum())
    run = pd.DataFrame([{"category": "__run__", "n_sibling_clusters": n_siblings,
                         "sibling_clusters": [],
                         "observed_distance": float(np.average(
                             out.observed_distance, weights=out.n_sibling_clusters)),
                         "null_mean": float(np.average(out.null_mean,
                                                       weights=out.n_sibling_clusters)),
                         "excess": float(np.average(out.excess,
                                                    weights=out.n_sibling_clusters)),
                         "p": float(out.p.min()), "q": float(out.q.min()),
                         "n_explained_clusters": n_explained,
                         "fraction_explained": n_explained / max(n_siblings, 1)}])
    return pd.concat([out, run], ignore_index=True)


def _touches(layer: str, used_columns) -> bool:
    """Was this layer, or the experiment behind it, among the inputs that built the map?

    Name-prefix matching rather than exact: `lopit_prob_er` built the map and `compartment` is what
    hyperLOPIT calls its answer, so an exact test would call that pair independent. The families are
    the ones `search.SAME_QUANTITY` already names.
    """
    from .search import SAME_QUANTITY
    used = {str(c) for c in used_columns}
    if layer in used:
        return True
    for family in SAME_QUANTITY.values():
        members = {str(m) for m in family}
        if layer in members or any(layer.startswith(str(m)) for m in members):
            if used & members or any(u.startswith(tuple(members)) for u in used):
                return True
    return False


# --------------------------------------------------------------------------- one number, for a search
def yield_score(findings: pd.DataFrame, alpha: float = 0.05) -> float:
    """How much a structure found, as one number an optimiser can climb.

    The product of strength and reach, summed over surviving findings and logged: a configuration
    that produces one overwhelming claim about four genes should not beat one that produces six
    solid claims about ninety. Findings that fail correction contribute nothing rather than a
    little, because a claim that does not survive multiplicity is not weak evidence -- it is a
    finding the run did not make.

    Circular findings are excluded outright. Their strength is arithmetic, and rewarding it would
    drive the optimiser straight into the configuration that puts a layer into the map and reads it
    back out, which is the failure this project has already made once.
    """
    if findings is None or findings.empty:
        return 0.0
    keep = findings[(findings.get("q", 1.0) <= alpha) & (~findings.get("circular", False))]
    if keep.empty:
        return 0.0
    strength = np.log1p(np.clip(keep["lift"].to_numpy(dtype=float), 0.0, 50.0))
    reach = np.log1p(keep["n_predicted"].to_numpy(dtype=float))
    return float(np.sum(strength * reach))
