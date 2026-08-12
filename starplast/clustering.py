#!/usr/bin/env python3
"""Clustering an embedding, and asking what — if anything — the clusters correspond to.

Two parts. `walk_*` searches hyperparameters and scores each setting. `battery` then asks, for every
feature, whether it separates the clusters, and reports it in the form a biologist actually wants:
*of all known tachyzoite genes, what fraction landed in cluster 1, and of cluster 1, what fraction are
tachyzoite?*

**The circularity guard is the point of this module.** If a feature fed the embedding, then finding that
it separates the resulting clusters is guaranteed and means nothing — the clusters were built to separate
it. Every result is therefore labelled `used` or `held_out`, and only `held_out` features are evidence.
Reporting a `used` feature as a discovery is the single easiest way to produce a confident artefact here,
so `battery` refuses to sort the two together.

Scores are chosen so features of different types can be ranked against each other:

* categorical feature x cluster -> **Cramer's V** (0-1) from the contingency table, plus per-cluster
  precision and recall and a Fisher exact test
* continuous feature x cluster -> **epsilon-squared** (0-1) from Kruskal-Wallis, plus per-cluster medians

Both are bounded, both mean "share of this feature's variation explained by cluster membership", and both
get Benjamini-Hochberg correction across the whole battery, which for ~90 features and k clusters is a
lot of tests.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

NOISE = -1


# --------------------------------------------------------------------------- hyperparameter walks
def _score(X: np.ndarray, labels: np.ndarray) -> dict:
    from sklearn.metrics import silhouette_score
    mask = labels != NOISE
    k = len(set(labels[mask]))
    out = {"n_clusters": k, "noise_frac": float((~mask).mean()),
           "largest_frac": 0.0, "silhouette": np.nan}
    if k:
        _, counts = np.unique(labels[mask], return_counts=True)
        out["largest_frac"] = float(counts.max() / len(labels))
    # Silhouette is computed on clustered points only; including noise makes every setting look bad.
    if 2 <= k <= 200 and mask.sum() > k + 10:
        try:
            out["silhouette"] = float(silhouette_score(X[mask], labels[mask]))
        except Exception:
            pass
    return out


def walk_dbscan(X: np.ndarray, eps_values=None, min_samples_values=(5, 10, 20, 50),
                log=print) -> pd.DataFrame:
    """Grid over eps x min_samples. eps defaults to a spread around the median 4-NN distance."""
    from sklearn.cluster import DBSCAN
    from sklearn.neighbors import NearestNeighbors
    if eps_values is None:
        d, _ = NearestNeighbors(n_neighbors=5).fit(X).kneighbors(X)
        base = float(np.median(d[:, -1]))
        eps_values = [round(base * m, 4) for m in (0.5, 0.75, 1.0, 1.5, 2.0, 3.0)]
        log(f"eps grid from median 5-NN distance {base:.3f}: {eps_values}")
    rows = []
    for eps in eps_values:
        for ms in min_samples_values:
            lab = DBSCAN(eps=eps, min_samples=ms).fit_predict(X)
            rows.append({"algorithm": "dbscan", "eps": eps, "min_samples": ms, **_score(X, lab)})
    return pd.DataFrame(rows).sort_values("silhouette", ascending=False, na_position="last")


def walk_hdbscan(X: np.ndarray, min_cluster_sizes=(10, 20, 50, 100),
                 min_samples_values=(None, 5, 15), log=print) -> pd.DataFrame:
    """HDBSCAN handles the variable density UMAP produces far better than a single global eps."""
    try:
        from sklearn.cluster import HDBSCAN
        def run(mcs, ms):
            return HDBSCAN(min_cluster_size=mcs, min_samples=ms).fit_predict(X)
    except ImportError:                                        # pragma: no cover
        try:
            import hdbscan as _h
            def run(mcs, ms):
                return _h.HDBSCAN(min_cluster_size=mcs, min_samples=ms).fit_predict(X)
        except ImportError:
            log("HDBSCAN unavailable (needs scikit-learn >= 1.3 or the hdbscan package)")
            return pd.DataFrame()
    rows = []
    for mcs in min_cluster_sizes:
        for ms in min_samples_values:
            lab = run(mcs, ms)
            rows.append({"algorithm": "hdbscan", "min_cluster_size": mcs, "min_samples": ms,
                         **_score(X, lab)})
    return pd.DataFrame(rows).sort_values("silhouette", ascending=False, na_position="last")


def cluster(X: np.ndarray, algorithm="hdbscan", **kw) -> np.ndarray:
    """Cluster an embedding, returning a label per point with -1 for noise.

    HDBSCAN by default: it allows clusters of differing density and calls the rest noise, which suits
    a proteome where most genes belong to no tight group. Calling that noise is the honest answer.
    """
    if algorithm == "dbscan":
        from sklearn.cluster import DBSCAN
        return DBSCAN(eps=kw.get("eps", 0.5), min_samples=kw.get("min_samples", 10)).fit_predict(X)
    try:
        from sklearn.cluster import HDBSCAN
        return HDBSCAN(min_cluster_size=kw.get("min_cluster_size", 25),
                       min_samples=kw.get("min_samples")).fit_predict(X)
    except ImportError:                                        # pragma: no cover
        import hdbscan as _h
        return _h.HDBSCAN(min_cluster_size=kw.get("min_cluster_size", 25),
                          min_samples=kw.get("min_samples")).fit_predict(X)


# --------------------------------------------------------------------------- statistics
def _bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg. The battery runs hundreds of tests; uncorrected p-values would be theatre."""
    p = np.asarray(p, dtype=float)
    ok = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    if not ok.any():
        return q
    v = p[ok]
    order = np.argsort(v)
    ranked = v[order]
    n = len(v)
    adj = ranked * n / np.arange(1, n + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


def _cramers_v(table: np.ndarray) -> float:
    from scipy.stats import chi2_contingency
    if table.shape[0] < 2 or table.shape[1] < 2 or table.sum() == 0:
        return np.nan
    try:
        chi2 = chi2_contingency(table)[0]
    except ValueError:
        return np.nan
    n = table.sum()
    return float(np.sqrt((chi2 / n) / max(min(table.shape) - 1, 1)))


def categorical_feature(labels: np.ndarray, values: pd.Series, min_count=5) -> tuple:
    """Return (summary_row, per_cluster_rows) for one categorical feature."""
    from scipy.stats import fisher_exact
    ok = (labels != NOISE) & values.notna().to_numpy()
    if ok.sum() < 20:
        return None, []
    lab, val = labels[ok], values[ok].astype(str).to_numpy()
    cats = [c for c, n in zip(*np.unique(val, return_counts=True)) if n >= min_count]
    clusters = sorted(set(lab))
    if len(cats) < 2 or len(clusters) < 2:
        return None, []
    table = np.array([[int(((lab == cl) & (val == c)).sum()) for c in cats] for cl in clusters])
    v = _cramers_v(table)
    rows = []
    for i, cl in enumerate(clusters):
        for j, c in enumerate(cats):
            a = table[i, j]                       # in cluster, is category
            b = table[i].sum() - a                # in cluster, other category
            cc = table[:, j].sum() - a            # other cluster, is category
            d = table.sum() - a - b - cc
            try:
                orr, p = fisher_exact([[a, b], [cc, d]], alternative="greater")
            except ValueError:
                orr, p = np.nan, np.nan
            rows.append({"cluster": int(cl), "category": c,
                         "n_in_cluster": int(a),
                         "precision": a / max(table[i].sum(), 1),   # of this cluster, how much is c
                         "recall": a / max(table[:, j].sum(), 1),   # of all c, how much is here
                         "odds_ratio": float(orr), "p": float(p)})
    return {"score": v, "score_type": "cramers_v", "n": int(ok.sum())}, rows


def continuous_feature(labels: np.ndarray, values: pd.Series) -> tuple:
    """Return (summary_row, per_cluster_rows) for one continuous feature."""
    from scipy.stats import kruskal
    v = pd.to_numeric(values, errors="coerce")
    ok = (labels != NOISE) & v.notna().to_numpy()
    if ok.sum() < 30:
        return None, []
    lab, val = labels[ok], v[ok].to_numpy()
    groups = [val[lab == cl] for cl in sorted(set(lab)) if (lab == cl).sum() >= 5]
    if len(groups) < 2:
        return None, []
    try:
        h, p = kruskal(*groups)
    except ValueError:
        return None, []
    n, k = len(val), len(groups)
    eps2 = float((h - k + 1) / max(n - k, 1))          # epsilon-squared, bounded effect size
    overall = np.median(val)
    rows = [{"cluster": int(cl), "n_in_cluster": int((lab == cl).sum()),
             "median": float(np.median(val[lab == cl])),
             "direction": "high" if np.median(val[lab == cl]) > overall else "low"}
            for cl in sorted(set(lab)) if (lab == cl).sum() >= 5]
    return ({"score": max(eps2, 0.0), "score_type": "epsilon_squared",
             "n": int(ok.sum()), "p": float(p)}, rows)


# --------------------------------------------------------------------------- the battery
def _association_with_inputs(nodes: pd.DataFrame, used: set, max_categories=30) -> dict:
    """Strongest association between each column and any embedding input (0-1).

    A hard threshold alone is not enough. `lopit_mcmc` is a different inference over the same
    experiment as the input `compartment` -- 0.74 associated, below any sensible cutoff, yet
    reporting it as a discovery is close to circular. Rather than tune the cutoff, the number is
    carried alongside every result so the reader can see how independent a finding really is.
    """
    return _assoc_impl(nodes, used, max_categories)


def _derived_from(nodes: pd.DataFrame, used: set, threshold=0.95, max_categories=30) -> set:
    """Columns that are near-perfect restatements of something the embedding already used.

    Naming the used columns is not enough. `compartment` fed the embedding, but `lopit_map` is its exact
    twin and `lopit_unified`, `compartment_best`, `compartment_source` are derivations of it -- all of
    which the first version of this battery happily reported as top held-out discoveries, at V = 0.96.
    A guard that can be defeated by renaming a column is not a guard, so association is measured rather
    than assumed.
    """
    a = _assoc_impl(nodes, used, max_categories)
    return {c for c, v in a.items() if v >= threshold}


def _correlation_ratio(cat: pd.Series, num: pd.Series, max_categories=30):
    """eta: how much of a continuous column's variance is explained by a categorical one, in [0, 1].

    The categorical-versus-continuous counterpart of Cramer's V, and on the same scale, so one
    threshold applies to both. Returns None where the comparison is not meaningful -- too many
    categories, or a constant column, whose variance ratio is 0/0 rather than 0.
    """
    g = num.groupby(cat.astype(str))
    k = g.ngroups
    if k < 2 or k > max_categories:
        return None
    total_var = float(np.var(num.to_numpy(), ddof=0))
    if not np.isfinite(total_var) or total_var <= 0:
        return None
    grand = float(num.mean())
    between = float(sum(len(v) * (float(v.mean()) - grand) ** 2 for _, v in g))
    eta_sq = between / (total_var * len(num))
    return float(np.sqrt(max(0.0, min(1.0, eta_sq))))


def _assoc_impl(nodes: pd.DataFrame, used: set, max_categories=30) -> dict:
    out = {}
    for u in used:
        if u not in nodes.columns:
            continue
        su = nodes[u]
        u_cat = not pd.api.types.is_numeric_dtype(su) or su.nunique(dropna=True) <= 12
        for c in nodes.columns:
            if c in used:
                continue
            sc = nodes[c]
            ok = su.notna() & sc.notna()
            if ok.sum() < 30:
                continue
            c_cat = not pd.api.types.is_numeric_dtype(sc) or sc.nunique(dropna=True) <= 12
            try:
                if u_cat and c_cat:
                    if su[ok].nunique() > max_categories or sc[ok].nunique() > max_categories:
                        continue
                    tab = pd.crosstab(su[ok].astype(str), sc[ok].astype(str)).to_numpy()
                    v = _cramers_v(tab)
                elif not u_cat and not c_cat:
                    v = abs(float(np.corrcoef(su[ok].astype(float), sc[ok].astype(float))[0, 1]))
                else:
                    # One categorical, one continuous. This branch used to `continue`, which made the
                    # guard structurally blind to every numeric column: a categorical label computed
                    # FROM continuous columns -- exactly what stage_enriched_derived is -- scored no
                    # association with any of its own sources, because the comparison was never made.
                    # The correlation ratio is the right measure here and is on the same 0-1 scale as
                    # Cramer's V, so the single reported number stays comparable across column kinds.
                    cat, num = (su, sc) if u_cat else (sc, su)
                    v = _correlation_ratio(cat[ok], num[ok].astype(float), max_categories)
                    if v is None:
                        continue
                if np.isfinite(v):
                    out[c] = max(out.get(c, 0.0), float(v))
            except Exception:
                continue
    return out


def battery(nodes: pd.DataFrame, labels: np.ndarray, used_features=(), features=None,
            max_categories=30, guard_derived=True, log=print) -> tuple:
    """Test every feature against the clustering.

    `used_features` are the node-table columns that fed the embedding. Those, and anything that restates
    them (see `_derived_from`), are reported but marked `used` / `derived`, because a feature the clusters
    were built from separates them by construction. Only `held_out` results are evidence.
    """
    used = set(used_features)
    assoc = _association_with_inputs(nodes, used, max_categories) if guard_derived else {}
    derived = {c for c, v in assoc.items() if v >= 0.95}
    if derived:
        log(f"battery: {len(derived)} column(s) are near-perfect restatements of an embedding input "
            f"and are marked derived, not held out: {', '.join(sorted(derived)[:6])}"
            + (" ..." if len(derived) > 6 else ""))
    cols = list(features) if features is not None else [
        c for c in nodes.columns if c not in ("gene_id", "sequence", "product")]

    summaries, details = [], []
    for c in cols:
        s = nodes[c]
        if pd.api.types.is_numeric_dtype(s) and s.nunique(dropna=True) > 12:
            sm, rows = continuous_feature(labels, s)
        else:
            if s.nunique(dropna=True) > max_categories:
                continue
            sm, rows = categorical_feature(labels, s)
        if sm is None:
            continue
        sm.update(feature=c, evidence="used" if c in used
                  else "derived" if c in derived else "held_out",
                  assoc_with_input=round(float(assoc.get(c, 0.0)), 3))
        summaries.append(sm)
        for r in rows:
            r.update(feature=c, evidence=sm["evidence"])
        details += rows

    S = pd.DataFrame(summaries)
    D = pd.DataFrame(details)
    if not S.empty and "p" in S:
        S["q"] = _bh(S["p"].to_numpy())
    if not D.empty and "p" in D:
        D["q"] = _bh(D["p"].to_numpy())
    if not S.empty:
        order = {"held_out": 0, "derived": 1, "used": 2}
        S = S.assign(_o=S.evidence.map(order)).sort_values(
            ["_o", "score"], ascending=[True, False]).drop(columns="_o")
        n = S.evidence.value_counts()
        log(f"battery: {len(S)} features tested -- {n.get('held_out',0)} held out "
            f"(evidence), {n.get('derived',0)} derived from an input, "
            f"{n.get('used',0)} used to build the embedding")
    return S, D


def per_category(detail: pd.DataFrame, evidence: str = "held_out",
                 min_in_cluster: int = 5) -> pd.DataFrame:
    """The best cluster for each CATEGORY of each held-out feature: precision, recall, F1.

    The battery reports one number per feature -- Cramer's V over the whole contingency table -- and
    a single number over 27 compartments hides the thing worth knowing. V = 0.2 is the same value
    whether every compartment is weakly smeared across every cluster or one falls out cleanly and
    the rest are noise, and those are entirely different findings. Only the per-category table
    separates them, and it uses the same "best single cluster per label" convention as
    `search.score_recovery`, so the two cannot disagree about what recovery means.

    F1 is computed here rather than read off the battery, whose rows are per (cluster, category)
    pair: the row that matters for a category is its best one.

    `min_in_cluster` is the singleton guard in its per-category form, and it is not optional.
    Ranked by lift without it, the top of this table was a cluster holding ONE gene of a rare
    category at 5x enrichment -- true, meaningless, and indistinguishable at a glance from a real
    result. The same floor `objectives.score` applies for the same reason.
    """
    cols = ["feature", "category", "cluster", "n_in_cluster", "precision", "prevalence", "lift",
            "recall", "f1", "q"]
    if detail is None or detail.empty or "category" not in detail.columns:
        return pd.DataFrame(columns=cols)
    d = detail[detail.evidence == evidence] if "evidence" in detail.columns else detail
    d = d[d.precision.notna() & d.recall.notna()]
    if d.empty:
        return pd.DataFrame(columns=cols)
    # Prevalence and lift are computed over every row, including the ones too small to be reported:
    # they describe how common the category is, which does not depend on which rows are shown.
    keep = d.n_in_cluster >= min_in_cluster
    denom = (d.precision + d.recall).to_numpy(dtype=float)
    f1 = np.divide(2 * d.precision.to_numpy(dtype=float) * d.recall.to_numpy(dtype=float),
                   denom, out=np.zeros(len(d)), where=denom > 0)
    # How common the category is across everything this feature scored, and how much more of it the
    # cluster holds than the map does. Without this the table's top row is whatever class dominates:
    # a cluster holding 90% of the genes "recovers" a 90%-prevalent label at F1 0.95 while telling
    # you nothing, which is the trivial partition this project has already been caught by once.
    # Lift 1.0 means the cluster is no more that category than the proteome is.
    total = d.groupby("feature").n_in_cluster.transform("sum").to_numpy(dtype=float)
    n_cat = d.groupby(["feature", "category"]).n_in_cluster.transform("sum").to_numpy(dtype=float)
    prevalence = np.divide(n_cat, total, out=np.full(len(d), np.nan), where=total > 0)
    lift = np.divide(d.precision.to_numpy(dtype=float), prevalence,
                     out=np.full(len(d), np.nan), where=prevalence > 0)
    d = d.assign(f1=f1, prevalence=prevalence, lift=lift)[keep.to_numpy()]
    if d.empty:
        return pd.DataFrame(columns=cols)
    best = (d.sort_values("f1", ascending=False)
             .groupby(["feature", "category"], as_index=False, sort=False).head(1))
    for c in cols:
        if c not in best.columns:
            best[c] = np.nan
    return best[cols].sort_values(["lift", "f1"], ascending=False).reset_index(drop=True)


def describe(summary: pd.DataFrame, detail: pd.DataFrame, top=6, min_score=0.15,
             q_max=0.05) -> list:
    """Turn the battery into the sentences a reader wants, held-out features only."""
    out = []
    if summary.empty:
        return out
    S = summary[(summary.evidence == "held_out") & (summary.score >= min_score)].head(top)
    for r in S.itertuples():
        if r.score_type == "cramers_v":
            d = detail[(detail.feature == r.feature) & (detail.evidence == "held_out")]
            d = d[(d.q <= q_max) & (d.n_in_cluster >= 5)].sort_values("precision", ascending=False)
            for x in d.head(2).itertuples():
                out.append(
                    f"cluster {x.cluster}: {x.precision:.0%} of its genes are "
                    f"{r.feature}={x.category} (that is {x.recall:.0%} of all such genes, "
                    f"OR {x.odds_ratio:.1f}, q={x.q:.1e}) [V={r.score:.2f}]")
        else:
            d = detail[detail.feature == r.feature].sort_values("median", ascending=False)
            if len(d) >= 2:
                hi, lo = d.iloc[0], d.iloc[-1]
                out.append(
                    f"{r.feature}: high in cluster {int(hi.cluster)} (median {hi['median']:.2f}), "
                    f"low in cluster {int(lo.cluster)} (median {lo['median']:.2f}) "
                    f"[eps2={r.score:.2f}, q={getattr(r, 'q', float('nan')):.1e}]")
    return out
