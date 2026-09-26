"""Is the slot grouping tight enough that holding one thing out holds out its information? Measured.

Every held-out analysis in this project rests on one promise: when a target is withheld, everything
that restates it is withheld too. `search.excluded_detail` keeps that promise with five mechanisms --
the target, measured association above a threshold, declared families, shared provenance, the same
quantity measured another way -- and this module checks them against the data rather than against
intentions.

Three questions, each answered with a number:

1. **How associated can two columns be when they measure DIFFERENT things?** The closure's cutoff
   (0.8) was chosen by hand. Here it is compared with the empirical distribution of association
   between columns from different axes of the slot catalogue -- expression against fitness,
   localization against sequence -- whose upper tail is the most that genuinely different
   measurements reach on this proteome. :func:`calibrate` reports that tail with a bootstrap
   interval, and a column pair above it is as associated as nothing measuring different biology gets.
2. **Does the closure's own measure see every near-copy?** It uses Pearson correlation between
   numeric columns, and a log-transformed copy of a column correlates perfectly on ranks and far less
   on values. :func:`association_matrix` reports both.
3. **After the closure, does anything left predict the target like a copy would?** For each target,
   :func:`residual_leaks` scores every PERMITTED slot on how well it alone predicts the target, out of
   fold, and flags the ones that stand out from the rest -- the undeclared derivations that no
   pairwise statistic can see, like a label that is the argmax of several columns.

Nothing here changes a closure; the audit reports, and fixes are made where the mechanisms are
declared (`datasets` provenance and `derived_from`, `search.SAME_QUANTITY`, slot families), so the
reason for each exclusion stays readable.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

#: Numeric columns with this many distinct values or fewer are treated as categories, as the
#: closure itself treats them.
CATEGORICAL_MAX_UNIQUE = 12
#: Pairs measured on fewer genes than this are not scored: an association on thirty genes is noise
#: with a large number attached.
MIN_OVERLAP = 50
#: Columns that are identifiers or free text, never measurements.
SKIP = ("gene_id", "sequence", "product", "alphafold_accession", "interpro_desc", "interpro_id",
        "pfam_id", "orthogroup", "ec_number")


def is_categorical(s: pd.Series) -> bool:
    """The closure's rule: non-numeric, boolean, or numeric with a dozen values or fewer."""
    if pd.api.types.is_bool_dtype(s) or not pd.api.types.is_numeric_dtype(s):
        return True
    return s.nunique(dropna=True) <= CATEGORICAL_MAX_UNIQUE


def _cramers_v(tab: np.ndarray) -> float:
    from .clustering import _cramers_v
    return float(_cramers_v(tab))


def _eta(cat: pd.Series, num: pd.Series) -> float:
    """Correlation ratio: the share of a numeric column's variance a category explains, sqrt'ed."""
    groups = [num[cat == k].to_numpy(dtype=float) for k in pd.unique(cat)]
    groups = [g for g in groups if len(g)]
    allv = num.to_numpy(dtype=float)
    total = ((allv - allv.mean()) ** 2).sum()
    if total <= 0:
        return float("nan")
    between = sum(len(g) * (g.mean() - allv.mean()) ** 2 for g in groups)
    return float(np.sqrt(between / total))


def measured_columns(nodes: pd.DataFrame, min_values: int = MIN_OVERLAP) -> list:
    """Every column that is a measurement: not an identifier, measured on enough genes, varying."""
    out = []
    for c in nodes.columns:
        if c in SKIP:
            continue
        s = nodes[c]
        known = s.dropna()
        if not pd.api.types.is_numeric_dtype(s) and not pd.api.types.is_bool_dtype(s):
            known = known.astype(str)
            known = known[~known.str.strip().str.lower().isin(
                ("", "nan", "none", "unassigned", "unknown"))]
        if len(known) >= min_values and known.nunique() >= 2:
            out.append(c)
    return out


def association_matrix(nodes: pd.DataFrame, columns=None, max_categories: int = 30) -> tuple:
    """(on ranks, on values): column-by-column association in [0, 1], pairwise complete.

    Numeric pairs are |Spearman| in the first matrix and |Pearson| -- the closure's own measure -- in
    the second; a pair involving a category is Cramer's V or the correlation ratio in both, as the
    closure computes them. Pairs measured on fewer than `MIN_OVERLAP` genes are NaN.
    """
    columns = list(columns) if columns is not None else measured_columns(nodes)
    kinds = {c: is_categorical(nodes[c]) for c in columns}
    num = [c for c in columns if not kinds[c]]
    cat = [c for c in columns if kinds[c]]
    ranks = pd.DataFrame(index=columns, columns=columns, dtype=float)
    values = pd.DataFrame(index=columns, columns=columns, dtype=float)
    if num:
        frame = nodes[num].apply(pd.to_numeric, errors="coerce").astype(float)
        sp = frame.rank().corr(min_periods=MIN_OVERLAP).abs()
        pe = frame.corr(min_periods=MIN_OVERLAP).abs()
        ranks.loc[num, num] = sp.to_numpy()
        values.loc[num, num] = pe.to_numpy()
    labels = {}
    for c in cat:
        s = nodes[c]
        text = s.astype("object").where(s.notna(), None)
        text = text.map(lambda v: None if v is None or str(v).strip().lower() in
                        ("", "nan", "none", "unassigned", "unknown") else str(v))
        labels[c] = text
    for i, a in enumerate(cat):
        la = labels[a]
        for b in cat[i + 1:]:
            lb = labels[b]
            ok = la.notna() & lb.notna()
            if ok.sum() < MIN_OVERLAP or la[ok].nunique() > max_categories \
                    or lb[ok].nunique() > max_categories:
                continue
            v = _cramers_v(pd.crosstab(la[ok], lb[ok]).to_numpy())
            ranks.loc[a, b] = ranks.loc[b, a] = values.loc[a, b] = values.loc[b, a] = v
        for b in num:
            nb = pd.to_numeric(nodes[b], errors="coerce")
            ok = la.notna() & nb.notna()
            if ok.sum() < MIN_OVERLAP or la[ok].nunique() > max_categories:
                continue
            ranks.loc[a, b] = ranks.loc[b, a] = _eta(la[ok], nb[ok].rank())
            values.loc[a, b] = values.loc[b, a] = _eta(la[ok], nb[ok])
    for c in columns:
        ranks.loc[c, c] = values.loc[c, c] = 1.0
    return ranks, values


def slot_of_columns(nodes: pd.DataFrame, organism: str) -> dict:
    """column -> the catalogue slot that declares it (the first, if several do)."""
    from . import slots
    out = {}
    for sl in slots.all_slots(organism):
        for c in slots.declared_columns(nodes, sl):
            out.setdefault(c, sl)
    return out


def relation(a: str, b: str, owner: dict, organism: str | None = None) -> str:
    """How the catalogue relates two columns: the same slot, family, quantity, experiment, or none.

    `different axis` pairs are the reference for what genuinely different measurements look like;
    `same axis` pairs measure related things (two expression series) and are reported, not used to
    calibrate.
    """
    from . import datasets, search
    for fam in search.SAME_QUANTITY.values():
        if search._in_family(a, fam) and search._in_family(b, fam):
            return "same quantity"
    da, db = datasets.provenance(a, organism), datasets.provenance(b, organism)
    if da is not None and db is not None and getattr(da, "key", None) == getattr(db, "key", 0):
        return "same experiment"
    sa, sb = owner.get(a), owner.get(b)
    if sa is None or sb is None:
        return "undeclared"
    if sa.key == sb.key:
        return "same slot"
    if getattr(sa, "target_family", None) and sa.target_family == getattr(sb, "target_family", None):
        return "same family"
    return "same axis" if sa.axis == sb.axis else "different axis"


def pair_table(ranks: pd.DataFrame, values: pd.DataFrame, owner: dict,
               organism: str | None = None) -> pd.DataFrame:
    """Every column pair once, with both associations, their slots and their declared relation."""
    cols = list(ranks.index)
    iu = np.triu_indices(len(cols), 1)
    r = ranks.to_numpy()[iu]
    v = values.to_numpy()[iu]
    keep = np.isfinite(r) | np.isfinite(v)
    a = np.array(cols)[iu[0]][keep]
    b = np.array(cols)[iu[1]][keep]
    out = pd.DataFrame({"a": a, "b": b, "rank_association": r[keep], "value_association": v[keep]})
    out["association"] = out[["rank_association", "value_association"]].max(axis=1)
    out["slot_a"] = [getattr(owner.get(x), "key", "") for x in out["a"]]
    out["slot_b"] = [getattr(owner.get(x), "key", "") for x in out["b"]]
    out["relation"] = [relation(x, y, owner, organism) for x, y in zip(out["a"], out["b"])]
    return out


def calibrate(pairs: pd.DataFrame, quantile: float = 0.999, n_boot: int = 1000,
              seed: int = 0) -> dict:
    """The most two columns from different axes are associated, with a bootstrap interval.

    The `quantile` of the association between columns whose slots sit on different axes of the
    catalogue: pairs measuring genuinely different things. Bootstrapped over SLOT PAIRS rather than
    column pairs, because the columns of one slot are not independent draws and resampling them
    as if they were would make the interval falsely narrow.
    """
    neg = pairs[pairs["relation"] == "different axis"]
    vals = neg["association"].to_numpy(dtype=float)
    vals = vals[np.isfinite(vals)]
    if not len(vals):
        return {"threshold": float("nan"), "low": float("nan"), "high": float("nan"), "n": 0}
    keys = (neg["slot_a"] + "|" + neg["slot_b"]).to_numpy()
    groups = {k: neg["association"].to_numpy()[keys == k] for k in np.unique(keys)}
    names = list(groups)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(int(n_boot)):
        pick = rng.choice(len(names), size=len(names), replace=True)
        sample = np.concatenate([groups[names[i]] for i in pick])
        sample = sample[np.isfinite(sample)]
        boots.append(np.quantile(sample, quantile))
    return {"threshold": float(np.quantile(vals, quantile)),
            "low": float(np.percentile(boots, 2.5)), "high": float(np.percentile(boots, 97.5)),
            "n": int(len(vals)), "slot_pairs": len(names), "quantile": quantile}


def closure_gaps(nodes: pd.DataFrame, pairs: pd.DataFrame, threshold: float,
                 targets=None) -> pd.DataFrame:
    """Pairs the closure lets through although they are more associated than `threshold`.

    For each target (every measured column by default), the columns `search.excluded_for` removes
    are compared with the columns associated with it at or above `threshold`; what is associated
    and NOT removed is a gap -- a copy of the target a held-out analysis would still be allowed to
    read.
    """
    from . import search
    strong = pairs[pairs["association"] >= threshold]
    partners: dict = {}
    for r in strong.itertuples():
        partners.setdefault(r.a, []).append((r.b, r.association, r.rank_association,
                                             r.value_association))
        partners.setdefault(r.b, []).append((r.a, r.association, r.rank_association,
                                             r.value_association))
    rows = []
    for t in (targets if targets is not None else sorted(partners)):
        if t not in partners:
            continue
        banned = search.excluded_for(nodes, t)
        for other, assoc, ra, va in partners[t]:
            if other not in banned:
                rows.append({"target": t, "column": other, "association": assoc,
                             "rank_association": ra, "value_association": va})
    return pd.DataFrame(rows).sort_values("association", ascending=False) if rows else \
        pd.DataFrame(columns=["target", "column", "association", "rank_association",
                              "value_association"])


def _predictability(y: pd.Series, X: np.ndarray, categorical: bool, seed: int = 0,
                    folds: int = 5) -> float:
    """Out-of-fold: Cohen's kappa for a category, correlation of prediction for a number."""
    from sklearn.linear_model import LogisticRegression, Ridge
    from sklearn.metrics import cohen_kappa_score
    known = np.flatnonzero(y.notna().to_numpy())
    if len(known) < 60 or X.shape[1] == 0:
        return float("nan")
    yk = y.iloc[known]
    if categorical:
        counts = yk.value_counts()
        yk = yk.where(yk.isin(counts[counts >= 10].index))
        keep = yk.notna().to_numpy()
        known, yk = known[keep], yk[keep]
        if yk.nunique() < 2:
            return float("nan")
    rng = np.random.default_rng(seed)
    fold = rng.integers(0, folds, size=len(known))
    pred = np.empty(len(known), dtype=object if categorical else float)
    for f in range(folds):
        tr, te = fold != f, fold == f
        if te.sum() == 0 or (categorical and pd.Series(yk.to_numpy()[tr]).nunique() < 2):
            continue
        if categorical:
            m = LogisticRegression(max_iter=300, C=1.0, class_weight="balanced")
            m.fit(X[known[tr]], yk.to_numpy()[tr].astype(str))
        else:
            m = Ridge(alpha=1.0).fit(X[known[tr]], yk.to_numpy(dtype=float)[tr])
        pred[te] = m.predict(X[known[te]])
    if categorical:
        return float(cohen_kappa_score(yk.to_numpy().astype(str), pred.astype(str)))
    ok = np.isfinite(pred.astype(float))
    if ok.sum() < 20:
        return float("nan")
    return float(np.corrcoef(pred.astype(float)[ok], yk.to_numpy(dtype=float)[ok])[0, 1])


def residual_leaks(nodes: pd.DataFrame, target: str, organism: str = "Tg", z: float = 4.0,
                   floor: float = 0.5, log=None) -> pd.DataFrame:
    """Every slot the closure PERMITS for `target`, scored on how well it alone predicts it.

    A slot that restates the target -- even as a joint function no pairwise statistic can see --
    predicts it far better than the other permitted slots do. Flagged when it is BOTH a robust
    outlier -- more than `z` median absolute deviations above the median slot -- AND at least
    `floor` (kappa or correlation): halfway to a copy. The outlier test alone flags ordinary
    biology whenever every other slot sits near zero -- expression predicts compartment at kappa
    0.1, six deviations above a median of 0.01 -- and the floor alone would flag any strong but
    genuine predictor. A flag says "look"; the table says what to look at.
    """
    from . import search, slots
    from .strategies import Context
    ctx = Context(nodes, graph={}, organism=organism)
    banned = search.excluded_for(nodes, target)
    y = nodes[target]
    categorical = is_categorical(y)
    if categorical:
        y = ctx.truth(target)
    else:
        y = pd.to_numeric(y, errors="coerce")
    rows = []
    for sl in slots.all_slots(organism):
        cols = [c for c in slots.declared_columns(nodes, sl, numeric_only=True)
                if c in nodes.columns and c not in SKIP]
        if not cols or set(cols) & banned:
            continue
        X = ctx.matrix([c for c in cols if pd.api.types.is_numeric_dtype(nodes[c])
                        or pd.api.types.is_bool_dtype(nodes[c])])
        if X.shape[1] == 0:
            continue
        if log:
            log(f"{target}: {sl.key}")
        rows.append({"target": target, "slot": sl.key, "axis": sl.axis, "columns": len(cols),
                     "predictability": _predictability(y, X, categorical)})
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    p = out["predictability"].to_numpy(dtype=float)
    med = np.nanmedian(p)
    mad = np.nanmedian(np.abs(p - med)) or 1e-9
    out["robust_z"] = (p - med) / (1.4826 * mad)
    out["flag"] = (out["robust_z"] > z) & (out["predictability"] >= floor)
    return out.sort_values("predictability", ascending=False).reset_index(drop=True)
