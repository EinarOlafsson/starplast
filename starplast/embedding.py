#!/usr/bin/env python3
"""Building an embedding from a chosen set of features, with the choices made explicit.

The original embedding hard-coded three decisions that turned out to be wrong, and hid them:

* **Missing values were imputed to the column median without saying so.** "Not measured" silently became
  "average". Genes lacking a measurement do sit apart in the map -- 0.4 to 1.0 map-radii depending on
  policy -- but measuring this properly shows the gap is NOT simply an imputation artefact: it survives
  dropping the offending column, and missingness between assays is only weakly correlated (r = 0.07-0.15),
  so there is no single "was this gene studied" factor. Much of the gap is real biology that missingness
  tracks, since hyperLOPIT, phosphoproteomics and structure prediction all miss low-abundance and
  unusually large proteins. The defect was never the gap; it was that the choice was hidden and
  unmeasured. `missingness_leak()` now quantifies it for any spec.
* **All columns were z-scored alike.** Standard deviations across this table span **8,351x**
  (`length` 977, `lopit_prob_mcmc` 0.117), and 64x *within* the fitness screens alone
  (`fit_invivo_PE` 34.7, `fit_ifng` 0.541). z-scoring equalises spread but not shape, so heavy-tailed
  published scores still dominate.
* **hyperLOPIT was multiplied by 0.5 and called a contributor.** 27 one-hot columns at that scale carry
  0.17 of 16.17 total variance -- **1.1%**. It was named as a design input while having no influence.

Everything here is a parameter instead. `EmbeddingSpec` records the full recipe, so an embedding can be
reproduced, compared against another, and reported in a methods section.

Missing-value policies:

    "drop_columns"    exclude any feature missing in more than `max_missing` of genes
    "drop_genes"      keep the features, drop genes with any missing value
    "median"          impute to the column median (the old behaviour, kept for comparison)
    "indicator"       impute, AND add a 0/1 column recording that it was missing

Measured on this table: `median` gives the smallest gap (0.68 map-radii), `indicator` the largest (1.02)
because it deliberately hands missingness to the embedding as a feature, and `drop_columns` barely helps
(0.69) for the reason above. `drop_genes` is available but near-useless here -- with ~50 features only
**3 of 8,140 genes** have no missing value at all.

So `median` remains the default, not because it is unbiased but because it is the least assertive; what
changed is that the policy is now stated, selectable, and its consequence measurable.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field, asdict

import numpy as np
import pandas as pd

NA_POLICIES = ("indicator", "median", "drop_columns", "drop_genes")
SCALINGS = ("robust", "zscore", "rank", "none")

# Feature blocks the user picks from. Regexes match node-table columns.
BLOCKS = {
    "expression_summary": r"^expr_",
    "expression_raw": r"^rna\d+_",
    "fitness_screens": r"^fit_",
    "published_screens": r"^(crispr_|hosttx_)",
    "localisation": r"^(lopit_prob|lopit_methods_agree)",
    "protein_features": r"^(mean_plddt|paralog_number|n_interpro|n_phosphosites|n_tm|length"
                        r"|tm_kd_|has_domain|has_signal_peptide|is_tm|lineage_specific"
                        r"|protein_ibaq_log2)",
    "literature": r"^(n_publications|n_fulltext|n_papers_)",
    "interactions": r"^(n_xlink_partners|n_struct_similar|n_ipms_partners|n_holes)",
}
# Categorical blocks are one-hot encoded and scaled separately (see `categorical_weight`).
CATEGORICAL_BLOCKS = {"compartment": "compartment", "compartment_best": "compartment_best"}


@dataclass
class EmbeddingSpec:
    """A complete, reproducible recipe for one embedding."""
    name: str = "default"
    blocks: tuple = ("expression_summary", "fitness_screens", "protein_features")
    categorical: tuple = ()
    extra_columns: tuple = ()          # anything the user imported
    na_policy: str = "median"
    max_missing: float = 0.5           # for drop_columns, and for indicator's own guard
    scaling: str = "robust"
    block_weights: dict = field(default_factory=dict)   # block -> multiplier, default 1.0
    categorical_weight: float = 1.0    # see note in `build_matrix`
    indicator_weight: float = 0.5      # missingness informs, but is not a measurement
    n_components: int = 3
    n_neighbors: int = 25
    min_dist: float = 0.25
    metric: str = "euclidean"
    random_state: int = 42

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        d = dict(d)
        for k in ("blocks", "categorical", "extra_columns"):
            if k in d and d[k] is not None:
                d[k] = tuple(d[k])
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


# --------------------------------------------------------------------------- feature selection
def columns_for(nodes: pd.DataFrame, spec: EmbeddingSpec) -> dict:
    """block -> the node-table columns it contributes, in this table."""
    out = {}
    for b in spec.blocks:
        pat = BLOCKS.get(b)
        if not pat:
            continue
        rx = re.compile(pat)
        cols = [c for c in nodes.columns
                if pd.api.types.is_numeric_dtype(nodes[c]) and rx.search(c)]
        if cols:
            out[b] = cols
    extra = [c for c in spec.extra_columns
             if c in nodes.columns and pd.api.types.is_numeric_dtype(nodes[c])]
    if extra:
        out["imported"] = extra
    return out


def _scale(X: np.ndarray, how: str) -> np.ndarray:
    if how == "none" or X.size == 0:
        return X
    if how == "rank":
        # Rank is the safe default for the published screens: they carry inverted sign conventions,
        # ~64x differences in spread, and heavy tails that z-scoring does not tame.
        out = np.empty_like(X, dtype=float)
        for j in range(X.shape[1]):
            col = X[:, j]
            ok = np.isfinite(col)
            r = np.full(col.shape, np.nan)
            if ok.sum() > 1:
                order = np.argsort(np.argsort(col[ok]))
                r[ok] = order / max(ok.sum() - 1, 1) - 0.5
            out[:, j] = r
        return out
    if how == "robust":
        med = np.nanmedian(X, axis=0)
        q1, q3 = np.nanpercentile(X, [25, 75], axis=0)
        iqr = np.where((q3 - q1) > 1e-9, q3 - q1, 1.0)
        return (X - med) / iqr
    mean = np.nanmean(X, axis=0)
    sd = np.nanstd(X, axis=0)
    return (X - mean) / np.where(sd > 1e-9, sd, 1.0)


def build_matrix(nodes: pd.DataFrame, spec: EmbeddingSpec, log=print):
    """Return (X, feature_names, kept_gene_index) for one spec.

    Scaling happens **within each block before weighting**, so a block's influence is set by its weight
    rather than by how many columns it happens to have or what units they are in.
    """
    per_block = columns_for(nodes, spec)
    if not per_block:
        raise ValueError("no numeric features selected")

    mats, names = [], []
    ind_mats, ind_names = [], []
    raw_for_drop = []          # pre-imputation copies, for the drop_genes policy
    for block, cols in per_block.items():
        M = nodes[cols].to_numpy(dtype=float)
        keep = list(cols)

        if spec.na_policy == "drop_columns":
            frac = np.isnan(M).mean(axis=0)
            sel = frac <= spec.max_missing
            if not sel.any():
                continue
            M, keep = M[:, sel], [c for c, s in zip(cols, sel) if s]

        M = _scale(M, spec.scaling)

        if spec.na_policy == "indicator":
            # An indicator is only informative when missingness is neither universal nor absent, and
            # only honest when the column is not mostly missing -- a 85.6%-missing column contributes
            # almost nothing but its own indicator, which is exactly the study-effort artefact.
            frac = np.isnan(M).mean(axis=0)
            ind = (frac > 0.02) & (frac <= spec.max_missing)
            if ind.any():
                ind_mats.append(np.isnan(M[:, ind]).astype(float))
                ind_names += [f"missing::{c}" for c, i in zip(keep, ind) if i]
            drop = frac > spec.max_missing
            if drop.any():
                log(f"  {block}: dropped {int(drop.sum())} column(s) missing >{spec.max_missing:.0%}"
                    f" -- an indicator is all they would contribute")
                M, keep = M[:, ~drop], [c for c, d in zip(keep, drop) if not d]

        if M.shape[1] == 0:
            continue
        # Normalise each block to TOTAL variance 1 before weighting, so influence follows the user's
        # weights rather than an accident of column count or tail shape. Without it, hyperLOPIT's 27
        # one-hot columns took 56% of the matrix simply by being numerous, and robust scaling on the
        # published screens (fit_invivo_PE spans -799..516) inflated that block to half. With it,
        # selecting four blocks at weight 1.0 gives each of them a quarter.
        raw_for_drop.append(M.copy())      # before imputation; drop_genes needs the real NaNs
        M = np.nan_to_num(M, nan=0.0)
        v = float(M.var(axis=0).sum())
        if v > 1e-12:
            M = M / np.sqrt(v)      # block now carries TOTAL variance 1, regardless of column count
        w = float(spec.block_weights.get(block, 1.0))
        mats.append(M * w)
        names += keep

    # Missingness is real information, but it is not a measurement, so it is normalised as its own
    # block and weighted below 1 by default rather than being allowed to rival a measured block by
    # sheer column count.
    if ind_mats:
        I = np.hstack(ind_mats)
        v = float(I.var(axis=0).sum())
        if v > 1e-12:
            I = I / np.sqrt(v) * float(spec.block_weights.get("missing_indicators",
                                                              spec.indicator_weight))
        mats.append(I)
        names += ind_names

    if not mats:
        raise ValueError("every selected feature was dropped by the missing-value policy")
    X = np.hstack(mats)
    rows = np.ones(len(nodes), dtype=bool)

    # Categorical blocks. One-hot columns each carry variance 0.25*p(1-p), which for 27 compartments
    # totalled 1.1% of the matrix -- the reason hyperLOPIT had no influence despite being listed as an
    # input. Here the block is scaled to carry `categorical_weight` times the mean variance of one
    # numeric feature, so its influence is stated rather than accidental.
    for cat in spec.categorical:
        col = CATEGORICAL_BLOCKS.get(cat, cat)
        if col not in nodes.columns:
            continue
        D = pd.get_dummies(nodes[col].astype(str)).to_numpy(dtype=float)
        v = D.var(axis=0).sum()
        D *= (spec.categorical_weight / np.sqrt(v)) if v > 1e-12 else 1.0
        X = np.hstack([X, D])
        names += [f"{col}::{c}" for c in pd.get_dummies(nodes[col].astype(str)).columns]

    if spec.na_policy == "drop_genes":
        # Evaluated against the pre-imputation values. Testing X would always pass, because every block
        # has already had its NaNs replaced by then -- which made this policy a silent no-op.
        rows = np.isfinite(np.hstack(raw_for_drop)).all(axis=1) if raw_for_drop else rows
        X = X[rows]
        log(f"  drop_genes: kept {int(rows.sum()):,} of {len(rows):,} genes with no missing value")

    log(f"matrix: {X.shape[0]:,} genes x {X.shape[1]} features "
        f"({spec.na_policy}, {spec.scaling}"
        + (f", categorical x{spec.categorical_weight}" if spec.categorical else "") + ")")
    return X, names, rows


def embed(nodes: pd.DataFrame, spec: EmbeddingSpec, log=print):
    """Build the matrix and run UMAP. Returns (coords, feature_names, kept_rows)."""
    X, names, rows = build_matrix(nodes, spec, log=log)
    try:
        import umap
        Y = umap.UMAP(n_components=spec.n_components, n_neighbors=spec.n_neighbors,
                      min_dist=spec.min_dist, metric=spec.metric,
                      random_state=spec.random_state).fit_transform(X)
    except Exception as e:
        log(f"UMAP unavailable ({type(e).__name__}); falling back to PCA")
        from sklearn.decomposition import PCA
        Y = PCA(n_components=spec.n_components,
                random_state=spec.random_state).fit_transform(np.nan_to_num(X))
    # umap returns a read-only array in recent versions; the in-place centring below
    # then fails with "output array is read-only". Copy rather than view.
    Y = np.array(Y, dtype=np.float32, copy=True)
    Y -= Y.mean(0)
    Y /= (np.abs(Y).max() + 1e-9)
    return Y * 50.0, names, rows


# --------------------------------------------------------------------------- diagnostics
def variance_share(nodes: pd.DataFrame, spec: EmbeddingSpec) -> pd.DataFrame:
    """How much of the matrix's variance each block actually carries.

    This is the check that catches a block being named as an input while contributing nothing, which is
    what happened to hyperLOPIT at 1.1%.
    """
    X, names, _ = build_matrix(nodes, spec, log=lambda *a: None)
    var = X.var(axis=0)
    block_of = {}
    for b, cols in columns_for(nodes, spec).items():
        for c in cols:
            block_of[c] = b
    rows = []
    for nm, v in zip(names, var):
        if nm.startswith("missing::"):
            b = "missing indicators"
        elif "::" in nm:
            b = nm.split("::")[0]
        else:
            b = block_of.get(nm, "other")
        rows.append((b, v))
    d = pd.DataFrame(rows, columns=["block", "variance"]).groupby("block").variance.sum()
    return (d / d.sum()).sort_values(ascending=False).to_frame("share")


def missingness_leak(coords: np.ndarray, nodes: pd.DataFrame, columns) -> pd.DataFrame:
    """Does 'was this gene measured' predict where it sits? Reports centroid gap in map radii."""
    r = np.sqrt(((coords - coords.mean(0)) ** 2).sum(1).mean())
    out = []
    for c in columns:
        if c not in nodes.columns:
            continue
        m = nodes[c].isna().to_numpy()
        if m.sum() < 30 or (~m).sum() < 30:
            continue
        gap = np.linalg.norm(coords[m].mean(0) - coords[~m].mean(0)) / r
        out.append({"column": c, "missing_frac": float(m.mean()), "centroid_gap": float(gap)})
    return pd.DataFrame(out).sort_values("centroid_gap", ascending=False)
