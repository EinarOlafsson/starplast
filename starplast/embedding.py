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
    "median"          impute to the column median (the old behavior, kept for comparison)
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

#: How the matrix is taken down to a few dimensions. Three answers to different questions, and the
#: differences matter more than "which one looks nicer".
#:
#: umap -- keeps neighbourhoods and makes some attempt at the space between them. The default, and
#: what every published figure from this project used.
#: tsne -- keeps neighbourhoods and abandons everything else. Distances BETWEEN clusters in a t-SNE
#: are not interpretable at all, which is a real cost here because the map is read as a space; what
#: it buys is that within-cluster structure is usually cleaner, so a clustering run on it splits
#: groups UMAP merges. Worth having for exactly that.
#: pca -- keeps variance and nothing else. Nobody should read biology off a PCA of this matrix, and
#: it belongs here anyway: it is the honest baseline, and a claim that survives on a PCA is a claim
#: that did not need the embedding to be true.
METHODS = ("umap", "tsne", "pca")

# Feature blocks the user picks from. Regexes match node-table columns.
BLOCKS = {
    "structure_af3": r"^af3_",
    "expression_summary": r"^expr_",
    "expression_raw": r"^rna\d+_",
    # These six assay families already ship in nodes.parquet.  Keeping each biological question in
    # its own block is the first executable part of the slot model: an in-vivo brain transcriptome
    # is not averaged with an alkaline-stress series merely because both arrived through
    # expression.py.  The two invivo patterns are split again because tissue-culture tachyzoites and
    # infected brain are different conditions and therefore different slots.
    "transcription_tachyzoite_comparator": r"^invivo_TZ_",
    "transcription_in_vivo_brain": r"^invivo_(?:WholeBrain|BZ_)",
    "transcription_stress_conversion": r"^stress_",
    "transcription_tf_chromatin_perturbation": r"^morc_",
    "fitness_screens": r"^fit_",
    "published_screens": r"^(crispr_|hosttx_)",
    "localization": r"^(lopit_prob|lopit_methods_agree)",
    "protein_abundance_perturbation": r"^proteome_",
    "protein_abundance_oocyst": r"^oocyst_",
    "phosphorylation_quantitative": r"^phospho_",
    "protein_features": r"^(mean_plddt|paralog_number|n_interpro|n_phosphosites|n_tm|length"
                        r"|tm_kd_|has_domain|has_signal_peptide|is_tm|lineage_specific"
                        r"|protein_ibaq_log2)",
    "literature": r"^(n_publications|n_fulltext|n_papers_)",
    "interactions": r"^(n_xlink_partners|n_struct_similar|n_ipms_partners|n_holes)",
}
# Biological-question slots are the primary blocks for new maps. The regex blocks above remain as
# recipe compatibility for maps saved before v0.31; they are no longer what the optimiser offers.
from .slots import all_slots as _all_slots
SLOT_BLOCKS = {slot.key: slot for slot in _all_slots("Tg")
               if slot.unit == "gene" and slot.patterns and slot.role == "feature"}
BLOCKS.update({key: "" for key in SLOT_BLOCKS})
#: Blocks that have been renamed: the spelling a recipe may carry -> what it is called now. A recipe
#: is a promise that a run can be rebuilt, and the embeddings saved before the rename name their
#: blocks as they were spelled when they were computed. An unrecognised block does not fail loudly --
#: `columns_for` skips it -- so a stale recipe would quietly rebuild a DIFFERENT embedding under the
#: same name. Translating on the way in is what keeps "saved recipe" and "rebuildable" the same thing.
RENAMED_BLOCKS = {"localisation": "localization"}

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
    method: str = "umap"
    n_components: int = 3
    n_neighbors: int = 25
    min_dist: float = 0.25
    # t-SNE's own two. Perplexity is its n_neighbors and behaves like one; exaggeration decides how
    # hard clusters are pushed apart early on, which is why a t-SNE looks more separated than it is.
    perplexity: float = 30.0
    early_exaggeration: float = 12.0
    metric: str = "euclidean"
    random_state: int = 42

    def __post_init__(self):
        # Every route a spec arrives by passes through here: a stored recipe, a row of a results
        # table saved last week, a hand-written call. So this is where a renamed block is translated,
        # once, rather than in each of them.
        self.blocks = tuple(RENAMED_BLOCKS.get(b, b) for b in self.blocks)
        if self.block_weights:
            self.block_weights = {RENAMED_BLOCKS.get(k, k): v for k, v in self.block_weights.items()}

    def to_dict(self):
        """The recipe as a plain dict, for storing beside an embedding."""
        return asdict(self)

    @classmethod
    def from_dict(cls, d):
        """Rebuild a spec from a stored recipe."""
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
        if b in SLOT_BLOCKS:
            from .slots import source_columns
            cols = list(source_columns(nodes, SLOT_BLOCKS[b]))
            if cols:
                out[b] = cols
            continue
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
        # On the CPU on purpose: the device version was six times slower on the real matrix and
        # disagreed by 1.4e-4 because float32 reorders near-ties. See the note in `gpu.py`.
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


def as_text(s):
    """A plain string Series with missing values as "", whatever pandas version is installed.

    `.astype(str)` keeps NA under pandas 3 and produces the literal "nan" under pandas 2, so one-hot
    encoding a column with missing values gives a different set of dummy columns on each. Missing
    becomes "", which is already an absence label everywhere else in this project.
    """
    return s.astype("object").where(s.notna(), "").astype(str)


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
        if block in SLOT_BLOCKS:
            from .slots import resolve
            resolved = resolve(nodes, SLOT_BLOCKS[block])
            M = resolved.values.to_numpy(dtype=float)
            keep = list(resolved.values.columns)
        else:
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
        # Normalize each block to TOTAL variance 1 before weighting, so influence follows the user's
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

    # Missingness is real information, but it is not a measurement, so it is normalized as its own
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
        dummies = pd.get_dummies(as_text(nodes[col]))
        # np.array, not to_numpy() alone. Under pandas 3 the array backing a one-hot frame comes
        # back READ-ONLY, so scaling it in place raises "output array is read-only" -- the third
        # place this project has hit that, after embedding.embed and build_graph.embed. It is
        # version-dependent, which is why it only appears on the newer pandas.
        D = np.array(dummies.to_numpy(dtype=float), dtype=float, copy=True)
        v = D.var(axis=0).sum()
        D *= (spec.categorical_weight / np.sqrt(v)) if v > 1e-12 else 1.0
        X = np.hstack([X, D])
        # Taken from the frame already built rather than rebuilt, which computed the same one-hot
        # encoding twice and could disagree with itself if the column had NA.
        names += [f"{col}::{c}" for c in dummies.columns]

    if spec.na_policy == "drop_genes":
        # Evaluated against the pre-imputation values. Testing X would always pass, because every block
        # has already had its NaNs replaced by then -- which made this policy a silent no-op.
        rows = np.isfinite(np.hstack(raw_for_drop)).all(axis=1) if raw_for_drop else rows
        X = X[rows]
        log(f"  drop_genes: kept {int(rows.sum()):,} of {len(rows):,} genes with no missing value")

    if X.shape[0] == 0:
        # The mirror of the empty-column error above, which already existed. One all-missing column is
        # enough for drop_genes to remove every gene, and the matrix then reported "0 genes x 4
        # features" and carried on -- so the failure surfaced from inside UMAP, far from its cause.
        raise ValueError("every gene was dropped by the missing-value policy; at least one selected "
                         "column is missing for every gene")

    log(f"matrix: {X.shape[0]:,} genes x {X.shape[1]} features "
        f"({spec.na_policy}, {spec.scaling}"
        + (f", categorical x{spec.categorical_weight}" if spec.categorical else "") + ")")
    return X, names, rows


def normalize(Y, scale: float = 50.0) -> np.ndarray:
    """Center an embedding and scale it to a fixed extent.

    Every map arrives at the same size, which is what lets one replace another in the view without
    the camera having to be re-framed, and what makes two thumbnails comparable. The scaling is
    uniform across axes rather than per-axis, so it moves and resizes the cloud without distorting
    it: neighborhoods, distance ranks and any clustering computed on the result are unchanged.

    Centerd in double precision and cast to float32 only at the end. Done the other way round -- the
    cast first, as this did -- a cloud whose spread is small next to its offset from the origin loses
    that spread to cancellation: 900 ± 0.001 in float32 has about three digits left to subtract with,
    and the map comes back quantised into bands. It is the cast that has to be last, not the
    subtraction. The GL widget still receives float32, which is what it wants.

    Building a new array rather than scaling in place also sidesteps the read-only input: UMAP
    returns one in recent versions, and in-place centering fails there with "output array is
    read-only" -- the fourth place this project has hit that.
    """
    Y = np.asarray(Y, dtype=np.float64)
    Y = Y - Y.mean(0)
    Y = Y / (np.abs(Y).max() + 1e-9)
    return (Y * scale).astype(np.float32)


def embed(nodes: pd.DataFrame, spec: EmbeddingSpec, log=print, return_matrix: bool = False):
    """Build and project a matrix.

    Returns ``(coords, feature_names, kept_rows)`` by default. With ``return_matrix=True`` the
    exact pre-projection matrix is appended, allowing callers to compute trustworthiness without
    rebuilding and potentially drifting from the matrix that actually produced the map.
    """
    from .logging_util import get_logger
    # The full recipe, at INFO. A map with no record of what produced it cannot be compared with
    # another or reported in a methods section, and the recipe is small next to the run.
    get_logger(__name__).info("embedding: %s", spec.to_dict())
    X, names, rows = build_matrix(nodes, spec, log=log)
    if spec.method == "pca":
        from sklearn.decomposition import PCA
        # Clamped to what the matrix can actually give. A block set with two columns and a request
        # for three components is an error sklearn raises, and raising it in the middle of a search
        # ends the whole run over one configuration that was never going to work.
        k = int(max(min(spec.n_components, min(np.shape(X))), 1))
        log(f"PCA: {k} components -- the baseline, not a map to read biology off")
        out = (normalize(PCA(n_components=k, random_state=spec.random_state)
                         .fit_transform(np.nan_to_num(X))), names, rows)
        return (*out, X) if return_matrix else out
    if spec.method == "tsne":
        # Perplexity has to stay under a third of the sample or the neighbourhoods it builds cover
        # the whole set; sklearn raises rather than clamping, which would end a sweep mid-run.
        perplexity = float(min(spec.perplexity, max((len(X) - 1) / 3.0, 2.0)))
        from . import gpu
        on_gpu = gpu.tsne_class()
        if on_gpu is not None and len(X) >= 1000 and spec.n_components == 2:
            try:
                log(f"t-SNE: {gpu.backend()['tsne']} on the GPU, perplexity {perplexity:g} -- a "
                    "different map from scikit-learn, not the same map faster")
                Y = on_gpu(n_components=2, perplexity=perplexity,
                           early_exaggeration=spec.early_exaggeration,
                           random_state=spec.random_state).fit_transform(X)
                out = (normalize(np.asarray(Y)), names, rows)
                return (*out, X) if return_matrix else out
            except Exception as exc:
                log(f"cuml t-SNE failed, using scikit-learn: {type(exc).__name__}: {exc}")
        elif on_gpu is not None and spec.n_components != 2:
            log(f"t-SNE: cuml supports 2 components here; using scikit-learn for "
                f"{spec.n_components} components")
        from sklearn.manifold import TSNE
        log(f"t-SNE: scikit-learn, perplexity {perplexity:g} -- neighbourhoods are meaningful, "
            f"distances between clusters are NOT")
        Y = TSNE(n_components=spec.n_components, perplexity=perplexity,
                 early_exaggeration=spec.early_exaggeration, metric=spec.metric,
                 init="pca", random_state=spec.random_state).fit_transform(np.nan_to_num(X))
        out = (normalize(np.asarray(Y)), names, rows)
        return (*out, X) if return_matrix else out
    try:
        from . import gpu
        on_gpu = gpu.umap_class()
        if on_gpu is not None:
            # cuml's UMAP is NOT the reference implementation, so this map is not identical to the
            # one the CPU builds -- it is a different map of the same data. Said out loud, because a
            # walk whose rows came from two implementations would be a comparison of the libraries.
            log(f"UMAP: {gpu.backend()['umap']} on the GPU -- a different map from the CPU path, "
                f"not the same map faster")
            Y = on_gpu(n_components=spec.n_components, n_neighbors=spec.n_neighbors,
                       min_dist=spec.min_dist, random_state=spec.random_state).fit_transform(X)
            out = (normalize(np.asarray(Y)), names, rows)
            return (*out, X) if return_matrix else out
        import umap
        # Named on the CPU path too. "No message" is not an answer to "which library ran": a user
        # who has just turned the switch on needs to see that it did nothing here and why.
        log(f"UMAP: umap-learn {umap.__version__} on the CPU"
            + ("" if gpu.available()["cuml"] else " (cuml not installed)"))
        Y = umap.UMAP(n_components=spec.n_components, n_neighbors=spec.n_neighbors,
                      min_dist=spec.min_dist, metric=spec.metric,
                      random_state=spec.random_state).fit_transform(X)
    except Exception as e:
        log(f"UMAP unavailable ({type(e).__name__}); falling back to PCA")
        from sklearn.decomposition import PCA
        Y = PCA(n_components=spec.n_components,
                random_state=spec.random_state).fit_transform(np.nan_to_num(X))
    out = (normalize(Y), names, rows)
    return (*out, X) if return_matrix else out


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
    cols = ["column", "missing_frac", "centroid_gap"]
    if not out:
        # An empty frame built from [] has no columns at all, so sorting raised KeyError -- meaning the
        # leak check crashed in exactly the case it should have reported as clean.
        return pd.DataFrame(columns=cols)
    return pd.DataFrame(out, columns=cols).sort_values("centroid_gap", ascending=False)
