#!/usr/bin/env python3
"""Embedding construction: the NA policies, the scalings, and the leak checks.

The missing-value policy is not a detail here. Roughly half this proteome has never been measured by
most assays, so "missing" is the most common value in the table, and how it is handled decides what the
map is a map OF. `drop_genes` must be evaluated before imputation or it silently becomes a no-op;
`indicator` deliberately makes missingness a feature so its contribution can be measured rather than
assumed away.

`missingness_leak` exists because of a claim that turned out to be overstated: median imputation was
said to encode study effort. Measured, missingness correlates with the embedding at r = 0.07-0.15
across assays -- much of the apparent structure is real biology. The function is what settles that
argument for any given embedding instead of re-arguing it.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast.embedding import (BLOCKS, EmbeddingSpec, NA_POLICIES, SCALINGS, build_matrix,  # noqa: E402
                                 columns_for, embed, missingness_leak, variance_share)


def _nodes(n=120, seed=0):
    rng = np.random.default_rng(seed)
    d = pd.DataFrame({f"fit_{i}": rng.normal(size=n) for i in range(4)})
    for i in range(3):
        d[f"expr_{i}"] = rng.normal(size=n)
    d.insert(0, "gene_id", [f"TGME49_{200000+i}" for i in range(n)])
    return d


# --------------------------------------------------------------------------- the spec
def test_a_spec_round_trips_through_a_dict():
    """Stored beside every saved embedding, so a result can be reopened with the settings that made it."""
    spec = EmbeddingSpec(blocks=("fitness_screens",), na_policy="median", scaling="rank")
    assert EmbeddingSpec.from_dict(spec.to_dict()) == spec


def test_unknown_keys_in_a_stored_spec_are_ignored():
    """An embedding saved by a later version must still open in an earlier one."""
    d = EmbeddingSpec(blocks=("fitness_screens",)).to_dict()
    d["a_setting_from_the_future"] = 7
    assert EmbeddingSpec.from_dict(d).blocks == ("fitness_screens",)


# --------------------------------------------------------------------------- column selection
def test_a_block_with_no_matching_columns_is_absent():
    out = columns_for(_nodes(), EmbeddingSpec(blocks=("fitness_screens", "literature")))
    assert "fitness_screens" in out and "literature" not in out


def test_an_unknown_block_name_is_ignored():
    assert columns_for(_nodes(), EmbeddingSpec(blocks=("no_such_block",))) == {}


def test_imported_columns_arrive_as_their_own_block():
    """A user's own CSV is kept separate so its contribution can be weighted and reported on its own."""
    d = _nodes()
    d["my_column"] = 1.0
    out = columns_for(d, EmbeddingSpec(blocks=(), extra_columns=("my_column",)))
    assert out == {"imported": ["my_column"]}


def test_a_non_numeric_extra_column_is_refused():
    d = _nodes()
    d["label"] = "text"
    assert columns_for(d, EmbeddingSpec(blocks=(), extra_columns=("label",))) == {}


# --------------------------------------------------------------------------- building the matrix
def test_selecting_no_features_is_an_explicit_error():
    """Silently returning an empty matrix would produce a UMAP of nothing that looks like a result."""
    with pytest.raises(ValueError, match="no numeric features"):
        build_matrix(_nodes(), EmbeddingSpec(blocks=("literature",)), log=lambda *_: None)


@pytest.mark.parametrize("scaling", SCALINGS)
def test_every_declared_scaling_produces_a_finite_matrix(scaling):
    X, names, keep = build_matrix(_nodes(), EmbeddingSpec(blocks=("fitness_screens",),
                                                          scaling=scaling), log=lambda *_: None)
    assert np.isfinite(X).all()
    assert X.shape[1] == len(names)


@pytest.mark.parametrize("policy", NA_POLICIES)
def test_every_declared_na_policy_produces_a_finite_matrix(policy):
    d = _nodes()
    d.loc[:20, "fit_0"] = np.nan
    X, names, keep = build_matrix(d, EmbeddingSpec(blocks=("fitness_screens",),
                                                   na_policy=policy), log=lambda *_: None)
    assert np.isfinite(X).all()
    assert keep.sum() == X.shape[0]


def test_drop_genes_actually_drops_genes():
    """Evaluated before imputation. Applied after, every gene has a value and the policy is a no-op --
    which is exactly what it silently was."""
    d = _nodes()
    d.loc[:20, "fit_0"] = np.nan
    X, _, keep = build_matrix(d, EmbeddingSpec(blocks=("fitness_screens",), na_policy="drop_genes"),
                              log=lambda *_: None)
    assert keep.sum() == len(d) - 21


def test_indicator_adds_a_column_that_says_what_was_missing():
    """Missingness is made a feature on purpose, so its contribution can be measured rather than
    assumed away."""
    d = _nodes()
    d.loc[:20, "fit_0"] = np.nan
    X, names, _ = build_matrix(d, EmbeddingSpec(blocks=("fitness_screens",), na_policy="indicator"),
                               log=lambda *_: None)
    assert any("missing" in n for n in names)


def test_drop_columns_removes_the_column_not_the_gene():
    d = _nodes()
    d["fit_0"] = np.nan
    X, names, keep = build_matrix(d, EmbeddingSpec(blocks=("fitness_screens",),
                                                   na_policy="drop_columns"), log=lambda *_: None)
    assert "fit_0" not in names
    assert keep.all()


def test_a_column_missing_beyond_the_threshold_is_dropped_and_reported():
    d = _nodes()
    d.loc[:100, "fit_0"] = np.nan          # >80% missing
    msgs = []
    X, names, _ = build_matrix(d, EmbeddingSpec(blocks=("fitness_screens",),
                                                na_policy="indicator", max_missing=0.5),
                               log=msgs.append)
    assert "fit_0" not in names
    assert any("dropped" in m for m in msgs), "a dropped column must be announced, not silently gone"


def test_dropping_every_column_is_an_explicit_error():
    """An empty matrix would produce a UMAP of nothing that looks exactly like a result."""
    d = _nodes()
    for c in [c for c in d.columns if c.startswith("fit_")]:
        d[c] = np.nan
    with pytest.raises(ValueError, match="dropped by the missing-value policy"):
        build_matrix(d, EmbeddingSpec(blocks=("fitness_screens",), na_policy="indicator",
                                      max_missing=0.5), log=lambda *_: None)


def test_dropping_every_gene_is_an_explicit_error():
    """One all-missing column is enough for drop_genes to remove every gene. It used to report
    "0 genes x 4 features" and carry on, so the failure surfaced from inside UMAP, far from its cause."""
    d = _nodes()
    d.loc[:, "fit_0"] = np.nan
    with pytest.raises(ValueError, match="every gene was dropped"):
        build_matrix(d, EmbeddingSpec(blocks=("fitness_screens",), na_policy="drop_genes"),
                     log=lambda *_: None)


def test_a_constant_column_does_not_produce_nan_under_scaling():
    """Zero spread divides by zero in both robust and z-score scaling."""
    d = _nodes()
    d["fit_0"] = 5.0
    for scaling in ("robust", "zscore"):
        X, _, _ = build_matrix(d, EmbeddingSpec(blocks=("fitness_screens",), scaling=scaling),
                               log=lambda *_: None)
        assert np.isfinite(X).all()


def test_blocks_are_normalized_to_equal_variance_before_weighting():
    """Otherwise a block with 70 columns drowns one with 3, and the map is about column counts."""
    share = variance_share(_nodes(), EmbeddingSpec(blocks=("fitness_screens", "expression_summary")))
    assert len(share) == 2
    assert share.share.sum() == pytest.approx(1.0, abs=1e-6)
    assert share.share.min() > 0.3, "two blocks should carry comparable variance by default"


def test_a_weight_shifts_the_variance_share():
    """Blocks are normalized to equal variance first, so a weight means what it says rather than
    depending on how many columns a block happens to have."""
    spec = EmbeddingSpec(blocks=("fitness_screens", "expression_summary"),
                         block_weights={"fitness_screens": 3.0})
    share = variance_share(_nodes(), spec)
    assert share.loc["fitness_screens", "share"] > share.loc["expression_summary", "share"]


def test_variance_share_names_a_block_contributing_almost_nothing():
    """The check that caught hyperLOPIT being named as an input while carrying 1.1% of the variance."""
    spec = EmbeddingSpec(blocks=("fitness_screens", "expression_summary"),
                         block_weights={"expression_summary": 0.01})
    share = variance_share(_nodes(), spec)
    assert share.loc["expression_summary", "share"] < 0.05


# --------------------------------------------------------------------------- embedding
def test_embedding_is_reproducible_under_a_fixed_seed():
    spec = EmbeddingSpec(blocks=("fitness_screens",), random_state=42)
    a, *_ = embed(_nodes(), spec, log=lambda *_: None)
    b, *_ = embed(_nodes(), spec, log=lambda *_: None)
    assert np.allclose(a, b)


def test_the_embedding_array_is_writable():
    """umap returns a read-only array in recent versions, and the in-place centring downstream raises
    "output array is read-only"."""
    Y, *_ = embed(_nodes(), EmbeddingSpec(blocks=("fitness_screens",)), log=lambda *_: None)
    Y[0, 0] = 1.0


# --------------------------------------------------------------------------- the leak check
def test_missingness_leak_is_large_when_the_layout_is_driven_by_missingness():
    """The check that settles the argument rather than re-arguing it. Reported as a centroid gap in map
    radii, so it is comparable across embeddings of different scale."""
    d = _nodes(200)
    d.loc[:99, "fit_0"] = np.nan
    coords = np.zeros((200, 3))
    coords[:100, 0] = 10.0                  # position determined entirely by whether fit_0 is missing
    out = missingness_leak(coords, d, ["fit_0"])
    assert out.iloc[0].centroid_gap > 1.5


def test_missingness_leak_is_small_when_the_layout_is_independent_of_it():
    """Measured across the real assays it sits at 0.07-0.15, so much of the apparent structure is real
    biology rather than study effort."""
    rng = np.random.default_rng(0)
    d = _nodes(200)
    d.loc[:99, "fit_0"] = np.nan
    out = missingness_leak(rng.normal(size=(200, 3)), d, ["fit_0"])
    assert out.iloc[0].centroid_gap < 0.5


@pytest.mark.parametrize("make", [
    lambda d: d,                                    # never missing
    lambda d: d.assign(fit_0=np.nan),               # always missing
])
def test_a_column_with_no_contrast_is_not_reported(make):
    """Nothing to measure either way, and a row of zeros is noise in the table."""
    rng = np.random.default_rng(0)
    out = missingness_leak(rng.normal(size=(200, 3)), make(_nodes(200)), ["fit_0"])
    assert out.empty


def test_a_column_that_is_not_in_the_table_is_skipped():
    rng = np.random.default_rng(0)
    assert missingness_leak(rng.normal(size=(200, 3)), _nodes(200), ["not_a_column"]).empty


# --------------------------------------------------------------------------- categorical blocks
def test_a_categorical_block_is_scaled_to_a_stated_influence():
    """27 one-hot compartment columns each carry variance 0.25*p(1-p), which totalled 1.1% of the
    matrix -- so hyperLOPIT had no influence at all while being listed as an input. The block is now
    scaled so its weight is stated rather than accidental."""
    d = _nodes(200)
    d["compartment"] = ["a", "b", "c", "d"] * 50
    light = build_matrix(d, EmbeddingSpec(blocks=("fitness_screens",), categorical=("compartment",),
                                          categorical_weight=0.1), log=lambda *_: None)[0]
    heavy = build_matrix(d, EmbeddingSpec(blocks=("fitness_screens",), categorical=("compartment",),
                                          categorical_weight=5.0), log=lambda *_: None)[0]
    cat_var_light = light[:, -4:].var(axis=0).sum()
    cat_var_heavy = heavy[:, -4:].var(axis=0).sum()
    assert cat_var_heavy > cat_var_light * 10


def test_categorical_columns_are_named_so_their_block_is_identifiable():
    d = _nodes(200)
    d["compartment"] = ["a", "b"] * 100
    _, names, _ = build_matrix(d, EmbeddingSpec(blocks=("fitness_screens",),
                                                categorical=("compartment",)), log=lambda *_: None)
    assert any(n.startswith("compartment::") for n in names)


def test_a_categorical_column_that_is_not_in_the_table_is_skipped():
    d = _nodes(200)
    X, names, _ = build_matrix(d, EmbeddingSpec(blocks=("fitness_screens",),
                                                categorical=("not_a_column",)), log=lambda *_: None)
    assert not any("::" in n for n in names)


def test_a_constant_categorical_does_not_divide_by_zero():
    d = _nodes(200)
    d["compartment"] = "same"
    X, _, _ = build_matrix(d, EmbeddingSpec(blocks=("fitness_screens",),
                                            categorical=("compartment",)), log=lambda *_: None)
    assert np.isfinite(X).all()


def test_variance_share_separates_categorical_and_indicator_contributions():
    """So "the map is mostly missingness indicators" is a visible statement rather than a suspicion."""
    d = _nodes(200)
    d.loc[:60, "fit_0"] = np.nan
    d["compartment"] = ["a", "b"] * 100
    share = variance_share(d, EmbeddingSpec(blocks=("fitness_screens",),
                                            categorical=("compartment",), na_policy="indicator"))
    assert "missing indicators" in share.index
    assert "compartment" in share.index


# --------------------------------------------------------------------------- the fallback
def test_without_umap_the_embedding_falls_back_to_pca_and_says_so(monkeypatch):
    """A heavy optional dependency. Its absence should cost the layout quality, not the session."""
    import builtins
    real = builtins.__import__

    def no_umap(name, *a, **k):
        if name == "umap":
            raise ImportError("no umap here")
        return real(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", no_umap)
    msgs = []
    Y, names, keep = embed(_nodes(), EmbeddingSpec(blocks=("fitness_screens",)), log=msgs.append)
    assert Y.shape[1] == 3
    assert any("falling back to PCA" in m for m in msgs)


def test_a_block_whose_columns_are_all_too_sparse_is_skipped_under_drop_columns():
    d = _nodes()
    for c in [c for c in d.columns if c.startswith("fit_")]:
        d.loc[:100, c] = np.nan
    X, names, _ = build_matrix(d, EmbeddingSpec(blocks=("fitness_screens", "expression_summary"),
                                                na_policy="drop_columns", max_missing=0.5),
                               log=lambda *_: None)
    assert not any(n.startswith("fit_") for n in names)
    assert any(n.startswith("expr_") for n in names)


def test_a_read_only_one_hot_array_is_copied_before_scaling():
    """The third home of "output array is read-only", and the one that hit a user.

    pandas 3 hands back a READ-ONLY array from a one-hot frame, so scaling it in place raises. It is
    version-dependent, so this builds the read-only array explicitly rather than relying on whichever
    pandas is installed -- otherwise the test passes on pandas 2 while the application fails on 3.
    """
    import numpy as np
    import pandas as pd
    from starplast.embedding import EmbeddingSpec, build_matrix

    real = pd.get_dummies

    def read_only_dummies(*a, **kw):
        out = real(*a, **kw)
        arr = out.to_numpy()
        arr.flags.writeable = False
        return pd.DataFrame(arr, columns=out.columns, index=out.index)

    nodes = pd.DataFrame({
        "gene_id": [f"g{i}" for i in range(40)],
        "compartment": ["a", "b"] * 20,
        "expr_tachy": np.linspace(0, 1, 40),
        "expr_cyst": np.linspace(1, 0, 40),
    })
    spec = EmbeddingSpec(blocks=("expression_summary",), categorical=("compartment",))
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(pd, "get_dummies", read_only_dummies)
        X, names, rows = build_matrix(nodes, spec, log=lambda *a: None)
    assert X.shape[0] == 40
    assert any("compartment::" in n for n in names)


def test_one_hot_column_names_come_from_the_frame_that_was_encoded():
    """They were built by encoding the column a second time, which is the same work done twice and
    two chances to disagree if the column has missing values."""
    import numpy as np
    import pandas as pd
    from starplast.embedding import EmbeddingSpec, build_matrix

    nodes = pd.DataFrame({
        "gene_id": [f"g{i}" for i in range(30)],
        "compartment": ["a", "b", None] * 10,
        "expr_tachy": np.linspace(0, 1, 30),
        "expr_cyst": np.linspace(1, 0, 30),
    })
    X, names, rows = build_matrix(
        nodes, EmbeddingSpec(blocks=("expression_summary",), categorical=("compartment",)),
        log=lambda *a: None)
    onehot = [n for n in names if n.startswith("compartment::")]
    assert len(onehot) == X.shape[1] - 2, "the names do not match the columns that were added"


# --------------------------------------------------------------------------- normalizing a map
def test_a_map_arrives_at_a_known_size_wherever_it_came_from():
    """Every embedding is put at the same extent, which is what lets one replace another in the view
    without the camera having to be re-framed, and what makes two thumbnails comparable."""
    from starplast.embedding import normalize
    Y = np.random.default_rng(0).normal(size=(50, 3)) * 0.001 + 900.0
    out = normalize(Y)
    assert np.allclose(out.mean(0), 0.0, atol=1e-3)
    assert np.isclose(np.abs(out).max(), 50.0, atol=1e-3)


def test_normalizing_moves_and_resizes_a_map_without_distorting_it():
    """It is applied AFTER trustworthiness and the clustering are computed, and the guarantee that it
    could not have changed them is that it scales both axes by one factor. Per-axis scaling would
    stretch a genuinely elongated map into a round one -- the difference between two configurations
    that a walk exists to show."""
    from scipy.spatial.distance import pdist
    from starplast.embedding import normalize
    Y = np.random.default_rng(1).normal(size=(40, 3)) * np.array([10.0, 1.0, 0.5])
    d0, d1 = pdist(Y), pdist(normalize(Y))
    assert np.allclose(d1 / d0, (d1 / d0)[0])


def test_a_read_only_embedding_is_copied_rather_than_scaled_in_place():
    """umap returns a read-only array in recent versions and the centring then fails with "output
    array is read-only" -- the fourth place this project has hit that, and it appears only on the
    newer library."""
    from starplast.embedding import normalize
    Y = np.random.default_rng(2).normal(size=(20, 3))
    Y.flags.writeable = False
    out = normalize(Y)
    assert out.shape == (20, 3) and not Y.flags.writeable


# --------------------------------------------------------------------------- a renamed block
def _lopit_nodes(n=120):
    """A table carrying the columns the renamed block selects."""
    d = _nodes(n)
    rng = np.random.default_rng(3)
    d["lopit_prob"] = rng.random(n)
    d["lopit_methods_agree"] = rng.integers(0, 2, n).astype(float)
    return d


def test_a_recipe_saved_under_the_old_block_name_still_selects_its_columns():
    """The one thing a rename can break in silence. `columns_for` skips a block it does not
    recognise, so a recipe carrying the pre-rename spelling would not fail -- it would build an
    embedding with a whole feature block missing, under the name of the run it is supposed to
    reproduce. Every embedding saved before the rename carries that spelling."""
    spec = EmbeddingSpec.from_dict({"blocks": ["localisation", "fitness_screens"],
                                    "block_weights": {"localisation": 2.0}})
    assert spec.blocks == ("localization", "fitness_screens")
    assert spec.block_weights == {"localization": 2.0}, "a weight must follow its block's new name"
    assert "localization" in columns_for(_lopit_nodes(), spec)


def test_a_results_row_saved_under_the_old_block_name_rebuilds_its_map():
    """The same translation by the other route. A results table saved before the rename holds its
    blocks as one `+`-joined string, and clicking that row is supposed to rebuild the map it
    describes -- from the same features, or it is a different map wearing the row's numbers."""
    from starplast.search import rebuild
    coords, genes, labels, features = rebuild(
        _lopit_nodes(80),
        {"blocks": "localisation", "na_policy": "median", "scaling": "rank", "n_neighbors": 5,
         "min_dist": 0.1, "seed": 0, "min_cluster_size": 5, "sample_size": 60},
        log=lambda *_: None)
    assert len(coords) == int(np.sum(genes))
    assert [f for f in features if f.startswith("lopit_")], "the renamed block contributed nothing"
