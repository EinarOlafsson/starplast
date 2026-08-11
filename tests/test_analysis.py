#!/usr/bin/env python3
"""Regression tests for the embedding, clustering and tuning layers.

As with the literature layer, each test pins something that was wrong at some point. The two that matter
most are the block-weight semantics (weights that do not control influence are worse than no weights) and
the circularity guard (a guard defeated by renaming a column is not a guard).
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import clustering as C          # noqa: E402
from starplast import tuning as T              # noqa: E402
from starplast.embedding import (              # noqa: E402
    EmbeddingSpec, build_matrix, columns_for, variance_share, missingness_leak)


@pytest.fixture
def toy():
    """Three well-separated groups, plus columns that restate the grouping."""
    rng = np.random.default_rng(0)
    n = 240
    grp = np.repeat([0, 1, 2], n // 3)
    df = pd.DataFrame({
        "gene_id": [f"TGME49_{200000 + i}" for i in range(n)],
        "expr_tachy": grp * 5.0 + rng.normal(0, .3, n),
        "expr_cyst": grp * 5.0 + rng.normal(0, .3, n),
        "expr_max": grp * 5.0 + rng.normal(0, .3, n),
        "fit_a": rng.normal(0, 1, n),
        "fit_b": rng.normal(0, 1, n),
        "mean_plddt": rng.normal(70, 10, n),
        "paralog_number": rng.integers(0, 5, n).astype(float),
        "compartment": pd.Series(grp).map({0: "ER", 1: "Golgi", 2: "nucleus"}),
    })
    df["compartment_copy"] = df.compartment              # exact restatement
    # Deliberately noisy: correlates ~0.85 with the expression inputs, so it separates the groups
    # strongly while staying below the derived-restatement threshold. An earlier version used low noise,
    # correlated 0.99 with expr_tachy, and was correctly flagged derived -- the fixture was the bug.
    df["held_out_signal"] = grp * 3.0 + rng.normal(0, 1.5, n)
    df["noise"] = rng.normal(0, 1, n)
    df.loc[df.index[:60], "mean_plddt"] = np.nan          # a measured-vs-not pattern
    return df, grp


# --------------------------------------------------------------------------- embedding
def test_block_weights_control_influence(toy):
    """A block's share must follow its weight, not its column count or tail shape.

    Before this, hyperLOPIT's 27 one-hot columns took 56% of the matrix by being numerous, and robust
    scaling on the published screens inflated that block to half.
    """
    df, _ = toy
    spec = EmbeddingSpec(blocks=("expression_summary", "fitness_screens", "protein_features"),
                         na_policy="median", scaling="zscore")
    v = variance_share(df, spec)["share"]
    assert v.max() - v.min() < 0.05, f"equal weights should give equal shares, got {v.to_dict()}"

    spec2 = EmbeddingSpec(blocks=("expression_summary", "fitness_screens", "protein_features"),
                          na_policy="median", scaling="zscore",
                          block_weights={"fitness_screens": 3.0})
    v2 = variance_share(df, spec2)["share"]
    assert v2["fitness_screens"] > 0.5


def test_categorical_weight_is_not_swamped_by_column_count(toy):
    """hyperLOPIT was named as a design input while carrying 1.1% of the variance."""
    df, _ = toy
    spec = EmbeddingSpec(blocks=("expression_summary",), categorical=("compartment",),
                         na_policy="median", scaling="zscore")
    v = variance_share(df, spec)["share"]
    assert 0.3 < v["compartment"] < 0.7, f"categorical block should be comparable, got {v.to_dict()}"


def test_drop_genes_actually_drops_genes(toy):
    """It tested for NaN after imputation had already run, so it silently dropped nothing."""
    df, _ = toy
    spec = EmbeddingSpec(blocks=("protein_features",), na_policy="drop_genes", scaling="zscore")
    _, _, rows = build_matrix(df, spec, log=lambda *a: None)
    assert rows.sum() < len(df)
    assert rows.sum() == int(df.mean_plddt.notna().sum())


def test_na_policies_all_produce_a_matrix(toy):
    df, _ = toy
    for pol in ("median", "indicator", "drop_columns"):
        spec = EmbeddingSpec(blocks=("expression_summary", "protein_features"),
                             na_policy=pol, scaling="rank")
        X, names, _ = build_matrix(df, spec, log=lambda *a: None)
        assert X.shape[0] == len(df) and X.shape[1] >= 1
        assert np.isfinite(X).all(), f"{pol} left non-finite values"


def test_indicator_policy_adds_missingness_columns(toy):
    df, _ = toy
    spec = EmbeddingSpec(blocks=("protein_features",), na_policy="indicator", scaling="zscore")
    _, names, _ = build_matrix(df, spec, log=lambda *a: None)
    assert any(n.startswith("missing::") for n in names)


def test_missingness_leak_is_measurable(toy):
    """The defect was never the gap itself; it was that the gap was hidden and unmeasured."""
    df, _ = toy
    spec = EmbeddingSpec(blocks=("expression_summary", "protein_features"),
                         na_policy="median", scaling="zscore")
    X, _, _ = build_matrix(df, spec, log=lambda *a: None)
    d = missingness_leak(X[:, :3], df, ["mean_plddt"])
    assert list(d.columns) == ["column", "missing_frac", "centroid_gap"]
    assert np.isfinite(d.centroid_gap).all()


# --------------------------------------------------------------------------- clustering
def test_walks_return_scored_grids(toy):
    df, grp = toy
    X = df[["expr_tachy", "expr_cyst", "expr_max"]].to_numpy()
    d = C.walk_dbscan(X, min_samples_values=(5, 10), log=lambda *a: None)
    assert {"eps", "min_samples", "n_clusters", "silhouette"} <= set(d.columns)
    assert (d.n_clusters > 0).any()


def test_battery_finds_a_held_out_signal(toy):
    df, grp = toy
    S, D = C.battery(df, grp, used_features=["expr_tachy", "expr_cyst", "expr_max"],
                     log=lambda *a: None)
    held = S[S.evidence == "held_out"]
    assert "held_out_signal" in set(held.feature)
    assert held.set_index("feature").loc["held_out_signal", "score"] > 0.5
    noise = S.set_index("feature").loc["noise", "score"]
    assert noise < 0.2, "pure noise should not score as separating"


def test_circularity_guard_is_not_defeated_by_renaming(toy):
    """`compartment` fed the embedding; `compartment_copy` is the same column under another name.

    The first version reported such aliases as the top held-out discoveries at V = 0.96.
    """
    df, grp = toy
    S, _ = C.battery(df, grp, used_features=["compartment"], log=lambda *a: None)
    ev = S.set_index("feature").evidence
    assert ev["compartment"] == "used"
    assert ev["compartment_copy"] == "derived", "an exact restatement must not count as evidence"


def test_every_result_reports_its_association_with_the_inputs(toy):
    """So a borderline case like a second inference over the same experiment is visible, not hidden."""
    df, grp = toy
    S, _ = C.battery(df, grp, used_features=["compartment"], log=lambda *a: None)
    assert "assoc_with_input" in S.columns
    assert S.set_index("feature").loc["compartment_copy", "assoc_with_input"] > 0.9
    assert S.set_index("feature").loc["noise", "assoc_with_input"] < 0.5


def test_precision_and_recall_are_reported_separately(toy):
    """'98% of cluster 1 is tachyzoite' and 'cluster 1 holds 98% of tachyzoite genes' differ."""
    df, grp = toy
    _, D = C.battery(df, grp, used_features=[], features=["compartment"], log=lambda *a: None)
    assert {"precision", "recall", "odds_ratio", "q"} <= set(D.columns)
    r = D[(D.feature == "compartment") & (D.category == "ER")].sort_values("precision").iloc[-1]
    assert r.precision > 0.9 and r.recall > 0.9


def test_multiple_testing_is_corrected(toy):
    df, grp = toy
    S, D = C.battery(df, grp, used_features=[], log=lambda *a: None)
    assert "q" in D.columns
    ok = D.p.notna() & D.q.notna()          # continuous rows carry no per-cluster p
    assert ok.any(), "no testable rows"
    assert (D.q[ok] >= D.p[ok] - 1e-12).all(), "q must never be smaller than p"


# --------------------------------------------------------------------------- tuning
def test_umap_walk_is_seeded_and_records_it(toy):
    pytest.importorskip("umap")
    df, _ = toy
    spec = EmbeddingSpec(blocks=("expression_summary",), na_policy="median", scaling="zscore",
                         n_components=2)
    w = T.walk_umap(df, spec, n_neighbors_values=(10,), min_dist_values=(0.1,),
                    sample_size=120, seed=42, cluster_check=False, log=lambda *a: None)
    assert len(w) == 1 and int(w.iloc[0].seed) == 42 and int(w.iloc[0].sample_size) == 120


def test_embedding_store_round_trips_the_recipe(tmp_path, toy):
    """Coordinates without their spec cannot be compared or reported."""
    df, _ = toy
    spec = EmbeddingSpec(name="t", blocks=("expression_summary",), scaling="rank", n_neighbors=11)
    s = T.EmbeddingStore(str(tmp_path))
    s.save("t", np.zeros((len(df), 3), dtype=np.float32), spec, gene_ids=df.gene_id)
    xyz, spec2, meta = s.load("t")
    assert xyz.shape == (len(df), 3)
    assert spec2.scaling == "rank" and spec2.n_neighbors == 11
    assert list(s.list().name) == ["t"]


def test_import_repairs_a_malformed_gene_column():
    df = pd.DataFrame({"id": ["TgME49.208830", "TgME49.205250", "junk"], "v": [1.0, 2.0, 3.0]})
    col, frac, needs = T.suggest_gene_column(df)
    assert col == "id" and needs is True, "a malformed column is exactly when a suggestion is needed"
    pat, rep = T.suggest_repair(df.id)
    out = T.import_table(df, col, pattern=pat, replacement=rep, log=lambda *a: None)
    assert list(out.index) == ["TGME49_205250", "TGME49_208830"]
    assert "imported_v" in out.columns


def test_import_routes_through_the_identity_resolver():
    """Superseded accessions must resolve, or a table silently contributes nothing."""
    df = pd.DataFrame({"id": ["TGME49_008830"], "v": [1.0]})
    out = T.import_table(df, "id", resolve=lambda a: "TGME49_208830", log=lambda *a: None)
    assert list(out.index) == ["TGME49_208830"]


# --------------------------------------------------------------------------- sources
def test_quantification_is_inferred_from_range_not_filename():
    """Both directions of the same mistake, from the same GEO series.

    Real FPKM reaching 16,520 must be logged; the already-logged derivative of it, still named FPKM,
    must not be logged twice. An earlier version returned "log_intensity" for anything non-integer and
    so silently skipped the log on real FPKM.
    """
    from starplast import sources as S
    rng = np.random.default_rng(0)
    real_fpkm = pd.Series(rng.gamma(1.2, 400, 4000))            # wide, positive, non-integer
    already_log = np.log2(real_fpkm + 1)
    assert S.infer_quant(real_fpkm, "GSE108740_FPKM") == "fpkm"
    assert S.infer_quant(already_log, "rna108740_Tachyzoites_FPKM") == "log_intensity"
    assert S.infer_quant(pd.Series(rng.normal(0, 2, 500)), "LFCs") == "lfc"


def test_normalise_never_logs_twice():
    from starplast import sources as S
    df = pd.DataFrame({"a": [0.0, 10.0, 16520.0, 3.0] * 25})
    once = S.normalise(df, "fpkm", log=lambda *a: None)
    twice = S.normalise(once, "log_intensity", log=lambda *a: None)
    assert once.max().max() < 20, "linear intensity should be logged"
    # centring is idempotent-ish; the point is that no second log is applied
    assert abs(twice.max().max() - once.max().max()) < 1e-9


def test_ratios_are_left_alone():
    """Centring a log-ratio moves its zero, which is the reference condition."""
    from starplast import sources as S
    df = pd.DataFrame({"lfc": [-3.0, 0.0, 2.5, 1.0] * 25})
    out = S.normalise(df, "lfc", log=lambda *a: None)
    assert (out.lfc == df.lfc).all()


def test_rank_normalise_puts_incomparable_units_on_one_axis():
    from starplast import sources as S
    df = pd.DataFrame({"fpkm": [1.0, 10.0, 100.0, 1000.0], "ibaq": [30.0, 25.0, 20.0, 15.0]})
    r = S.rank_normalise(df)
    assert r.min().min() >= -0.5 and r.max().max() <= 0.5
    assert r.fpkm.corr(r.ibaq) < 0, "opposite orderings must stay opposite"
