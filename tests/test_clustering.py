#!/usr/bin/env python3
"""Clustering walks and the hypothesis battery.

The battery runs hundreds of tests against one clustering, so two things decide whether its output
means anything: the multiple-testing correction, and the evidence tier. A feature the embedding was
built from separates the clusters by construction, and so does anything that restates it -- which is
how `compartment` produced its own twin `lopit_map` as the top "discovery" at V = 0.96. Only `held_out`
results are evidence, and the tier is assigned by measuring association rather than by trusting names.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import clustering as CL  # noqa: E402
from starplast.clustering import NOISE  # noqa: E402


def _blobs(n=180, sep=6.0, seed=0):
    rng = np.random.default_rng(seed)
    half = n // 2
    return np.vstack([rng.normal(-sep / 2, 0.5, (half, 3)), rng.normal(sep / 2, 0.5, (n - half, 3))])


def _labels(n=180):
    return np.array([0] * (n // 2) + [1] * (n - n // 2))


# --------------------------------------------------------------------------- scoring a clustering
def test_a_clustering_is_scored_on_clustered_points_only():
    """Including noise in the silhouette makes every setting look bad, so the number would rank
    settings by how much they refused to classify."""
    X = _blobs()
    lab = _labels()
    lab[:10] = NOISE
    out = CL._score(X, lab)
    assert out["n_clusters"] == 2
    assert out["noise_frac"] == pytest.approx(10 / len(lab))
    assert out["silhouette"] > 0.5


def test_a_single_cluster_has_no_silhouette():
    """It is undefined with one group, and reporting zero would rank it against real partitions."""
    out = CL._score(_blobs(), np.zeros(180, dtype=int))
    assert out["n_clusters"] == 1
    assert "silhouette" not in out or np.isnan(out.get("silhouette", np.nan))


def test_everything_being_noise_is_scored_without_crashing():
    out = CL._score(_blobs(), np.full(180, NOISE))
    assert out["n_clusters"] == 0
    assert out["noise_frac"] == pytest.approx(1.0)


def test_the_largest_cluster_share_is_reported():
    """One cluster holding 95% of the proteome is not a structure, and the silhouette will not say so."""
    lab = np.array([0] * 170 + [1] * 10)
    assert CL._score(_blobs(), lab)["largest_frac"] == pytest.approx(170 / 180)


# --------------------------------------------------------------------------- the walks
def test_the_dbscan_grid_defaults_to_a_spread_around_the_data_scale():
    """A fixed eps is meaningless across embeddings whose scale changes with the feature set."""
    msgs = []
    out = CL.walk_dbscan(_blobs(), min_samples_values=(5,), log=msgs.append)
    assert len(out) == 6
    assert any("eps grid from median" in m for m in msgs)


def test_an_explicit_eps_grid_is_used_as_given():
    out = CL.walk_dbscan(_blobs(), eps_values=[0.5, 1.0], min_samples_values=(5,),
                         log=lambda *_: None)
    assert sorted(out.eps.unique()) == [0.5, 1.0]


def test_walks_are_ranked_by_silhouette_with_unscorable_settings_last():
    out = CL.walk_hdbscan(_blobs(), min_cluster_sizes=(10, 20), min_samples_values=(None,),
                          log=lambda *_: None)
    assert len(out) == 2
    sil = out.silhouette.dropna()
    assert sil.is_monotonic_decreasing


def test_hdbscan_absence_is_reported_rather_than_raised(monkeypatch):
    import builtins
    real = builtins.__import__

    def no_hdbscan(name, *a, **k):
        if name == "hdbscan" or (name == "sklearn.cluster" and "HDBSCAN" in (a[2] or ())):
            raise ImportError("nope")
        return real(name, *a, **k)

    monkeypatch.setattr(builtins, "__import__", no_hdbscan)
    msgs = []
    out = CL.walk_hdbscan(_blobs(), log=msgs.append)
    assert out.empty
    assert any("HDBSCAN unavailable" in m for m in msgs)


def test_cluster_dispatches_to_dbscan_when_asked():
    lab = CL.cluster(_blobs(), algorithm="dbscan", eps=1.0, min_samples=5)
    assert len(set(lab[lab != NOISE])) >= 2


def test_cluster_defaults_to_hdbscan():
    lab = CL.cluster(_blobs(), min_cluster_size=20)
    assert len(set(lab[lab != NOISE])) >= 2


# --------------------------------------------------------------------------- multiple testing
def test_benjamini_hochberg_is_monotone_and_bounded():
    p = np.array([0.001, 0.01, 0.02, 0.5, 0.9])
    q = CL._bh(p)
    assert np.all(q >= p - 1e-12)
    assert np.all((q >= 0) & (q <= 1))
    assert np.all(np.diff(q[np.argsort(p)]) >= -1e-12)


def test_missing_p_values_stay_missing():
    q = CL._bh(np.array([0.01, np.nan, 0.2]))
    assert np.isnan(q[1]) and np.isfinite(q[0])


def test_all_missing_p_values_return_all_missing():
    assert np.isnan(CL._bh(np.array([np.nan, np.nan]))).all()


def test_a_single_p_value_is_unchanged():
    assert CL._bh(np.array([0.03]))[0] == pytest.approx(0.03)


# --------------------------------------------------------------------------- association measures
def test_cramers_v_is_near_one_for_a_perfect_association():
    """Not exactly one: chi2_contingency applies Yates' continuity correction to a 2x2, which pulls it
    slightly down. Worth knowing, since the derived-column guard thresholds at 0.95."""
    assert CL._cramers_v(np.array([[50, 0], [0, 50]])) > 0.95


def test_cramers_v_is_near_zero_for_independence():
    assert CL._cramers_v(np.array([[25, 25], [25, 25]])) == pytest.approx(0.0, abs=1e-6)


def test_a_degenerate_table_has_no_association():
    assert np.isnan(CL._cramers_v(np.array([[1, 2]])))
    assert np.isnan(CL._cramers_v(np.zeros((2, 2))))


def test_the_correlation_ratio_detects_a_categorical_driving_a_continuous_column():
    """The branch that used to be skipped entirely, which made the circularity guard blind to every
    numeric column."""
    cat = pd.Series(["a"] * 50 + ["b"] * 50)
    num = pd.Series([0.0] * 50 + [10.0] * 50)
    assert CL._correlation_ratio(cat, num) == pytest.approx(1.0, abs=1e-6)


def test_the_correlation_ratio_is_near_zero_when_the_groups_do_not_differ():
    rng = np.random.default_rng(0)
    cat = pd.Series(["a"] * 200 + ["b"] * 200)
    num = pd.Series(rng.normal(size=400))
    assert CL._correlation_ratio(cat, num) < 0.2


def test_the_correlation_ratio_declines_degenerate_comparisons():
    """A constant column has zero variance, so the ratio is 0/0 rather than 0."""
    assert CL._correlation_ratio(pd.Series(["a", "b"]), pd.Series([1.0, 1.0])) is None
    assert CL._correlation_ratio(pd.Series(["a", "a"]), pd.Series([1.0, 2.0])) is None
    many = pd.Series([f"c{i}" for i in range(50)])
    assert CL._correlation_ratio(many, pd.Series(np.arange(50.0)), max_categories=30) is None


# --------------------------------------------------------------------------- per-feature tests
def test_a_categorical_feature_that_matches_the_clustering_scores_high():
    lab = _labels()
    vals = pd.Series(["in"] * 90 + ["out"] * 90)
    summary, rows = CL.categorical_feature(lab, vals)
    assert summary["score"] > 0.95
    assert summary["score_type"] == "cramers_v"
    best = max(rows, key=lambda r: r["precision"])
    assert best["precision"] == pytest.approx(1.0)


def test_a_categorical_feature_with_too_little_data_is_not_scored():
    lab = np.array([0] * 5 + [1] * 5)
    assert CL.categorical_feature(lab, pd.Series(["a"] * 5 + ["b"] * 5)) == (None, [])


def test_a_feature_with_one_category_is_not_scored():
    """Nothing to separate: every gene has the same value."""
    lab = _labels()
    assert CL.categorical_feature(lab, pd.Series(["same"] * 180)) == (None, [])


def test_rare_categories_are_excluded_from_the_table():
    lab = _labels()
    vals = pd.Series(["in"] * 88 + ["rare"] * 2 + ["out"] * 90)
    _, rows = CL.categorical_feature(lab, vals, min_count=5)
    assert "rare" not in {r["category"] for r in rows}


def test_precision_and_recall_are_reported_separately_per_cluster():
    """"This cluster is 90% rhoptry" and "90% of rhoptry proteins are here" are different claims."""
    lab = _labels()
    vals = pd.Series(["in"] * 90 + ["out"] * 90)
    _, rows = CL.categorical_feature(lab, vals)
    r = [x for x in rows if x["cluster"] == 0 and x["category"] == "in"][0]
    assert r["precision"] == pytest.approx(1.0)
    assert r["recall"] == pytest.approx(1.0)


def test_a_continuous_feature_that_differs_between_clusters_is_detected():
    lab = _labels()
    vals = pd.Series([0.0] * 90 + [10.0] * 90)
    summary, rows = CL.continuous_feature(lab, vals)
    assert summary["score"] > 0.5
    assert summary["score_type"] == "epsilon_squared"
    assert len(rows) == 2


def test_a_continuous_feature_with_too_little_data_is_not_scored():
    lab = np.array([0] * 10 + [1] * 10)
    assert CL.continuous_feature(lab, pd.Series(np.arange(20.0))) == (None, [])


def test_a_continuous_feature_that_is_constant_is_not_scored():
    """Kruskal-Wallis on identical groups is undefined, not uninteresting."""
    lab = _labels()
    assert CL.continuous_feature(lab, pd.Series([1.0] * 180)) == (None, [])


def test_a_single_usable_cluster_is_not_scored():
    lab = np.array([0] * 176 + [1] * 4)     # a cluster of 4 is below the 5-member minimum
    assert CL.continuous_feature(lab, pd.Series(np.arange(180.0))) == (None, [])


def test_noise_points_are_excluded_from_both_kinds_of_test():
    lab = _labels()
    lab[:40] = NOISE
    summary, _ = CL.categorical_feature(lab, pd.Series(["in"] * 90 + ["out"] * 90))
    assert summary["n"] == 140


# --------------------------------------------------------------------------- evidence tiers
def _battery_nodes():
    truth = ["in"] * 90 + ["out"] * 90
    return pd.DataFrame({
        "gene_id": [f"g{i}" for i in range(180)],
        "compartment": truth,
        "compartment_best": truth,                    # an exact restatement under another name
        "unrelated": np.random.default_rng(1).normal(size=180),
    })


def test_a_feature_the_embedding_used_is_marked_used_not_discovered():
    nodes = _battery_nodes()
    res, _ = CL.battery(nodes, _labels(), used_features=("compartment",), log=lambda *_: None)
    assert res.set_index("feature").loc["compartment", "evidence"] == "used"


def test_a_renamed_copy_of_a_used_feature_is_marked_derived():
    """A guard that can be defeated by renaming a column is not a guard. This is the exact failure that
    reported lopit_map as a top discovery at V = 0.96."""
    nodes = _battery_nodes()
    res, _ = CL.battery(nodes, _labels(), used_features=("compartment",), log=lambda *_: None)
    assert res.set_index("feature").loc["compartment_best", "evidence"] == "derived"


def test_an_independent_feature_is_held_out_and_is_the_only_evidence():
    nodes = _battery_nodes()
    res, _ = CL.battery(nodes, _labels(), used_features=("compartment",), log=lambda *_: None)
    assert res.set_index("feature").loc["unrelated", "evidence"] == "held_out"


def test_every_result_carries_its_association_with_the_inputs():
    """Carried alongside rather than only used as a cutoff: lopit_mcmc sits at 0.74, below any sensible
    threshold and still close to circular, so the reader needs the number."""
    nodes = _battery_nodes()
    res, _ = CL.battery(nodes, _labels(), used_features=("compartment",), log=lambda *_: None)
    assert "assoc_with_input" in res.columns


def test_p_values_are_corrected_because_the_battery_runs_hundreds_of_tests():
    nodes = _battery_nodes()
    res, _ = CL.battery(nodes, _labels(), used_features=(), log=lambda *_: None)
    assert "q" in res.columns
    assert (res.q.dropna() >= res.p.dropna()).all()


def test_the_guard_can_be_turned_off_explicitly():
    nodes = _battery_nodes()
    res, _ = CL.battery(nodes, _labels(), used_features=("compartment",),
                        guard_derived=False, log=lambda *_: None)
    assert "derived" not in set(res.evidence)


def test_identifier_columns_are_never_tested_as_features():
    """gene_id separates every cluster perfectly and means nothing."""
    nodes = _battery_nodes()
    res, _ = CL.battery(nodes, _labels(), log=lambda *_: None)
    assert "gene_id" not in set(res.feature)


def test_the_battery_reports_what_it_marked_derived():
    nodes = _battery_nodes()
    msgs = []
    CL.battery(nodes, _labels(), used_features=("compartment",), log=msgs.append)
    assert any("derived" in m for m in msgs)


def test_describe_summarises_a_result_table():
    nodes = _battery_nodes()
    res, per = CL.battery(nodes, _labels(), used_features=("compartment",), log=lambda *_: None)
    lines = CL.describe(res, per)
    assert isinstance(lines, list)
    assert all(isinstance(l, str) for l in lines)


# --------------------------------------------------------------------------- the awkward edges
def test_a_silhouette_that_cannot_be_computed_is_omitted_not_faked(monkeypatch):
    """sklearn raises on degenerate input, and a fabricated zero would rank that setting against real
    partitions. Patched on sklearn itself because _score imports it inside the function."""
    import sklearn.metrics

    def boom(*a, **k):
        raise ValueError("degenerate")

    monkeypatch.setattr(sklearn.metrics, "silhouette_score", boom)
    assert np.isnan(CL._score(_blobs(), _labels())["silhouette"])


def test_cramers_v_declines_a_table_scipy_rejects(monkeypatch):
    import starplast.clustering as mod

    def boom(*a, **k):
        raise ValueError("bad table")

    monkeypatch.setattr("scipy.stats.chi2_contingency", boom)
    assert np.isnan(mod._cramers_v(np.array([[5, 5], [5, 5]])))


def test_a_fisher_test_that_cannot_run_leaves_the_row_without_a_p_value(monkeypatch):
    """The row is still reported -- the counts are real -- but with no p rather than a made-up one."""
    monkeypatch.setattr("scipy.stats.fisher_exact",
                        lambda *a, **k: (_ for _ in ()).throw(ValueError("bad")))
    _, rows = CL.categorical_feature(_labels(), pd.Series(["in"] * 90 + ["out"] * 90))
    assert rows and all(np.isnan(r["p"]) for r in rows)


def test_association_is_measured_both_ways_round_for_continuous_inputs():
    """A continuous input against a continuous column uses correlation; the guard has to see it whether
    the input or the candidate is the numeric one."""
    rng = np.random.default_rng(3)
    x = rng.normal(size=200)
    nodes = pd.DataFrame({"used_num": x, "copy_num": x * 2 + 1,
                          "used_cat": ["a"] * 100 + ["b"] * 100})
    a = CL._association_with_inputs(nodes, {"used_num"})
    assert a["copy_num"] == pytest.approx(1.0, abs=1e-6)


def test_a_column_with_too_few_shared_observations_is_not_compared():
    nodes = pd.DataFrame({"used": [1.0] * 5 + [np.nan] * 195,
                          "other": np.arange(200.0)})
    assert "other" not in CL._association_with_inputs(nodes, {"used"})


def test_a_used_column_that_is_not_in_the_table_is_ignored():
    nodes = pd.DataFrame({"a": np.arange(50.0)})
    assert CL._association_with_inputs(nodes, {"not_a_column"}) == {}


def test_too_many_categories_are_skipped_on_both_sides():
    """A 4,000-level column is an identifier, and crosstabbing it would build a huge useless table."""
    nodes = pd.DataFrame({"used": [f"c{i}" for i in range(100)],
                          "other": [f"d{i}" for i in range(100)]})
    assert CL._association_with_inputs(nodes, {"used"}, max_categories=30) == {}


def test_an_association_that_raises_is_skipped_rather_than_fatal(monkeypatch):
    monkeypatch.setattr("numpy.corrcoef",
                        lambda *a, **k: (_ for _ in ()).throw(RuntimeError("boom")))
    nodes = pd.DataFrame({"used": np.arange(100.0), "other": np.arange(100.0)})
    assert CL._association_with_inputs(nodes, {"used"}) == {}


def test_derived_from_names_near_perfect_restatements():
    x = np.arange(200.0)
    nodes = pd.DataFrame({"used": x, "twin": x, "unrelated": np.random.default_rng(0).normal(size=200)})
    out = CL._derived_from(nodes, {"used"})
    assert "twin" in out and "unrelated" not in out


def test_a_feature_with_more_categories_than_allowed_is_skipped_by_the_battery():
    nodes = pd.DataFrame({"gene_id": [f"g{i}" for i in range(180)],
                          "many": [f"c{i}" for i in range(180)],
                          "real": ["in"] * 90 + ["out"] * 90})
    res, _ = CL.battery(nodes, _labels(), max_categories=30, log=lambda *_: None)
    assert "many" not in set(res.feature)
    assert "real" in set(res.feature)


def test_a_battery_over_nothing_returns_empty_frames():
    nodes = pd.DataFrame({"gene_id": [f"g{i}" for i in range(180)]})
    res, per = CL.battery(nodes, _labels(), log=lambda *_: None)
    assert res.empty and per.empty


def test_describe_says_nothing_when_nothing_clears_the_bar():
    """Silence is the correct output for a clustering that organises nothing, and inventing a sentence
    about the least-bad feature would be worse than saying nothing."""
    nodes = pd.DataFrame({"gene_id": [f"g{i}" for i in range(180)],
                          "noise": np.random.default_rng(0).normal(size=180)})
    res, per = CL.battery(nodes, _labels(), log=lambda *_: None)
    assert CL.describe(res, per, min_score=0.99) == []


def test_describe_reports_a_categorical_and_a_continuous_finding():
    rng = np.random.default_rng(0)
    nodes = pd.DataFrame({"gene_id": [f"g{i}" for i in range(180)],
                          "cat": ["in"] * 90 + ["out"] * 90,
                          "num": np.concatenate([rng.normal(0, 1, 90), rng.normal(8, 1, 90)])})
    res, per = CL.battery(nodes, _labels(), log=lambda *_: None)
    lines = CL.describe(res, per, min_score=0.0, q_max=1.0)
    assert any("cat" in l for l in lines)
    assert any("num" in l for l in lines)


def test_a_correlation_ratio_that_declines_leaves_no_association_recorded():
    """A categorical input with more levels than allowed: the ratio is undefined rather than zero, and
    recording zero would claim the comparison was made and came back clean."""
    nodes = pd.DataFrame({"used_cat": [f"c{i % 50}" for i in range(200)],
                          "num": np.arange(200.0)})
    assert "num" not in CL._association_with_inputs(nodes, {"used_cat"}, max_categories=30)


def test_a_feature_that_cannot_be_scored_at_all_is_absent_from_the_battery():
    """One category everywhere: nothing to separate, so no row rather than a row of zeros."""
    nodes = pd.DataFrame({"gene_id": [f"g{i}" for i in range(180)],
                          "same": ["one_value"] * 180,
                          "real": ["in"] * 90 + ["out"] * 90})
    res, _ = CL.battery(nodes, _labels(), log=lambda *_: None)
    assert "same" not in set(res.feature)


def test_describe_over_an_empty_battery_returns_nothing():
    assert CL.describe(pd.DataFrame(), pd.DataFrame()) == []


# --------------------------------------------------------------------------- per category
def test_the_per_category_table_takes_each_category_s_best_cluster():
    """One number per feature hides the thing worth knowing: V = 0.2 describes 27 compartments
    weakly smeared across every cluster and one compartment falling out cleanly, and those are
    entirely different findings."""
    detail = pd.DataFrame([
        {"feature": "compartment", "category": "apicoplast", "cluster": 0, "n_in_cluster": 5,
         "precision": 0.2, "recall": 0.2, "evidence": "held_out", "q": 0.5},
        {"feature": "compartment", "category": "apicoplast", "cluster": 1, "n_in_cluster": 40,
         "precision": 0.9, "recall": 0.8, "evidence": "held_out", "q": 1e-9},
        {"feature": "compartment", "category": "nucleus", "cluster": 0, "n_in_cluster": 10,
         "precision": 0.1, "recall": 0.1, "evidence": "held_out", "q": 0.4},
    ])
    out = CL.per_category(detail, min_in_cluster=1)
    assert list(out.category) == ["apicoplast", "nucleus"], "not ranked by how well each is matched"
    best = out.iloc[0]
    assert best.cluster == 1 and best.f1 == pytest.approx(2 * 0.9 * 0.8 / 1.7)


def test_a_category_matched_by_nothing_scores_zero_rather_than_dividing_by_zero():
    detail = pd.DataFrame([{"feature": "f", "category": "x", "cluster": 0, "n_in_cluster": 0,
                            "precision": 0.0, "recall": 0.0, "evidence": "held_out", "q": 1.0}])
    assert CL.per_category(detail, min_in_cluster=0).f1.iloc[0] == 0.0


def test_only_held_out_features_are_offered_as_evidence():
    """A feature the map was built from separates the clusters by construction, and its per-category
    scores would be the most flattering rows in the table."""
    detail = pd.DataFrame([
        {"feature": "used", "category": "x", "cluster": 0, "n_in_cluster": 9, "precision": 1.0,
         "recall": 1.0, "evidence": "used", "q": 0.0},
        {"feature": "free", "category": "y", "cluster": 1, "n_in_cluster": 9, "precision": 0.5,
         "recall": 0.5, "evidence": "held_out", "q": 0.01},
    ])
    assert list(CL.per_category(detail, min_in_cluster=1).feature) == ["free"]


def test_a_battery_with_no_categorical_features_gives_an_empty_table_with_columns():
    """An empty frame with no columns at all breaks the caller that fills a table from it."""
    for d in (pd.DataFrame(), None,
              pd.DataFrame([{"feature": "n", "cluster": 0, "median": 1.0, "evidence": "held_out"}])):
        out = CL.per_category(d)
        assert list(out.columns)[:3] == ["feature", "category", "cluster"]
        assert out.empty


def test_the_per_category_table_comes_from_a_real_battery():
    """The columns are produced by categorical_feature, so a rename there must not silently empty
    this table."""
    nodes = pd.DataFrame({"gene_id": [f"g{i}" for i in range(180)],
                          "cat": (["in"] * 60 + ["out"] * 120)})
    S, D = CL.battery(nodes, _labels(), log=lambda *_: None)
    out = CL.per_category(D)
    assert set(out.category) == {"in", "out"}
    assert (out.precision.between(0, 1)).all() and (out.recall.between(0, 1)).all()


def test_categories_with_no_precision_at_all_give_an_empty_table_not_a_row_of_nan():
    """Fisher's test declines on a degenerate table, and a row of nan reads as a measured zero."""
    detail = pd.DataFrame([{"feature": "f", "category": "x", "cluster": 0, "n_in_cluster": 3,
                            "precision": np.nan, "recall": np.nan, "evidence": "held_out",
                            "q": np.nan}])
    out = CL.per_category(detail, min_in_cluster=1)
    assert out.empty and list(out.columns)[:2] == ["feature", "category"]


def test_a_battery_without_q_values_still_produces_the_table():
    """`q` is added only when the battery had p-values to correct, so the column can be absent --
    and a table that vanishes because one column is missing loses the whole result."""
    detail = pd.DataFrame([{"feature": "f", "category": "x", "cluster": 1, "n_in_cluster": 9,
                            "precision": 0.9, "recall": 0.5, "evidence": "held_out"}])
    out = CL.per_category(detail, min_in_cluster=1)
    assert len(out) == 1 and np.isnan(out.q.iloc[0])
    assert out.f1.iloc[0] == pytest.approx(2 * 0.9 * 0.5 / 1.4)


def test_a_category_that_is_most_of_the_data_is_not_reported_as_a_discovery():
    """One cluster holding nearly everything "recovers" a dominant class at F1 0.95 while telling
    you nothing. Lift is what says so: 1.0 means the cluster is no more that category than the map
    is, and the table is ranked by it."""
    detail = pd.DataFrame([
        {"feature": "f", "category": "common", "cluster": 0, "n_in_cluster": 90, "precision": 0.9,
         "recall": 1.0, "evidence": "held_out", "q": 1e-12},
        {"feature": "f", "category": "rare", "cluster": 1, "n_in_cluster": 8, "precision": 0.8,
         "recall": 0.8, "evidence": "held_out", "q": 1e-6},
        {"feature": "f", "category": "rare", "cluster": 0, "n_in_cluster": 2, "precision": 0.02,
         "recall": 0.2, "evidence": "held_out", "q": 0.9},
    ])
    out = CL.per_category(detail)
    assert list(out.category) == ["rare", "common"], "the majority class outranked a real one"
    common = out[out.category == "common"].iloc[0]
    assert common.prevalence == pytest.approx(0.9) and common.lift == pytest.approx(1.0)
    rare = out[out.category == "rare"].iloc[0]
    assert rare.lift > 5


def test_a_feature_with_no_genes_at_all_has_no_prevalence_rather_than_zero():
    detail = pd.DataFrame([{"feature": "f", "category": "x", "cluster": 0, "n_in_cluster": 0,
                            "precision": 0.5, "recall": 0.5, "evidence": "held_out", "q": 0.1}])
    out = CL.per_category(detail, min_in_cluster=0)
    assert np.isnan(out.prevalence.iloc[0]) and np.isnan(out.lift.iloc[0])


def test_a_cluster_holding_one_gene_of_a_category_is_not_a_finding():
    """Rendered and looked at: ranked by lift, the top of this table was a cluster holding ONE
    apicoplast protein at 5x enrichment -- true, meaningless, and indistinguishable at a glance from
    a real result. It is the singleton exploit in its per-category form."""
    detail = pd.DataFrame([
        {"feature": "f", "category": "rare", "cluster": 0, "n_in_cluster": 1, "precision": 0.17,
         "recall": 0.07, "evidence": "held_out", "q": 0.2},
        {"feature": "f", "category": "common", "cluster": 0, "n_in_cluster": 60, "precision": 0.6,
         "recall": 0.9, "evidence": "held_out", "q": 1e-9},
    ])
    assert list(CL.per_category(detail).category) == ["common"]
    assert list(CL.per_category(detail, min_in_cluster=1).category) == ["rare", "common"]


def test_everything_below_the_floor_gives_an_empty_table_with_its_columns():
    detail = pd.DataFrame([{"feature": "f", "category": "x", "cluster": 0, "n_in_cluster": 2,
                            "precision": 1.0, "recall": 1.0, "evidence": "held_out", "q": 0.1}])
    out = CL.per_category(detail, min_in_cluster=5)
    assert out.empty and "lift" in out.columns
