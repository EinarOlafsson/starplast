"""Strategies at the edges: tables missing what a strategy needs, and the guards that say so.

Every strategy is written against the shipped tables, where a crosslink layer, a second organism
and a literature layer all exist. A user's own table, or the other arm, may have none of them, and
a strategy must then say what is missing rather than fail with a traceback or -- worse -- return
an empty result that reads as "nothing found". Each test here removes one thing and checks the
answer is a refusal, an explanation or an inconclusive verdict.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import strategies as S  # noqa: E402
from starplast import strategy_catalog as C  # noqa: E402


@pytest.fixture(scope="module")
def planted():
    return S.planted_context()


@pytest.fixture(scope="module")
def planted_big():
    return S.planted_context(n=1500)


def _bare(n=60, seed=0, **extra):
    """A table with a label and a few measurements, and no graph and no second organism."""
    rng = np.random.default_rng(seed)
    frame = {"gene_id": [f"G{i:04d}" for i in range(n)],
             "label": rng.choice(["a", "b", "c"], size=n),
             "m1": rng.normal(size=n), "m2": rng.normal(size=n), "m3": rng.normal(size=n)}
    frame.update(extra)
    nodes = pd.DataFrame(frame)
    return S.Context(nodes, graph={}, other=pd.DataFrame({"gene_id": ["X"], "orthogroup": ["O"]}),
                     organism="Tg")


# --------------------------------------------------------------------------- the context
def test_a_table_without_products_orthogroups_or_columns_still_answers():
    ctx = _bare()
    assert ctx.product([0, 1]) == ["", ""]
    assert len(set(ctx.groups())) == ctx.n
    assert ctx.matrix([]).shape == (ctx.n, 0)
    pos, missing = ctx.resolve_genes(",,G0001,,")
    assert list(pos) == [1] and missing == []


def test_a_closure_that_cannot_be_computed_bans_the_label_and_says_so(monkeypatch):
    from starplast import search
    ctx = _bare()
    said = []
    ctx.log = said.append
    monkeypatch.setattr(search, "excluded_for", lambda *a, **k: 1 / 0)
    monkeypatch.setattr(search, "excluded_edges", lambda *a, **k: 1 / 0)
    assert ctx.banned("label") == {"label"} and said
    assert ctx.banned_layers("label") == set()


def test_a_table_the_catalogue_cannot_describe_is_grouped_by_prefix(monkeypatch):
    from starplast import embedding, slots
    monkeypatch.setattr(embedding, "default_spec", lambda *a, **k: 1 / 0)
    monkeypatch.setattr(slots, "all_slots", lambda *a, **k: 1 / 0)
    ctx = _bare()
    assert set(ctx.blocks()) == {"m1", "m2", "m3"}
    assert ctx.family_of("m1") == "m1" and ctx.family_of("Tg_fitness_x") == "fitness"


def test_the_analytic_null_reports_its_spread():
    import time
    r = S.judge("k", "m", 0.5, [], min_effect=0.1, n_hidden=50, hidden="h", null_kind="n",
                t0=time.monotonic(), analytic=(0.2, 0.03))
    assert r.null_sd == 0.03


def test_too_little_to_hide_is_inconclusive_in_every_pattern(planted):
    t = S.label_transfer_test(planted, "k", "compartment", lambda vis: vis,
                              restrict=np.zeros(planted.n, bool))
    assert t.verdict == "INCONCLUSIVE"
    t = S.set_expansion_test(planted, "k", [1, 2, 3], lambda q: np.zeros(planted.n))
    assert t.verdict == "INCONCLUSIVE"
    t = S.value_test(planted, "k", pd.Series([1.0, 2.0] + [np.nan] * (planted.n - 2)),
                     lambda v: np.zeros(planted.n))
    assert t.verdict == "INCONCLUSIVE"


# --------------------------------------------------------------------------- defaults and fallbacks
def test_defaults_fall_back_when_the_preferred_columns_are_absent():
    ctx = _bare(extra1=np.arange(60) % 7 * 1.0)
    assert C._layer_default(["xlms"])(ctx) is None
    assert C._literature_layer_default(ctx) is None
    assert C._second_label(ctx, "label") is None
    assert C._second_label(_bare(other_label=np.array(["x", "y"] * 30)), "label") == "other_label"
    assert C._family_default(ctx) in ctx.families()
    assert C._condition_default("condition")(ctx) in ctx.numeric_targets()
    assert C._condition_default("baseline")(ctx) in ctx.numeric_targets()
    assert C._transfer_target_default(ctx) == S.default_numeric(ctx)
    assert C._source_default(ctx) is None                         # the other has no measurement
    lonely = _bare()
    lonely._other = False
    assert C._source_default(lonely) is None
    rich = _bare()
    rich._other = S.Context(pd.DataFrame({"gene_id": ["X", "Y", "Z"] * 20,
                                          "orthogroup": ["O"] * 60,
                                          "other_value": np.arange(60) * 1.0}), graph={})
    assert C._source_default(rich) == "other_value"


def test_a_list_strategy_finds_an_example_in_another_label_or_says_there_is_none():
    ctx = _bare(n=1000, big=np.where(np.arange(1000) < 60, "set", "rest"))
    members, what, target = C._list_or_example(ctx, {"genes": "", "exclude": "label"})
    assert target == "big" and len(members) == 60
    members, what, target = C._list_or_example(_bare(), {"genes": ""})
    assert len(members) == 0 and what == "no known set"
    assert C._list_or_example(ctx, {"genes": ""})[2] == "big"
    r = S.get("geneset_hunt").test(_bare(), genes="")
    assert r.verdict == "INCONCLUSIVE"


# --------------------------------------------------------------------------- small helpers
def test_the_hidden_f1_of_a_mapping_with_nothing_placed_is_undefined():
    per = pd.DataFrame({"label": ["a"], "cluster": [0]})
    assert np.isnan(C._hidden_f1(np.array([0, 0]), [0, 1], per, pd.Series(["a", "a"]), [5]))
    assert C._hidden_f1(np.array([0, 1]), [0, 1], per.iloc[:0], pd.Series(["a", "b"]), [0, 1]) == 0


def test_an_atlas_skips_classes_too_small_to_score():
    coords = np.random.default_rng(0).normal(size=(20, 3))
    lab = np.array(["a"] * 17 + ["b"] * 3, dtype=object)
    table = C._atlas(coords, lab, np.arange(20), 5)
    assert list(table["category"]) == ["a"]


def test_the_atlas_test_is_inconclusive_with_too_few_labels():
    ctx = _bare(n=80)
    r = S.get("recoverability_atlas").test(ctx, target="label", sample=0)
    assert r.verdict in ("INCONCLUSIVE", "FAIL")


def test_the_battery_skips_what_it_cannot_test():
    labels = np.array([0] * 10 + [1] * 10)
    frame = pd.DataFrame({"few": [np.nan] * 15 + ["x"] * 5, "flat": [1.0] * 20,
                          "one_class": ["x"] * 20})
    assert C._associations(labels, frame, ["few", "one_class"], ["flat"]).empty


def test_the_battery_refuses_a_family_that_does_not_exist(planted):
    with pytest.raises(ValueError, match="no measurement family"):
        S.get("blind_battery").run(planted, map_from="nope")


def test_cross_validation_needs_enough_labels():
    assert np.isnan(C._cv_knn(np.zeros((4, 1)), pd.Series(["a", "b", np.nan, np.nan]),
                              np.arange(4)))


def test_ablation_can_score_every_block(planted):
    r = S.get("block_ablation").run(planted, unit="block")
    assert set(r.tables["evidence"]["unit"]) == set(planted.blocks("compartment"))


def test_enrichment_ignores_clusters_with_too_few_labelled_members():
    pred, table = C._enriched(np.array([0, 0, 1, 1, 1, 1]),
                              pd.Series([np.nan, "a", "a", "a", "a", "b"]), 1.0)
    assert table["cluster"].tolist() == [1] or table.empty or 0 not in table["cluster"].tolist()


def test_outliers_need_labels_to_corrupt():
    r = S.get("label_outliers").test(_bare(n=30), target="label")
    assert r.verdict == "INCONCLUSIVE"


def test_ec_numbers_are_cut_and_other_labels_kept():
    t = pd.Series(["3.1.3.16 (x); 2.7.1.1", "kinase", np.nan, "4.2"], dtype=object)
    assert C._ec_level(t, 2).tolist()[:2] == ["3.1", "kinase"]
    assert C._ec_level(t, 0).tolist()[0] == "3.1.3.16 (x)"
    assert pd.isna(C._ec_level(t, 1).iloc[2])


# --------------------------------------------------------------------------- missing layers
def test_network_strategies_say_what_is_missing():
    ctx = _bare()
    r = S.get("physical_partners").run(ctx, target="label")
    assert not r.ok and "crosslink" in r.summary
    assert S.get("physical_partners").test(ctx, target="label").verdict == "INCONCLUSIVE"
    with pytest.raises(ValueError, match="two permitted"):
        S.get("multiplex_modules").run(ctx, target="label")
    with pytest.raises(ValueError, match="no edge layer"):
        S.get("link_prediction").run(ctx, layer="xlms")
    with pytest.raises(ValueError, match="no edge layer"):
        S.get("link_prediction").test(ctx, layer="xlms")
    with pytest.raises(ValueError, match="co-mention"):
        S.get("attention_correction").run(ctx, target="label", layer=None)
    r = S.get("unwritten_links").run(ctx)
    assert not r.ok and "measurement layers" in r.summary
    assert S.get("unwritten_links").test(ctx).verdict == "INCONCLUSIVE"
    assert not S.get("paralog_divergence").run(ctx).ok
    assert S.get("paralog_divergence").test(ctx, target="label").verdict == "INCONCLUSIVE"


def test_link_prediction_caps_its_candidates_and_says_when_there_are_none(planted, monkeypatch):
    monkeypatch.setattr(C, "LINK_CANDIDATES", 5)
    r = S.get("link_prediction").run(planted, layer="xlms", top=3)
    assert len(r.tables["predicted links"]) == 3
    graph = {"only__a": np.array([0]), "only__b": np.array([1]), "only__w": np.array([1.0])}
    lonely = S.Context(_bare().nodes, graph=graph)
    assert not S.get("link_prediction").run(lonely, layer="only").ok
    assert S.get("link_prediction").test(lonely, layer="only").verdict == "INCONCLUSIVE"


def test_attention_needs_enough_labelled_pairs():
    graph = {"comention__a": np.array([0, 1]), "comention__b": np.array([1, 2]),
             "comention__w": np.array([1.0, 2.0])}
    ctx = S.Context(_bare().nodes, graph=graph)
    r = S.get("attention_correction").test(ctx, target="label", layer="comention")
    assert r.verdict == "INCONCLUSIVE"


def test_unwritten_links_needs_the_literature_too():
    graph = {"coexpression__a": np.array([0, 1]), "coexpression__b": np.array([1, 2]),
             "coexpression__w": np.array([1.0, 1.0])}
    ctx = S.Context(_bare().nodes, graph=graph)
    assert S.get("unwritten_links").test(ctx).verdict == "INCONCLUSIVE"


# --------------------------------------------------------------------------- learning
def test_a_classifier_with_one_class_calls_nothing():
    pred, prob, model = C._logistic(np.zeros((4, 1)), pd.Series(["a", "a", np.nan, np.nan]), [2, 3])
    assert model is None and pred.isna().all()


def test_a_regression_with_too_few_values_predicts_nothing():
    assert np.isnan(C._fit_predict(np.zeros((5, 1)), pd.Series([1.0] * 5), "ridge")).all()


def test_ridge_is_offered_beside_boosting(planted):
    r = S.get("trait_regression").run(planted, model="ridge")
    assert r.ok and np.isfinite(r.numbers["oof_spearman"])


def test_a_shift_needs_genes_measured_in_both():
    ctx = _bare(c=[np.nan] * 60)
    with pytest.raises(ValueError, match="fewer than 20"):
        C._shift(ctx, "c", "m1")


def test_a_profile_skips_what_it_cannot_test(planted):
    few = planted.gene_ids[:2]
    prof = C._profile(planted, planted.resolve_genes(list(few))[0], "compartment")
    assert "compartment" not in set(prof.get("feature", []))
    graph = {"empty__a": np.array([], dtype=int), "empty__b": np.array([], dtype=int)}
    ctx = S.Context(_bare().nodes, graph=graph)
    assert "empty" not in set(C._profile(ctx, [0, 1, 2], None).get("feature", []))


# --------------------------------------------------------------------------- contrasts
def test_contrasts_refuse_the_same_column_twice_and_numbers_for_a_conjunction(planted):
    with pytest.raises(ValueError, match="two different"):
        S.get("split_clusters").run(planted, a="compartment", b="compartment")
    with pytest.raises(ValueError, match="two categorical"):
        S.get("conjunctions").run(planted, a="compartment", b="fit_invitro_hff")
    with pytest.raises(ValueError, match="two categorical"):
        S.get("conjunctions").test(planted, a="compartment", b="fit_invitro_hff")


def test_a_split_on_a_measurement_is_replicated_as_two_modes_wider_than_chance():
    rng = np.random.default_rng(0)
    labels = np.r_[np.zeros(40, int), np.ones(80, int)]
    a = np.array(["x"] * 120, dtype=object)
    finding = {"cluster": 0, "category": "x"}
    bimodal = np.r_[rng.normal(0, 0.3, 20), rng.normal(10, 0.3, 20), rng.normal(5, 3.0, 80)]
    assert C.split_replicates(finding, labels, a, bimodal, np.arange(120), True)
    assert not C.split_replicates(finding, labels, a, bimodal, np.arange(5), True)
    # A cluster whose values are just a draw of the rest splits in two as well -- two-means always
    # does -- and that is exactly what must NOT count as replicating.
    noise = rng.normal(size=120)
    assert not C.split_replicates(finding, labels, a, noise, np.arange(120), True)


def test_a_split_on_a_measurement_runs_end_to_end(planted_big):
    r = S.get("split_clusters").run(planted_big, b="fit_invitro_hff")
    assert "split fit_invitro_hff" in " ".join(r.tables)


def test_paralogs_come_from_the_orthogroup_column_when_there_is_no_layer():
    ogs = ["O1", "O1", "O2", "O2", "O2"] + [f"S{i}" for i in range(55)]
    ctx = _bare(orthogroup=ogs)
    assert len(C._paralog_pairs(ctx)) == 1 + 3


# --------------------------------------------------------------------------- across species
def test_orthologs_need_the_other_table_its_column_and_orthogroups():
    ctx = _bare()
    ctx._other = False
    with pytest.raises(ValueError, match="not available"):
        C._through_orthologs(ctx, "x")
    ctx = _bare(orthogroup=["O"] * 60)
    with pytest.raises(ValueError, match="no column"):
        C._through_orthologs(ctx, "nope")
    ctx._other = S.Context(pd.DataFrame({"gene_id": ["X"], "v": [1.0]}), graph={})
    with pytest.raises(ValueError, match="orthogroup"):
        C._through_orthologs(ctx, "v")


def test_labels_cross_species_in_every_combination(planted):
    r = S.get("ortholog_transfer").test(planted, target="compartment",
                                        source="localization_other")
    assert r.verdict == "PASS", r.summary()
    src = pd.Series([1.0, 2.0, 3.0, np.nan, 5.0, 6.0])
    cat = pd.Series(["a", "a", "b", "b", np.nan, "a"], dtype=object)
    num = pd.Series([1.0, 1.5, 3.0, 3.5, np.nan, 1.2])
    assert C._transfer_map(cat, num, False, True)(cat).notna().sum() >= 4
    assert C._transfer_map(src, cat, True, False)(src).notna().sum() >= 4
    assert C._transfer_map(cat, cat, False, False)(cat).notna().sum() >= 4
    r = S.get("ortholog_transfer").run(planted, target="compartment",
                                       source="localization_other")
    assert r.ok


# --------------------------------------------------------------------------- strata
def test_strata_are_read_from_whatever_the_table_carries():
    ctx = _bare(ortholog_number=[0] * 10 + [5] * 50, product=["hypothetical protein"] * 5
                + ["kinase"] * 55, attention_depth=["focal"] * 20 + [""] * 40)
    assert C.stratum_mask(ctx, "lineage-specific").sum() == 10
    assert C.stratum_mask(ctx, "conserved").sum() == 50
    assert C.stratum_mask(ctx, "hypothetical protein").sum() == 5
    assert C.stratum_mask(ctx, "understudied").sum() == 40
    bare = _bare()
    assert C.stratum_mask(bare, "lineage-specific").sum() == 0
    assert C.stratum_mask(bare, "understudied").all()
    with pytest.raises(ValueError, match="unknown stratum"):
        C.stratum_mask(bare, "famous")


def test_the_family_walk_tries_each_kind_of_evidence(planted):
    r = S.get("holdout_search").run(planted, features="families", n_neighbors="15",
                                    min_dist="0.1", min_cluster_size="20", selection="eom")
    assert len(set(r.tables["configurations"]["features"])) > 1


def test_a_map_needs_something_left_to_build_from():
    ctx = _bare(copy=["a", "b", "c"] * 20)
    ctx._cache[("banned", "label")] = {"label", "m1", "m2", "m3"}
    with pytest.raises(ValueError, match="nothing is left"):
        C._blocksets(ctx, "label", "all")


# --------------------------------------------------------------------------- the last guards
def test_an_unknown_way_of_joining_communities_is_refused():
    from starplast import methods
    with pytest.raises(ValueError, match="join must be"):
        methods.multiplex_communities(["a"], 3, graph={}, join="union")


def test_the_atlas_needs_labels_to_learn_from_and_two_categories_to_rank():
    few = _bare(n=40)
    assert S.get("recoverability_atlas").test(few, target="label", sample=0).verdict \
        == "INCONCLUSIVE"
    lopsided = _bare(n=120, label=np.array(["a"] * 112 + ["b"] * 8))
    assert S.get("recoverability_atlas").test(lopsided, target="label",
                                              sample=0).verdict == "INCONCLUSIVE"


def test_link_prediction_with_nothing_left_to_predict():
    n = 120
    matching = {"m__a": np.arange(0, n, 2), "m__b": np.arange(1, n, 2),
                "m__w": np.ones(n // 2)}
    ctx = S.Context(_bare(n=n).nodes, graph=matching)
    r = S.get("link_prediction").run(ctx, layer="m")
    assert not r.ok and "No candidate" in r.summary
    k = 12
    a, b = np.triu_indices(k, 1)
    clique = {"c__a": a, "c__b": b, "c__w": np.ones(len(a))}
    ctx = S.Context(_bare(n=k).nodes, graph=clique)
    r = S.get("link_prediction").run(ctx, layer="c")
    assert not r.ok and "already an edge" in r.summary


def test_a_profile_of_one_gene_skips_the_network_test():
    graph = {"g__a": np.array([0, 1, 2]), "g__b": np.array([1, 2, 3]), "g__w": np.ones(3)}
    ctx = S.Context(_bare().nodes, graph=graph)
    prof = C._profile(ctx, [0], None)
    assert "g" not in set(prof.get("feature", []))


@pytest.fixture(scope="module")
def app():
    # Held by the fixture: an application created and not referenced is collected at once, and
    # the next widget built without one aborts the interpreter.
    from PyQt6 import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


def test_the_panel_offers_nothing_from_an_organism_that_is_not_there(app, monkeypatch):
    from starplast.strategy_panel import StrategyPanel
    monkeypatch.setattr(S.Context, "shipped", classmethod(lambda cls, *a, **k: 1 / 0))
    p = StrategyPanel(_bare().nodes, graph={}, organism="Tg")
    p.select("ortholog_transfer")
    _param, combo = p.inputs["source"]
    assert combo.count() == 0
    p.deleteLater()
