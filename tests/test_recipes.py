#!/usr/bin/env python3
"""Recipes: a question, its inputs, its holdout, and the control that says whether to believe it.

These tests are written as ATTACKS wherever a guard exists. A recipe is a hand-picked set of inputs,
which is precisely the situation the level sweep protected us from automatically, so the interesting
question is never "does it run" -- it is "can I get it to answer a question with the answer".
"""
import json
import os

import numpy as np
import pandas as pd
import pytest

from starplast import recipes as R


# --------------------------------------------------------------------------- fixtures
@pytest.fixture(scope="module")
def nodes():
    """A real slice of the shipped table. Real, because a recipe addresses the real hierarchy."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(root, "starplast", "data", "nodes.parquet")
    if not os.path.exists(path):
        pytest.skip("built node table not present")
    return pd.read_parquet(path).sample(n=400, random_state=0).reset_index(drop=True)


@pytest.fixture(scope="module")
def full_nodes():
    """The WHOLE table. The graph's edge endpoints are positions in it, so anything that walks the
    graph has to be given it entire; every other test here uses the 400-row sample for speed."""
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(root, "starplast", "data", "nodes.parquet")
    if not os.path.exists(path):
        pytest.skip("built node table not present")
    return pd.read_parquet(path)


RNA = ("biology", ("gene expression", "RNA abundance"))
FITNESS = ("biology", ("parasite phenotype", "fitness and essentiality"))
LOCALIZATION = ("biology", ("cell organization", "localization and topology"))


# --------------------------------------------------------------------------- the object
def test_an_address_survives_the_round_trip_through_json():
    """JSON has no tuples. A recipe reloaded from disk must address the same branch it was saved
    with, or a saved question quietly becomes a different question."""
    r = R.Recipe(question="q", inputs=[RNA, "a_block"], holdout="compartment")
    back = R.Recipe.from_dict(json.loads(json.dumps(r.to_dict())))
    assert back.inputs == r.inputs == (RNA, "a_block")
    assert back.to_dict() == r.to_dict()


def test_a_one_part_path_may_be_written_as_a_bare_string():
    assert R.Recipe(question="q", inputs=[("biology", "metabolism")]).inputs == \
           (("biology", ("metabolism",)),)


def test_an_unknown_field_in_a_saved_recipe_is_ignored_rather_than_fatal():
    """A recipe saved by a later version must still load, or an upgrade destroys saved work."""
    r = R.Recipe.from_dict({"question": "q", "inputs": [], "invented_later": 7})
    assert r.question == "q" and not hasattr(r, "invented_later")


def test_the_question_is_the_name():
    assert R.Recipe(question="Where does GRA16 go?").name == "Where does GRA16 go?"


def test_an_address_reads_as_one_line_and_resolves_to_slots():
    assert R.address_label(RNA) == "biology: gene expression > RNA abundance"
    assert R.address_label("a_block") == "a_block"
    assert len(R.address_blocks(RNA)) > 1
    assert R.address_blocks("a_block") == ("a_block",)


# --------------------------------------------------------------------------- labels
def test_a_quantity_is_unusable_as_a_holdout_until_it_is_binned():
    """The failure this prevents is silent: scored raw, a continuous column produces one class per
    gene, nothing clears the floor, and the recipe reports that nothing recovered -- for a question
    the data could have answered."""
    frame = pd.DataFrame({"x": np.linspace(0, 1, 200)})
    assert R.label_series(frame, "x").nunique() == 200
    binned = R.label_series(frame, "x", bins=4)
    assert binned.nunique() == 4 and (binned.value_counts() >= 15).all()


def test_unassigned_is_a_gene_waiting_to_be_named_not_a_class():
    """hyperLOPIT writes `unassigned` where it could not place a protein. Left as a string it becomes
    the commonest label in the table, clusters get "recovered" as unassigned, and the genes this
    module exists to name are counted as already labelled and never predicted."""
    frame = pd.DataFrame({"compartment": ["nucleus", "unassigned", "", None, "cytosol"]})
    got = R.label_series(frame, "compartment")
    assert list(got.notna()) == [True, False, False, False, True]


def test_a_binned_quantitys_missing_values_do_not_become_a_bin():
    frame = pd.DataFrame({"x": [1.0, 2.0, 3.0, 4.0, np.nan]})
    assert R.label_series(frame, "x", bins=2).isna().sum() == 1


# --------------------------------------------------------------------------- closure
def test_a_recipe_that_predicts_a_label_from_itself_is_refused_not_trimmed(nodes):
    """The sharpest attack on this module. Asked to predict compartment FROM localisation, a builder
    that trimmed silently would return a smaller recipe answering a question nobody asked."""
    c = R.close(nodes, R.Recipe(question="q", inputs=[LOCALIZATION], holdout="compartment"))
    assert not c.ok
    assert "restate a holdout" in c.refusal and "localization and topology" in c.refusal


def test_class_scope_is_what_closes_the_sequence_predictors(nodes):
    """Measured, not assumed, and it is why `scope` does not default to the sweep's `target_family`.

    `target_family` closes the family estimating the same quantity and deliberately KEEPS sequence
    predictors -- there is a test in test_search.py asserting exactly that. A localisation recipe
    scoped that way accepts membrane topology as an input and recovers hyperLOPIT by predicting it
    from signal peptide and TM count, which is the failure this module exists to prevent.
    """
    ask = dict(question="q", inputs=[LOCALIZATION], holdout="compartment")
    assert R.close(nodes, R.Recipe(**ask, scope="target_family")).ok
    assert not R.close(nodes, R.Recipe(**ask, scope="biology")).ok


def test_a_control_that_is_a_copy_of_the_primary_is_refused(nodes):
    """Six columns on this table are views of one hyperLOPIT experiment. Any pair of them would look
    like a question with an independent control."""
    c = R.close(nodes, R.Recipe(question="q", inputs=[RNA], holdout="compartment",
                                validation_holdout="lopit_map"))
    assert not c.ok and "controls nothing" in c.refusal


def test_a_holdout_that_is_not_a_column_is_refused_before_anything_is_built(nodes):
    assert "not a column" in R.close(nodes, R.Recipe(
        question="q", inputs=[RNA], holdout="no_such_column")).refusal
    assert "not a column" in R.close(nodes, R.Recipe(
        question="q", inputs=[RNA], holdout="")).refusal


def test_a_control_that_is_not_a_column_is_refused(nodes):
    c = R.close(nodes, R.Recipe(question="q", inputs=[RNA], holdout="compartment",
                                validation_holdout="no_such_column"))
    assert not c.ok and "validation holdout" in c.refusal


def test_a_recipe_with_no_inputs_is_refused(nodes):
    assert "names no inputs" in R.close(
        nodes, R.Recipe(question="q", inputs=[], holdout="compartment")).refusal


def test_a_block_with_no_columns_here_is_reported_rather_than_silently_missing(nodes):
    c = R.close(nodes, R.Recipe(question="q", inputs=[RNA, "not_a_real_block"],
                                holdout="compartment"))
    assert c.dropped_blocks["not_a_real_block"] == "no columns in this cache"


def test_closure_says_which_mechanism_removed_each_column(nodes):
    """"Excluded" is not a report. Which of the five mechanisms caught it is the fact a reader
    checks, and it is what distinguishes a copy of the target from a same-experiment sibling."""
    c = R.close(nodes, R.Recipe(question="q", inputs=[FITNESS], holdout="compartment",
                                validation_holdout="cellcycle_phase"))
    assert c.ok, c.refusal
    assert any(w.startswith("primary:") for w in c.excluded.values())
    assert any(w.startswith("control:") for w in c.excluded.values())
    assert "excluded" in c.report() and str(len(c.blocks)) in c.report()


def test_a_refused_closure_reports_the_refusal_and_nothing_else(nodes):
    assert R.close(nodes, R.Recipe(question="q", inputs=[LOCALIZATION],
                                   holdout="compartment")).report().startswith("REFUSED:")


def test_the_same_block_named_by_two_addresses_is_fed_in_once(nodes):
    """A block weighted twice is a block the map was told to care about twice."""
    c = R.close(nodes, R.Recipe(question="q", inputs=[RNA, RNA], holdout="compartment"))
    assert len(c.blocks) == len(set(c.blocks))


def test_a_cell_cycle_question_may_not_be_asked_of_rna(nodes):
    """A real trap, found by this test suite rather than reasoned about, and it constrains which
    questions instruction 45 may ship: `cellcycle_phase` comes from single-parasite RNA sequencing
    and lives under `gene expression > RNA abundance`, so a recipe feeding transcription while
    holding out cell-cycle phase is asking RNA to predict a label derived from RNA."""
    c = R.close(nodes, R.Recipe(question="q", inputs=[RNA], holdout="cellcycle_phase"))
    assert not c.ok and "RNA abundance" in c.refusal


def test_an_address_that_names_nothing_is_refused_rather_than_passed_over(nodes):
    """A misspelt address resolved to nothing and skipped quietly would build a map from the OTHER
    inputs and answer a question nobody asked."""
    c = R.close(nodes, R.Recipe(question="q", inputs=[RNA, ("biology", ("no such branch",))],
                                holdout="compartment"))
    assert not c.ok and "name nothing in this catalogue" in c.refusal


def test_an_input_absent_from_this_build_is_not_reported_as_circular(nodes):
    """Emptied by the guard and emptied by an empty cache are different facts. Reporting the second
    as the first tells a user their question is circular when the truth is that the data is missing."""
    c = R.close(nodes, R.Recipe(question="q", inputs=["not_a_real_block"], holdout="compartment"))
    assert not c.ok and "no columns in this build" in c.refusal
    assert "restate" not in c.refusal


def test_a_banned_column_reaching_the_matrix_is_caught_even_if_the_lookups_disagree(nodes, monkeypatch):
    """The belt-and-braces guard, attacked by making the thing it guards against actually happen.

    Blocks are dropped whole, so a banned column can only survive if the per-block lookup and the
    combined lookup disagree. That should be impossible -- which is exactly why it is worth one test
    rather than a comment claiming it cannot happen.
    """
    real = R.columns_for

    def disagreeing(frame, spec):
        if len(spec.blocks) == 1:
            return real(frame, spec)
        return {"pretend": ["compartment"]}
    monkeypatch.setattr(R, "columns_for", disagreeing)
    c = R.close(nodes, R.Recipe(question="q", inputs=[RNA], holdout="compartment"))
    assert not c.ok and "survived closure" in c.refusal


# --------------------------------------------------------------------------- the control
def _labels(*sizes):
    return np.concatenate([np.full(n, i) for i, n in enumerate(sizes)])


def test_a_control_enriched_over_its_base_rate_agrees():
    labels = _labels(40, 40)
    control = pd.Series(["host"] * 35 + ["other"] * 5 + ["other"] * 40)
    got = R.dominant_by_cluster(labels, control)
    assert got.loc[0, "agrees"] and got.loc[0, "enrichment"] > 1.9
    assert not got.loc[1, "agrees"]


def test_a_control_sitting_at_its_base_rate_is_not_agreement():
    """The failure mode this exists for: a label carried by 60% of everything is not evidence when
    it is carried by 60% of a cluster."""
    labels = _labels(40, 40)
    got = R.dominant_by_cluster(labels, pd.Series(["host"] * 24 + ["x"] * 16 +
                                                  ["host"] * 24 + ["x"] * 16))
    assert not got.agrees.any() and "sits at" in got.loc[0, "why"]


def test_too_few_control_labels_is_reported_differently_from_disagreement():
    """A cluster nobody could evaluate must not read like a cluster the control refuted."""
    labels = _labels(40, 40)
    control = pd.Series([None] * 36 + ["host"] * 4 + [None] * 40, dtype=object)
    got = R.dominant_by_cluster(labels, control)
    assert not got.loc[0, "agrees"] and "the floor is 15" in got.loc[0, "why"]
    assert got.loc[1, "n_labelled"] == 0
    assert "no gene in this cluster carries a label" in got.loc[1, "why"]


def test_the_base_rate_a_cluster_is_measured_against_includes_that_cluster():
    """There is no zero-denominator guard because there cannot be a zero denominator, and this is the
    test that says so rather than a comment claiming it."""
    labels = _labels(30, 30)
    got = R.dominant_by_cluster(labels, pd.Series([None] * 30 + ["only_here"] * 30, dtype=object))
    assert got.loc[1, "base_rate"] == 1.0 and got.loc[1, "enrichment"] == 1.0
    assert np.isfinite(got.enrichment.dropna()).all()


def test_noise_genes_are_not_a_cluster():
    labels = np.array([-1] * 20 + [0] * 20)
    assert list(R.dominant_by_cluster(labels, pd.Series(["a"] * 40)).cluster) == [0]


# --------------------------------------------------------------------------- per cluster
def test_a_label_split_across_clusters_is_visible_per_cluster_not_only_at_its_best():
    """The reading failure this table exists for. A label cleanly split across three clusters is a
    real finding about sub-structure, and `score_recovery` reports it as one mediocre best-cluster
    F1 -- which looks like a label that did not recover."""
    labels = _labels(30, 30, 30)
    truth = pd.Series(["split"] * 90)
    got = R.recovery_by_cluster(labels, truth)
    assert len(got) == 3 and (got.precision == 1.0).all()
    assert got.recall.sum() == pytest.approx(1.0)


def test_a_cluster_holding_none_of_a_label_is_a_zero_not_a_missing_row():
    """A missing row reads as "not computed"; a zero says the cluster was checked and holds none."""
    labels = _labels(30, 30)
    truth = pd.Series(["a"] * 30 + ["b"] * 30)
    got = R.recovery_by_cluster(labels, truth)
    assert len(got) == 4 and (got.f1 == 0).sum() == 2


def test_a_label_under_the_floor_is_not_scored_per_cluster():
    labels = _labels(30, 30)
    truth = pd.Series(["a"] * 50 + ["rare"] * 10)
    assert set(R.recovery_by_cluster(labels, truth).label) == {"a"}


def test_an_unlabelled_table_scores_nothing_per_cluster():
    assert R.recovery_by_cluster(_labels(20, 20), pd.Series([None] * 40, dtype=object)).empty


# --------------------------------------------------------------------------- running one
def test_a_refused_recipe_returns_a_result_rather_than_raising(nodes):
    """A caller that must catch an exception to learn "this cannot be asked" will sooner or later
    report it as a failure of the program instead of a finding about the data."""
    res = R.run(nodes, R.Recipe(question="q", inputs=[LOCALIZATION], holdout="compartment"),
                log=lambda *a: None)
    assert not res.ok and "restate a holdout" in res.stopped_because
    assert res.recovery.empty and res.inference.empty


def test_a_degenerate_map_returns_the_verdict_and_no_score(nodes, monkeypatch):
    """A recovery number computed on a bisection reads exactly like a real one."""
    monkeypatch.setattr("starplast.clustering.cluster",
                        lambda X, **k: np.zeros(len(X), dtype=int))
    res = R.run(nodes, R.Recipe(question="q", inputs=[RNA, FITNESS], holdout="compartment"),
                tune=False, log=lambda *a: None)
    assert not res.ok and "cannot support a claim" in res.stopped_because
    assert res.summary == {} and res.recovery.empty
    assert res.labels.size and res.genes.any()


def test_a_recipe_no_umap_setting_can_map_says_so(nodes, monkeypatch):
    monkeypatch.setattr(R, "tune_umap", lambda *a, **k: pd.DataFrame({"usable": [False]}))
    res = R.run(nodes, R.Recipe(question="q", inputs=[RNA], holdout="compartment"),
                log=lambda *a: None)
    assert "no UMAP setting" in res.stopped_because


def test_an_empty_tuning_table_is_not_an_exception(nodes, monkeypatch):
    monkeypatch.setattr(R, "tune_umap", lambda *a, **k: pd.DataFrame())
    assert "no UMAP setting" in R.run(nodes, R.Recipe(question="q", inputs=[RNA],
                                                     holdout="compartment"),
                                      log=lambda *a: None).stopped_because


def _fake_clusters(n):
    """Three clusters with real noise -- a partition with no noise is degenerate and scores nothing."""
    lab = np.arange(n) % 3
    lab[::7] = -1
    return lab


def test_a_run_names_genes_and_puts_the_controls_verdict_on_the_same_rows(nodes, monkeypatch):
    """The whole deliverable in one assertion: the inference is genes, and the control that qualifies
    it is on the inference's own rows rather than in a table beside it that nobody reads."""
    lab = _fake_clusters(len(nodes))
    monkeypatch.setattr("starplast.clustering.cluster", lambda X, **k: _fake_clusters(len(X)))
    truth = np.array([None] * len(nodes), dtype=object)
    control = np.array([None] * len(nodes), dtype=object)
    for cl, name in ((0, "pvm"), (1, "cytosol")):
        truth[np.where(lab == cl)[0][:40]] = name      # pure, and well above the floor
    # The control carries evidence about DIFFERENT genes -- the ones the primary leaves unlabelled,
    # which is the real situation and also the only one the independence guard permits: a control
    # built from the primary's own rows is measurably a copy of it and is refused, as it should be.
    control[np.where(lab == 0)[0][40:75]] = "host"
    control[np.where(lab == 1)[0][40:55]] = "not host"
    control[np.where(lab == 2)[0][:40]] = "not host"
    frame = nodes.copy()
    frame["compartment"] = truth
    frame["a_control"] = control
    res = R.run(frame, R.Recipe(question="q", inputs=[RNA, FITNESS], holdout="compartment",
                                validation_holdout="a_control"), tune=False, log=lambda *a: None)
    assert res.ok, res.stopped_because
    assert len(res.inference), "a pure cluster with 40 labelled genes named nobody"
    assert {"control_agrees", "control_says", "enrichment"} <= set(res.inference.columns)
    assert set(res.inference.predicted) <= {"pvm", "cytosol"}
    assert res.inference.control_agrees.any(), "a control confined to one cluster did not agree"
    assert set(res.tables()) == {"recovery", "per_cluster", "inference", "control"}
    assert res.control_summary
    assert len(res.per_cluster) == len(set(res.per_cluster.label)) * len(set(res.per_cluster.cluster))


def test_a_run_without_a_control_still_infers(nodes, monkeypatch):
    monkeypatch.setattr("starplast.clustering.cluster", lambda X, **k: _fake_clusters(len(X)))
    res = R.run(nodes, R.Recipe(question="q", inputs=[RNA, FITNESS], holdout="compartment"),
                tune=False, log=lambda *a: None)
    assert res.ok and res.control.empty and res.control_summary == {}


def test_the_cluster_floor_comes_from_coverage_and_is_stated(nodes, monkeypatch, capsys):
    """The finding this exists for: tuning maximises evenness, which on the real catalogue chose 121
    clusters of ~40 genes -- an even, usable, entirely useless map that named ZERO genes, because no
    cluster held fifteen LABELLED ones at 80% purity. A cluster must be big enough that the floor is
    reachable, and how big depends on what fraction of genes carry a label at all."""
    seen = {}
    monkeypatch.setattr("starplast.clustering.cluster", lambda X, **k: _fake_clusters(len(X)))
    monkeypatch.setattr(R, "best_clustering",
                        lambda coords, grid=None, log=None: seen.update(grid or {}) or {})
    frame = nodes.copy()
    frame["compartment"] = np.where(np.arange(len(frame)) % 4 == 0, "pvm", None)   # 25% labelled
    R.run(frame, R.Recipe(question="q", inputs=[RNA], holdout="compartment"),
          tune=False, log=lambda *a: None)
    assert min(seen["min_cluster_size"]) >= R.MIN_LABEL / 0.25
    assert all(v >= 60 for v in seen["min_cluster_size"])


def test_tuning_chooses_the_settings_and_records_them(nodes, monkeypatch):
    """A run whose settings are not recorded cannot be rebuilt, and an unrebuildable hit is not a
    result."""
    monkeypatch.setattr(R, "tune_umap", lambda *a, **k: pd.DataFrame(
        [{"n_neighbors": 7, "min_dist": 0.3, "usable": True, "score": 1.0, "stage": "full"},
         {"n_neighbors": 9, "min_dist": 0.9, "usable": True, "score": 0.1, "stage": "full"}]))
    monkeypatch.setattr("starplast.clustering.cluster",
                        lambda X, **k: np.repeat(np.arange(3), len(X) // 3 + 1)[:len(X)])
    res = R.run(nodes, R.Recipe(question="q", inputs=[RNA], holdout="compartment"),
                log=lambda *a: None)
    assert res.settings["n_neighbors"] == 7 and res.settings["min_dist"] == 0.3


# --------------------------------------------------------------------------- keeping them
def test_a_saved_recipe_reloads_as_the_same_question(tmp_path):
    store = R.RecipeStore(str(tmp_path))
    store.add(R.Recipe(question="Where does GRA16 go?", inputs=[RNA], holdout="compartment"))
    back = R.RecipeStore(str(tmp_path)).load_all()
    assert len(back) == 1 and back[0].inputs == (RNA,)


def test_asking_the_same_question_twice_overwrites_rather_than_accumulates(tmp_path):
    store = R.RecipeStore(str(tmp_path))
    store.add(R.Recipe(question="same", inputs=[RNA]))
    store.add(R.Recipe(question="same", inputs=[FITNESS]))
    assert len(store.recipes) == 1 and store.recipes[0].inputs == (FITNESS,)
    assert len(os.listdir(tmp_path)) == 1


def test_a_question_with_punctuation_still_gets_a_filename(tmp_path):
    store = R.RecipeStore(str(tmp_path))
    store.add(R.Recipe(question="Which proteins face the host? (PVM/dense granule)"))
    name = os.listdir(tmp_path)[0]
    assert name.endswith(".json") and "/" not in name and "?" not in name


def test_a_file_that_is_not_a_recipe_is_skipped_rather_than_fatal(tmp_path):
    (tmp_path / "notes.txt").write_text("not a recipe")
    (tmp_path / "broken.json").write_text("{ not json")
    (tmp_path / "wrong.json").write_text(json.dumps(["a list, not a recipe"]))
    assert R.RecipeStore(str(tmp_path)).load_all() == []


def test_loading_twice_does_not_duplicate(tmp_path):
    store = R.RecipeStore(str(tmp_path))
    store.add(R.Recipe(question="q"))
    assert len(store.load_all()) == 1 and len(store.load_all()) == 1


def test_a_store_with_no_root_keeps_recipes_in_memory_only():
    store = R.RecipeStore()
    store.add(R.Recipe(question="q"))
    assert store.get("q") is not None and store.load_all() == store.recipes
    assert store.remove("q") and store.get("q") is None


def test_removing_a_question_removes_its_file(tmp_path):
    store = R.RecipeStore(str(tmp_path))
    store.add(R.Recipe(question="q"))
    assert store.remove("q") and os.listdir(tmp_path) == []
    assert not store.remove("q")


def test_a_store_whose_directory_disappeared_returns_what_it_has(tmp_path):
    store = R.RecipeStore(str(tmp_path / "gone"))
    os.rmdir(store.root)
    assert store.load_all() == []


def test_the_same_recipe_run_twice_gives_the_same_answer(nodes):
    """"Rebuilt from its seed" is an acceptance criterion, so it is asserted by rebuilding rather
    than by observing that a seed is stored. A hit nobody can reproduce is not a result -- and every
    stage here is a place a rebuild can silently diverge: the matrix, the UMAP, the clustering.

    Run without tuning, because tuning subsamples and this asks about the map, not the search.
    """
    recipe = R.Recipe(question="q", inputs=[RNA, FITNESS], holdout="compartment", seed=7)
    first = R.run(nodes, recipe, tune=False, log=lambda *a: None)
    second = R.run(nodes, recipe, tune=False, log=lambda *a: None)
    assert np.array_equal(first.labels, second.labels)
    assert np.array_equal(first.genes, second.genes)
    assert first.settings == second.settings
    assert first.summary == second.summary
    pd.testing.assert_frame_equal(first.inference, second.inference)


def test_a_different_seed_is_a_different_map(nodes):
    """The other half of the same claim: if the seed did not reach the embedding, two seeds would
    agree and the reproducibility above would be vacuous."""
    ask = dict(question="q", inputs=[RNA, FITNESS], holdout="compartment")
    a = R.run(nodes, R.Recipe(**ask, seed=7), tune=False, log=lambda *a: None)
    b = R.run(nodes, R.Recipe(**ask, seed=8), tune=False, log=lambda *a: None)
    assert not np.array_equal(a.labels, b.labels)


# --------------------------------------------------------------------------- the edge guard (47)
def test_a_graph_method_may_not_traverse_the_layer_its_holdout_was_built_from(nodes):
    """The largest leak in this project, and until instruction 47 nothing could see it.

    `compartment` is the most-used holdout here AND a 118,712-edge layer made by joining genes that
    share a compartment. Propagating a label along that layer recovers it perfectly, and every
    column-based guard passes, because none of them can see an edge.
    """
    c = R.close(nodes, R.Recipe(question="q", inputs=[FITNESS], holdout="compartment"))
    assert c.ok
    assert "compartment" in c.excluded_layers
    assert "compartment" not in c.layers
    assert "Tg_shared_compartment" in c.excluded_layers["compartment"]


def test_the_edge_guard_holds_at_class_scope_which_is_what_recipes_use(nodes):
    """Measured, and it inverted the assumption this was built on. Asked at CLASS scope -- what every
    recipe uses -- the slot selection bans nothing for `compartment`, because `Tg_shared_compartment`
    sits under `molecular relationships` while the label sits under `cell organization`. The class of
    the target does not contain the layer built from the target. The family is therefore closed
    unconditionally, whatever scope is asked for.
    """
    from starplast import slots
    from starplast.search import excluded_edges, scoped_slots
    class_layers = {layer for slot in scoped_slots("compartment", "biology")
                    for layer in slots.edge_types(slot)}
    assert "compartment" not in class_layers, \
        "the biology class of compartment now contains its own edge layer; this test is stale"
    for scope in ("direct", "target_family", "biology", "evidence", "context"):
        assert "compartment" in excluded_edges("compartment", scope), scope


def test_a_holdout_no_edge_layer_is_built_from_bans_no_layer(nodes):
    """The guard must not simply ban everything: a graph method is useless with no graph."""
    c = R.close(nodes, R.Recipe(question="q", inputs=[FITNESS], holdout="cellcycle_phase"))
    assert c.ok and not c.excluded_layers
    assert len(c.layers) == len(R.all_edge_layers("Tg"))


def test_the_controls_layers_are_banned_as_well_as_the_primarys(nodes):
    """A control is a holdout too, and a map that traversed the control's own layer would corroborate
    itself."""
    # Fitness inputs, not RNA: `cellcycle_phase` IS RNA-derived, so an RNA-fed recipe is refused
    # before it reaches the edge guard -- as this file's own earlier test asserts.
    c = R.close(nodes, R.Recipe(question="q", inputs=[FITNESS], holdout="cellcycle_phase",
                                validation_holdout="compartment"))
    assert c.ok, c.refusal
    assert c.excluded_layers["compartment"].startswith("control:")


def test_the_closure_report_says_which_layers_may_be_traversed(nodes):
    report = R.close(nodes, R.Recipe(question="q", inputs=[FITNESS],
                                     holdout="compartment")).report()
    assert "edge layer" in report and "compartment" in report


def test_alternative_maps_are_built_and_are_actually_different(nodes, monkeypatch):
    """The report shows what the winner was chosen against, so the alternatives have to be real
    builds rather than the same map listed twice. The walk scores every configuration once on a
    sample and again in full, so the same settings arrive twice and would otherwise become "the
    alternative" to themselves."""
    monkeypatch.setattr(R, "tune_umap", lambda *a, **k: pd.DataFrame([
        {"n_neighbors": 7, "min_dist": 0.0, "usable": True, "score": 1.0, "stage": "full"},
        {"n_neighbors": 7, "min_dist": 0.0, "usable": True, "score": 0.9, "stage": "sample"},
        {"n_neighbors": 40, "min_dist": 0.3, "usable": True, "score": 0.5, "stage": "full"},
    ]))
    monkeypatch.setattr("starplast.clustering.cluster", lambda X, **k: _fake_clusters(len(X)))
    res = R.run(nodes, R.Recipe(question="q", inputs=[RNA], holdout="compartment"),
                alternatives=1, log=lambda *a: None)
    assert len(res.alternatives) == 1
    assert res.settings["n_neighbors"] == 7
    assert res.alternatives[0].settings["n_neighbors"] == 40, "the duplicate was taken as the alternative"


def test_no_alternatives_are_built_unless_asked_for(nodes, monkeypatch):
    """They cost a full embedding each."""
    monkeypatch.setattr("starplast.clustering.cluster", lambda X, **k: _fake_clusters(len(X)))
    res = R.run(nodes, R.Recipe(question="q", inputs=[RNA], holdout="compartment"),
                tune=False, log=lambda *a: None)
    assert res.alternatives == []


# --------------------------------------------------------------------------- method dispatch (47)
def test_propagation_without_a_named_layer_is_refused_with_the_choices(nodes):
    """There is no sensible merged default: the layers mean different things and two of them are
    attention-biased, so a walk over the merged graph carries the literature's popularity contest
    into every answer."""
    res = R.run(nodes, R.Recipe(question="q", inputs=[FITNESS], holdout="compartment",
                                method="propagation"), tune=False, log=lambda *a: None)
    assert not res.ok and "name the layer" in res.stopped_because
    assert "xlms" in res.stopped_because


def test_propagation_may_not_walk_the_layer_its_holdout_was_built_from(nodes):
    """The edge guard enforced where it matters rather than merely reported. Diffusing a label across
    the layer built out of that label recovers it perfectly and means nothing."""
    res = R.run(nodes, R.Recipe(question="q", inputs=[FITNESS], holdout="compartment",
                                method="propagation:compartment"), tune=False, log=lambda *a: None)
    assert not res.ok and "excluded for this holdout" in res.stopped_because


def test_an_unknown_layer_is_refused_and_the_real_ones_are_listed(nodes):
    res = R.run(nodes, R.Recipe(question="q", inputs=[FITNESS], holdout="compartment",
                                method="propagation:invented"), tune=False, log=lambda *a: None)
    assert not res.ok and "no layer named" in res.stopped_because


def test_a_classifier_answers_the_same_recipe_and_says_what_it_used(nodes):
    """Every method produces a partition, so everything downstream is one code path. The classifier
    additionally reports which measurements carried the answer, which a clustering cannot."""
    res = R.run(nodes, R.Recipe(question="q", inputs=[FITNESS], holdout="compartment",
                                method="logistic"), tune=False, log=lambda *a: None)
    assert res.ok, res.stopped_because
    assert res.summary.get("n_labels_scored", 0) > 0
    assert len(res.coefficients), "the linear method reported no coefficients"
    assert set(res.coefficients.columns) == {"label", "feature", "coefficient"}
    if len(res.inference):
        assert "predicted_class" in res.inference.columns


def test_a_holdout_with_too_few_labels_to_train_on_says_so(nodes):
    frame = nodes.copy()
    frame["sparse_label"] = None
    frame.loc[frame.index[:8], "sparse_label"] = "a"
    res = R.run(frame, R.Recipe(question="q", inputs=[FITNESS], holdout="sparse_label",
                                method="logistic"), tune=False, log=lambda *a: None)
    assert not res.ok and "too few labelled genes" in res.stopped_because


def test_a_walk_over_a_subset_is_refused_because_edge_ids_are_positions(nodes):
    """Found by a test that used the 400-row fixture: the graph's endpoints are positions in the FULL
    table, so a subset walks between the wrong genes -- and stays in range once the subset is large
    enough, which is the silent version of the same bug."""
    res = R.run(nodes, R.Recipe(question="q", inputs=[FITNESS], holdout="compartment",
                                method="propagation:xlms"), tune=False, log=lambda *a: None)
    assert not res.ok and "full node table" in res.stopped_because


def test_propagation_answers_a_recipe_and_the_control_judges_it_the_same_way(full_nodes):
    """The unification, asserted end to end: a walk over an edge layer produces a partition, and from
    there it is the clustering path's own code -- the same scoring, the same naming, and the same
    independent control on the same groups."""
    res = R.run(full_nodes, R.Recipe(question="q", inputs=[FITNESS], holdout="compartment",
                                     validation_holdout="cellcycle_phase",
                                     method="propagation:xlms"), tune=False, log=lambda *a: None)
    assert res.ok, res.stopped_because
    assert res.settings["method"] == "propagation:xlms"
    assert res.summary.get("n_labels_scored", 0) > 0
    assert not res.control.empty, "the control was not put on the walk's own groups"
    assert res.coefficients.empty, "a walk has no per-column weights to report"


def test_propagation_answers_a_holdout_too_sparse_for_any_clustering(full_nodes):
    """The constraint propagation exists to remove. A clustering needs MIN_LABEL labelled genes
    INSIDE one cluster before it can say anything; a label carried by 40 genes spread over 8,140
    cannot reach that in any partition, and the whole question is unanswerable by that route.
    Seeded with the same 40, a walk scores every gene in the graph."""
    frame = full_nodes.copy()
    rng = np.random.default_rng(0)
    picked = rng.choice(len(frame), 40, replace=False)
    sparse = np.array([None] * len(frame), dtype=object)
    sparse[picked[:20]] = "near"
    sparse[picked[20:]] = "far"
    frame["sparse_label"] = sparse
    res = R.run(frame, R.Recipe(question="q", inputs=[FITNESS], holdout="sparse_label",
                                method="propagation:xlms"), tune=False, log=lambda *a: None)
    assert res.ok, res.stopped_because
    assert (res.labels >= 0).sum() == len(res.labels), "a walk scores every gene, labelled or not"


# --------------------------------------------------------------------------- relevance (43)
def _named(n=12, enrichment=5.0, agrees=True):
    return pd.DataFrame({"gene_id": [f"TGME49_{200000+i}" for i in range(n)],
                         "predicted": ["pvm"] * n, "cluster": [1] * n,
                         "cluster_precision": [0.6] * n, "enrichment": [enrichment] * n,
                         "base_rate": [0.12] * n, "n_labelled_in_cluster": [30] * n,
                         "control_agrees": [agrees] * n})


def test_relevance_keeps_its_terms_as_columns_rather_than_folding_them_away(nodes):
    """Instruction 43 required the score never override the raw numbers it summarises. A reader who
    disagrees with the weighting re-sorts on the term they care about."""
    out = R.relevance(_named(), nodes)
    assert {"strength", "reach", "novelty", "corroborated", "relevance"} <= set(out.columns)
    assert (out.enrichment == 5.0).all(), "the raw numbers were replaced by the score"


def test_a_corroborated_claim_outranks_an_identical_uncorroborated_one(nodes):
    """The term no other ranking in this project has: the difference between "the map found
    structure" and "the map found the structure I asked about"."""
    backed = R.relevance(_named(agrees=True), nodes).relevance.iloc[0]
    alone = R.relevance(_named(agrees=False), nodes).relevance.iloc[0]
    assert backed > alone == pytest.approx(backed / 2, rel=1e-6)


def test_an_uncorroborated_claim_is_weakened_rather_than_voided(nodes):
    """Several shipped questions have controls too sparse to corroborate anything at all; scoring
    those to zero would rank a question's answers by whether its control happened to be dense."""
    assert R.relevance(_named(agrees=False), nodes).relevance.iloc[0] > 0


def test_strength_saturates_so_a_rare_label_does_not_win_by_being_rare(nodes):
    modest = R.relevance(_named(enrichment=10.0), nodes).strength.iloc[0]
    extreme = R.relevance(_named(enrichment=300.0), nodes).strength.iloc[0]
    assert modest == pytest.approx(1.0, abs=0.01) and extreme == pytest.approx(1.0, abs=0.01)


def test_reach_saturates_so_a_cluster_that_swallowed_a_compartment_does_not_win(nodes):
    small = R.relevance(_named(n=40), nodes).reach.iloc[0]
    huge = R.relevance(_named(n=400), nodes).reach.iloc[0]
    assert small == pytest.approx(1.0, abs=0.02) and huge == pytest.approx(1.0, abs=0.02)


def test_nothing_named_ranks_to_an_empty_table_rather_than_an_error(nodes):
    assert R.relevance(pd.DataFrame(), nodes).empty
    assert R.relevance(None, nodes).empty


def test_relevance_is_sorted_so_the_first_row_is_the_one_to_read(nodes):
    mixed = pd.concat([_named(n=4, enrichment=2.0, agrees=False),
                       _named(n=4, enrichment=9.0, agrees=True)], ignore_index=True)
    out = R.relevance(mixed, nodes)
    assert out.relevance.is_monotonic_decreasing and out.enrichment.iloc[0] == 9.0
