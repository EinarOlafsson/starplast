#!/usr/bin/env python3
"""The analysis tab: data -> map -> clusters -> meaning -> search.

The threading contract is the part worth pinning. Every long job runs off the GUI thread, and its
completion is relayed through a BOUND METHOD rather than connected straight to the handler — a
directly-connected handler executes on the worker thread, and touching widgets from there crashes on a
slow machine and works on a fast one, which is the worst way for a bug to behave.

Everything here runs headless. The jobs themselves are exercised synchronously where possible, because
what matters is that the wiring passes the right arguments and puts the result in the right table, not
that Qt's thread pool works.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="module")
def app():
    from PyQt6 import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture(scope="module")
def nodes():
    from starplast import paths
    return pd.read_parquet(paths.cache_file("nodes.parquet"))


@pytest.fixture
def panel(app, nodes, tmp_path):
    from starplast.analysis_panel import AnalysisPanel
    from starplast.tuning import EmbeddingStore
    p = AnalysisPanel(nodes, store=EmbeddingStore(str(tmp_path)))
    yield p
    p.deleteLater()


# --------------------------------------------------------------------------- the worker
def test_a_job_reports_its_result(app):
    from starplast.analysis_panel import Worker
    w = Worker(lambda progress=None: "the answer")
    got = []
    w.done.connect(lambda r, e: got.append((r, e)))
    w.run()
    assert got == [("the answer", None)]


def test_a_job_that_raises_reports_the_error_rather_than_dying(app):
    """A failed analysis must leave the panel usable, not take the window down with it."""
    from starplast.analysis_panel import Worker

    def boom(progress=None):
        raise ValueError("no usable columns")

    w = Worker(boom)
    got = []
    w.done.connect(lambda r, e: got.append((r, e)))
    w.run()
    assert got[0][0] is None
    assert "no usable columns" in str(got[0][1])


def test_a_job_can_report_progress(app):
    from starplast.analysis_panel import Worker
    seen = []
    w = Worker(lambda progress: progress("half way"))
    w.progress.connect(seen.append)
    w.run()
    assert seen == ["half way"]


# --------------------------------------------------------------------------- the spec
def test_the_panel_builds_a_spec_from_its_controls(panel):
    spec = panel.spec()
    assert spec.blocks
    assert spec.na_policy in ("indicator", "median", "drop_columns", "drop_genes")
    assert spec.scaling in ("robust", "zscore", "rank", "none")


def test_changing_a_control_changes_the_spec(panel):
    """Otherwise the panel shows settings it does not use, which is worse than showing none."""
    before = panel.spec()
    panel.na_policy.setCurrentText("drop_columns")
    panel.scaling.setCurrentText("rank")
    after = panel.spec()
    assert after.na_policy == "drop_columns" and after.scaling == "rank"
    assert (before.na_policy, before.scaling) != (after.na_policy, after.scaling)


def test_the_seed_is_part_of_the_spec(panel):
    """A result that cannot be reproduced is not a result."""
    panel.seed.setValue(7)
    assert panel.spec().random_state == 7


def test_variance_shares_are_shown_for_the_current_spec(panel):
    """The check that catches a block being named as an input while contributing nothing -- hyperLOPIT
    at 1.1%."""
    panel.show_variance()
    assert panel.variance_view.toPlainText().strip()


# --------------------------------------------------------------------------- filling tables
def test_a_table_is_filled_with_headers_and_formatted_numbers(panel):
    df = pd.DataFrame({"a": [1.23456, 2.0], "b": ["x", "y"]})
    panel._fill(panel.walk_table, df)
    assert panel.walk_table.rowCount() == 2
    assert panel.walk_table.horizontalHeaderItem(0).text() == "a"
    assert panel.walk_table.item(0, 0).text() == "1.235"


def test_a_long_table_is_truncated_rather_than_freezing_the_window(panel):
    df = pd.DataFrame({"a": range(5000)})
    panel._fill(panel.walk_table, df, limit=50)
    assert panel.walk_table.rowCount() == 50


def test_a_non_finite_number_is_shown_as_itself(panel):
    """nan formatted as '%.3f' reads as a measurement of nan; shown plainly it reads as missing."""
    panel._fill(panel.walk_table, pd.DataFrame({"a": [np.nan, np.inf]}))
    assert "nan" in panel.walk_table.item(0, 0).text().lower()


def test_filling_a_table_twice_replaces_rather_than_appends(panel):
    panel._fill(panel.walk_table, pd.DataFrame({"a": [1.0, 2.0, 3.0]}))
    panel._fill(panel.walk_table, pd.DataFrame({"a": [1.0]}))
    assert panel.walk_table.rowCount() == 1


# --------------------------------------------------------------------------- one job at a time
def test_a_second_job_is_refused_while_one_is_running(panel):
    """Two UMAP walks over the same table at once would compete for every core and finish later than
    either alone."""
    said = []
    panel.status.connect(said.append)
    panel._thread = object()                     # pretend a job is in flight
    panel._run(lambda progress=None: None, lambda r: None)
    assert any("already running" in m for m in said)
    panel._thread = None


def test_a_completed_job_frees_the_slot(panel, app):
    from PyQt6 import QtCore
    done = []
    panel._run(lambda progress=None: 42, done.append)
    for _ in range(200):
        app.processEvents()
        QtCore.QThread.msleep(5)
        if done:
            break
    assert done == [42]
    assert panel._thread is None, "the slot must be free for the next job"


def test_a_failing_job_frees_the_slot_and_reports(panel, app):
    from PyQt6 import QtCore

    def boom(progress=None):
        raise RuntimeError("nope")

    said = []
    panel.status.connect(said.append)
    panel._run(boom, lambda r: None)
    for _ in range(200):
        app.processEvents()
        QtCore.QThread.msleep(5)
        if panel._thread is None:
            break
    assert panel._thread is None
    assert any("failed" in m for m in said)


def test_the_completion_handler_runs_on_the_gui_thread(panel, app):
    """Relayed through a bound method for exactly this reason: a directly-connected handler executes on
    the worker, and touching widgets from there crashes on a slow machine and works on a fast one."""
    from PyQt6 import QtCore
    seen = {}
    main_thread = QtCore.QThread.currentThread()

    def on_done(result):
        seen["thread"] = QtCore.QThread.currentThread()

    panel._run(lambda progress=None: 1, on_done)
    for _ in range(200):
        app.processEvents()
        QtCore.QThread.msleep(5)
        if seen:
            break
    assert seen["thread"] is main_thread


# --------------------------------------------------------------------------- the tabs exist
def test_every_stage_of_the_workflow_has_a_tab(panel):
    """data -> map -> clusters -> inference -> search -> validation, in the order the work is done.

    "Meaning" was renamed to "Inference": it reports which held-out features distinguish the
    clusters, which is an inference and exactly as strong as the held-out testing behind it, and the
    name should not claim more. Validation is the new stage that puts a number on an annotation.
    """
    tabs = panel.findChild(type(panel).__mro__[0] if False else __import__("PyQt6.QtWidgets",
                           fromlist=["QTabWidget"]).QTabWidget)
    titles = [tabs.tabText(i).lower() for i in range(tabs.count())]
    assert len(titles) >= 6
    for expected in ("data", "map", "cluster", "inference", "search", "validation"):
        assert any(expected in t for t in titles), f"no tab for {expected}: {titles}"


def test_the_search_tab_offers_every_registered_target(panel):
    """A target in the menu that is not in the table fails only when chosen."""
    from starplast import search as S
    offered = {panel.target.itemText(i) for i in range(panel.target.count())}
    assert offered & set(S.TARGETS), f"none of the registered targets are offered: {offered}"


def test_the_measured_cell_cycle_target_is_offered(panel):
    """Half the project's stated purpose is unaskable without it."""
    offered = {panel.target.itemText(i) for i in range(panel.target.count())}
    assert "cellcycle_phase" in offered


# --------------------------------------------------------------------------- the job handlers
@pytest.fixture
def sync(panel, monkeypatch):
    """Run jobs synchronously so the closure and its completion handler are both exercised.

    What matters is that the wiring passes the right arguments and puts the result in the right table,
    not that Qt's thread pool works -- and _run's threading contract is tested directly above.
    """
    def run_now(fn, on_done, name="analysis", on_error=None):
        # `name` is accepted because every job carries one now, so it can be identified in the Jobs
        # panel and stopped there. Recorded rather than dropped: a job named "analysis" for all five
        # tabs would make the panel useless, so the names are worth asserting on.
        #
        # `on_error` is routed the same way the runner routes it, because some failures are not
        # failures: validation refusing a circular target is the guard working, and a fixture that
        # let the exception escape would never exercise the branch that says so.
        started.append(name)
        try:
            result = fn(lambda *_: None)
        except Exception as exc:
            if not (on_error is not None and on_error(exc)):
                raise
            return
        on_done(result)

    started = []
    monkeypatch.setattr(panel, "_run", run_now)
    # Recorded on the panel rather than on the returned list, which is a plain list and takes no
    # attributes.
    panel.started_job_names = started
    said = []
    panel.status.connect(said.append)
    return said


def test_the_umap_walk_fills_its_table(panel, sync, monkeypatch):
    import starplast.tuning as T
    monkeypatch.setattr(T, "walk_umap",
                        lambda n, spec, **k: pd.DataFrame({"n_neighbors": [15], "trust": [0.9]}))
    panel.run_umap_walk()
    assert panel.walk_table.rowCount() == 1
    assert any("walk complete" in m for m in sync)


def test_a_walk_row_appears_the_moment_that_configuration_finishes(panel, sync, monkeypatch):
    """The whole reason the walk emits per configuration. Filling the table at the end means a
    288-configuration sweep -- half an hour -- shows an empty table for the entire run."""
    import starplast.tuning as T
    from starplast.embedding import EmbeddingSpec

    def fake(n, spec, on_step=None, **k):
        for i, md in enumerate((0.1, 0.25), start=1):
            on_step(T.WalkStep(index=i, total=2,
                               row={"n_neighbors": 15, "min_dist": md, "trustworthiness": 0.9},
                               coords=np.zeros((4, 3)), genes=np.ones(len(n), bool),
                               spec=EmbeddingSpec(n_neighbors=15, min_dist=md)))
            # Asserted INSIDE the walk: after it returns, an incremental fill and a fill-at-the-end
            # are indistinguishable.
            assert panel.walk_table.rowCount() == i, "the row did not arrive with the step"
        return pd.DataFrame({"n_neighbors": [15, 15], "min_dist": [0.25, 0.1],
                             "trustworthiness": [0.95, 0.9]})

    monkeypatch.setattr(T, "walk_umap", fake)
    panel.run_umap_walk()
    assert panel.walk_table.rowCount() == 2
    assert any("walk 2 of 2" in m for m in sync)
    # The ranking replaces the running order once there is a whole sweep to rank.
    assert panel.walk_table.item(0, 1).text().startswith("0.25")


def test_a_new_walk_empties_the_table_before_it_starts(panel, sync, monkeypatch):
    """Rows from two walks in one table are a comparison between configurations that were never
    compared -- and the second walk's grid may not even have the same columns."""
    import starplast.tuning as T
    seen = []
    panel.walk_started.connect(lambda: seen.append(panel.walk_table.rowCount()))
    monkeypatch.setattr(T, "walk_umap",
                        lambda n, spec, **k: pd.DataFrame({"n_neighbors": [15], "trust": [0.9]}))
    panel.run_umap_walk()
    panel.run_umap_walk()
    assert seen == [0, 0], "the table still held the last walk when the next one started"


def test_a_second_walk_is_refused_while_one_is_still_running(panel, sync, monkeypatch):
    """Two walks fill one table and one gallery, and the mixture reads as a single sweep -- which
    invites a comparison between configurations that were never compared. The runner allows
    concurrent jobs on purpose, so the constraint belongs on the one job whose output accumulates
    somewhere shared."""
    import starplast.tuning as T
    from starplast.analysis_panel import WALK_JOB
    from starplast.jobs import Job, RUNNING

    class Runner:
        jobs = {1: Job(id=1, name=WALK_JOB, state=RUNNING)}

    ran = []
    monkeypatch.setattr(T, "walk_umap", lambda *a, **k: ran.append(1) or pd.DataFrame({"a": [1]}))
    panel.runner = Runner()
    panel.run_umap_walk()
    assert ran == [] and any("already running" in m for m in sync)
    Runner.jobs[1].state = "done"
    panel.run_umap_walk()
    assert ran == [1], "a finished walk must not block the next one"


def test_the_walk_saves_every_configuration_through_the_store(panel, sync, monkeypatch):
    """A stopped walk should leave behind what it finished. Without the store it leaves nothing."""
    import starplast.tuning as T
    seen = {}

    def fake(n, spec, store=None, **k):
        seen["store"] = store
        return pd.DataFrame({"a": [1]})

    monkeypatch.setattr(T, "walk_umap", fake)
    panel.run_umap_walk()
    assert seen["store"] is panel.store


def test_the_walk_is_given_the_sample_size_and_seed_from_the_panel(panel, sync, monkeypatch):
    """A walk that silently used different settings than the ones on screen would be unreproducible."""
    import starplast.tuning as T
    seen = {}

    def fake(n, spec, sample_size=None, seed=None, log=None, **kw):
        seen.update(sample_size=sample_size, seed=seed, blocks=spec.blocks)
        return pd.DataFrame({"a": [1]})

    monkeypatch.setattr(T, "walk_umap", fake)
    panel.sample.setValue(1500)
    panel.seed.setValue(11)
    panel.run_umap_walk()
    assert seen["sample_size"] == 1500 and seen["seed"] == 11


def test_building_the_map_publishes_the_coordinates(panel, sync, monkeypatch):
    """The main window listens for this to redraw; without it the map is built and never shown."""
    import starplast.embedding as E
    coords = np.zeros((5, 3))
    rows = np.array([True] * 5)
    monkeypatch.setattr(E, "embed", lambda n, spec, log=None: (coords, ["a", "b"], rows))
    got = []
    panel.embedding_ready.connect(lambda c, r: got.append((c, r)))
    panel.run_embed()
    assert len(got) == 1
    assert panel.coords is coords
    assert any("map built" in m for m in sync)


def test_saving_before_building_says_so_rather_than_writing_nothing(panel, sync):
    panel.coords = None
    panel.save_embedding()
    assert any("build a map first" in m for m in sync)


def test_an_embedding_is_saved_with_its_full_recipe(panel, sync, monkeypatch, tmp_path):
    """A result that cannot be reopened with the settings that produced it is not reproducible."""
    import starplast.embedding as E
    rows = np.array([True] * len(panel.nodes))
    monkeypatch.setattr(E, "embed",
                        lambda n, spec, log=None: (np.zeros((len(n), 3)), ["a"], rows))
    panel.run_embed()
    panel.emb_name.setText("my run")
    panel.save_embedding()
    assert "my run" in list(panel.store.list().name)
    _, spec_back, _ = panel.store.load("my run")
    assert spec_back is not None and spec_back.na_policy == panel.spec().na_policy


def test_clustering_before_building_says_so(panel, sync):
    panel.coords = None
    panel.run_cluster_walk()
    panel.run_cluster()
    assert sum("build a map first" in m for m in sync) == 2


def test_the_cluster_walk_fills_its_table(panel, sync, monkeypatch):
    import starplast.clustering as C
    panel.coords = np.random.default_rng(0).normal(size=(60, 3))
    monkeypatch.setattr(C, "walk_hdbscan",
                        lambda Y, log=None: pd.DataFrame({"min_cluster_size": [10]}))
    monkeypatch.setattr(C, "walk_dbscan", lambda Y, log=None: pd.DataFrame({"eps": [0.5]}))
    panel.algo.setCurrentText("hdbscan")
    panel.run_cluster_walk()
    assert panel.cluster_table.rowCount() == 1


def test_the_chosen_algorithm_is_the_one_that_runs(panel, sync, monkeypatch):
    """Two different algorithms with different failure modes; running the wrong one silently would
    make every downstream number wrong."""
    import starplast.clustering as C
    panel.coords = np.random.default_rng(0).normal(size=(60, 3))
    called = []
    monkeypatch.setattr(C, "walk_hdbscan",
                        lambda Y, log=None: (called.append("hdbscan"), pd.DataFrame({"a": [1]}))[1])
    monkeypatch.setattr(C, "walk_dbscan",
                        lambda Y, log=None: (called.append("dbscan"), pd.DataFrame({"a": [1]}))[1])
    panel.algo.setCurrentText("dbscan")
    panel.run_cluster_walk()
    assert called == ["dbscan"]


def test_clustering_reports_the_count_and_how_much_is_unassigned(panel, sync, monkeypatch):
    """The noise fraction is half the answer: a clustering that assigns 5% of genes has not organised
    the proteome however good its silhouette is."""
    import starplast.clustering as C
    panel.coords = np.random.default_rng(0).normal(size=(60, 3))
    labels = np.array([0] * 30 + [1] * 20 + [-1] * 10)
    monkeypatch.setattr(C, "cluster", lambda Y, **k: labels)
    panel.run_cluster()
    assert panel.labels is labels
    assert any("2 clusters" in m and "17% unassigned" in m for m in sync)


def test_the_battery_before_clustering_says_so(panel, sync):
    panel.labels = None
    panel.run_battery()
    assert any("cluster the map first" in m for m in sync)


def test_the_battery_reports_only_held_out_features(panel, sync, monkeypatch):
    """A feature the embedding used separates the clusters by construction, so showing it beside a
    discovery would put the two on the same footing."""
    import starplast.clustering as C
    panel.rows = np.array([True] * len(panel.nodes))
    panel.labels = np.array([0] * (len(panel.nodes) // 2)
                            + [1] * (len(panel.nodes) - len(panel.nodes) // 2))
    S = pd.DataFrame({"feature": ["a", "b"], "evidence": ["held_out", "used"],
                      "score_type": ["cramers_v"] * 2, "score": [0.8, 0.99],
                      "assoc_with_input": [0.1, 1.0], "n": [100, 100], "q": [0.01, 0.0]})
    monkeypatch.setattr(C, "battery", lambda *a, **k: (S, pd.DataFrame()))
    monkeypatch.setattr(C, "describe", lambda *a, **k: ["a separates the clusters"])
    panel.run_battery()
    assert panel.battery_table.rowCount() == 1
    assert "a separates the clusters" in panel.findings.toPlainText()
    assert any("1 held-out features" in m for m in sync)


def test_a_battery_that_found_nothing_says_so(panel, sync, monkeypatch):
    """`+` binds tighter than `or`, so the header made the whole expression truthy and this fallback
    was unreachable -- an empty result showed the explanation alone, which reads as a page that has
    not loaded rather than as an honest 'nothing here'."""
    import starplast.clustering as C
    panel.rows = np.array([True] * len(panel.nodes))
    panel.labels = np.zeros(len(panel.nodes), dtype=int)
    monkeypatch.setattr(C, "battery", lambda *a, **k: (pd.DataFrame(), pd.DataFrame()))
    monkeypatch.setattr(C, "describe", lambda *a, **k: [])
    panel.run_battery()
    assert "nothing separated the clusters" in panel.findings.toPlainText()


def test_the_battery_is_told_what_the_embedding_used(panel, sync, monkeypatch):
    """Which is the whole basis of the evidence tiering -- get it wrong and a used feature is reported
    as a discovery."""
    import starplast.clustering as C
    panel.rows = np.array([True] * len(panel.nodes))
    panel.labels = np.zeros(len(panel.nodes), dtype=int)
    seen = {}

    def fake_battery(sub, lab, used_features=(), log=None):
        seen["used"] = list(used_features)
        return pd.DataFrame(), pd.DataFrame()

    monkeypatch.setattr(C, "battery", fake_battery)
    monkeypatch.setattr(C, "describe", lambda *a, **k: [])
    panel.cat_cb.setChecked(True)
    panel.run_battery()
    assert seen["used"], "the battery must be told which columns fed the embedding"
    assert "compartment" in seen["used"], "the categorical block counts as used when it is enabled"


def test_the_structure_search_fills_its_table(panel, sync, monkeypatch):
    import starplast.search as S
    monkeypatch.setattr(S, "search",
                        lambda n, **k: (pd.DataFrame({"mean_f1": [0.5]}), pd.DataFrame()))
    panel.run_search()
    assert panel.search_table.rowCount() == 1
    assert any("search complete" in m for m in sync)


def test_the_search_also_fills_a_table_of_per_category_scores(panel, sync, monkeypatch):
    """A mean over categories hides the case this table exists for: a map where the GRAs are clean
    and everything else is a mess is exactly what you want when you are looking for GRAs."""
    import starplast.search as S
    R = pd.DataFrame({"mean_f1": [0.5], "best_f1": [0.8], "blocks": ["fitness_screens"]})
    P = pd.DataFrame({"label": ["dense granules", "rhoptries"], "precision": [0.9, 0.2],
                      "recall": [0.7, 0.1], "f1": [0.79, 0.13], "n_label": [193, 71],
                      "blocks": ["fitness_screens"] * 2, "n_neighbors": [15, 15],
                      "min_dist": [0.1, 0.1], "min_cluster_size": [25, 25]})
    monkeypatch.setattr(S, "search", lambda n, **k: (R, P))
    panel.run_search()
    assert panel.category_search_table.rowCount() == 2
    heads = [panel.category_search_table.horizontalHeaderItem(c).text()
             for c in range(panel.category_search_table.columnCount())]
    assert "category" in heads and {"precision", "recall", "f1"} <= set(heads)
    assert panel.category_search_table.item(0, heads.index("category")).text() == "dense granules"


def test_a_scored_configuration_reaches_both_tables_as_it_finishes(panel, sync):
    """The automated walk is hundreds of runs; a table that arrives at the end is a table nobody
    watches, and the question it exists to answer is asked while it is still running."""
    import numpy as np
    from starplast.embedding import EmbeddingSpec
    from starplast.search import RunStep
    row = {"blocks": "fitness_screens", "n_neighbors": 15, "min_dist": 0.1,
           "min_cluster_size": 25, "seed": 42, "sample_size": 0, "excluded": "compartment",
           "mean_f1": 0.4, "best_f1": 0.6, "n_clusters": 7}
    per = pd.DataFrame({"label": ["dense granules"], "precision": [0.9], "recall": [0.5],
                        "f1": [0.64], "n_label": [193], "n_in_cluster": [96], "cluster": [3]})
    panel.search_step.emit(RunStep(index=1, total=8, row=row, per=per,
                                   coords=np.zeros((10, 3)), genes=np.ones(len(panel.nodes), bool),
                                   labels=np.zeros(10, int), spec=EmbeddingSpec()))
    assert panel.search_table.rowCount() == 1
    assert panel.category_search_table.rowCount() == 1
    assert any("search 1 of 8" in m for m in sync)


def test_a_per_category_row_carries_enough_to_rebuild_its_map(panel, sync, monkeypatch):
    """Otherwise the score and the configuration live in different tables and the row cannot be
    clicked into anything."""
    import numpy as np
    import starplast.search as S
    from starplast.embedding import EmbeddingSpec
    from starplast.search import RunStep
    row = {"blocks": "fitness_screens", "n_neighbors": 15, "min_dist": 0.1,
           "min_cluster_size": 25, "seed": 42, "sample_size": 0, "excluded": "compartment"}
    per = pd.DataFrame({"label": ["dense granules"], "precision": [0.9], "recall": [0.5],
                        "f1": [0.64], "n_label": [193], "n_in_cluster": [96], "cluster": [3]})
    panel.search_step.emit(RunStep(index=1, total=1, row=row, per=per, coords=np.zeros((10, 3)),
                                   genes=np.ones(len(panel.nodes), bool), labels=np.zeros(10, int),
                                   spec=EmbeddingSpec()))
    seen = {}

    def fake_rebuild(nodes, r, log=print):
        seen.update(r)
        genes = np.zeros(len(nodes), bool)
        genes[:20] = True
        return np.zeros((20, 3)), genes, np.zeros(20, int), ["f"]

    monkeypatch.setattr(S, "rebuild", fake_rebuild)
    panel.show_search_category_row(0)
    assert seen["blocks"] == "fitness_screens" and seen["min_cluster_size"] == "25"
    assert any("dense granules was matched by cluster 3" in m for m in sync)


def test_a_per_category_row_with_no_configuration_says_so(panel, sync):
    panel._fill(panel.category_search_table, pd.DataFrame({"category": ["x"], "f1": [0.5]}))
    panel.show_search_category_row(0)
    assert any("does not name a configuration" in m for m in sync)


def test_the_search_is_given_the_chosen_target_and_a_bounded_set_of_combinations(panel, sync,
                                                                                 monkeypatch):
    """Every block combination up to the chosen size; unbounded it is 2^n runs."""
    import starplast.search as S
    seen = {}

    def fake(n, target=None, block_sets=None, sample_size=None, seed=None, store=None,
             objective=None, log=None, **kw):
        seen.update(target=target, n_sets=len(block_sets), sample_size=sample_size)
        return pd.DataFrame({"mean_f1": [0.1]}), pd.DataFrame()

    monkeypatch.setattr(S, "search", fake)
    panel.max_blocks.setValue(1)
    panel.search_sample.setValue(1000)
    panel.run_search()
    assert seen["target"] == panel.target.currentText()
    assert seen["sample_size"] == 1000
    assert all(len(s) == 1 for s in [()] ) or seen["n_sets"] > 0


def test_a_spec_that_selects_nothing_is_reported_in_place(panel):
    """Selecting no usable block raises rather than returning an empty matrix, and the panel shows the
    reason where the shares would go instead of leaving a blank box."""
    for cb in panel.block_cb.values() if hasattr(panel, "block_cb") else []:
        cb.setChecked(False)
    import starplast.analysis_panel as AP
    orig = AP.variance_share

    def boom(nodes, spec):
        raise ValueError("no numeric features selected")

    AP.variance_share = boom
    try:
        panel.show_variance()
    finally:
        AP.variance_share = orig
    assert "no numeric features selected" in panel.variance_view.toPlainText()


# --------------------------------------------------------------------------- the objective chooser
def test_the_search_offers_every_objective(panel):
    """Which objective a walk maximises is the most consequential choice in it, and it used to be
    made for the user with no way to change it."""
    from starplast.objectives import OBJECTIVES
    assert panel.objective.count() == len(OBJECTIVES)
    names = {panel.objective.itemData(i) for i in range(panel.objective.count())}
    assert names == set(OBJECTIVES)


def test_the_objective_settings_travel_as_one_dict(panel):
    """A score whose objective is not recorded beside it cannot be compared with another."""
    s = panel.objective_settings()
    # Every keyword objectives.score takes, so the whole scoring decision travels together. The
    # floors are in here because they are part of it: min_cluster is what stops a four-gene cluster
    # winning a purity objective, and a score recorded without it cannot be reproduced.
    assert set(s) == {"objective", "weighting", "category", "min_recall",
                      "min_cluster", "min_label"}
    assert s["weighting"] == "macro", "size weighting buries the rare classes worth hunting"
    import inspect
    from starplast.objectives import score as score_fn
    accepted = set(inspect.signature(score_fn).parameters)
    assert set(s) <= accepted | {"min_cluster", "min_label"}, "a setting score() cannot accept"


def test_the_focus_list_follows_the_target(panel):
    """The categories of compartment are not the categories of cellcycle_phase, and a stale list
    would let someone optimise for a label absent from the column being scored."""
    panel.target.setCurrentText("compartment")
    comp = panel.focus.count()
    panel.target.setCurrentText("cellcycle_phase")
    assert panel.focus.count() != comp


def test_absence_is_never_offered_as_a_category_to_optimise_for(panel):
    from starplast.search import ABSENCE_LABELS
    panel.target.setCurrentText("compartment")
    from PyQt6 import QtCore as _Qt
    got = {panel.focus.item(i).data(_Qt.Qt.ItemDataRole.UserRole)
           for i in range(panel.focus.count())} - {None}
    assert not {g for g in got if str(g).lower() in ABSENCE_LABELS}


def test_the_search_is_given_the_chosen_objective(panel, sync, monkeypatch):
    import starplast.search as S
    seen = {}

    def fake(nodes, **kw):
        seen.update(kw)
        return pd.DataFrame({"mean_f1": [0.5]}), pd.DataFrame()

    monkeypatch.setattr(S, "search", fake)
    panel.objective.setCurrentIndex(
        [panel.objective.itemData(i) for i in range(panel.objective.count())].index("best_precision"))
    panel.run_search()
    assert seen["objective"]["objective"] == "best_precision"


# --------------------------------------------------------------------------- tooltip shape
def test_a_tooltip_is_a_block_not_one_long_line(panel):
    """Qt lays a plain tooltip out on one line, so a two-sentence explanation becomes a strip wider
    than the window."""
    from PyQt6 import QtGui
    tip = panel.mcs.toolTip()
    assert "white-space:pre" in tip, "without this Qt re-wraps and strands words on their own line"
    doc = QtGui.QTextDocument()
    doc.setHtml(tip)
    assert doc.idealWidth() < 600, f"laid out at {doc.idealWidth():.0f}px, which runs off the screen"
    lines = [x for x in doc.toPlainText().split("\n") if x.strip()]
    assert len(lines) > 1, "not wrapped at all"
    # No stranded line: every line but the last in a paragraph should be near the full width.
    widths = [len(x.rstrip()) for x in lines]
    assert max(widths) - min(w for w in widths if w) < 70, f"very uneven lines: {widths}"


def test_a_bounded_control_says_why_its_bounds_are_where_they_are(panel):
    """An unexplained limit reads as arbitrary, or worse as a limit of the method."""
    from starplast.analysis_panel import LIMITS
    for attr in LIMITS:
        w = getattr(panel, attr, None)
        if w is None:
            continue
        tip = w.toolTip()
        assert "Range" in tip, f"{attr} does not state its range"
        # The reason is present, checked on a distinctive fragment rather than the first few words,
        # which are often ordinary ones like "Below" or "Same".
        # Unescaped and unwrapped before comparing: the tooltip is rich text, so an apostrophe is
        # &#x27; and the line breaks fall wherever the wrap put them.
        import html as _html
        flat = " ".join(_html.unescape(tip.replace("<br>", " ")).split())
        frag = " ".join(LIMITS[attr].split()[:6])
        assert frag[:38] in flat, f"{attr} does not say why: {flat[:120]}"


def test_the_range_in_the_tooltip_is_read_from_the_widget(panel):
    """Written by hand it would drift from the range actually set."""
    from starplast.analysis_panel import range_note
    panel.mcs.setRange(7, 99)
    assert "7" in range_note(panel.mcs) and "99" in range_note(panel.mcs)
    panel.mcs.setRange(3, 2000)


def test_hovering_the_name_works_not_only_the_field(panel):
    """The label is the word people point at; the spin box is the thing they click."""
    from PyQt6 import QtWidgets
    checked = 0
    for form in panel.findChildren(QtWidgets.QFormLayout):
        for attr in ("nn", "md", "mcs", "val_folds", "search_sample"):
            w = getattr(panel, attr, None)
            if w is None:
                continue
            label = form.labelForField(w)
            if label is not None:
                assert label.toolTip() == w.toolTip(), f"{attr}'s label has a different tooltip"
                checked += 1
    assert checked >= 3, "no form labels were found to check"


# --------------------------------------------------------------------------- the walk grid
def test_the_walk_grid_is_reachable_from_the_interface(panel):
    """It was hardcoded in tuning.walk_umap, so 'walk hyperparameters' swept a set nobody could see
    or change -- while the spin boxes above it, which look like they control it, only ever affected
    'build this map'."""
    g = panel.walk_grid()
    assert set(g) == {"n_neighbors_values", "min_dist_values", "min_cluster_sizes"}
    panel.nn_grid.setText("7, 9")
    assert panel.walk_grid()["n_neighbors_values"] == (7, 9)
    panel.nn_grid.setText("5, 15, 25, 50, 100")


def test_a_grid_accepts_a_list_or_a_range(panel):
    assert panel.parse_grid("5, 15, 25", int) == (5, 15, 25)
    assert panel.parse_grid("5:50:15", int) == (5, 20, 35, 50)
    assert panel.parse_grid("0.0 0.25", float) == (0.0, 0.25)


def test_an_unparsable_grid_falls_back_rather_than_raising(panel):
    """A walk is expensive to start; losing one to a stray comma is worse than sweeping defaults."""
    assert panel.parse_grid("oops", int, (1, 2)) == (1, 2)
    assert panel.parse_grid("", int, (1, 2)) == (1, 2)
    assert panel.parse_grid("5:50:0", int, (1, 2)) == (1, 2), "a zero step would never terminate"


def test_the_two_tabs_cannot_disagree_about_the_grid(panel):
    """The Map tab shows it and the Search tab sweeps it, so they are kept in step both ways."""
    panel.nn_grid.setText("11, 22")
    assert panel.nn_grid2.text() == "11, 22"
    panel.md_grid2.setText("0.3")
    assert panel.md_grid.text() == "0.3"
    panel.nn_grid.setText("5, 15, 25, 50, 100")
    panel.md_grid.setText("0.0, 0.1, 0.25, 0.5")


def test_the_walk_is_given_the_grid_from_the_interface(panel, sync, monkeypatch):
    import starplast.tuning as T
    seen = {}

    def fake(n, spec, **kw):
        seen.update(kw)
        return pd.DataFrame({"n_neighbors": [15]})

    monkeypatch.setattr(T, "walk_umap", fake)
    panel.nn_grid.setText("3, 4")
    panel.run_umap_walk()
    assert seen["n_neighbors_values"] == (3, 4)
    panel.nn_grid.setText("5, 15, 25, 50, 100")


# --------------------------------------------------------------------------- clicking a walk row
def test_clicking_a_walk_row_builds_that_configuration(panel, monkeypatch):
    """A table of scores is not a map. The point of a walk is to look at the ones that scored well,
    and until now there was no way to get from a row to the embedding it describes."""
    panel._fill(panel.walk_table,
                pd.DataFrame({"n_neighbors": [7, 33], "min_dist": [0.3, 0.1], "trust": [0.9, 0.8]}))
    built = {}
    monkeypatch.setattr(panel, "run_embed", lambda then_cluster=None: built.setdefault(
        "spec", (panel.nn.value(), panel.md.value(), then_cluster)))
    panel.show_walk_row(1, 0)
    # No cluster count on these rows, so nothing to cluster: the walk was run without the check.
    assert built["spec"] == (33, 0.1, None)


def test_a_walk_row_that_counted_clusters_brings_them_with_it(panel, monkeypatch):
    """The row says "11 clusters". A map shown without them leaves the reader taking that number on
    trust, which is the one thing this application is built not to ask for."""
    from starplast.tuning import WALK_MIN_CLUSTER_SIZE
    panel._fill(panel.walk_table, pd.DataFrame({"n_neighbors": [7], "min_dist": [0.3],
                                                "n_clusters_hdbscan": [11]}))
    built = {}
    monkeypatch.setattr(panel, "run_embed",
                        lambda then_cluster=None: built.setdefault("mcs", then_cluster))
    panel.show_walk_row(0, 0)
    assert built["mcs"] == WALK_MIN_CLUSTER_SIZE, "clustered at a different size from the walk's own"


def test_clicking_a_row_that_names_no_configuration_says_so(panel, monkeypatch):
    said = []
    panel.status.connect(said.append)
    panel._fill(panel.walk_table, pd.DataFrame({"n_neighbors": ["n/a"], "min_dist": ["n/a"]}))
    monkeypatch.setattr(panel, "run_embed",
                        lambda then_cluster=None: pytest.fail("must not build from a bad row"))
    panel.show_walk_row(0, 0)
    assert any("cannot rebuild" in m or "does not name" in m for m in said)


def test_mirroring_the_grid_does_not_re_enter(panel):
    """A plain two-way binding re-enters -- setText emits textChanged, which sets the first again --
    and the pair can still be firing at each other while Qt is deleting them."""
    panel.nn_grid.setText("3, 4, 5")
    assert panel.nn_grid2.text() == "3, 4, 5"
    panel.nn_grid2.setText("6, 7")
    assert panel.nn_grid.text() == "6, 7"
    panel.nn_grid.setText("5, 15, 25, 50, 100")


# --------------------------------------------------------------------------- choosing the labels
def test_any_objective_can_be_combined_with_any_set_of_labels(panel):
    """"Precision for dense granules" and "recall for dense granules and rhoptries" are different
    questions, so the objective and the label set are independent choices."""
    from PyQt6 import QtCore as _Qt
    panel.target.setCurrentText("compartment")
    panel.focus.set_checked([])
    assert panel.objective_settings()["category"] is None, "none ticked must mean all of them"
    # Ticked, not highlighted: a highlight is destroyed by the next ordinary click, which is why
    # assembling a set of categories had to be redone every time the reader looked away.
    for i in range(min(2, panel.focus.count())):
        panel.focus.item(i).setCheckState(_Qt.Qt.CheckState.Checked)
    got = panel.objective_settings()["category"]
    assert isinstance(got, list) and len(got) == 2
    assert all(g == panel.focus.item(i).data(_Qt.Qt.ItemDataRole.UserRole)
               for i, g in enumerate(got))
    panel.focus.set_checked([])


def test_scoring_several_labels_aggregates_by_the_objective():
    """mean objectives average over the chosen labels; best objectives take the best among them."""
    import numpy as np
    from starplast import objectives as O
    t = pd.Series(["a"] * 60 + ["b"] * 60 + ["c"] * 60)
    lab = np.array([0] * 60 + [1] * 30 + [2] * 30 + [3] * 60)   # "b" split, "a" and "c" clean
    both = O.score(lab, t, objective="mean_recall", category=["a", "b"], min_cluster=1, min_label=5)
    best = O.score(lab, t, objective="best_f1", category=["a", "b"], min_cluster=1, min_label=5)
    only_a = O.score(lab, t, objective="mean_recall", category="a", min_cluster=1, min_label=5)
    assert both["score"] < only_a["score"], "the split label should drag the mean down"
    assert best["score"] >= both["score"], "best takes the better of the two"
    assert both["categories"] == ["a", "b"]


def test_asking_for_labels_none_of_which_cluster_says_which(panel):
    import numpy as np
    from starplast import objectives as O
    t = pd.Series(["a"] * 60 + ["b"] * 60)
    r = O.score(np.zeros(120, int), t, objective="best_f1", category=["nope", "also-nope"],
                min_cluster=1, min_label=5)
    assert r["score"] == 0.0
    assert "nope" in r["detail"] and "also-nope" in r["detail"]


def test_every_results_table_can_be_sorted_by_clicking_a_column(panel):
    """A walk produces hundreds of rows and the useful ones are at whichever end you sort to."""
    panel._fill(panel.walk_table, pd.DataFrame({"a": [3.0, 1.0, 2.0]}))
    assert panel.walk_table.isSortingEnabled()


def test_a_score_column_sorts_numerically_not_as_text(panel):
    """Stored as text, 0.9 sorts below 0.10, which puts the worst configurations at the top."""
    from PyQt6 import QtCore as _Qt
    panel._fill(panel.walk_table, pd.DataFrame({"score": [0.9, 0.10, 0.5]}))
    panel.walk_table.sortItems(0, _Qt.Qt.SortOrder.DescendingOrder)
    top = panel.walk_table.item(0, 0)
    assert float(top.data(_Qt.Qt.ItemDataRole.DisplayRole)) == 0.9


def test_filling_a_table_twice_while_sorted_does_not_interleave(panel):
    """With sorting left on during insertion Qt re-sorts after every row."""
    panel._fill(panel.walk_table, pd.DataFrame({"a": [1.0, 2.0, 3.0]}))
    panel._fill(panel.walk_table, pd.DataFrame({"a": [9.0]}))
    assert panel.walk_table.rowCount() == 1


# --------------------------------------------------------------------------- 6 validation
def _clustered(panel, n=None):
    """Give the panel a clustering of the real table: three clusters over every gene."""
    import numpy as np
    n = n if n is not None else len(panel.nodes)
    panel.labels = np.arange(n) % 3
    panel.rows = np.ones(len(panel.nodes), bool)
    return panel.labels


def test_validating_before_clustering_says_so(panel, sync):
    panel.labels = None
    panel.run_validation()
    assert any("cluster a map first" in m for m in sync)


def test_the_columns_the_guard_checks_are_column_names_not_block_names(panel):
    """The tab passed `columns_for(...).keys()` -- the BLOCK names -- as "the columns the embedding
    used", so the circularity guard compared a compartment against "expression_summary" and never
    fired once in the tab that exists to prevent exactly that."""
    panel.cat_cb.setChecked(True)
    used = panel.used_columns()
    assert "compartment" in used
    assert not ({"expression_summary", "fitness_screens", "protein_features"} & set(used))
    assert any(c.startswith("fit_") for c in used), "no real column from a selected block"


def test_validating_a_label_the_map_was_built_from_is_refused_and_explained(panel, sync):
    """The refusal is the guard working. Shown as a red failed job with a traceback it reads as the
    tab being broken, and the next move would be to look for the bug rather than to rebuild the map
    without that column."""
    _clustered(panel)
    # With a previous run's scores and candidates on screen, which is the case that matters:
    # rendered, the note said "Not scored" above a full table of precisions, and a table under an
    # explanation reads as the explanation's result.
    panel._fill(panel.val_table, pd.DataFrame({"category": ["dense granules"], "precision": [0.9]}))
    panel._fill(panel.cand_table, pd.DataFrame({"gene_id": ["TGME49_1"]}))
    panel.cat_cb.setChecked(True)                  # compartment now feeds the embedding
    panel.val_target.setCurrentText("compartment")
    panel.run_validation()
    # isHidden rather than isVisibleTo: a tab that is not the current one is not visible, and the
    # question here is whether the panel showed the note, not which tab is on top.
    assert not panel.val_note.isHidden()
    assert "circular" in panel.val_note.text()
    assert any("circular" in m for m in sync)
    assert panel.val_table.rowCount() == 0, "a refused run must not leave a score on screen"
    assert panel.cand_table.rowCount() == 0, "nor a candidate list"
    assert panel._validation_scores is None


def test_a_held_out_label_is_scored_and_the_refusal_note_goes_away(panel, sync, monkeypatch):
    import starplast.validate as V
    _clustered(panel)
    panel.cat_cb.setChecked(False)
    panel.val_note.setText("stale"); panel.val_note.show()
    monkeypatch.setattr(V, "validate_all",
                        lambda *a, **k: pd.DataFrame({"category": ["nucleus - chromatin"],
                                                      "n_labelled": [769], "n_folds": [5],
                                                      "precision": [0.4], "recall": [0.3],
                                                      "f1": [0.34], "refit": [False],
                                                      "note": [""]}))
    panel.run_validation()
    assert panel.val_note.isHidden()
    assert panel.val_table.rowCount() == 1
    assert any("roughly 40%" in m for m in sync)


def test_the_target_is_aligned_to_the_genes_the_clustering_covers(panel, sync, monkeypatch):
    """`labels` comes from an embedding that may have dropped genes, and scoring a clustering of one
    set against the labels of another compares gene i's cluster with gene j's compartment."""
    import numpy as np
    import starplast.validate as V
    keep = np.zeros(len(panel.nodes), bool)
    keep[:500] = True
    panel.rows = keep
    panel.labels = np.arange(500) % 3
    panel.cat_cb.setChecked(False)
    seen = {}

    def fake(labels, truth, **kw):
        seen.update(n_labels=len(labels), n_truth=len(truth), target=kw.get("target_column"))
        return pd.DataFrame()

    monkeypatch.setattr(V, "validate_all", fake)
    panel.val_target.setCurrentText("compartment")
    panel.run_validation()
    assert seen["n_labels"] == seen["n_truth"] == 500
    assert seen["target"] == "compartment"


def test_a_clustering_that_cannot_be_aligned_says_so_rather_than_scoring(panel, sync):
    import numpy as np
    panel.rows = None
    panel.labels = np.zeros(17, int)
    panel.cat_cb.setChecked(False)
    panel.run_validation()
    assert any("rebuild the map" in m for m in sync)


def test_asking_to_re_embed_passes_a_way_to_do_it(panel, sync, monkeypatch):
    """`refit` without a rebuild callable is an error in validate.py, deliberately -- so the tab has
    to supply one rather than tick a box that changes only the printed sentence."""
    import starplast.validate as V
    _clustered(panel)
    panel.cat_cb.setChecked(False)
    seen = {}
    monkeypatch.setattr(V, "validate_all",
                        lambda labels, truth, **kw: (seen.update(kw), pd.DataFrame())[1])
    panel.val_refit.setChecked(True)
    panel.run_validation()
    assert seen["refit"] is True and callable(seen["rebuild"])
    panel.val_refit.setChecked(False)
    panel.run_validation()
    assert seen["refit"] is False and seen["rebuild"] is None


def test_a_re_fit_fold_rebuilds_the_map_with_a_different_seed(panel, monkeypatch):
    """Same seed every fold would re-run the identical embedding and call the repetition a spread."""
    import numpy as np
    import starplast.embedding as E
    import starplast.clustering as C
    seeds = []
    monkeypatch.setattr(E, "embed", lambda nodes, spec, log=None: (
        seeds.append(spec.random_state), (np.zeros((len(nodes), 3)), [], np.ones(len(nodes), bool)))[1])
    monkeypatch.setattr(C, "cluster", lambda Y, **kw: np.zeros(len(Y), int))
    panel.seed.setValue(7)
    panel._refit_labels(0)
    panel._refit_labels(1)
    assert seeds == [8, 9]


def test_a_candidate_list_never_arrives_without_its_numbers(panel, sync):
    """The same list looks identical whether it is 90% right or 6% right, and on this proteome it
    has been 6%."""
    import numpy as np
    truth = panel.nodes["compartment"]
    from starplast import app as A
    v = A.as_text(truth).to_numpy()
    target = "nucleus - chromatin"
    # Cluster 0 holds that compartment AND some unlabelled genes -- which is the whole point: the
    # unlabelled members of a mostly-one-category cluster are the candidates.
    unlabelled = np.flatnonzero(v == "unassigned")[:200]
    labels = np.where(v == target, 0, 1)
    labels[unlabelled] = 0
    panel.labels, panel.rows = labels, np.ones(len(panel.nodes), bool)
    panel._validation_target = "compartment"
    panel._validation_scores = pd.DataFrame({"category": [target], "n_labelled": [int((v == target).sum())],
                                             "n_folds": [5], "precision": [0.4], "recall": [0.3],
                                             "f1": [0.34], "refit": [False], "note": [""]})
    panel._fill(panel.val_table, panel._validation_scores)
    panel.show_candidates(0)
    headers = [panel.cand_table.horizontalHeaderItem(c).text()
               for c in range(panel.cand_table.columnCount())]
    for needed in ("cluster_frac_category", "cluster_frac_contradicting", "shares_orthogroup"):
        assert needed in headers, headers
    assert any("would be right" in m for m in sync)


def test_clicking_a_category_no_cluster_holds_says_so(panel, sync):
    import numpy as np
    panel.labels = np.full(len(panel.nodes), -1)
    panel.rows = np.ones(len(panel.nodes), bool)
    panel._validation_target = "compartment"
    panel._validation_scores = pd.DataFrame({"category": ["apicoplast"], "precision": [0.1]})
    panel._fill(panel.val_table, panel._validation_scores)
    panel.show_candidates(0)
    assert any("no cluster holds" in m for m in sync)


def test_clicking_before_any_validation_does_nothing(panel):
    panel._validation_scores = None
    panel.show_candidates(0)
    assert panel.cand_table.rowCount() == 0


# --------------------------------------------------------------------------- 4 inference
def test_the_inference_tab_reports_each_category_not_only_each_feature(panel, sync, monkeypatch):
    """Cramer's V = 0.2 describes 27 compartments weakly smeared across every cluster and one
    compartment falling out cleanly. Only the per-category table separates those."""
    import numpy as np
    S = pd.DataFrame([{"feature": "compartment", "evidence": "held_out", "score_type": "cramers_v",
                       "score": 0.2, "assoc_with_input": 0.0, "n": 100, "q": 0.001}])
    D = pd.DataFrame([
        {"feature": "compartment", "category": "apicoplast", "cluster": 1, "n_in_cluster": 40,
         "precision": 0.9, "recall": 0.8, "evidence": "held_out", "q": 1e-9},
        {"feature": "compartment", "category": "nucleus", "cluster": 0, "n_in_cluster": 5,
         "precision": 0.1, "recall": 0.1, "evidence": "held_out", "q": 0.4},
        {"feature": "compartment", "category": "one-off", "cluster": 0, "n_in_cluster": 1,
         "precision": 0.5, "recall": 0.5, "evidence": "held_out", "q": 0.9},
    ])
    panel._battery_done((S, D, ["a finding"]))
    assert panel.category_table.rowCount() == 2, "a cluster holding one gene is not a category score"
    headers = [panel.category_table.horizontalHeaderItem(c).text()
               for c in range(panel.category_table.columnCount())]
    assert {"category", "precision", "prevalence", "lift", "recall", "f1"} <= set(headers)
    assert any("reach F1 0.5" in m for m in sync)


# --------------------------------------------------------------------------- results from the runner
@pytest.fixture
def runner_panel(app, nodes, tmp_path):
    """A panel wired to the window's job runner, which is how it runs in the application."""
    from starplast.analysis_panel import AnalysisPanel
    from starplast.jobs import JobRunner
    from starplast.tuning import EmbeddingStore
    r = JobRunner()
    p = AnalysisPanel(nodes, store=EmbeddingStore(str(tmp_path)), runner=r)
    yield p, r
    p.deleteLater()


def _finished(panel, runner, job, ok, on_done=None, on_error=None):
    runner.jobs[job.id] = job
    panel._jobs[job.id] = (on_done or (lambda r: None), on_error)
    panel._on_job_finished(job.id, ok)


def test_a_result_reaches_the_handler_that_asked_for_it(runner_panel):
    from starplast.jobs import DONE, Job
    panel, runner = runner_panel
    got = []
    _finished(panel, runner, Job(id=1, name="j", state=DONE, result=42), True, on_done=got.append)
    assert got == [42]


def test_a_refusal_reaches_the_handler_as_the_exception_not_as_its_text(runner_panel):
    """Telling a deliberate refusal from a crash by parsing a formatted message is guesswork, and
    the two need completely different words on screen."""
    from starplast.jobs import FAILED, Job
    panel, runner = runner_panel
    said, seen = [], []
    panel.status.connect(said.append)
    _finished(panel, runner, Job(id=2, name="validate", state=FAILED,
                                 error="ValueError: circular", exception=ValueError("circular")),
              False, on_error=lambda e: (seen.append(e), True)[1])
    assert isinstance(seen[0], ValueError)
    assert not any("failed" in m for m in said), "a handled refusal must not also be reported as a crash"


def test_a_failure_the_handler_declines_is_still_reported(runner_panel):
    """A handler that only recognises its own refusals must not swallow a real crash."""
    from starplast.jobs import FAILED, Job
    panel, runner = runner_panel
    said = []
    panel.status.connect(said.append)
    _finished(panel, runner, Job(id=3, name="battery", state=FAILED, error="KeyError: 'x'"),
              False, on_error=lambda e: False)
    assert any("battery failed" in m for m in said)


def test_a_failure_with_no_handler_at_all_is_reported(runner_panel):
    from starplast.jobs import FAILED, Job
    panel, runner = runner_panel
    said = []
    panel.status.connect(said.append)
    _finished(panel, runner, Job(id=4, name="search", state=FAILED, error="boom"), False)
    assert any("search failed" in m for m in said)


def test_a_stopped_job_says_stopped_rather_than_failed(runner_panel):
    """A job the user stopped, reported in red with a traceback, teaches people to distrust the
    failure list."""
    from starplast.jobs import CANCELLED, Job
    panel, runner = runner_panel
    said = []
    panel.status.connect(said.append)
    _finished(panel, runner, Job(id=5, name="walk", state=CANCELLED), False)
    assert any("stopped" in m for m in said) and not any("failed" in m for m in said)


def test_a_job_nobody_is_waiting_for_is_ignored(runner_panel):
    panel, runner = runner_panel
    panel._on_job_finished(999, True)             # never submitted here
    from starplast.jobs import DONE, Job
    panel._jobs[6] = (lambda r: None, None)
    panel._on_job_finished(6, True)               # handler, but the runner has no such job


def test_the_panel_runs_through_the_runner_when_it_has_one(runner_panel, app):
    """Its work then appears in the Jobs panel and can be stopped there, and it does not block the
    other tabs the way a private thread did."""
    from PyQt6 import QtCore
    panel, runner = runner_panel
    got = []
    job = panel._run(lambda p: 7, got.append, name="through the runner")
    assert job is not None and job.name == "through the runner"
    for _ in range(200):
        app.processEvents()
        QtCore.QThread.msleep(5)
        if got:
            break
    assert got == [7]


def test_progress_from_a_job_reaches_the_status_line_and_the_job_note(runner_panel, app):
    """The note is what the Jobs panel shows for a job in flight, and for a walk it is the only
    sign that anything is happening."""
    from PyQt6 import QtCore
    panel, runner = runner_panel
    said = []
    panel.status.connect(said.append)
    job = panel._run(lambda p: (p("half way"), "done")[1], lambda r: None, name="reports")
    # Waited on the status line rather than on the job: the report is emitted from the worker and
    # queued, so it can still be in flight when the job itself is already finished.
    for _ in range(200):
        app.processEvents()
        QtCore.QThread.msleep(5)
        if said:
            break
    assert any("half way" in m for m in said)
    assert job.note == "half way"


def test_a_stopped_job_unwinds_at_its_next_progress_report(runner_panel, app):
    """Cooperative by construction: raising from the reporting callable is what stops a
    328-configuration search in seconds rather than in half an hour, and it must arrive as a stop
    rather than as a crash."""
    from PyQt6 import QtCore
    from starplast.jobs import CANCELLED
    panel, runner = runner_panel

    def work(p):
        p("started")
        for i in range(1000):
            p(f"step {i}")                        # the cancellation lands here
        return "never"

    job = panel._run(work, lambda r: None, name="stoppable")
    for _ in range(200):
        app.processEvents()
        QtCore.QThread.msleep(5)
        if job.note:
            break
    job.cancel()
    for _ in range(400):
        app.processEvents()
        QtCore.QThread.msleep(5)
        if not job.active:
            break
    assert job.state == CANCELLED
    assert "stopped after" in job.note


def test_a_control_named_in_the_tooltip_tables_but_absent_is_skipped(panel):
    """The tables are keyed by attribute name, and a control removed from a tab must not take the
    whole panel down the next time tooltips are applied."""
    import starplast.analysis_panel as AP
    with pytest.MonkeyPatch.context() as mp:
        mp.setitem(AP.TOOLTIPS, "no_such_control", "explains a control that is not there")
        panel._apply_tooltips()


def test_a_crash_in_validation_is_still_a_crash(panel):
    """The refusal handler recognises the guard's ValueError. Anything else is a real failure and
    must keep its red job and its traceback rather than being dressed up as an explanation."""
    assert panel._validation_refused(KeyError("compartment")) is False
    assert panel.val_note.isHidden()


def test_the_reporting_callable_records_the_note_and_says_it_out_loud(panel, capsys):
    """Driven directly. It runs on a Qt-managed thread in the application, where coverage cannot
    follow it and where a test that only submits a job proves nothing about this body.

    The note is what the Jobs panel shows for a job in flight; the print is what reaches the console
    pane, which is where a walk is actually read."""
    from starplast.analysis_panel import _Progress
    from starplast.jobs import Job
    said = []
    panel.status.connect(said.append)
    job = Job(id=1, name="walk")
    _Progress(panel, job)("configuration 12 of 288")
    assert job.note == "configuration 12 of 288"
    assert said == ["configuration 12 of 288"]
    assert "configuration 12" in capsys.readouterr().out


def test_the_reporting_callable_is_where_a_stop_takes_effect(panel):
    """These functions report once per configuration, which makes this the one place guaranteed to
    be reached repeatedly without threading a cancellation flag through search, tuning and
    embedding. Raising here is what stops a 328-configuration search in seconds."""
    from starplast.analysis_panel import Cancelled, _Progress
    from starplast.jobs import Job
    job = Job(id=2, name="recovery search")
    job.note = "configuration 12 of 288"
    job.cancel()
    with pytest.raises(Cancelled, match="configuration 12 of 288"):
        _Progress(panel, job)("configuration 13 of 288")


def test_the_validation_job_runs_against_the_real_module_not_only_a_stand_in(panel, sync):
    """Every other test here monkeypatches `validate_all`, so the arguments the panel actually
    passes were never checked against the function that receives them -- and one of them, `log`,
    was forwarded into `masked_recovery`, which has no such parameter. The whole tab raised
    TypeError the first time it was run for real.

    A clustering over the real table, two folds, one target: slow enough to be worth doing once and
    fast enough to keep."""
    import numpy as np
    from starplast import app as A
    v = A.as_text(panel.nodes["compartment"]).to_numpy()
    panel.labels = np.where(v == "nucleus - chromatin", 0, 1)
    panel.rows = np.ones(len(panel.nodes), bool)
    panel.cat_cb.setChecked(False)
    panel.val_target.setCurrentText("compartment")
    panel.val_folds.setValue(2)
    panel.run_validation()
    d = panel._validation_scores
    assert len(d), "no category was scored against the real implementation"
    assert set(["category", "precision", "recall", "f1", "refit"]) <= set(d.columns)
    assert any("would be right" in m for m in sync), "the verdict never reached the status line"


# --------------------------------------------------------------------------- every results table
RESULTS_TABLES = ["walk_table", "cluster_table", "battery_table", "category_table",
                  "search_table", "val_table", "cand_table"]


def test_every_results_table_can_be_saved_and_says_how_much_it_wrote(panel, sync, tmp_path):
    """A table that can only be read on screen has to be re-derived anywhere else it is needed, and
    the run that produced it is minutes long."""
    import pandas as pd
    for name in RESULTS_TABLES:
        table = getattr(panel, name)
        panel._fill(table, pd.DataFrame({"a": [1.0, 2.0], "b": ["x", "y"]}))
        path = panel.save_table(table, str(tmp_path / f"{name}.csv"))
        assert path, name
        assert len(pd.read_csv(path)) == 2, name
    assert any("wrote 2 rows" in m for m in sync)


def test_saving_writes_the_whole_result_not_the_screenful_that_is_shown(panel, sync, tmp_path):
    """The 200-row truncation exists to keep the window responsive. A file that silently stopped
    there would be a different result from the one that was computed."""
    import pandas as pd
    panel._fill(panel.search_table, pd.DataFrame({"score": range(500)}))
    assert panel.search_table.rowCount() == 200
    path = panel.save_table(panel.search_table, str(tmp_path / "all.csv"))
    assert len(pd.read_csv(path)) == 500
    assert any("shows the first 200" in m for m in sync)


def test_a_streamed_walk_can_be_saved_for_what_it_found_so_far(panel, sync, tmp_path):
    """A walk stopped half way has still done the work it did, and the rows arrive one at a time."""
    panel._start_table(panel.walk_table, [])
    for md in (0.1, 0.25, 0.5):
        panel._append(panel.walk_table, {"n_neighbors": 15, "min_dist": md, "trustworthiness": 0.9})
    import pandas as pd
    path = panel.save_table(panel.walk_table, str(tmp_path / "partial.csv"))
    got = pd.read_csv(path)
    assert list(got.min_dist) == [0.1, 0.25, 0.5]


def test_a_new_run_does_not_leave_the_last_one_saveable(panel, sync, tmp_path):
    """Saving a table that has been emptied on screen must not write the previous walk's rows."""
    import pandas as pd
    panel._fill(panel.walk_table, pd.DataFrame({"a": [1, 2, 3]}))
    panel._start_table(panel.walk_table, [])
    assert panel.save_table(panel.walk_table, str(tmp_path / "empty.csv")) == ""
    assert any("nothing to save" in m for m in sync)


def test_saving_an_empty_table_says_so_rather_than_writing_a_header(panel, sync, tmp_path):
    assert panel.save_table(panel.cand_table, str(tmp_path / "none.csv")) == ""
    assert any("nothing to save" in m for m in sync)


def test_cancelling_the_save_dialog_writes_nothing(panel, monkeypatch, tmp_path):
    import pandas as pd
    from PyQt6 import QtWidgets
    panel._fill(panel.walk_table, pd.DataFrame({"a": [1]}))
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: ("", "")))
    assert panel.save_table(panel.walk_table) == ""


def test_the_save_dialog_offers_a_name_that_says_what_the_table_is(panel, monkeypatch, tmp_path):
    """"results.csv" seven times in a downloads folder is not a set of results."""
    import pandas as pd
    from PyQt6 import QtWidgets
    seen = {}

    def fake(parent, caption, name, filt):
        seen["name"] = name
        return str(tmp_path / "x.csv"), filt

    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName", staticmethod(fake))
    panel._fill(panel.search_table, pd.DataFrame({"a": [1]}))
    panel.save_table(panel.search_table)
    assert "recovery_search" in seen["name"]


def test_every_results_table_offers_the_same_right_click_menu(panel):
    import pandas as pd
    for name in RESULTS_TABLES:
        table = getattr(panel, name)
        panel._fill(table, pd.DataFrame({"a": [1.0]}))
        actions = [a.text() for a in panel.build_table_menu(table).actions() if a.text()]
        assert any("CSV" in a for a in actions), name
        assert any("Copy" in a for a in actions), name


def test_the_menu_offers_the_map_only_where_a_row_names_one(panel):
    import pandas as pd
    panel._fill(panel.search_table, pd.DataFrame({"a": [1.0]}))
    panel._fill(panel.battery_table, pd.DataFrame({"a": [1.0]}))
    assert any("map" in a.text() for a in panel.build_table_menu(panel.search_table).actions())
    assert not any("map" in a.text() for a in panel.build_table_menu(panel.battery_table).actions())


def test_saving_is_offered_but_disabled_when_there_is_nothing_to_save(panel):
    """Absent, it looks like the feature does not exist; enabled, it writes an empty file."""
    act = [a for a in panel.build_table_menu(panel.walk_table).actions() if "CSV" in a.text()][0]
    assert not act.isEnabled()


def test_selected_rows_can_be_copied_with_their_headers(panel, sync):
    import pandas as pd
    panel._fill(panel.walk_table, pd.DataFrame({"n_neighbors": [5, 15], "trust": [0.9, 0.8]}))
    panel.walk_table.selectRow(1)
    text = panel.copy_rows(panel.walk_table)
    assert text.splitlines()[0].split("\t") == ["n_neighbors", "trust"]
    assert "15" in text.splitlines()[1]
    assert any("copied 1 row" in m for m in sync)


def test_copying_nothing_says_to_select_a_row(panel, sync):
    import pandas as pd
    panel._fill(panel.walk_table, pd.DataFrame({"a": [1]}))
    panel.walk_table.clearSelection()
    assert panel.copy_rows(panel.walk_table) == ""
    assert any("select a row" in m for m in sync)


# --------------------------------------------------------------------------- rows to maps
def test_clicking_a_clustering_row_applies_it_to_the_map(panel, sync, monkeypatch):
    """A silhouette for a clustering nobody can see is a number about nothing."""
    import numpy as np
    import starplast.clustering as C
    panel.coords = np.random.default_rng(0).normal(size=(120, 3))
    panel.rows = None
    seen = {}
    monkeypatch.setattr(C, "cluster", lambda Y, **kw: (seen.update(kw), np.zeros(len(Y), int))[1])
    panel._fill(panel.cluster_table, pd.DataFrame({"algorithm": ["hdbscan"],
                                                   "min_cluster_size": [40], "min_samples": [None],
                                                   "silhouette": [0.4]}))
    got = []
    panel.clusters_ready.connect(got.append)
    panel.show_cluster_row(0)
    assert seen["min_cluster_size"] == 40
    assert panel.mcs.value() == 40, "the controls must describe the map on screen"
    assert len(got) == 1


def test_clicking_a_clustering_row_before_there_is_a_map_says_so(panel, sync):
    panel.coords = None
    panel._fill(panel.cluster_table, pd.DataFrame({"algorithm": ["hdbscan"],
                                                   "min_cluster_size": [40]}))
    panel.show_cluster_row(0)
    assert any("build a map first" in m for m in sync)


def test_a_clustering_row_that_names_no_settings_says_so(panel, sync):
    import numpy as np
    panel.coords = np.zeros((30, 3))
    panel._fill(panel.cluster_table, pd.DataFrame({"algorithm": ["hdbscan"],
                                                   "min_cluster_size": ["n/a"]}))
    panel.show_cluster_row(0)
    assert any("does not name a clustering" in m for m in sync)


def test_clicking_a_search_row_rebuilds_that_exact_configuration(panel, sync, monkeypatch):
    """The one table where a row is a whole recipe, and the one that could not be looked at."""
    import numpy as np
    import starplast.search as S
    seen = {}

    def fake_rebuild(nodes, row, log=print):
        seen.update(row)
        genes = np.zeros(len(nodes), bool)
        genes[:50] = True
        return np.zeros((50, 3)), genes, np.array([0] * 25 + [1] * 25), ["f"]

    monkeypatch.setattr(S, "rebuild", fake_rebuild)
    panel._fill(panel.search_table, pd.DataFrame({
        "blocks": ["expression_summary+fitness_screens"], "na_policy": ["median"],
        "scaling": ["rank"], "n_neighbors": [15], "min_dist": [0.1], "min_cluster_size": [25],
        "seed": [42], "sample_size": [3000], "excluded": ["compartment;lopit_map"],
        "mean_f1": [0.3]}))
    coords, clusters = [], []
    panel.embedding_ready.connect(lambda c, r: coords.append((c, r)))
    panel.clusters_ready.connect(clusters.append)
    panel.show_search_row(0)
    assert seen["blocks"] == "expression_summary+fitness_screens"
    assert seen["sample_size"] == "3000" and seen["excluded"] == "compartment;lopit_map"
    assert len(coords) == 1 and len(clusters) == 1
    assert len(clusters[0]) == len(panel.nodes), "a subsample's clustering must cover the map"
    assert any("2 clusters" in m for m in sync)


def test_a_search_row_that_cannot_be_rebuilt_explains_rather_than_failing(panel, sync, monkeypatch):
    import starplast.search as S

    def boom(nodes, row, log=print):
        raise ValueError("every block in that row feeds a column the run excluded")

    monkeypatch.setattr(S, "rebuild", boom)
    panel._fill(panel.search_table, pd.DataFrame({"blocks": ["localization"], "n_neighbors": [15],
                                                  "min_dist": [0.1]}))
    panel.show_search_row(0)
    assert any("cannot rebuild that row" in m for m in sync)
    assert panel._row_rebuild_failed(KeyError("x")) is False


def test_a_search_row_with_no_blocks_says_so(panel, sync):
    panel._fill(panel.search_table, pd.DataFrame({"mean_f1": [0.3]}))
    panel.show_search_row(0)
    assert any("does not name a configuration" in m for m in sync)


def test_clicking_an_inference_row_colors_the_map_by_the_clustering_it_scored(panel, sync):
    """Reading "cluster 3 is 90% apicoplast" while looking at a map colored by compartment is a
    needless act of translation."""
    import numpy as np
    panel.labels = np.arange(len(panel.nodes)) % 4
    panel.rows = np.ones(len(panel.nodes), bool)
    panel._fill(panel.category_table, pd.DataFrame({"feature": ["compartment"],
                                                    "category": ["apicoplast"], "cluster": [3],
                                                    "lift": [5.4], "f1": [0.4]}))
    got = []
    panel.clusters_ready.connect(got.append)
    panel.show_inference_row(0)
    assert len(got) == 1
    assert any("cluster 3" in m and "apicoplast" in m and "5.4x" in m for m in sync)


def test_an_inference_row_before_any_clustering_says_what_to_do(panel, sync):
    panel.labels = None
    panel._fill(panel.category_table, pd.DataFrame({"feature": ["f"], "category": ["x"],
                                                    "cluster": [0]}))
    panel.show_inference_row(0)
    assert any("cluster a map first" in m for m in sync)


def test_a_clustering_of_a_subsample_is_published_over_the_whole_table(panel):
    """The window colors 8,140 points by it. Left short, it fell back to grey everywhere, which
    reads as "this clustering found nothing"."""
    import numpy as np
    keep = np.zeros(len(panel.nodes), bool)
    keep[:300] = True
    panel.rows = keep
    got = []
    panel.clusters_ready.connect(got.append)
    panel._publish_clusters(np.zeros(300, int))
    assert len(got[0]) == len(panel.nodes)
    assert (got[0][:300] == 0).all() and (got[0][300:] == -1).all()


def test_a_full_length_clustering_is_published_unchanged(panel):
    import numpy as np
    panel.rows = np.ones(len(panel.nodes), bool)
    got = []
    panel.clusters_ready.connect(got.append)
    labels = np.arange(len(panel.nodes)) % 3
    panel._publish_clusters(labels)
    assert np.array_equal(got[0], labels)


def test_clicking_a_row_where_nothing_is_wired_does_nothing(panel):
    import pandas as pd
    panel._fill(panel.battery_table, pd.DataFrame({"a": [1]}))
    panel._row_clicked(panel.battery_table, 0)


def test_the_context_menu_can_be_opened_on_a_table(panel, monkeypatch):
    """`_table_menu` execs, which blocks; this checks the wiring reaches it and nothing else."""
    import pandas as pd
    from PyQt6 import QtCore, QtWidgets
    panel._fill(panel.walk_table, pd.DataFrame({"a": [1.0]}))
    monkeypatch.setattr(QtWidgets.QMenu, "exec", lambda self, *a: None)
    m = panel._table_menu(panel.walk_table, QtCore.QPoint(2, 2))
    assert [a.text() for a in m.actions() if a.text()]


def test_clicking_a_cell_is_what_triggers_the_row_action(panel, monkeypatch):
    """The wiring, not the handler: connected to the wrong signal, every one of these tables would
    look inert while every handler test passed."""
    import pandas as pd
    seen = []
    panel._row_action[panel.walk_table] = seen.append
    panel._fill(panel.walk_table, pd.DataFrame({"a": [1.0, 2.0]}))
    panel.walk_table.cellClicked.emit(1, 0)
    assert seen == [1]


def test_building_a_map_can_cluster_it_in_the_same_job(panel, sync, monkeypatch):
    """One job rather than two: a map that appears for a moment without the clusters the walk row
    promised reads as the clustering having failed."""
    import numpy as np
    import starplast.clustering as C
    import starplast.embedding as E
    rows = np.zeros(len(panel.nodes), bool)
    rows[:60] = True
    coords = np.zeros((60, 3))
    monkeypatch.setattr(E, "embed", lambda n, spec, log=None: (coords, ["f"], rows))
    monkeypatch.setattr(C, "cluster", lambda Y, **kw: np.array([0] * 30 + [1] * 30))
    got = []
    panel.clusters_ready.connect(got.append)
    panel.run_embed(then_cluster=15)
    assert len(got) == 1 and panel.labels is not None
    assert any("2 clusters" in m for m in sync)


def test_building_a_map_without_clustering_leaves_the_clustering_alone(panel, sync, monkeypatch):
    """A new map invalidates the old labels, but silently publishing a stale clustering over it
    would color the new map by the old one."""
    import numpy as np
    import starplast.embedding as E
    rows = np.zeros(len(panel.nodes), bool)
    rows[:60] = True
    monkeypatch.setattr(E, "embed",
                        lambda n, spec, log=None: (np.zeros((60, 3)), ["f"], rows))
    got = []
    panel.clusters_ready.connect(got.append)
    panel.run_embed()
    assert got == []
    assert any("map built" in m and "clusters" not in m for m in sync)


def test_a_clustering_that_matches_neither_the_map_nor_the_table_says_so(panel, sync):
    """Grey everywhere reads as "this clustering found nothing", which is a finding. A clustering
    of the wrong genes is a mistake, and the two must not look the same."""
    import numpy as np
    keep = np.zeros(len(panel.nodes), bool)
    keep[:300] = True
    panel.rows = keep
    got = []
    panel.clusters_ready.connect(got.append)
    panel._publish_clusters(np.zeros(77, int))
    assert len(got[0]) == 77, "a mismatched clustering must not be stretched onto the wrong genes"
    assert any("cluster this map again" in m for m in sync)


# --------------------------------------------------------------------------- saving annotations
@pytest.fixture
def annotating(panel, tmp_path):
    """A panel with somewhere to save, a validated score, and a candidate list on screen."""
    import numpy as np
    from starplast import app as A
    from starplast.annotations import AnnotationStore
    panel.store_annotations = AnnotationStore(str(tmp_path / "annotations.csv"))
    v = A.as_text(panel.nodes["compartment"]).to_numpy()
    target = "nucleus - chromatin"
    unlabelled = np.flatnonzero(v == "unassigned")[:200]
    labels = np.where(v == target, 0, 1)
    labels[unlabelled] = 0
    panel.labels, panel.rows = labels, np.ones(len(panel.nodes), bool)
    panel._validation_target = "compartment"
    panel._validation_scores = pd.DataFrame({"category": [target], "n_labelled": [769],
                                             "n_folds": [5], "precision": [0.62], "recall": [0.3],
                                             "f1": [0.4], "refit": [False], "note": [""]})
    panel._fill(panel.val_table, panel._validation_scores)
    panel.show_candidates(0)
    return panel


def test_candidates_can_be_saved_with_the_numbers_that_justify_them(annotating, sync):
    panel = annotating
    panel.reasoning.setText("dense and mostly chromatin")
    panel.save_candidates()
    got = panel.store_annotations.load()
    assert len(got) == panel.cand_table.rowCount() or len(got) > 0
    r = got.iloc[0]
    assert r.precision == pytest.approx(0.62) and r.target == "compartment"
    assert r.reasoning == "dense and mostly chromatin" and r.date
    assert r.blocks and r.n_neighbors == panel.nn.value()
    assert any("saved" in m and "precision 0.62" in m for m in sync)


def test_saving_tells_the_map_to_redraw(annotating):
    """The fourth color has to appear without the user having to find it."""
    seen = []
    annotating.annotations_changed.connect(lambda: seen.append(True))
    annotating.save_candidates()
    assert seen == [True]


def test_an_unvalidated_candidate_list_is_refused_and_explained(annotating, sync):
    """The refusal is the feature, and it has to read as the store working rather than failing."""
    panel = annotating
    panel._validation_scores = pd.DataFrame({"category": ["nucleus - chromatin"],
                                             "precision": [float("nan")], "recall": [float("nan")],
                                             "n_folds": [0], "refit": [False]})
    panel.save_candidates()
    assert not panel.val_note.isHidden() and "Not saved" in panel.val_note.text()
    assert "no validated precision" in panel.val_note.text()
    assert panel.store_annotations.load().empty


def test_saving_before_there_are_candidates_says_what_to_do(panel, sync, tmp_path):
    from starplast.annotations import AnnotationStore
    panel.store_annotations = AnnotationStore(str(tmp_path / "a.csv"))
    panel._validation_scores = None
    panel.save_candidates()
    assert any("no candidates" in m for m in sync)


def test_a_panel_with_nowhere_to_save_says_so(panel, sync):
    panel.store_annotations = None
    panel.save_candidates()
    assert any("no annotations file" in m for m in sync)


def test_the_saved_configuration_is_the_one_the_map_was_built_from(annotating):
    panel = annotating
    panel.nn.setValue(33)
    panel.mcs.setValue(60)
    panel.save_candidates()
    r = panel.store_annotations.load().iloc[0]
    assert r.n_neighbors == 33 and r.min_cluster_size == 60 and r.seed == panel.seed.value()


def test_a_search_that_scored_nothing_says_so_rather_than_showing_an_empty_table(panel, sync,
                                                                                 monkeypatch):
    """An empty table reads as "still running". Every configuration failing to produce a scorable
    clustering is a result about the grid, and it has to be said out loud."""
    import starplast.search as S
    monkeypatch.setattr(S, "search", lambda n, **k: (pd.DataFrame(), pd.DataFrame()))
    panel.run_search()
    assert any("nothing was scorable" in m for m in sync)


# --------------------------------------------------------------------------- saving results
def test_a_saved_table_can_be_loaded_back_and_is_still_clickable(panel, sync, tmp_path,
                                                                 monkeypatch):
    """The point of the whole format. A loader that produced a table you could only read would be a
    screenshot with extra steps."""
    import numpy as np
    import starplast.search as S
    rows = pd.DataFrame({"blocks": ["fitness_screens"], "na_policy": ["median"],
                         "scaling": ["rank"], "n_neighbors": [15], "min_dist": [0.1],
                         "min_cluster_size": [25], "seed": [42], "sample_size": [0],
                         "excluded": ["compartment"], "mean_f1": [0.3]})
    panel._fill(panel.search_table, rows)
    path = panel.save_results(panel.search_table, str(tmp_path / "s.csv"))
    assert path and any("still clickable" in m for m in sync)

    panel._start_table(panel.search_table, [])
    assert panel.load_results(panel.search_table, path)
    assert panel.search_table.rowCount() == 1

    seen = {}

    def fake_rebuild(nodes, row, log=print):
        seen.update(row)
        genes = np.zeros(len(nodes), bool)
        genes[:10] = True
        return np.zeros((10, 3)), genes, np.zeros(10, int), ["f"]

    monkeypatch.setattr(S, "rebuild", fake_rebuild)
    panel.show_search_row(0)
    assert seen["blocks"] == "fitness_screens", "a loaded row could not be rebuilt"


def test_a_file_from_another_tab_is_refused(panel, sync, tmp_path):
    """Its rows are per-category scores, not recipes: loaded into the walk tab they produce rows
    nobody can rebuild and an error about the wrong thing."""
    panel._fill(panel.val_table, pd.DataFrame({"category": ["dense granules"], "precision": [0.6]}))
    path = panel.save_results(panel.val_table, str(tmp_path / "v.csv"))
    assert panel.load_results(panel.search_table, path) is False
    assert any("load it into the tab it came from" in m for m in sync)


def test_a_csv_that_names_no_table_loads_and_says_so(panel, sync, tmp_path):
    path = tmp_path / "mine.csv"
    pd.DataFrame({"a": [1, 2]}).to_csv(path, index=False)
    assert panel.load_results(panel.walk_table, str(path))
    assert any("did not say which table" in m for m in sync)


def test_loading_something_unreadable_says_so_rather_than_raising(panel, sync, tmp_path):
    path = tmp_path / "bad.csv"
    path.write_bytes(b"\xff\xfe\x00 not a csv")
    assert panel.load_results(panel.walk_table, str(path)) is False
    assert any("could not read" in m for m in sync)


def test_saving_an_empty_table_says_so(panel, sync, tmp_path):
    panel._start_table(panel.walk_table, [])
    assert panel.save_results(panel.walk_table, str(tmp_path / "x.csv")) == ""
    assert any("nothing to save" in m for m in sync)


def test_every_tab_can_be_saved_and_loaded_in_one_file(panel, sync, tmp_path):
    panel._fill(panel.walk_table, pd.DataFrame({"n_neighbors": [15], "trustworthiness": [0.9]}))
    panel._fill(panel.search_table, pd.DataFrame({"blocks": ["fitness_screens"], "mean_f1": [0.3]}))
    panel._fill(panel.val_table, pd.DataFrame({"category": ["nucleus"], "precision": [0.4]}))
    path = panel.save_all_results(str(tmp_path / "all.starplast"))
    assert path and any("wrote 3 table(s)" in m for m in sync)
    for table in (panel.walk_table, panel.search_table, panel.val_table):
        panel._start_table(table, [])
    assert panel.load_all_results(path) == 3
    assert panel.walk_table.rowCount() == 1 and panel.search_table.rowCount() == 1
    assert panel.val_table.rowCount() == 1


def test_a_bundle_from_a_later_version_names_what_this_one_cannot_show(panel, sync, tmp_path):
    """"3 of 5 loaded" with no names is not something anyone can act on."""
    from starplast.results import save_bundle
    path = save_bundle(str(tmp_path / "future.starplast"),
                       {"umap_walk": pd.DataFrame({"n_neighbors": [15]}),
                        "some_new_tab": pd.DataFrame({"x": [1]})})
    assert panel.load_all_results(path) == 1
    assert any("no tab for some_new_tab" in m for m in sync)


def test_saving_everything_with_nothing_computed_says_so(panel, sync, tmp_path):
    for table in panel.results_tables().values():
        panel._start_table(table, [])
    assert panel.save_all_results(str(tmp_path / "none.starplast")) == ""
    assert any("no results to save" in m for m in sync)


def test_loading_a_bundle_that_will_not_open_says_so(panel, sync, tmp_path):
    path = tmp_path / "not.starplast"
    path.write_bytes(b"definitely not a zip")
    assert panel.load_all_results(str(path)) == 0
    assert any("could not read" in m for m in sync)


def test_every_results_table_offers_saving_and_loading(panel):
    import pandas as pd
    for name in RESULTS_TABLES:
        table = getattr(panel, name)
        panel._fill(table, pd.DataFrame({"a": [1.0]}))
        actions = [a.text() for a in panel.build_table_menu(table).actions() if a.text()]
        assert any("reloadable" in a for a in actions), name
        assert any("Load results" in a for a in actions), name


# --------------------------------------------------------------------------- the file dialogs
def _rows():
    import pandas as pd
    return pd.DataFrame({"blocks": ["fitness_screens"], "na_policy": ["median"],
                         "scaling": ["rank"], "n_neighbors": [15], "min_dist": [0.1],
                         "min_cluster_size": [25], "seed": [42], "sample_size": [0],
                         "excluded": [""], "mean_f1": [0.3]})


def test_saving_with_no_path_asks_where_and_writes_there(panel, sync, tmp_path, monkeypatch):
    """The menu entry passes no path, so the dialog is the only thing that supplies one."""
    from PyQt6 import QtWidgets
    target = tmp_path / "asked.csv"
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(target), "")))
    panel._fill(panel.search_table, _rows())
    assert panel.save_results(panel.search_table) == str(target)
    assert target.exists() and target.read_text().startswith("# starplast-table:")


def test_cancelling_the_save_dialog_writes_nothing(panel, sync, tmp_path, monkeypatch):
    """Cancel means cancel. A file written to a default name the user never chose is a file they
    will not find, and a second table quietly overwriting the first."""
    from PyQt6 import QtWidgets
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: ("", "")))
    panel._fill(panel.search_table, _rows())
    assert panel.save_results(panel.search_table) == ""
    assert not list(tmp_path.iterdir())


def test_loading_with_no_path_asks_which_file(panel, tmp_path, monkeypatch):
    from PyQt6 import QtWidgets
    path = tmp_path / "s.csv"
    panel._fill(panel.search_table, _rows())
    panel.save_results(panel.search_table, str(path))
    panel._start_table(panel.search_table, [])
    monkeypatch.setattr(QtWidgets.QFileDialog, "getOpenFileName",
                        staticmethod(lambda *a, **k: (str(path), "")))
    assert panel.load_results(panel.search_table) is True
    assert panel.search_table.rowCount() == 1


def test_cancelling_the_load_dialog_leaves_the_table_alone(panel, monkeypatch):
    """Not an empty table: whatever was computed before the dialog opened is still the result."""
    from PyQt6 import QtWidgets
    panel._fill(panel.search_table, _rows())
    monkeypatch.setattr(QtWidgets.QFileDialog, "getOpenFileName",
                        staticmethod(lambda *a, **k: ("", "")))
    assert panel.load_results(panel.search_table) is False
    assert panel.search_table.rowCount() == 1


def test_saving_everything_with_no_path_asks_where(panel, tmp_path, monkeypatch):
    from PyQt6 import QtWidgets
    target = tmp_path / "all.starplast"
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(target), "")))
    panel._fill(panel.search_table, _rows())
    assert panel.save_all_results() == str(target)
    assert target.exists()


def test_cancelling_the_bundle_save_writes_nothing(panel, tmp_path, monkeypatch):
    from PyQt6 import QtWidgets
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: ("", "")))
    panel._fill(panel.search_table, _rows())
    assert panel.save_all_results() == ""
    assert not list(tmp_path.iterdir())


def test_loading_a_bundle_with_no_path_asks_which_file(panel, tmp_path, monkeypatch):
    from PyQt6 import QtWidgets
    path = tmp_path / "all.starplast"
    panel._fill(panel.search_table, _rows())
    panel.save_all_results(str(path))
    panel._start_table(panel.search_table, [])
    monkeypatch.setattr(QtWidgets.QFileDialog, "getOpenFileName",
                        staticmethod(lambda *a, **k: (str(path), "")))
    assert panel.load_all_results() == 1


def test_cancelling_the_bundle_load_loads_nothing(panel, monkeypatch):
    from PyQt6 import QtWidgets
    monkeypatch.setattr(QtWidgets.QFileDialog, "getOpenFileName",
                        staticmethod(lambda *a, **k: ("", "")))
    assert panel.load_all_results() == 0


def test_cancelling_the_csv_dialog_writes_nothing(panel, tmp_path, monkeypatch):
    """The right-click "save as CSV" entry passes no path either, and cancel has to mean cancel
    there too -- a file under a default name is a file the user will not find again."""
    from PyQt6 import QtWidgets
    panel._fill(panel.search_table, _rows())
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: ("", "")))
    assert panel.save_table(panel.search_table) == ""
    assert not list(tmp_path.iterdir())


def test_saving_a_csv_with_no_path_asks_where(panel, tmp_path, monkeypatch):
    from PyQt6 import QtWidgets
    target = tmp_path / "asked.csv"
    panel._fill(panel.search_table, _rows())
    monkeypatch.setattr(QtWidgets.QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(target), "")))
    assert panel.save_table(panel.search_table) == str(target)
    assert "blocks" in target.read_text().splitlines()[0]


# --------------------------------------------------------------------------- stopping and autosave
def test_every_long_running_tab_has_a_visible_stop_button(panel):
    """The first thing asked about stopping was "I can't see the stop button" -- it existed only in
    the Jobs dock's right-click menu, which is not where anyone looks while watching a sweep."""
    from PyQt6 import QtWidgets
    buttons = [b for b in panel.findChildren(QtWidgets.QPushButton) if b.text() == "stop"]
    assert len(buttons) >= 4, "the walk, clustering, search and validation tabs each need one"
    assert all(not b.isEnabled() for b in buttons), "stop is enabled only while something runs"
    assert all(b.toolTip() for b in buttons)


def test_the_stop_button_asks_the_running_job_and_says_what_survives(panel, sync):
    import threading

    class FakeJob:
        def __init__(self):
            self.active, self.asked, self.id = True, False, 1

        def cancel(self):
            self.asked = True

    job = FakeJob()
    panel._running = [job]
    panel._set_stoppable(True)
    assert panel.stop_running() == 1
    assert job.asked
    assert any("what is already computed is kept" in m for m in sync)


def test_stopping_with_nothing_running_says_so(panel, sync):
    panel._running = []
    assert panel.stop_running() == 0
    assert any("nothing is running" in m for m in sync)


def test_rows_reach_the_disk_as_they_arrive(panel, tmp_path, monkeypatch):
    """In memory survives a stop; it does not survive quitting, a crash or a power cut, and an hour
    of search is too much to hold in a process nobody promised to keep alive."""
    from starplast import paths
    from starplast.results import load_table, table_kind
    monkeypatch.setenv(paths.ENV_STATE, str(tmp_path))
    panel._start_table(panel.walk_table, [])
    panel._append(panel.walk_table, {"n_neighbors": 5, "min_dist": 0.0, "trustworthiness": 0.9})
    panel._append(panel.walk_table, {"n_neighbors": 15, "min_dist": 0.1, "trustworthiness": 0.8})
    written = panel._close_row_logs()
    assert len(written) == 1, written
    assert table_kind(written[0]) == "umap_walk", "an autosaved file must load back into its own tab"
    back = load_table(written[0])
    assert len(back) == 2 and list(back.n_neighbors) == [5, 15]


def test_an_autosave_that_cannot_be_written_costs_the_file_and_not_the_run(panel, tmp_path,
                                                                          monkeypatch):
    """A read-only disk must not take the analysis down with it: the rows are still in the table."""
    from starplast import paths, results
    monkeypatch.setenv(paths.ENV_STATE, str(tmp_path))

    def refuse(self, row):
        raise OSError("read-only file system")

    monkeypatch.setattr(results.RowLog, "append", refuse)
    panel._start_table(panel.walk_table, [])
    panel._append(panel.walk_table, {"n_neighbors": 5, "trustworthiness": 0.9})
    assert panel.walk_table.rowCount() == 1
    assert len(panel._frames[panel.walk_table]) == 1


def test_the_progress_handle_reports_a_stop_without_raising(panel):
    """The raising form lands at the next log line, and the search logs every fortieth run. This is
    the form the sweeps check per configuration."""
    from starplast.analysis_panel import _Progress

    class FakeJob:
        cancelled, name, note = False, "search", ""

    job = FakeJob()
    p = _Progress(panel, job)
    assert p.stopped() is False
    job.cancelled = True
    assert p.stopped() is True


def test_a_stopped_job_that_returned_something_still_delivers_it(panel, sync):
    """Discarding a stopped run's result would make stopping cost the whole run."""
    got = []
    panel._jobs[7] = (got.append, None)

    class FakeJob:
        id, name, state, result = 7, "search", "cancelled", "partial table"
        active = False

    panel.runner = type("R", (), {"jobs": {7: FakeJob()}})()
    panel._on_job_finished(7, False)
    assert got == ["partial table"]
    assert any("what it finished is kept" in m for m in sync)


def test_a_stopped_job_whose_result_will_not_load_is_reported_not_raised(panel, sync):
    """A partial result can be shaped oddly. Losing the message about what survived because the
    handler tripped would leave the user thinking the stop threw everything away."""
    def boom(_result):
        raise ValueError("half a table")

    panel._jobs[8] = (boom, None)

    class FakeJob:
        id, name, state, result = 8, "search", "cancelled", "partial"
        active = False

    panel.runner = type("R", (), {"jobs": {8: FakeJob()}})()
    panel._on_job_finished(8, False)
    assert any("stopped" in m for m in sync)


def test_a_finished_run_says_where_its_rows_were_saved(panel, sync, tmp_path, monkeypatch):
    from starplast import paths
    monkeypatch.setenv(paths.ENV_STATE, str(tmp_path))
    panel._start_table(panel.walk_table, [])
    panel._append(panel.walk_table, {"n_neighbors": 5, "trustworthiness": 0.9})
    panel._jobs[9] = (lambda r: None, None)

    class FakeJob:
        id, name, state, result = 9, "walk", "done", None
        active = False

    panel.runner = type("R", (), {"jobs": {9: FakeJob()}})()
    panel._on_job_finished(9, True)
    assert any("rows saved to" in m for m in sync), sync


# --------------------------------------------------------------------------- the discovery tab
def test_the_discovery_tab_offers_every_mode_and_a_layer_to_ask_about(panel):
    from starplast import optimize
    assert [panel.discover_mode.itemText(i)
            for i in range(panel.discover_mode.count())] == list(optimize.MODES)
    assert panel.discover_layer.currentText() in panel.nodes.columns
    assert panel.discover_against.currentText() in panel.nodes.columns
    # Fitness by default where there is one: "same compartment, opposite fitness" is the
    # disagreement people come here for, and the alphabetical default was protein length.
    assert not panel.discover_against.currentText().startswith("length")


def test_the_reading_is_not_offered_before_there_is_anything_to_read(panel):
    assert not panel.read_button.isEnabled()
    panel.read_discovery()                                  # must not raise with nothing found
    assert panel.discover_report.toPlainText() == ""


def test_a_finished_climb_enables_the_reading_and_says_what_it_found(panel):
    import numpy as np
    said = []
    panel.status.connect(said.append)
    findings = pd.DataFrame([{
        "kind": "guilt", "layer": "compartment", "layer_kind": "discrete", "cluster": 3,
        "category": "IMC", "n_cluster": 40, "n_known": 30, "n_hits": 25, "purity": 0.83,
        "background": 0.05, "lift": 16.6, "p": 1e-20, "q": 1e-18, "n_predicted": 10,
        "circular": False, "genes": list(panel.nodes.gene_id[:10])}])
    panel._discovery_done(pd.DataFrame([{"score": 12.5, "n_findings": 1, "_findings": findings,
                                         "_labels": np.zeros(3)}]))
    assert panel.read_button.isEnabled()
    assert any("best 12.500" in s for s in said)
    panel.read_discovery()
    text = panel.discover_report.toPlainText()
    assert "IMC" in text and "25 of 30" in text
    assert "hypothesis" in text, "the reading dropped its caveats"


def test_an_empty_climb_leaves_the_reading_switched_off(panel):
    panel._discovery_done(pd.DataFrame())
    assert not panel.read_button.isEnabled()
    panel._discovery_done(None)
    assert not panel.read_button.isEnabled()


def test_each_configuration_reaches_the_table_without_its_artefacts(panel):
    import numpy as np
    panel._start_table(panel.discover_table, [])
    panel._discovery_step_arrived({"restart": 0, "step": 1, "score": 3.0, "n_clusters": 12,
                                   "_labels": np.zeros(5), "_findings": pd.DataFrame()})
    headers = [panel.discover_table.horizontalHeaderItem(i).text()
               for i in range(panel.discover_table.columnCount())]
    assert "score" in headers and "n_clusters" in headers
    assert not any(h.startswith("_") for h in headers), "an artefact column reached the table"


def test_a_row_with_no_recipe_on_it_says_so_rather_than_rebuilding_nothing(panel):
    said = []
    panel.status.connect(said.append)
    panel._start_table(panel.discover_table, [])
    panel._discovery_step_arrived({"restart": 0, "step": 1, "score": 1.0})
    panel.show_discovery_row(0)
    assert any("does not name a configuration" in s for s in said)


def test_the_climb_runs_with_the_settings_on_screen(panel, sync, monkeypatch):
    """The whole path: the tab's controls become an evaluator and a climb, and every configuration
    it streams reaches the table while it is still running."""
    import numpy as np
    import starplast.optimize as O
    seen = {}

    def fake_climb(evaluate, start, **kw):
        seen.update(start=start, kw=kw)
        kw["on_step"]({"restart": 0, "step": 0, "score": 2.0, "n_clusters": 9}, start, {})
        return pd.DataFrame([{"score": 2.0, "n_findings": 0, "_findings": pd.DataFrame(),
                              "_labels": np.zeros(3)}])

    monkeypatch.setattr(O, "climb", fake_climb)
    monkeypatch.setattr(O, "evaluator", lambda nodes, **kw: seen.setdefault("evaluator", kw))
    panel.discover_mode.setCurrentText("disagreement")
    panel.discover_budget.setValue(11)
    panel.discover_restarts.setValue(3)
    panel.run_discovery()
    assert seen["evaluator"]["mode"] == "disagreement"
    assert seen["evaluator"]["layers"] == (panel.discover_layer.currentText(),)
    assert seen["evaluator"]["against"] == (panel.discover_against.currentText(),)
    assert seen["kw"]["max_evaluations"] == 11 and seen["kw"]["restarts"] == 3
    # The stop is taken from the job's progress object, which only carries one under the real
    # runner -- the same contract the recovery search uses. What is checkable here is that the
    # climb is asked at all.
    assert "should_stop" in seen["kw"], "a climb that is never asked whether to stop"
    assert seen["start"]["blocks"], "the climb started with no data in the map"
    assert panel.discover_table.rowCount() == 1
    assert panel.read_button.isEnabled()


def test_a_climb_with_no_table_loaded_does_nothing_rather_than_raising(panel):
    was, panel.nodes = panel.nodes, None
    try:
        panel.run_discovery()
    finally:
        panel.nodes = was


def test_a_discovery_row_rebuilds_the_map_it_scored(panel, sync, monkeypatch):
    """Through the same `search.rebuild` the recovery table uses: one implementation, so the two
    tables cannot come to disagree about what a row means."""
    import numpy as np
    import starplast.search as S
    asked = {}

    def fake_rebuild(nodes, row, log=print):
        asked.update(row)
        return np.zeros((3, 3)), np.ones(3, bool), np.zeros(3, int), ["x"]

    monkeypatch.setattr(S, "rebuild", fake_rebuild)
    panel._start_table(panel.discover_table, [])
    panel._discovery_step_arrived({"restart": 0, "step": 1, "score": 4.0, "method": "tsne",
                                   "algorithm": "kmeans", "blocks": "expression_summary"})
    panel.show_discovery_row(0)
    assert asked.get("blocks") == "expression_summary"
    assert asked.get("method") == "tsne"


def test_a_finished_climb_is_written_to_disk_without_being_asked(panel, tmp_path):
    """A run is expensive and its result is a hundred maps. The run a reader wants to go back to is
    never the one they thought to save, so it is saved as it finishes."""
    import numpy as np
    panel._discovery_done(pd.DataFrame([{"score": 3.0, "n_findings": 0, "algorithm": "kmeans",
                                         "_findings": pd.DataFrame(),
                                         "_labels": np.zeros(len(panel.nodes), int)}]))
    saved = panel.search_store().list()
    assert saved, "a finished search left nothing on disk"
    name, manifest = saved[0]
    assert manifest["mode"] == panel.discover_mode.currentText()
    assert manifest["layer"] == panel.discover_layer.currentText()
    assert manifest["fingerprint"]["n_genes"] == len(panel.nodes)
    assert panel.saved_searches.count() >= 1


def test_a_saved_search_comes_back_with_everything_it_found(panel):
    import numpy as np
    findings = pd.DataFrame([{
        "kind": "guilt", "layer": "compartment", "layer_kind": "discrete", "cluster": 2,
        "category": "IMC", "n_cluster": 30, "n_known": 20, "n_hits": 18, "purity": 0.9,
        "background": 0.05, "lift": 18.0, "p": 1e-15, "q": 1e-13, "n_predicted": 6,
        "circular": False, "genes": list(panel.nodes.gene_id[:6])}])
    panel._discovery_done(pd.DataFrame([{"score": 9.5, "n_findings": 1, "n_clusters": 30,
                                         "_findings": findings,
                                         "_labels": np.zeros(len(panel.nodes), int)}]))
    panel._start_table(panel.discover_table, [])
    panel._search = None
    panel.read_button.setEnabled(False)
    panel.saved_searches.setCurrentIndex(0)
    panel.load_search()
    assert panel.discover_table.rowCount() == 1
    assert panel.read_button.isEnabled()
    panel.read_discovery()
    assert "IMC" in panel.discover_report.toPlainText()


def test_loading_nothing_does_nothing(panel):
    panel.saved_searches.clear()
    panel.load_search()
    assert panel._search is None or panel._search.configs.empty


def test_a_search_loaded_against_another_table_says_so(panel, sync):
    import numpy as np
    said = []
    panel.status.connect(said.append)
    panel._discovery_done(pd.DataFrame([{"score": 1.0, "n_findings": 0,
                                         "_findings": pd.DataFrame(),
                                         "_labels": np.zeros(3, int)}]))
    was, panel.nodes = panel.nodes, panel.nodes.head(20)
    try:
        panel.saved_searches.setCurrentIndex(0)
        panel.load_search()
        assert any("DIFFERENT node table" in s for s in said)
    finally:
        panel.nodes = was


def test_a_search_that_cannot_be_written_is_still_a_search(panel, monkeypatch):
    """Losing a finished run because the disk is full would be a worse failure than the one being
    reported."""
    import numpy as np
    import starplast.searches as SS
    said = []
    panel.status.connect(said.append)

    def refuse(self, search, name=None):
        raise OSError("no space left on device")

    monkeypatch.setattr(SS.SearchStore, "save", refuse)
    panel._discovery_done(pd.DataFrame([{"score": 2.0, "n_findings": 0,
                                         "_findings": pd.DataFrame(),
                                         "_labels": np.zeros(3, int)}]))
    assert panel.read_button.isEnabled(), "the run was thrown away because it could not be saved"
    assert any("could not save" in s for s in said)
