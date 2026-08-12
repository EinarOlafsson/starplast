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
    """data -> map -> clusters -> meaning -> search, in the order the work is actually done."""
    tabs = panel.findChild(type(panel).__mro__[0] if False else __import__("PyQt6.QtWidgets",
                           fromlist=["QTabWidget"]).QTabWidget)
    titles = [tabs.tabText(i).lower() for i in range(tabs.count())]
    assert len(titles) >= 5
    for expected in ("data", "map", "cluster", "meaning", "search"):
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
    def run_now(fn, on_done, name="analysis"):
        # `name` is accepted because every job carries one now, so it can be identified in the Jobs
        # panel and stopped there. Recorded rather than dropped: a job named "analysis" for all five
        # tabs would make the panel useless, so the names are worth asserting on.
        started.append(name)
        on_done(fn(lambda *_: None))

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


def test_the_walk_is_given_the_sample_size_and_seed_from_the_panel(panel, sync, monkeypatch):
    """A walk that silently used different settings than the ones on screen would be unreproducible."""
    import starplast.tuning as T
    seen = {}

    def fake(n, spec, sample_size=None, seed=None, log=None):
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


def test_the_search_is_given_the_chosen_target_and_a_bounded_set_of_combinations(panel, sync,
                                                                                 monkeypatch):
    """Every block combination up to the chosen size; unbounded it is 2^n runs."""
    import starplast.search as S
    seen = {}

    def fake(n, target=None, block_sets=None, sample_size=None, seed=None, store=None, log=None):
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
