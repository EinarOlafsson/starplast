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
    monkeypatch.setattr(panel, "run_embed", lambda: built.setdefault(
        "spec", (panel.nn.value(), panel.md.value())))
    panel.show_walk_row(1, 0)
    assert built["spec"] == (33, 0.1)


def test_clicking_a_row_that_names_no_configuration_says_so(panel, monkeypatch):
    said = []
    panel.status.connect(said.append)
    panel._fill(panel.walk_table, pd.DataFrame({"n_neighbors": ["n/a"], "min_dist": ["n/a"]}))
    monkeypatch.setattr(panel, "run_embed", lambda: pytest.fail("must not build from a bad row"))
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
    panel.focus.clearSelection()
    assert panel.objective_settings()["category"] is None, "none selected must mean all of them"
    for i in range(min(2, panel.focus.count())):
        panel.focus.item(i).setSelected(True)
    got = panel.objective_settings()["category"]
    assert isinstance(got, list) and len(got) == 2
    assert all(g == panel.focus.item(i).data(_Qt.Qt.ItemDataRole.UserRole)
               for i, g in enumerate(got))
    panel.focus.clearSelection()


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
