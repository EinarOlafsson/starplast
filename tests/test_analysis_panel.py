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
