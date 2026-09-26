"""The Strategies tab: the list, the guide, the settings, the jobs and the results.

Built on the planted organism, so every strategy can actually run here and the tests check that
what a person presses reaches the strategy and what it returns reaches the screen -- not the
science, which `test_strategies` checks.
"""
from __future__ import annotations

import json
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from starplast import strategies as S  # noqa: E402


@pytest.fixture(scope="module")
def app():
    from PyQt6 import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture(scope="module")
def planted():
    return S.planted_context()


@pytest.fixture
def panel(app, planted):
    from starplast.strategy_panel import StrategyPanel
    p = StrategyPanel(planted.nodes, graph=planted.graph, other=planted.other(), organism="Tg")
    yield p
    p.deleteLater()


def _leaves(panel):
    out = []
    for i in range(panel.tree.topLevelItemCount()):
        fam = panel.tree.topLevelItem(i)
        out += [fam.child(j) for j in range(fam.childCount())]
    return out


# --------------------------------------------------------------------------- the list and guide
def test_every_strategy_is_listed_under_its_family_with_a_tooltip(panel):
    leaves = _leaves(panel)
    assert len(leaves) == len(S.catalog()) >= 30
    assert panel.tree.topLevelItemCount() == len(S.families())
    for leaf in leaves:
        assert leaf.toolTip(0) and leaf.toolTip(1)


def test_the_first_strategy_is_selected_and_explained(panel):
    assert panel.current.key == "holdout_search"
    text = panel.guide.toPlainText()
    for part in ("Walkthrough", "How it is tested", "Settings", "Hold out a category"):
        assert part in text


def test_selecting_a_strategy_rebuilds_its_form(panel):
    panel.select("feature_knn")
    assert panel.form.rowCount() == len(S.get("feature_knn").params)
    assert set(panel.settings()) == {p.name for p in S.get("feature_knn").params}
    panel.select("seed_expansion")
    assert "genes" in panel.settings()


def test_the_filter_hides_what_does_not_match(panel):
    panel.filter.setText("gene list")
    shown = [leaf for leaf in _leaves(panel) if not leaf.isHidden()]
    assert shown and all("list" in (S.get(l.data(0, 256)).title + S.get(l.data(0, 256)).question
                                     + S.get(l.data(0, 256)).tooltip).lower() for l in shown)
    panel.filter.setText("zzzz-no-match")
    assert all(panel.tree.topLevelItem(i).isHidden() for i in range(panel.tree.topLevelItemCount()))
    panel.filter.setText("")
    assert not any(leaf.isHidden() for leaf in _leaves(panel))


def test_every_control_in_the_panel_explains_itself(panel):
    """The application-wide rule, applied to every strategy's form in turn."""
    from PyQt6 import QtWidgets
    missing = []
    for s in S.catalog():
        panel.select(s.key)
        for cls in (QtWidgets.QComboBox, QtWidgets.QSpinBox, QtWidgets.QDoubleSpinBox,
                    QtWidgets.QPushButton, QtWidgets.QLineEdit, QtWidgets.QPlainTextEdit):
            for w in panel.findChildren(cls):
                # A spin box's and a combo box's own text fields are parts of those controls,
                # which carry the tooltip; they are not controls of their own.
                if isinstance(w.parent(), (QtWidgets.QAbstractSpinBox, QtWidgets.QComboBox)):
                    continue
                if not w.toolTip().strip():
                    missing.append(f"{s.key}: {cls.__name__} "
                                   f"{getattr(w, 'text', lambda: '')()!r}")
    assert not missing, missing


def test_optional_parameters_offer_none(panel):
    panel.select("geneset_hunt")
    _p, combo = panel.inputs["exclude"]
    assert combo.itemText(0) == "(none)" and panel.settings()["exclude"] is None


def test_a_column_from_the_other_organism_is_offered_there(panel):
    panel.select("ortholog_transfer")
    _p, combo = panel.inputs["source"]
    assert "piggybac_mis" in [combo.itemText(i) for i in range(combo.count())]


def test_settings_can_be_set_and_an_impossible_one_is_refused(panel):
    panel.select("feature_knn")
    panel.set_setting("k", 9)
    panel.set_setting("min_share", 0.4)
    panel.set_setting("target", "cellcycle_phase")
    assert panel.settings() == {"target": "cellcycle_phase", "k": 9, "min_share": 0.4}
    with pytest.raises(ValueError):
        panel.set_setting("target", "not_a_column")
    panel.select("holdout_search")
    panel.set_setting("n_neighbors", "10, 20")
    assert panel.settings()["n_neighbors"] == "10, 20"
    panel.select("seed_expansion")
    panel.set_setting("genes", ["A", "B"])
    assert panel.settings()["genes"] == "A\nB"


# --------------------------------------------------------------------------- gene lists
def test_a_gene_list_can_be_loaded_from_a_file_an_example_or_the_gate(panel, planted, tmp_path):
    panel.select("seed_expansion")
    _p, host = panel.inputs["genes"]
    path = tmp_path / "list.csv"
    path.write_text("gene,score\n" + "\n".join(f"{g},1" for g in planted.gene_ids[:5]) + "\n")
    assert panel.load_genes_file(host.box, str(path)) == 6
    assert planted.gene_ids[0] in host.box.toPlainText()
    assert panel.fill_example(host.box) >= 30
    panel.gated = lambda: list(planted.gene_ids[:3])
    assert panel.fill_from_gate(host.box) == 3
    panel.gated = ["x"]
    assert panel.fill_from_gate(host.box) == 1


def test_a_cancelled_file_dialog_changes_nothing(panel, monkeypatch):
    from PyQt6 import QtWidgets
    monkeypatch.setattr(QtWidgets.QFileDialog, "getOpenFileName", lambda *a, **k: ("", ""))
    panel.select("seed_expansion")
    _p, host = panel.inputs["genes"]
    host.box.setPlainText("kept")
    assert panel.load_genes_file(host.box) == 0 and host.box.toPlainText() == "kept"


def test_an_example_with_no_usable_category_says_so(app):
    from starplast.strategy_panel import StrategyPanel
    tiny = pd.DataFrame({"gene_id": list("abcd"), "x": [1.0, 2.0, 3.0, 4.0]})
    p = StrategyPanel(tiny, graph={}, other=tiny, organism="Tg")
    said = []
    p.status.connect(said.append)
    from PyQt6 import QtWidgets
    assert p.fill_example(QtWidgets.QPlainTextEdit()) == 0
    assert said and "no category" in said[-1]
    p.deleteLater()


# --------------------------------------------------------------------------- running
def test_a_run_puts_its_summary_tables_and_map_on_screen(panel):
    got = {}
    panel.result_ready.connect(lambda r: got.setdefault("result", r))
    panel.embedding_ready.connect(lambda c, m: got.setdefault("map", (c, m)))
    panel.clusters_ready.connect(lambda l: got.setdefault("clusters", l))
    panel.select("map_neighbours")
    result = panel.run_current()
    assert got["result"] is result and panel.summary.text() == result.summary
    assert panel.result_tabs.count() >= 1 and panel.map_btn.isEnabled()
    assert panel.show_on_map()
    coords, mask = got["map"]
    assert mask.sum() == len(coords) and len(got["clusters"]) == panel.ctx.n
    assert len(panel.frame(0)) == len(result.tables["calls"])


def test_a_clicked_row_selects_its_gene(panel):
    picked = []
    panel.gene_selected.connect(picked.append)
    panel.select("feature_knn")
    panel.run_current()
    panel._row(panel.tables[0], 0)
    assert picked and picked[0] in set(panel.ctx.gene_ids)
    panel.select("link_prediction")
    panel.run_current()
    panel._row(panel.tables[0], 0)
    assert len(picked) == 2
    from PyQt6 import QtWidgets
    empty = QtWidgets.QTableWidget(1, 1)
    panel._row(empty, 0)                   # a table with no gene column does nothing
    assert len(picked) == 2


def test_a_self_test_shows_its_verdict_once(panel):
    got = []
    panel.test_ready.connect(got.append)
    panel.select("physical_partners")
    test = panel.test_current()
    assert got == [test] and test.verdict == "PASS"
    assert panel.verdict.text().count("PASS") == 1
    assert "Hidden:" in panel.verdict.text()


def test_a_run_that_fails_says_why_instead_of_crashing(panel):
    panel.select("layer_propagation")
    panel.set_setting("layer", "compartment")          # built from the label: refused
    assert panel.run_current() is None
    assert "Could not run" in panel.summary.text() and "built from" in panel.summary.text()


def test_showing_a_map_before_anything_ran_does_nothing(panel):
    panel.last_result = None
    assert panel.show_on_map() is False
    panel.current = None
    assert panel._submit("run") is None


def test_jobs_run_on_the_runner_and_can_be_stopped(app, planted):
    from PyQt6 import QtWidgets
    from starplast.jobs import JobRunner
    from starplast.strategy_panel import StrategyPanel
    runner = JobRunner()
    p = StrategyPanel(planted.nodes, runner=runner, graph=planted.graph, other=planted.other(),
                      organism="Tg")
    got = []
    p.result_ready.connect(got.append)
    p.select("feature_knn")
    job = p.run_current()
    assert job is not None and p.stop_btn.isEnabled()
    runner.wait(120000)
    for _ in range(50):
        QtWidgets.QApplication.processEvents()
        if got:
            break
    assert got and got[0].strategy == "feature_knn" and not p.stop_btn.isEnabled()
    # A job asked to stop before it starts is reported as stopped, not delivered.
    p.select("masked_imputation")
    job = p.run_current()
    assert p.stop_running() >= 0
    runner.wait(120000)
    for _ in range(50):
        QtWidgets.QApplication.processEvents()
    # A failing job reaches the summary rather than a traceback.
    p.select("layer_propagation")
    p.set_setting("layer", "compartment")
    p.run_current()
    runner.wait(120000)
    for _ in range(50):
        QtWidgets.QApplication.processEvents()
    assert "Could not run" in p.summary.text()
    p._on_job_finished(10 ** 6, True)                    # someone else's job is ignored
    p.deleteLater()


def test_progress_reports_and_honours_a_stop(app, panel):
    from starplast.jobs import Job, Stopped
    from starplast.strategy_panel import _Progress
    job = Job(id=1, name="x")
    said = []
    panel.status.connect(said.append)
    progress = _Progress(panel, job)
    progress("working")
    assert job.note == "working" and said[-1] == "working" and not progress.stopped()
    job.cancel()
    with pytest.raises(Stopped):
        progress("more")


# --------------------------------------------------------------------------- with the analysis panel
def test_results_tables_are_wired_through_the_analysis_panel(app, planted, tmp_path):
    from starplast.analysis_panel import AnalysisPanel
    from starplast.strategy_panel import StrategyPanel
    from starplast.tuning import EmbeddingStore
    analysis = AnalysisPanel(planted.nodes, store=EmbeddingStore(str(tmp_path)))
    p = StrategyPanel(planted.nodes, analysis=analysis, graph=planted.graph,
                      other=planted.other(), organism="Tg")
    assert "strategy result 1" in analysis.results_tables()
    p.select("feature_knn")
    result = p.run_current()
    assert len(p.frame(0)) == len(result.tables["calls"])
    assert len(analysis._frames[p.tables[0]]) == len(result.tables["calls"])
    p.deleteLater()
    analysis.deleteLater()


# --------------------------------------------------------------------------- measured verdicts
def test_measured_verdicts_are_read_and_shown(app, planted, tmp_path, monkeypatch):
    from starplast import strategy_panel as SP
    path = tmp_path / "m.json"
    path.write_text(json.dumps({"Tg": {
        "feature_knn": {"verdict": "PASS", "metric": "m", "observed": 0.5, "null_mean": 0.1,
                        "null_kind": "n", "n_hidden": 100},
        "layer_vote": {"verdict": "INCONCLUSIVE", "note": "too few"}}}))
    monkeypatch.setattr(SP, "MEASURED", str(path))
    # No calibration, so the list falls back to the single self-test verdict it is testing here.
    monkeypatch.setattr(SP.CAL, "load", lambda p=None: {"meta": {}, "organisms": {}})
    p = SP.StrategyPanel(planted.nodes, graph=planted.graph, other=planted.other(), organism="Tg")
    assert p.items["feature_knn"].text(1) == "PASS"
    assert "INCONCLUSIVE" in p._measured_line("layer_vote")
    p.select("feature_knn")
    assert "PASS on the shipped data" in p.guide.toPlainText()
    monkeypatch.setattr(SP, "MEASURED", str(tmp_path / "absent.json"))
    assert SP.load_measured("Tg") == {}
    p.deleteLater()


def test_the_shipped_verdicts_cover_every_strategy():
    from starplast.strategy_panel import MEASURED
    with open(MEASURED) as fh:
        data = json.load(fh)
    for org in ("Tg", "Pf"):
        assert set(data[org]) == {s.key for s in S.catalog()}, org
        for r in data[org].values():
            assert r["verdict"] in ("PASS", "FAIL", "INCONCLUSIVE", "NOT RUN")


# --------------------------------------------------------------------------- in the window
def test_the_window_docks_strategies_beside_evidence_and_analysis(app):
    from starplast.app import Window
    w = Window()
    try:
        tabbed = {d.windowTitle() for d in w.tabifiedDockWidgets(w.analysis_dock)}
        assert {"evidence", "strategies"} <= tabbed
        assert w.strategy_panel.ctx.organism == "Tg"
        w.gated = np.array([0, 1])
        assert len(w.strategy_panel.gated()) == 2
        w.gated = None
        assert w.strategy_panel.gated() == []
        picked = []
        w.search.textChanged.connect(picked.append)
        w.strategy_panel.gene_selected.emit(str(w.nodes["gene_id"].iloc[0]))
        assert picked
    finally:
        w.close()
        w.deleteLater()


def test_the_window_opens_without_the_strategies_panel(app, monkeypatch):
    """A panel that cannot be imported costs its tab, never the window."""
    from starplast.app import Window
    monkeypatch.setitem(sys.modules, "starplast.strategy_panel", None)
    w = Window()
    try:
        from PyQt6 import QtWidgets
        assert not hasattr(w, "strategy_panel")
        titles = {d.windowTitle() for d in w.findChildren(QtWidgets.QDockWidget)}
        assert "analysis" in titles and "strategies" not in titles
    finally:
        w.close()
        w.deleteLater()


def test_a_calibrated_strategy_shows_its_grade_and_offers_its_tuned_setting(app, planted, tmp_path,
                                                                         monkeypatch):
    """Where the sweep measured a strategy, the list shows the grade rather than the one verdict --
    a grade rests on hundreds of held-out tests -- and the tuned setting is one click away."""
    from starplast import calibration as C
    from starplast import strategy_panel as SP
    rows = []
    for target in ("a", "b"):
        for seed in (1, 2, 3, 4, 5):
            rows.append({"organism": "Tg", "strategy": "feature_knn",
                         "settings": {"target": target, "k": 5}, "seed": seed, "verdict": "PASS",
                         "metric": "correct calls", "observed": 0.8, "null_mean": 0.2,
                         "note": "", "wall_seconds": 1.0})
    path = C.write(C.summarise(C.runs_frame(rows)), str(tmp_path / "cal.json"),
                   meta={"date": "2026-09-26", "runs": len(rows)})
    # Every reader goes through `load`, so pointing it at the planted file isolates all of them.
    planted_calibration = C.load(path)
    monkeypatch.setattr(C, "load", lambda p=None: planted_calibration)
    panel = SP.StrategyPanel(planted.nodes, graph=planted.graph, other=planted.other(),
                             organism="Tg")
    assert panel.items["feature_knn"].text(1) == "reliable"
    panel.select("feature_knn")
    assert "RELIABLE" in panel.guide.toPlainText()
    assert panel.tuned_btn.isEnabled()
    applied = panel.use_tuned_settings()
    assert applied == {"k": 5} and panel.settings()["k"] == 5
    # A strategy the sweep never measured offers nothing to apply.
    panel.select("layer_vote")
    assert not panel.tuned_btn.isEnabled()
    assert panel.use_tuned_settings() == {}
    panel.deleteLater()
