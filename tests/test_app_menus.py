"""The reworked control surface: menus, exports, the generic category filter, and the busy spinner.

Everything here used to be a widget stacked in the left panel, which put most of the application's
settings in a column beside the view. These tests hold the behaviour that moved, not the layout: that
a menu action still changes what is drawn, that an export writes what is on screen, and that spinning
to show work in progress puts the camera back where it found it.
"""
import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")

from PyQt6 import QtCore, QtWidgets  # noqa: E402

import starplast.app as A  # noqa: E402
from starplast import app as APP  # noqa: E402
from starplast.app import category_columns as category_columns_for  # noqa: E402


@pytest.fixture(scope="module")
def win():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    w = A.Window()
    yield w
    w.console.remove()


# --------------------------------------------------------------------------- menus
def test_every_setting_is_reachable_from_the_menu_bar(win):
    """The point of the rework: no setting is only reachable from a panel."""
    titles = [a.menu().title() for a in win.menuBar().actions() if a.menu()]
    assert titles == ["&File", "&View", "&Edges", "&Tools", "&Help"]


def test_the_edge_menu_holds_every_type_plus_the_two_global_switches(win):
    """Items 13 and 15: the twelve types and the two things that govern all of them, together."""
    assert len(win.edge_act) == len(A.EDGE_TYPES)
    labels = [a.text() for a in win.edge_menu.actions() if not a.isSeparator()]
    assert "Draw all active edges" in labels
    assert "Attention-corrected co-mention" in labels


def test_toggling_an_edge_type_changes_what_is_drawn(win):
    """The complaint was that these appeared to do nothing -- because picking was broken, so no gene
    could be selected, and unselected the layers draw nothing unless draw-all is on."""
    win.reset()
    win.all_edges_act.setChecked(True)
    for k in win.edge_act:
        win.set_edge(k, False)
    win.redraw()
    assert win.edge_items == []
    key = next(k for k, _ in A.EDGE_TYPES if k in win.edges)
    win.set_edge(key, True)
    assert win.edge_items, "an enabled edge type drew nothing with draw-all on"
    win.all_edges_act.setChecked(False)
    for k in win.edge_act:
        win.set_edge(k, False)


def test_selecting_a_gene_draws_its_edges_without_draw_all(win):
    """The behaviour the user could never reach while picking was broken."""
    win.all_edges_act.setChecked(False)
    key = next(k for k, _ in A.EDGE_TYPES if k in win.edges and len(win.edges[k]["a"]))
    for k in win.edge_act:
        win.set_edge(k, k == key)
    win.sel = int(win.edges[key]["a"][0])
    win.redraw()
    assert win.edge_items, "a selected gene drew none of its own edges"
    win.sel = None
    win.set_edge(key, False)


def test_the_status_line_says_how_to_see_edges_when_none_are_drawn(win):
    """Silence was the whole problem: nothing drawn and nothing said."""
    win.reset()
    for k in win.edge_act:
        win.set_edge(k, False)
    win.redraw()
    assert "Draw all active edges" in win.status.currentMessage()




def test_point_size_changes_the_drawn_size(win):
    """Item 10. Automatic follows the point style; the rest are absolute."""
    win.set_point_size(None)
    auto = float(np.max(win.scatter.size))
    win.set_point_size(16.0)
    big = float(np.max(win.scatter.size))
    assert big > auto
    win.set_point_size(2.0)
    assert float(np.max(win.scatter.size)) < big
    win.set_point_size(None)


def test_the_explanations_are_offered_and_are_not_jargon(win):
    """Item 12: "never merged" read as jargon, so the menu now explains it in words."""
    assert "borrow the credibility" in A.EDGE_EXPLANATION
    assert "crosslink" in A.EDGE_EXPLANATION.lower()
    assert "negative control" in A.MAP_EXPLANATION.lower()


# --------------------------------------------------------------------------- category filter
def test_any_categorical_column_can_drive_the_filter(win):
    """Item 9. Compartment was hardcoded, and it is the worst-recovered property in the map."""
    assert len(win.categories) > 5
    assert "compartment" in win.categories and "cellcycle_phase" in win.categories


def test_categories_are_found_whatever_dtype_pandas_uses_for_strings():
    """pandas 3 gives string columns dtype `str` where pandas 2 gave `object`.

    Testing for `object` found nothing but a boolean column, so every filter category vanished under
    a newer pandas while the suite went on passing against the older one -- the bug was invisible in
    the test environment and total in the user's. Built explicitly here rather than read off the
    shipped table, so this holds regardless of which pandas is installed.
    """
    import pandas as pd
    d = pd.DataFrame({
        "gene_id": [f"g{i}" for i in range(6)],
        "compartment": pd.Series(["a", "b", "a", "b", "c", "c"], dtype="str"),
        "phase": pd.Series(["x", "y", "x", "y", "x", "y"], dtype="object"),
        "grade": pd.Series(["lo", "hi", "lo", "hi", "lo", "hi"], dtype="category"),
        "flag": pd.Series([True, False, True, False, True, False]),
        "count": pd.Series([1, 2, 3, 4, 5, 6]),           # a quantity, not a class
        "score": pd.Series([0.5, 1.5, 2.5, 3.5, 4.5, 5.5]),
        "when": pd.to_datetime(["2020-01-01"] * 6),
    })
    got = category_columns_for(d)
    assert {"compartment", "phase", "grade", "flag"} <= set(got)
    assert "count" not in got and "score" not in got, "a numeric quantity is not a class"
    assert "when" not in got and "gene_id" not in got


def test_missing_values_become_absence_not_the_string_nan_or_a_float():
    """The second pandas-3 failure, and the one that actually crashed the application.

    Under pandas 3, `.astype(str)` on a string column keeps NA, so the result mixes `str` with `nan`:
    sorting raises TypeError and the NaNs compare unequal to everything, including themselves. Under
    pandas 2 the same cast produced the literal "nan". Missing becomes "", which is already an
    absence label, so the value lists and the held-out search agree about what missing means.
    """
    import pandas as pd
    for dtype in ("str", "object"):
        s = pd.Series(["a", None, "b"], dtype=dtype)
        out = A.as_text(s)
        assert list(out) == ["a", "", "b"]
        sorted(out.unique())                       # must not raise
        assert "" in set(out) and "nan" not in set(out)


def test_an_identifier_column_is_not_offered_as_a_category(win):
    """orthogroup has 7,331 values. A list that long is not a filter."""
    assert "orthogroup" not in win.categories
    assert "gene_id" not in win.categories
    assert "product" not in win.categories


def test_switching_category_rebuilds_the_list_and_the_filter(win):
    win.on_category_changed("cellcycle_phase")
    assert win.category == "cellcycle_phase"
    # Through the application's own text helper, not `.astype(str)`: under pandas 3 that cast keeps
    # NA, so the count excluded the unlabelled genes the list shows as "" and this read 5 against 6.
    # `as_text` is what the window itself uses, and a test that asks the question a different way is
    # testing a different question.
    assert win.comp_list.count() == APP.as_text(win.nodes.cellcycle_phase).nunique()
    win.comp_list.item(0).setCheckState(QtCore.Qt.CheckState.Checked)
    vis = win.visible_mask()
    assert 0 < vis.sum() < win.n, "the filter is not restricting anything"
    win.comp_list.set_checked([])
    win.on_category_changed("compartment")


def test_absence_values_sink_to_the_bottom_of_the_list(win):
    """"unassigned" is the largest class in this proteome; sorted by size it would head every list."""
    win.on_category_changed("compartment")
    labels = [win.comp_list.item(i).data(QtCore.Qt.ItemDataRole.UserRole)
              for i in range(win.comp_list.count())]
    absent = [i for i, v in enumerate(labels) if str(v).lower() in A.ABSENCE]
    if absent:
        assert min(absent) == len(labels) - len(absent), "an absence value is not at the end"


def test_absence_values_are_marked_as_absence_not_as_a_class(win):
    win.on_category_changed("compartment")
    for i in range(win.comp_list.count()):
        it = win.comp_list.item(i)
        if str(it.data(QtCore.Qt.ItemDataRole.UserRole)).lower() in A.ABSENCE:
            assert "not the same as measuring zero" in it.toolTip()


# --------------------------------------------------------------------------- galaxy tier






# --------------------------------------------------------------------------- exports
def test_export_image_writes_a_png(win, tmp_path):
    """Skipped where there is no real GL context.

    Under the offscreen platform pyqtgraph warns that QOpenGLWidget is unsupported, and grabbing a
    framebuffer there is undefined -- it segfaulted the interpreter rather than failing. The export
    now refuses instead of crashing, and this checks whichever of the two applies.
    """
    p = tmp_path / "shot.png"
    got = win.export_image(str(p))
    if got is None:
        assert "OpenGL context" in win.status.currentMessage()
        pytest.skip("no usable GL context on this platform")
    assert got == str(p)
    assert p.exists() and p.stat().st_size > 0


def test_export_visible_writes_only_the_visible_genes_with_coordinates(win, tmp_path):
    """A gene list without coordinates cannot reproduce what was on screen."""
    import pandas as pd
    win.on_category_changed("compartment")
    win.comp_list.clearSelection()
    win.comp_list.item(0).setSelected(True)
    p = tmp_path / "genes.csv"
    # Columns passed explicitly: left to None the picker opens, and exec() blocks a test forever.
    win.export_visible(str(p), columns=["product", "compartment"])
    d = pd.read_csv(p)
    assert len(d) == int(win.visible_mask().sum())
    assert {"gene_id", "x", "y", "z"} <= set(d.columns)
    win.comp_list.clearSelection()


def test_export_graphml_keeps_the_edge_type_on_every_edge(win, tmp_path):
    """A merged edge in an exported file outlives the session, which is where it does most damage."""
    key = next(k for k, _ in A.EDGE_TYPES if k in win.edges and len(win.edges[k]["a"]))
    for k in win.edge_act:
        win.set_edge(k, k == key)
    p = tmp_path / "g.graphml"
    win.export_graphml(str(p))
    text = p.read_text()
    assert text.startswith("<?xml")
    assert '<data key="etype">' in text
    assert f">{key}<" in text
    import xml.etree.ElementTree as ET
    ET.fromstring(text)              # must be well-formed, not merely written
    win.set_edge(key, False)


def test_relationships_export_carries_both_endpoints_and_the_edge_type(win, tmp_path):
    """"Export the relationships I can see" means a spreadsheet, not a GraphML for Cytoscape.

    Both endpoints carry the chosen columns so the file reads without joining it back to anything,
    and the edge type is never collapsed: a row saying two genes are related at 0.8 has lost whether
    they were measured touching or merely mentioned together.
    """
    import pandas as pd
    key = next(k for k, _ in A.EDGE_TYPES if k in win.edges and len(win.edges[k]["a"]))
    for k in win.edge_act:
        win.set_edge(k, k == key)
    win.all_edges_act.setChecked(True)
    p = tmp_path / "edges.csv"
    win.export_relationships(str(p), columns=["product", "compartment"])
    d = pd.read_csv(p)
    assert {"gene_a", "gene_b", "edge_type", "weight", "weight_kind"} <= set(d.columns)
    assert {"product_a", "product_b", "compartment_a", "compartment_b"} <= set(d.columns)
    assert set(d.edge_type) == {key}
    win.all_edges_act.setChecked(False)
    win.set_edge(key, False)


def test_relationships_export_matches_what_is_drawn(win, tmp_path):
    """Exporting a whole layer while the screen shows one gene's neighbourhood would hand back a
    different graph from the one being looked at."""
    import pandas as pd
    key = next(k for k, _ in A.EDGE_TYPES if k in win.edges and len(win.edges[k]["a"]))
    for k in win.edge_act:
        win.set_edge(k, k == key)
    win.all_edges_act.setChecked(False)
    win.sel = int(win.edges[key]["a"][0])
    sel_id = str(win.nodes.gene_id.iloc[win.sel])
    p = tmp_path / "one.csv"
    win.export_relationships(str(p), columns=[])
    one = pd.read_csv(p)
    assert ((one.gene_a == sel_id) | (one.gene_b == sel_id)).all(), \
        "an edge not touching the selected gene was exported"
    win.all_edges_act.setChecked(True)
    p2 = tmp_path / "all.csv"
    win.export_relationships(str(p2), columns=[])
    assert len(pd.read_csv(p2)) > len(one)
    win.sel = None
    win.all_edges_act.setChecked(False)
    win.set_edge(key, False)


def test_exporting_relationships_with_none_visible_says_what_to_do(win, tmp_path):
    """Writing an empty file would look like "there are no relationships", which is a different claim."""
    win.reset()
    for k in win.edge_act:
        win.set_edge(k, False)
    assert win.export_relationships(str(tmp_path / "none.csv"), columns=[]) is None
    assert "no edges visible" in win.status.currentMessage()


def test_the_column_picker_defaults_to_what_is_on_screen(win):
    """The columns worth exporting are usually the ones being looked at."""
    win.on_category_changed("cellcycle_phase")
    win.set_color_mode("publications")
    got = win.default_export_columns()
    assert got[0] == "gene_id"
    assert "cellcycle_phase" in got, "the active filter column is not offered"
    assert "n_publications" in got, "the active coloring is not offered"
    win.on_category_changed("compartment")
    win.set_color_mode("compartment")


def test_the_column_picker_offers_every_column_and_reads_back_the_ticks(win):
    d = win.choose_export_columns(preselect=["gene_id", "product"])
    assert d.column_list.count() == len(win.nodes.columns)
    assert set(win._ticked(d)) == {"gene_id", "product"}


def test_a_cancelled_column_picker_writes_nothing(win, tmp_path, monkeypatch):
    """Cancel must abandon the export, not fall through to writing every column."""
    monkeypatch.setattr(A.Window, "choose_export_columns",
                        lambda self, preselect=None: _RejectedDialog())
    p = tmp_path / "never.csv"
    assert win.export_visible(str(p)) is None
    assert win.export_relationships(str(p)) is None
    assert not p.exists()


class _RejectedDialog:
    """A dialog whose exec() reports Cancel, so the modal loop is never entered."""

    column_list = None

    def exec(self):
        return 0


def test_a_cancelled_save_dialog_writes_nothing(win, monkeypatch):
    monkeypatch.setattr(A.QtWidgets.QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: ("", "")))
    assert win.export_image() is None
    assert win.export_visible() is None
    assert win.export_relationships() is None
    assert win.export_graphml() is None


def test_a_gene_id_containing_markup_is_escaped_in_graphml(win, tmp_path, monkeypatch):
    """An unescaped identifier would produce a file no parser will read."""
    nodes = win.nodes.copy()
    nodes.loc[nodes.index[0], "gene_id"] = "A&B<C>"
    monkeypatch.setattr(win, "nodes", nodes)
    p = tmp_path / "esc.graphml"
    win.export_graphml(str(p))
    import xml.etree.ElementTree as ET
    ET.fromstring(p.read_text())
    assert "A&amp;B&lt;C&gt;" in p.read_text()


# --------------------------------------------------------------------------- jobs and spinning
def test_a_running_job_shows_the_progress_bar(win):
    assert win.progress_bar.isHidden()
    win._on_busy(True)
    assert not win.progress_bar.isHidden()
    win._on_busy(False)
    assert win.progress_bar.isHidden()


def test_spinning_for_a_job_returns_the_camera_to_where_it_started(win):
    """Item 16. A spinner that leaves the map at a random angle has destroyed the view you set up."""
    win.spin_busy_act.setChecked(True)
    win.spin_act.setChecked(False)
    win.view.setCameraPosition(distance=55, elevation=21, azimuth=137)
    before = win._camera_state()
    win._on_busy(True)
    assert win.spin_act.isChecked()
    for _ in range(20):
        win._spin_step()
    assert win._camera_state() != before, "it never actually span"
    win._on_busy(False)
    assert not win.spin_act.isChecked()
    assert win._camera_state() == pytest.approx(before), "the camera was left where the spin stopped"


def test_a_spin_the_user_started_is_not_stopped_by_a_job_ending(win):
    """Only the spin the job started belongs to the job."""
    win.spin_act.setChecked(True)
    win._on_busy(True)
    win._on_busy(False)
    assert win.spin_act.isChecked(), "a job ending stopped a spin the user had asked for"
    win.spin_act.setChecked(False)


def test_the_spinner_can_be_turned_off(win):
    win.spin_busy_act.setChecked(False)
    win.spin_act.setChecked(False)
    win._on_busy(True)
    assert not win.spin_act.isChecked()
    win._on_busy(False)
    win.spin_busy_act.setChecked(True)


def test_a_failed_job_is_listed_with_its_error(win):
    """Item 18, and the reason finished jobs are kept."""
    def boom():
        raise RuntimeError("deliberate")
    job = win.run_job(boom, "failing job")
    win.jobs.wait(4000)
    QtWidgets.QApplication.processEvents()
    win._refresh_jobs()
    rows = [win.jobs_view.topLevelItem(i) for i in range(win.jobs_view.topLevelItemCount())]
    row = next(r for r in rows if r.text(0) == "failing job")
    assert row.text(1) == A.FAILED
    assert "deliberate" in row.text(2)
    assert "RuntimeError" in row.toolTip(2)
    assert job.state == A.FAILED


# --------------------------------------------------------------------------- context menu
def test_the_right_click_menu_carries_the_map_actions(win):
    """Item 14: spin moved off the panel, and item 11: export lives here.

    Built rather than shown: exec() would block on a modal loop until a human closed the menu.
    """
    menu = win.build_context_menu()
    labels = [a.text() for a in menu.actions() if a.text()]
    assert "Export image…" in labels
    assert "Export graph (GraphML)…" in labels
    assert "Reset view / clear filters" in labels
    assert any("Spin" in x for x in labels)
    menu.deleteLater()


def test_right_clicking_the_map_shows_the_menu(win, monkeypatch):
    """The shown path, which build_context_menu exists to keep out of the tests: exec() blocks."""
    shown = []
    monkeypatch.setattr(A.QtWidgets.QMenu, "exec", lambda self, *a: shown.append(self))
    menu = win._context_menu(QtCore.QPoint(3, 4))
    assert shown == [menu]


def test_the_explanations_are_shown_when_asked_for(win):
    """Item 12 again: the wording is only worth having if the menu entry actually shows it.

    Shown in a glass window beside the map rather than a modal box, so each call returns the window
    it put up; what is checked is still that it is on screen and says what it should."""
    shown = []
    for explain in (win.explain_edges, win.explain_map):
        d = explain()
        assert d.isVisible()
        shown.append((d.windowTitle(), d.text))
        d.close()
    assert "never" in shown[0][0].lower() or "separate" in shown[0][0].lower()
    assert "borrow the credibility" in shown[0][1]
    assert "negative control" in shown[1][1].lower()


def test_restoring_a_camera_with_nothing_recorded_does_nothing(win):
    """_on_busy(False) can arrive with no stored angle if spinning was never started for a job."""
    before = win._camera_state()
    win._restore_camera(None)
    assert win._camera_state() == before


def test_the_column_filter_hides_columns_without_unticking_them(win):
    """Filtering must not silently discard ticks the user can no longer see."""
    d = win.choose_export_columns(preselect=["gene_id", "product"])
    filt = d.findChild(QtWidgets.QLineEdit)
    filt.setText("product")
    lst = d.column_list
    visible = [lst.item(i).text() for i in range(lst.count()) if not lst.item(i).isHidden()]
    assert visible and all("product" in v for v in visible)
    # gene_id is hidden by the filter but must still be ticked underneath.
    assert set(win._ticked(d)) == {"gene_id", "product"}
    filt.setText("")
    assert all(not lst.item(i).isHidden() for i in range(lst.count()))


def test_an_edge_type_with_nothing_visible_is_skipped_not_written_empty(win, tmp_path):
    """Two layers on, one with no surviving edge for the selected gene: the empty one contributes no
    rows rather than a block of nulls."""
    import numpy as np
    import pandas as pd
    present = [k for k, _ in A.EDGE_TYPES if k in win.edges and len(win.edges[k]["a"])]
    if len(present) < 2:
        pytest.skip("this build has fewer than two populated edge types")
    # A gene in the first layer and in none of the others, so the others yield nothing at all.
    a, b = present[0], present[1]
    in_a = set(win.edges[a]["a"].tolist()) | set(win.edges[a]["b"].tolist())
    in_b = set(win.edges[b]["a"].tolist()) | set(win.edges[b]["b"].tolist())
    only_a = sorted(in_a - in_b)
    if not only_a:
        pytest.skip("no gene appears in one layer and not the other")
    win.reset()
    for k in win.edge_act:
        win.set_edge(k, k in (a, b))
    win.all_edges_act.setChecked(False)
    win.sel = int(only_a[0])
    p = tmp_path / "sparse.csv"
    win.export_relationships(str(p), columns=[])
    d = pd.read_csv(p)
    assert set(d.edge_type) == {a}, "the layer with no surviving edge still wrote rows"
    assert not d.empty
    win.sel = None
    for k in win.edge_act:
        win.set_edge(k, False)


class _AcceptedDialog:
    """A dialog whose exec() reports OK, carrying a fixed set of ticks."""

    def __init__(self, columns):
        self._columns = list(columns)

        class _Item:
            def __init__(self, text):
                self._t = text

            def text(self):
                return self._t

            def checkState(self):
                return QtCore.Qt.CheckState.Checked

        class _List:
            def __init__(self, cols):
                self._items = [_Item(c) for c in cols]

            def count(self):
                return len(self._items)

            def item(self, i):
                return self._items[i]

        self.column_list = _List(self._columns)

    def exec(self):
        return 1


def test_an_accepted_column_picker_is_what_gets_written(win, tmp_path, monkeypatch):
    """The path a user actually takes: open the picker, tick, press OK."""
    import pandas as pd
    monkeypatch.setattr(A.Window, "choose_export_columns",
                        lambda self, preselect=None: _AcceptedDialog(["gene_id", "compartment"]))
    p = tmp_path / "picked.csv"
    win.export_visible(str(p))
    d = pd.read_csv(p)
    assert set(d.columns) == {"gene_id", "compartment", "x", "y", "z"}

    key = next(k for k, _ in A.EDGE_TYPES if k in win.edges and len(win.edges[k]["a"]))
    for k in win.edge_act:
        win.set_edge(k, k == key)
    win.all_edges_act.setChecked(True)
    p2 = tmp_path / "picked_edges.csv"
    win.export_relationships(str(p2))
    e = pd.read_csv(p2)
    assert {"compartment_a", "compartment_b"} <= set(e.columns)
    win.all_edges_act.setChecked(False)
    win.set_edge(key, False)


def test_the_assistant_is_told_what_is_on_screen(win):
    win.sel = 0
    state = win.describe_state()
    assert str(win.nodes.gene_id.iloc[0]) in state
    assert "genes visible" in state
    win.sel = None
    assert "No gene is selected" in win.describe_state()


# --------------------------------------------------------------------------- stopping and resources
def test_an_analysis_shares_the_windows_job_runner(win):
    """The bug behind three complaints at once: the analysis panel had a private thread, so its work
    never appeared in the Jobs panel, could not be stopped, and blocked every other tab while a
    search ran -- which reads as four broken tabs."""
    panel = win.analysis_dock.widget()
    assert panel.runner is win.jobs


def test_a_running_job_can_be_stopped_and_is_not_reported_as_a_failure(win):
    """A stop is a decision, not a crash. Reported in red with a traceback, it teaches people to
    ignore the failure list."""
    import time
    from starplast.analysis_panel import _Progress

    class FakePanel:
        class _S:
            def emit(self, *a):
                pass
        status = _S()

    def work(job):
        p = _Progress(FakePanel(), job)
        for i in range(100000):
            p(f"run {i}")
            time.sleep(0.001)
        return "should never finish"

    job = win.jobs.submit(work, "stoppable walk")
    for _ in range(200):                       # wait for it to actually be running
        if job.note:
            break
        time.sleep(0.01)
    job.cancel()
    assert win.jobs.wait(8000)
    QtWidgets.QApplication.processEvents()
    assert job.state == "cancelled"
    assert not job.error, "a deliberate stop was recorded as an error"
    assert not job.traceback
    assert job.result is None


def test_stop_all_stops_everything_running(win):
    import time
    from starplast.analysis_panel import _Progress

    class FakePanel:
        class _S:
            def emit(self, *a):
                pass
        status = _S()

    def work(job):
        p = _Progress(FakePanel(), job)
        for i in range(100000):
            p(f"run {i}")
            time.sleep(0.001)

    jobs = [win.jobs.submit(work, f"walk {i}") for i in range(2)]
    for _ in range(200):
        if all(j.note or not j.active for j in jobs):
            break
        time.sleep(0.01)
    win.stop_all_jobs()
    assert win.jobs.wait(8000)
    QtWidgets.QApplication.processEvents()
    assert all(j.state == "cancelled" for j in jobs)


def test_stopping_with_nothing_running_says_so(win):
    win.jobs.cancel_all()
    win.stop_all_jobs()
    assert "nothing is running" in win.status.currentMessage()


def test_stopping_with_no_job_selected_says_so(win):
    win.jobs_view.clearSelection()
    win.stop_selected_job()
    assert "select a job" in win.status.currentMessage()


def test_a_finished_job_reports_that_rather_than_being_stopped(win):
    job = win.jobs.submit(lambda: 1, "quick")
    win.jobs.wait(4000)
    QtWidgets.QApplication.processEvents()
    win._refresh_jobs()
    row = next(win.jobs_view.topLevelItem(i) for i in range(win.jobs_view.topLevelItemCount())
               if win.jobs_view.topLevelItem(i).data(0, QtCore.Qt.ItemDataRole.UserRole) == job.id)
    row.setSelected(True)
    win.stop_selected_job()
    assert "already finished" in win.status.currentMessage()


def test_every_job_row_carries_its_id_so_it_can_be_acted_on(win):
    win.jobs.submit(lambda: 1, "identified")
    win.jobs.wait(4000)
    win._refresh_jobs()
    for i in range(win.jobs_view.topLevelItemCount()):
        assert win.jobs_view.topLevelItem(i).data(0, QtCore.Qt.ItemDataRole.UserRole) is not None


def test_the_resource_line_reports_what_it_can_and_never_raises(win):
    s = win.resource_summary()
    assert "CPUs" in s
    assert win.resources.text() or True


def test_freeing_memory_drops_the_rebuildable_cache_and_says_what_it_did(win):
    """"Clear RAM" that silently killed a running search would be a data-loss button wearing a
    housekeeping label, so it only drops what can be rebuilt for free."""
    msg = win.free_memory()
    assert "freed" in msg and "unreachable objects" in msg
    # And a job in flight is untouched.
    job = win.jobs.submit(lambda: "survived", "during cleanup")
    win.free_memory()
    win.jobs.wait(4000)
    assert job.result == "survived"


def test_the_job_context_menu_offers_stop_and_the_traceback(win):
    m = win.build_job_menu()
    labels = [a.text() for a in m.actions() if a.text()]
    assert "Stop this job" in labels
    assert "Stop all running jobs" in labels
    assert any("traceback" in x.lower() for x in labels)


# --------------------------------------------------------------------------- explanations
def test_no_control_anywhere_is_left_without_an_explanation(win):
    """An audit, not a spot check. 85 controls had no tooltip when this was first run, almost all of
    them in the analysis panel, which is the part where choosing wrongly costs the most."""
    missing = []
    # Ok/Cancel inside a QDialogButtonBox are excluded: they are universally understood, and Qt
    # supplies them rather than this application.
    standard = {b for box in win.findChildren(QtWidgets.QDialogButtonBox) for b in box.buttons()}
    for cls in (QtWidgets.QComboBox, QtWidgets.QSpinBox, QtWidgets.QDoubleSpinBox,
                QtWidgets.QCheckBox, QtWidgets.QPushButton):
        for wid in win.findChildren(cls):
            if wid in standard:
                continue
            if not wid.toolTip().strip():
                label = wid.text() if hasattr(wid, "text") else wid.objectName()
                missing.append(f"{cls.__name__}: {label!r}")
    assert not missing, f"controls with no tooltip: {missing}"


def test_the_scoring_explainer_is_reachable_and_carries_the_table(win):
    """The precision/recall table belongs in the application, not only in a commit message."""
    d = win.explain_scoring()
    text = d.findChild(QtWidgets.QPlainTextEdit).toPlainText()
    assert "precision = |c and l| / |c|" in text
    assert "all singletons" in text and "one big cluster" in text
    assert "1.000, 1.000 and\n1.000" in text or "1.000, 1.000" in text
    d.deleteLater()


def test_the_explainer_is_the_one_in_the_code_not_a_copy(win):
    """Two copies of an explanation drift, and the one on screen is the one people believe."""
    from starplast.objectives import EXPLANATION
    d = win.explain_scoring()
    assert d.findChild(QtWidgets.QPlainTextEdit).toPlainText() == EXPLANATION
    d.deleteLater()


def test_help_offers_both_explanations(win):
    helps = [a for a in win.menuBar().actions() if a.menu() and "Help" in a.menu().title()]
    labels = [x.text() for x in helps[0].menu().actions()]
    assert any("does and does not show" in x for x in labels)
    assert any("gamed" in x for x in labels)


# --------------------------------------------------------------------------- coloring by cluster
def test_a_clustering_can_be_seen_on_the_map(win):
    """The Clusters tab computed labels, printed how many there were, and threw them away. There was
    no cluster color mode at all, so the one thing this application is for -- looking at structure
    beside a held-out variable -- could not be done."""
    import numpy as np
    assert "clusters" in A.COLOR_MODES
    lab = np.random.default_rng(0).integers(-1, 5, size=win.n)
    win.use_clusters(lab)
    assert win.color_mode == "clusters", "clustering must show the clusters, not stay on compartment"
    cols = win.colors(np.ones(win.n, bool))
    assert len(set(map(tuple, cols[:, :3].round(4)))) >= 5, "clusters are not drawn apart"


def test_unclustered_genes_are_grey_like_everything_else_unknown(win):
    """HDBSCAN calling a gene unclustered is a finding about that gene, not a gap in the drawing."""
    import numpy as np
    lab = np.random.default_rng(1).integers(-1, 4, size=win.n)
    win.use_clusters(lab)
    cols = win.colors(np.ones(win.n, bool))
    grey = tuple(round(float(x), 4) for x in A.TH.unknown_color(win.theme)[:3])
    assert tuple(round(float(x), 4) for x in cols[lab < 0][0][:3]) == grey
    assert len(set(map(tuple, cols[lab < 0][:, :3]))) == 1, "noise must be one color"


def test_coloring_by_clusters_before_any_exist_is_not_a_crash(win):
    """Selecting the mode from the menu with nothing clustered yet must be grey, not an exception."""
    import numpy as np
    win.cluster_labels = None
    win.set_color_mode("clusters")
    cols = win.colors(np.ones(win.n, bool))
    assert len(set(map(tuple, cols[:, :3]))) == 1
    win.set_color_mode("compartment")


def test_a_stale_clustering_is_not_drawn_against_the_wrong_genes(win):
    """A clustering of a subsample has fewer labels than the map has genes, and lining them up by
    position would color genes by somebody else's cluster."""
    import numpy as np
    win.cluster_labels = np.zeros(7, int)
    win.set_color_mode("clusters")
    cols = win.colors(np.ones(win.n, bool))
    assert len(set(map(tuple, cols[:, :3]))) == 1, "a mismatched clustering must not be drawn"
    win.cluster_labels = None
    win.set_color_mode("compartment")


# --------------------------------------------------------------------------- acting on a job
def test_a_running_job_is_asked_to_stop(win):
    """The other half of the stop button: not "already finished", not "select a job", but a job that
    is actually running being told to stop."""
    import threading
    go = threading.Event()
    job = win.jobs.submit(lambda: go.wait(10), "long one")
    try:
        win._refresh_jobs()
        row = next(win.jobs_view.topLevelItem(i) for i in range(win.jobs_view.topLevelItemCount())
                   if win.jobs_view.topLevelItem(i).data(0, QtCore.Qt.ItemDataRole.UserRole)
                   == job.id)
        win.jobs_view.clearSelection()
        row.setSelected(True)
        win.stop_selected_job()
        assert "asked long one to stop" in win.status.currentMessage()
        assert job.cancelled
    finally:
        go.set()
        win.jobs.wait(4000)


def test_right_clicking_the_job_list_opens_the_menu(win, monkeypatch):
    """`exec` blocks on a real menu, so it is stubbed -- what is checked is that the handler builds a
    menu at the point clicked rather than raising when nothing is selected."""
    monkeypatch.setattr(QtWidgets.QMenu, "exec", lambda self, *a: None)
    win.jobs_view.clearSelection()
    m = win._job_menu(QtCore.QPoint(4, 4))
    assert [a.text() for a in m.actions() if a.text()][:1] == ["Stop this job"]


def test_copying_the_error_of_a_failed_job_puts_the_traceback_on_the_clipboard(win, monkeypatch):
    """A traceback that can only be read off a status bar is a traceback nobody will paste into a
    bug report."""
    held = {}

    class FakeClipboard:
        def setText(self, text):
            held["text"] = text

    def boom():
        raise ValueError("this one failed")

    job = win.jobs.submit(boom, "failing")
    win.jobs.wait(4000)
    win._refresh_jobs()
    row = next(win.jobs_view.topLevelItem(i) for i in range(win.jobs_view.topLevelItemCount())
               if win.jobs_view.topLevelItem(i).data(0, QtCore.Qt.ItemDataRole.UserRole) == job.id)
    win.jobs_view.clearSelection()
    row.setSelected(True)
    monkeypatch.setattr(QtWidgets.QApplication, "clipboard", staticmethod(lambda: FakeClipboard()))
    win.copy_job_error()
    assert "this one failed" in held["text"]
    assert "copied the traceback" in win.status.currentMessage()


def test_copying_an_error_with_no_clipboard_does_not_raise(win, monkeypatch):
    """There is no clipboard on some headless platforms, and a status line must not be able to take
    the application down."""
    monkeypatch.setattr(QtWidgets.QApplication, "clipboard", staticmethod(lambda: None))
    win.copy_job_error()


# --------------------------------------------------------------------------- the resource line
def test_the_resource_line_survives_a_machine_with_no_proc(win, monkeypatch):
    """Read from /proc where it exists, with no hard dependency on it. A status line that raises on a
    machine that reports memory differently is worse than one that says nothing about memory."""
    import builtins
    real = builtins.open

    def no_proc(path, *a, **k):
        if str(path).startswith("/proc"):
            raise OSError("no /proc here")
        return real(path, *a, **k)

    monkeypatch.setattr(builtins, "open", no_proc)
    s = win.resource_summary()
    assert "CPUs" in s and "this process" not in s and "system" not in s


class _FakeCuda:
    def __init__(self, up=True, reserved=2 * 2**30, raises=False):
        self._up, self._reserved, self._raises = up, reserved, raises
        self.emptied = False

    def is_initialized(self):
        if self._raises:
            raise RuntimeError("CUDA is in a bad way")
        return self._up

    def memory_reserved(self):
        return self._reserved

    def empty_cache(self):
        self.emptied = True


class _FakeTorch:
    def __init__(self, **kw):
        self.cuda = _FakeCuda(**kw)


def test_the_gpu_is_reported_only_when_torch_is_already_loaded(win, monkeypatch):
    """Never imported here. Importing torch initialises CUDA, and initialising CUDA inside a running
    OpenGL application segfaulted this program on a three-second timer."""
    import sys
    monkeypatch.setitem(sys.modules, "torch", _FakeTorch())
    assert "GPU 2.0 GB reserved" in win.resource_summary()


def test_a_torch_that_will_not_answer_costs_the_gpu_line_and_nothing_else(win, monkeypatch):
    import sys
    monkeypatch.setitem(sys.modules, "torch", _FakeTorch(raises=True))
    s = win.resource_summary()
    assert "CPUs" in s and "GPU" not in s


def test_freeing_memory_empties_the_gpu_cache_when_there_is_one(win, monkeypatch):
    import sys
    fake = _FakeTorch()
    monkeypatch.setitem(sys.modules, "torch", fake)
    msg = win.free_memory()
    assert fake.cuda.emptied and "GPU cache" in msg


def test_freeing_memory_survives_a_gpu_that_will_not_answer(win, monkeypatch):
    """The same rule as the status line: housekeeping must not be able to crash the session."""
    import sys
    monkeypatch.setitem(sys.modules, "torch", _FakeTorch(raises=True))
    assert "freed" in win.free_memory()


# --------------------------------------------------------------------------- exporting the view
def test_an_export_that_cannot_grab_the_view_says_so_rather_than_crashing(win, monkeypatch,
                                                                         tmp_path):
    """Grabbing a framebuffer without a real GL context is undefined -- it segfaulted rather than
    failing -- so the failure is caught and reported."""
    monkeypatch.setattr(win.view, "isValid", lambda: True)

    def boom():
        raise RuntimeError("no context")

    monkeypatch.setattr(win.view, "grabFramebuffer", boom)
    assert win.export_image(str(tmp_path / "x.png")) is None
    assert "could not capture the view" in win.status.currentMessage()


def test_an_export_with_a_real_frame_writes_it(win, monkeypatch, tmp_path):
    """The success path, on a platform that may have no GL: the frame is supplied, so what is tested
    is that a captured image reaches the disk under the name asked for."""
    from PyQt6 import QtGui
    img = QtGui.QImage(8, 8, QtGui.QImage.Format.Format_ARGB32)
    img.fill(QtGui.QColor("red"))
    monkeypatch.setattr(win.view, "isValid", lambda: True)
    monkeypatch.setattr(win.view, "grabFramebuffer", lambda: img)
    p = tmp_path / "frame.png"
    assert win.export_image(str(p)) == str(p)
    assert p.exists() and p.stat().st_size > 0
    assert f"wrote {p}" in win.status.currentMessage()


def test_an_export_that_cannot_be_written_says_which_file(win, monkeypatch, tmp_path):
    from PyQt6 import QtGui
    img = QtGui.QImage(8, 8, QtGui.QImage.Format.Format_ARGB32)
    img.fill(QtGui.QColor("red"))
    monkeypatch.setattr(win.view, "isValid", lambda: True)
    monkeypatch.setattr(win.view, "grabFramebuffer", lambda: img)
    bad = tmp_path / "no_such_directory" / "frame.png"
    assert win.export_image(str(bad)) is None
    assert "could not write" in win.status.currentMessage()


# --------------------------------------------------------------------------- display in the right-click
def test_the_right_click_menu_carries_every_display_choice(win):
    """The ask: the display settings belong where a reader already is when they want them, which is
    right-clicking the map -- and in Preferences too, not instead of it. Built from one table so the
    two cannot drift, and this asserts the table is what the menu shows."""
    from PyQt6 import QtWidgets
    m = win.build_context_menu()
    display = [a for a in m.actions() if a.text() == "Display"]
    assert display, "the right-click menu has no Display section"
    sub = display[0].menu()
    labels = [a.text() for a in sub.actions() if a.menu()]
    assert labels == [row[0] for row in win.display_choices()]
    for action in sub.actions():
        if action.menu() is None:
            continue
        options = [x.text() for x in action.menu().actions()]
        checked = [x.text() for x in action.menu().actions() if x.isChecked()]
        assert len(checked) == 1, f"{action.text()} has {len(checked)} options checked"
        assert len(options) == len(set(options)), f"{action.text()} repeats an option"


def test_every_display_menu_item_applies_the_option_it_names(win):
    """The classic failure for a menu built in a loop is every action carrying the LAST option,
    because the lambda closed over the loop variable. Driven for real: trigger each action and read
    the setting back.

    Every key it touches is restored afterwards. These settings are GLOBAL -- one QSettings scope for
    the whole run -- so a test that leaves `display/light_target_marker` on whatever it last clicked
    makes the restart test in test_display.py fail, and only when the two run in the same session.
    That is exactly what happened while this test was being written.
    """
    from PyQt6 import QtCore
    keys = ("display/lighting", "display/light_target_marker", "display/light_point_mode",
            "display/light_source")
    s = QtCore.QSettings("starplast", "starplast")
    before = {k: (s.value(k) if s.contains(k) else None) for k in keys}
    m = win.build_context_menu()
    sub = [a for a in m.actions() if a.text() == "Display"][0].menu()
    by_label = {a.text(): a for a in sub.actions() if a.menu()}
    try:
        for label in ("Light render mode", "Target marker", "Point render mode", "Light target"):
            for act in by_label[label].menu().actions():
                act.trigger()
                current = dict((row[0], row[2]) for row in win.display_choices())[label]
                assert str(current) == act.text(), f"{label}: chose {act.text()}, got {current}"
    finally:
        for key, value in before.items():
            s.remove(key) if value is None else s.setValue(key, value)
        s.sync()


def test_the_display_menu_toggles_drive_the_same_setters_as_preferences(win):
    """Two controls for one setting must not be two implementations of it."""
    m = win.build_context_menu()
    sub = [a for a in m.actions() if a.text() == "Display"][0].menu()
    toggles = {a.text(): a for a in sub.actions() if a.isCheckable() and a.menu() is None}
    assert set(toggles) == {label for label, _s, _f in win.display_toggles()}
    was = win.depth_cue
    toggles["Fade with distance"].trigger()
    assert win.depth_cue is not was
    toggles["Fade with distance"].trigger()
    assert win.depth_cue is was


def test_the_display_menu_offers_the_full_dialog_too(win):
    """The menu holds the choices; the sliders live in Preferences, so the menu has to point there."""
    m = win.build_context_menu()
    sub = [a for a in m.actions() if a.text() == "Display"][0].menu()
    assert any(a.text() == "All display settings…" for a in sub.actions())


# --------------------------------------------------------------------------- drag, then tick
def test_a_drag_and_a_right_click_tick_several_categories_at_once(win):
    """The ask: drag across the categories to select them, then right-click to check what is selected.
    Selection is how rows are chosen; ticks are what the map reads."""
    from PyQt6 import QtCore as Q
    win.on_category_changed("compartment")
    cl = win.comp_list
    cl.set_checked([])
    assert cl.checked() == [], "nothing is ticked to begin with, which means everything is shown"
    assert int(win.visible_mask().sum()) == win.n
    for i in range(3):                                   # what a drag across three rows leaves behind
        cl.item(i).setSelected(True)
    menu = cl.build_menu()
    labels = [a.text() for a in menu.actions() if a.text()]
    assert labels[0].startswith("Check selected (3)"), labels
    next(a for a in menu.actions() if a.text().startswith("Check selected")).trigger()
    assert len(cl.checked()) == 3
    assert 0 < int(win.visible_mask().sum()) < win.n, "ticking did not filter the map"
    cl.set_checked([])


def test_an_ordinary_click_no_longer_destroys_a_set_that_was_built_up(win):
    """The reason ticks exist rather than a highlight. With selection AS the choice, four ctrl-clicks
    of work were thrown away by one careless click, and there was no way to keep a set while clicking
    elsewhere to look at something."""
    from PyQt6 import QtCore as Q
    win.on_category_changed("compartment")
    cl = win.comp_list
    wanted = [cl.item(i).data(Q.Qt.ItemDataRole.UserRole) for i in range(3)]
    cl.set_checked(wanted)
    before = int(win.visible_mask().sum())
    cl.clearSelection()
    cl.item(5).setSelected(True)                        # the careless click
    assert cl.checked() == wanted
    assert int(win.visible_mask().sum()) == before
    cl.set_checked([])


def test_the_right_click_menu_offers_only_what_makes_sense(win):
    """"Check selected" with nothing selected is a menu item that cannot do anything."""
    cl = win.comp_list
    cl.clearSelection()
    menu = cl.build_menu()
    by_text = {a.text(): a for a in menu.actions() if a.text()}
    assert not by_text["Check selected (0)"].isEnabled()
    assert not by_text["Uncheck selected (0)"].isEnabled()
    assert by_text["Check all"].isEnabled() and by_text["Clear all"].isEnabled()


def test_check_only_selected_replaces_the_set_rather_than_adding_to_it(win):
    from PyQt6 import QtCore as Q
    win.on_category_changed("compartment")
    cl = win.comp_list
    cl.set_checked([cl.item(0).data(Q.Qt.ItemDataRole.UserRole)])
    cl.clearSelection()
    cl.item(4).setSelected(True)
    next(a for a in cl.build_menu().actions() if a.text() == "Check only selected").trigger()
    assert cl.checked() == [cl.item(4).data(Q.Qt.ItemDataRole.UserRole)]
    cl.set_checked([])


def test_a_bulk_tick_redraws_once_rather_than_once_per_category(win):
    """This list drives a redraw of 8,140 points. Ticking twenty-seven categories should cost one."""
    cl = win.comp_list
    cl.set_checked([])
    fired = []
    cl.checkedChanged.connect(lambda: fired.append(1))
    next(a for a in cl.build_menu().actions() if a.text() == "Check all").trigger()
    assert len(cl.checked()) == cl.count() > 1
    assert len(fired) == 1, f"{cl.count()} rows ticked emitted {len(fired)} signals"
    cl.set_checked([])


def test_space_toggles_whatever_is_selected(win):
    from PyQt6 import QtCore as Q
    from PyQt6.QtGui import QKeyEvent
    cl = win.comp_list
    cl.set_checked([])
    cl.clearSelection()
    for i in range(2):
        cl.item(i).setSelected(True)
    press = QKeyEvent(Q.QEvent.Type.KeyPress, Q.Qt.Key.Key_Space, Q.Qt.KeyboardModifier.NoModifier)
    cl.keyPressEvent(press)
    assert len(cl.checked()) == 2
    cl.keyPressEvent(press)                              # again, and it turns them back off
    assert cl.checked() == []


def test_clicking_the_diagram_ticks_rather_than_only_highlighting(win):
    """Both directions have to be the same operation: clicking an organelle in the drawing must do
    what ticking its row does, or one of the two halves stops filtering."""
    if win.diagram is None:
        pytest.skip("the artwork is not present in this checkout")
    from starplast.celldiagram import COMPARTMENT_SL
    win.on_category_changed("compartment")
    name = next(c for c in COMPARTMENT_SL if c in win.color_of)
    win.comp_list.set_checked([])
    win.select_compartment(name)
    assert win.comp_list.checked() == [name]
    win.comp_list.set_checked([])


def test_reset_clears_the_ticks_and_not_merely_the_highlight(win):
    """A reset that left the map filtered while saying it had reset it would be worse than no button."""
    from PyQt6 import QtCore as Q
    win.on_category_changed("compartment")
    win.comp_list.set_checked([win.comp_list.item(0).data(Q.Qt.ItemDataRole.UserRole)])
    assert int(win.visible_mask().sum()) < win.n
    win.reset()
    assert win.comp_list.checked() == []
    assert int(win.visible_mask().sum()) == win.n


# --------------------------------------------------------------------------- the second species
def test_both_species_can_be_opened(qapp):
    """Instruction 39's arm, made visible. Until this worked, `app.load` opened `nodes.parquet` by
    name and the whole Plasmodium arm -- node table, graph, host bridge, 41 filled slots -- was data
    the browser could not open."""
    from starplast.app import Window, available_species
    assert available_species() == ["Toxoplasma gondii", "Plasmodium falciparum"]
    for name in available_species():
        w = Window(species=name)
        try:
            assert w.species == name
            assert w.n > 1000, f"{name} loaded {w.n} genes"
            assert len(w.xyz) == w.n, "the layout does not cover the table"
            assert w.edges, f"{name} has no edge layers"
            assert name in w.windowTitle()
        finally:
            w.close()


def test_the_two_species_are_different_tables_and_not_one_merged_one(qapp):
    """Nothing is merged: one table per species. A union would put two identifier spaces in one index
    and make "cluster 5 is 71% IMC" a claim about a mixture of organisms."""
    from starplast.app import Window
    tg = Window(species="Toxoplasma gondii")
    pf = Window(species="Plasmodium falciparum")
    try:
        assert set(tg.nodes.gene_id) & set(pf.nodes.gene_id) == set(), "the two tables share genes"
        assert tg.nodes.gene_id.str.startswith("TGME49_").all()
        assert pf.nodes.gene_id.str.startswith("PF3D7_").all()
        assert tg.n != pf.n
    finally:
        tg.close()
        pf.close()


def test_a_layer_without_an_attention_residual_still_loads(qapp):
    """`r` is the attention-corrected residual and only the co-mention layers have one. Indexing it
    unconditionally worked on the Toxoplasma graph, whose builder writes it for every layer, and
    raised KeyError on the first graph built by anything else."""
    from starplast.app import load
    _nodes, _xyz, edges, _models = load("Plasmodium falciparum")
    assert "orthogroup" in edges
    for name, layer in edges.items():
        assert len(layer["r"]) == len(layer["w"]), name


def test_switching_species_opens_the_other_arm_and_remembers_it(qapp):
    """Restores the stored species afterwards. It is one global QSettings scope for the whole run, and
    leaving it on Plasmodium moved seventeen later tests onto a table whose genes they do not
    contain -- silently, because a Window built with no species argument reads this key."""
    from PyQt6 import QtCore
    from starplast.app import Window
    settings = QtCore.QSettings("starplast", "starplast")
    before = settings.value("data/species") if settings.contains("data/species") else None
    w = Window(species="Toxoplasma gondii")
    other = w.open_species("Plasmodium falciparum")
    try:
        assert other is not w
        assert other.species == "Plasmodium falciparum"
        assert str(other.settings().value("data/species")) == "Plasmodium falciparum"
        # Asking for the species already open is a no-op, not a second window.
        assert other.open_species("Plasmodium falciparum") is other
        assert other.open_species("Neospora caninum") is other
    finally:
        other.close()
        settings.remove("data/species") if before is None else \
            settings.setValue("data/species", before)
        settings.sync()


def test_the_species_menu_marks_the_one_that_is_open(qapp):
    from PyQt6 import QtWidgets
    from starplast.app import Window, available_species
    w = Window(species="Plasmodium falciparum")
    try:
        menu = next(m for m in w.menuBar().findChildren(QtWidgets.QMenu) if m.title() == "Species")
        assert [a.text() for a in menu.actions()] == available_species()
        assert [a.text() for a in menu.actions() if a.isChecked()] == ["Plasmodium falciparum"]
    finally:
        w.close()
