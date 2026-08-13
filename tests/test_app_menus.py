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


def test_level_of_detail_is_a_menu_group_with_one_choice_active(win):
    assert len(win.level_group.actions()) == 3
    for i, act in enumerate(win.level_group.actions()):
        act.trigger()
        assert win.level_idx == i
    win.set_level(2)


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
    assert win.comp_list.count() == win.nodes.cellcycle_phase.astype(str).nunique()
    win.comp_list.item(0).setSelected(True)
    vis = win.visible_mask()
    assert 0 < vis.sum() < win.n, "the filter is not restricting anything"
    win.comp_list.clearSelection()
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
def test_the_galaxy_tier_no_longer_collapses_to_the_centre(win):
    """Item 6. The bug: compartment centroids averaged genes spread over the whole map, so every one
    of them landed in the middle."""
    win.set_level(0)
    lab = win.galaxy_labels()
    from starplast import lod
    pos, num, spread = lod.centroids(win.xyz, lab)
    assert len(pos) >= 2
    centre = win.xyz.mean(0)
    R = float(np.linalg.norm(win.xyz - centre, axis=1).max())
    # Each centroid must be tighter than it is far from the middle, or it is a central blob again.
    assert max(np.linalg.norm(p - centre) for p in pos) > 0.25 * R
    assert float(np.median(spread)) < 0.25 * R
    win.set_level(2)


def test_the_galaxy_tier_says_what_the_blobs_are(win):
    """Five unexplained spheres is a mystery, not an abstraction."""
    win.set_level(0)
    win.redraw()
    assert win._galaxy_info
    msg = win.status.currentMessage()
    assert "genes (" in msg
    win.set_level(2)


def test_a_rebuilt_embedding_invalidates_the_cached_tier(win):
    """Left cached, the coarse tier would go on describing the map it replaced.

    The map is put back afterwards. It was not, and once a rebuilt embedding began HIDING the genes
    it does not cover, every export test after this one in the file was quietly exporting a
    100-gene map -- which is also what those tests were doing before the change, without the
    hiding to make it visible."""
    win.galaxy_labels()
    assert win._galaxies is not None
    before, placed_before = win.xyz.copy(), win.placed
    rows = np.zeros(win.n, bool)
    rows[:100] = True
    win.use_embedding(np.random.default_rng(0).normal(size=(100, 3)), rows)
    assert win._galaxies is None
    win.xyz = before
    win.view.xyz = before
    win.placed = placed_before
    win.view.pickable = placed_before
    win.redraw()


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


def test_the_explanations_are_shown_when_asked_for(win, monkeypatch):
    """Item 12 again: the wording is only worth having if the menu entry actually shows it."""
    shown = []
    monkeypatch.setattr(A.QtWidgets.QMessageBox, "information",
                        staticmethod(lambda parent, title, text: shown.append((title, text))))
    win.explain_edges()
    win.explain_map()
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
    win.galaxy_labels()
    assert win._galaxies is not None
    msg = win.free_memory()
    assert win._galaxies is None
    assert "freed" in msg and "level-of-detail" in msg
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
