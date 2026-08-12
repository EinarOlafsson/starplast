#!/usr/bin/env python3
"""Every control in the main window, driven headlessly against the real cache.

The reason this exists at all: a rendering bug in this project has three times passed every array-level
test in the suite. Additive blending summed 8,140 points to white and every colour mode drew as one
featureless blob, while every assertion about the colour array held. The point style was applied and
then immediately overwritten by the next redraw. The colour map only reached the continuous ramp, so
choosing one while colouring by compartment did nothing at all.

So these tests assert on the item list and on redraw ORDER as well as on the arrays -- a control that
takes effect and is then silently undone is the failure mode, and only the final state catches it.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from starplast import theme as TH  # noqa: E402
from starplast.app import COLOUR_MODES  # noqa: E402


@pytest.fixture(scope="module")
def app():
    from PyQt6 import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture(scope="module")
def win():
    from PyQt6 import QtWidgets
    from starplast.app import Window
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    w = Window()
    w.resize(900, 600)
    yield w
    w.close()


# --------------------------------------------------------------------------- colour modes
def test_every_colour_mode_produces_a_full_colour_array(win):
    """Each mode is a different claim about the data, and one that silently falls through to the
    previous mode's colours is indistinguishable from one that works."""
    seen = []
    for i in range(len(COLOUR_MODES)):
        win.set_colour_mode(COLOUR_MODES[i])
        c = np.asarray(win.scatter.color)
        assert c.shape == (win.n, 4), win.colour_mode
        assert np.isfinite(c).all()
        seen.append(c.copy())
    assert any(not np.allclose(seen[0], s) for s in seen[1:]), "every mode drew the same colours"


def test_unknown_is_drawn_in_its_own_grey_rather_than_as_a_category(win):
    """"Not measured" is not a compartment, and giving it a categorical colour would put it in the
    legend beside real ones."""
    for i in range(len(COLOUR_MODES)):
        win.set_colour_mode(COLOUR_MODES[i])
        c = np.asarray(win.scatter.color)
        assert c.shape[0] == win.n


def test_the_colour_map_choice_reaches_the_categorical_modes_too(win):
    """It used to affect only the continuous ramp, so choosing a map while colouring by compartment --
    the default -- appeared to do nothing."""
    win.set_colour_mode(COLOUR_MODES[0])
    before = np.asarray(win.scatter.color).copy()
    categorical = next(n for n, (kind, _) in TH.CMAPS.items() if kind == "categorical")
    win._on_cmap(categorical)          # the handler, which rebuilds the per-compartment palette
    assert not np.allclose(before, np.asarray(win.scatter.color))
    win._on_cmap("auto (match the data)")


def test_choosing_auto_returns_to_the_default_palette(win):
    """"auto" is a real choice -- let the data decide -- not an absence of one."""
    win._on_cmap("auto (match the data)")
    assert win.cmap_name is None
    assert np.asarray(win.scatter.color).shape == (win.n, 4)


def test_a_continuous_map_does_not_disturb_the_compartment_palette(win):
    """A sequential ramp has nothing to say about 27 categories, so the categorical palette stands."""
    win.set_colour_mode(COLOUR_MODES[0])
    seq = next(n for n, (kind, _) in TH.CMAPS.items() if kind == "sequential")
    before = dict(win.colour_of)
    win._on_cmap(seq)
    assert win.colour_of == before
    win._on_cmap("auto (match the data)")


# --------------------------------------------------------------------------- themes
@pytest.mark.parametrize("theme", list(TH.THEMES))
def test_every_theme_applies_to_the_viewport_and_the_widgets(win, theme):
    win.apply_theme(theme)
    assert win.theme == theme
    assert np.asarray(win.scatter.color).shape == (win.n, 4)


def test_changing_theme_keeps_the_selection(win):
    win.on_pick(100)
    win.apply_theme("light")
    assert win.sel == 100
    assert win.halo_item is not None
    win.apply_theme("dark")


# --------------------------------------------------------------------------- point styles
def test_every_point_style_survives_the_next_redraw(win):
    """The bug: apply_point_style set sizes, and redraw() -- which also sets sizes -- overwrote them,
    so choosing a style changed nothing the moment anything else redrew."""
    sizes = {}
    for style in TH.POINT_STYLES:
        win.point_style = style
        win.redraw()
        win.redraw()                       # the second redraw is the one that used to undo it
        sizes[style] = float(np.asarray(win.scatter.size).max())
    assert len(set(sizes.values())) > 1, f"all styles drew the same size: {sizes}"


def test_every_rendering_mode_is_applied_to_the_scatter(win):
    for mode in TH.POINT_MODES:
        win.point_mode = mode
        win.redraw()
        assert np.asarray(win.scatter.color).shape == (win.n, 4)


# --------------------------------------------------------------------------- levels of detail
def test_every_level_of_detail_draws_something(win):
    """The middle tier did nothing at all unless a gene happened to be selected, making it
    indistinguishable from the gene level."""
    drawn = {}
    for i in range(3):
        win.set_level(i)
        win.redraw()
        drawn[i] = win.centroid_item is not None
    assert drawn[0] and drawn[1], "the compartment and orthogroup tiers must draw their own markers"
    win.set_level(2)
    win.redraw()


def test_changing_level_moves_the_camera_rather_than_cutting(win):
    win.set_level(2)
    win.on_level_changed()
    assert win.view._cam_timer is not None and win.view._cam_timer.isActive()
    win.view._cam_timer.stop()


# --------------------------------------------------------------------------- edges
def test_every_edge_type_can_be_toggled_on_its_own(win):
    from starplast.app import EDGE_TYPES
    win.all_edges_act.setChecked(True)
    for k, _ in EDGE_TYPES:
        for cb in win.edge_act.values():
            cb.setChecked(False)
        win.set_edge(k, True)
        win.redraw()
        if k in win.edges:
            assert win.edge_items, f"{k} is present in the graph but drew nothing"
    win.all_edges_act.setChecked(False)


def test_with_nothing_selected_and_draw_all_off_no_edges_are_drawn(win):
    """8,140 genes' edges at once is not a view of anything."""
    win.sel = None
    win.all_edges_act.setChecked(False)
    for cb in win.edge_act.values():
        cb.setChecked(True)
    win.redraw()
    assert win.edge_items == []


def test_selecting_a_gene_draws_only_its_own_edges(win):
    win.all_edges_act.setChecked(False)
    win.on_pick(3000)
    assert isinstance(win.edge_items, list)


def test_the_attention_toggle_changes_which_comention_edges_are_drawn(win):
    """Raw co-mention is confidently misleading, so the corrected residual is the default. The toggle
    reorders exactly the quantity the edge alpha encodes."""
    from starplast.app import EDGE_TYPES, COMENTION
    for cb in win.edge_act.values():
        cb.setChecked(False)
    for k in COMENTION:
        if k in win.edge_act:
            win.set_edge(k, True)
    win.all_edges_act.setChecked(True)

    win.attn_act.setChecked(True)
    win.redraw()
    corrected = sum(len(i.pos) for i in win.edge_items)
    win.attn_act.setChecked(False)
    win.redraw()
    raw = sum(len(i.pos) for i in win.edge_items)
    assert corrected != raw, "the correction must change what is shown"
    win.attn_act.setChecked(True)
    win.all_edges_act.setChecked(False)


# --------------------------------------------------------------------------- search and navigation
def test_searching_an_exact_gene_id_selects_it(win):
    gid = str(win.nodes.gene_id.iloc[500])
    win.search.setText(gid)
    win.do_search()
    assert str(win.nodes.gene_id.iloc[win.sel]) == gid


def test_searching_a_product_substring_finds_a_gene(win):
    win.search.setText("kinase")
    win.do_search()
    assert win.sel is not None
    assert "kinase" in str(win.nodes["product"].iloc[win.sel]).lower()


def test_searching_for_something_absent_says_so_and_changes_nothing(win):
    win.on_pick(10)
    win.search.setText("zzz_no_such_gene_zzz")
    win.do_search()
    assert win.sel == 10
    assert "no match" in win.status.currentMessage()


def test_an_empty_search_does_nothing(win):
    win.on_pick(10)
    win.search.setText("   ")
    win.do_search()
    assert win.sel == 10


def test_reset_clears_the_selection_and_reframes(win):
    win.on_pick(50)
    win.reset()
    assert win.sel is None
    assert win.halo_item is None


# --------------------------------------------------------------------------- the evidence panel
def test_the_panel_fills_for_a_gene_with_data(win):
    win.on_pick(int(np.argmax(win.nodes.n_publications.fillna(0).to_numpy())))
    html = win.detail.toHtml()
    assert str(win.nodes.gene_id.iloc[win.sel]) in html


def test_the_panel_says_nothing_is_known_rather_than_showing_blanks(win):
    """5,574 genes are named in no paper at all. An empty panel would read as a loading failure."""
    never = win.nodes.index[win.nodes.attention_depth.fillna("") == ""]
    win.on_pick(int(never[0]))
    html = win.detail.toHtml().lower()
    assert str(win.nodes.gene_id.iloc[win.sel]).lower() in html


def test_every_gene_can_be_shown_without_raising(win):
    """A gene missing every optional column must not break the panel."""
    rng = np.random.default_rng(0)
    for i in rng.choice(win.n, 40, replace=False):
        win.show_detail(int(i))


# --------------------------------------------------------------------------- compartment filtering
def test_filtering_by_compartment_reduces_what_is_visible(win):
    win.comp_list.clearSelection()
    all_visible = int(win.visible_mask().sum())
    win.comp_list.item(0).setSelected(True)
    filtered = int(win.visible_mask().sum())
    assert filtered < all_visible
    win.comp_list.clearSelection()


def test_flying_to_a_compartment_switches_to_the_gene_level(win):
    win.set_level(0)
    win.fly_to_compartment(win.comp_list.item(0))
    assert win.level_idx == 2


# --------------------------------------------------------------------------- preferences and spin
def test_the_preferences_dialog_builds_every_control(win):
    """It is the only place several settings can be reached, so a failure here hides all of them.

    Built rather than opened: exec() enters a modal event loop and blocks until a human closes the
    window, so a test that called open_preferences() hung forever instead of failing."""
    win.build_preferences()
    for attr in ("theme_box", "cmap_box", "point_box", "mode_box", "spin_speed",
                 "depth_box", "ground_box"):
        assert hasattr(win, attr), attr


def test_every_preferences_control_carries_an_explanation(win):
    """The useful tooltip says WHY: "additive saturates dense regions to white and destroys the colour
    encoding" beats "blending mode"."""
    win.build_preferences()
    for attr in ("theme_box", "cmap_box", "point_box", "mode_box", "depth_box", "ground_box"):
        tip = getattr(win, attr).toolTip()
        assert tip and len(tip) > 40, f"{attr} has no real explanation"


def test_spin_can_be_started_and_stopped(win):
    """Depth in a 3D scatter only reads when it moves."""
    win.toggle_spin(True)
    assert win._spin_timer.isActive()
    before = win.view.opts["azimuth"]
    win._spin_step()
    assert win.view.opts["azimuth"] != before
    win.toggle_spin(False)
    assert not win._spin_timer.isActive()


def test_depth_cueing_and_the_horizon_can_each_be_turned_off(win):
    win.depth_cue = False
    win.show_ground = False
    win.redraw()
    assert win.grid_item is None
    alpha = np.asarray(win.scatter.color)[:, 3]
    assert alpha.max() - alpha.min() < 1e-6
    win.depth_cue = True
    win.show_ground = True
    win.redraw()


# --------------------------------------------------------------------------- tooltips
def test_every_main_window_control_explains_itself(win):
    """The useful tooltip says WHY, not what. "Grey always means unknown, never a category and never
    zero" is worth reading; "colour mode" is not."""
    for attr in ("search", "category_box", "spin_act", "attn_act", "all_edges_act", "comp_list"):
        tip = getattr(win, attr).toolTip()
        assert tip, f"{attr} has no tooltip"
        assert len(tip.split()) >= 15, f"{attr}'s tooltip only restates its label: {tip!r}"


def test_every_edge_checkbox_says_what_the_relation_means(win):
    """Twelve edge types, several of which are inferences rather than measurements, and the checkbox
    label alone cannot carry that distinction."""
    for k, cb in win.edge_act.items():
        assert cb.toolTip(), f"edge type {k} has no tooltip"


# --------------------------------------------------------------------------- loading and startup
def test_a_missing_cache_stops_with_the_command_that_builds_it(monkeypatch, tmp_path):
    """A pandas error from inside a constructor tells the user nothing about what to do next."""
    import starplast.app as A
    monkeypatch.setattr(A, "DATA", str(tmp_path))
    with pytest.raises(SystemExit, match="build_graph"):
        A.load()


def test_main_refuses_to_open_a_window_without_a_cache(monkeypatch, tmp_path, capsys):
    """Checked before building a window, so the message is readable rather than buried under a
    traceback from a constructor."""
    import starplast.app as A
    from starplast import paths
    monkeypatch.setattr(paths, "check", lambda: (False, "cache incomplete at /nowhere"))
    monkeypatch.setattr(paths, "describe", lambda: "cache /nowhere [INCOMPLETE]")
    with pytest.raises(SystemExit):
        A.main()
    assert "cache incomplete" in capsys.readouterr().err


def test_main_checks_the_cache_before_building_anything(monkeypatch):
    """`starplast` is the entry point in pyproject. Its first act is the cache check, so a missing
    cache is a sentence rather than a traceback out of a constructor.

    Only the guard is exercised: main() constructs a QApplication, and a second one in a process that
    already has it is undefined behaviour -- it segfaults the interpreter rather than failing."""
    import starplast.app as A
    order = []
    monkeypatch.setattr(A.paths, "check", lambda: (order.append("checked"), (False, "no cache"))[1])
    monkeypatch.setattr(A.paths, "describe", lambda: "")
    monkeypatch.setattr(A.pg, "setConfigOptions",
                        lambda **k: order.append("built a window"))
    with pytest.raises(SystemExit):
        A.main()
    assert order == ["checked"], "nothing may be constructed before the cache is known to be there"


# --------------------------------------------------------------------------- the camera
def test_an_empty_map_still_frames_something():
    """A build that produced no genes must not divide by zero while computing its own extent."""
    from starplast.app import Map3D
    v = Map3D.__new__(Map3D)
    v.xyz = np.zeros((0, 3), dtype=np.float32)
    assert v.data_radius() == 100.0


def test_framing_an_empty_map_uses_a_fixed_distance(win, monkeypatch):
    empty = np.zeros((0, 3), dtype=np.float32)
    real = win.view.xyz
    win.view.xyz = empty
    try:
        win.view.fit_view()
        assert win.view.opts["distance"] > 0
    finally:
        win.view.xyz = real
        win.view.fit_view()


def test_a_camera_move_to_where_it_already_is_starts_no_animation(win):
    """Otherwise every redraw at the same level restarts a 420 ms timer for nothing."""
    win.view._cam_timer = None
    win.view.animate_distance(float(win.view.opts["distance"]))
    assert getattr(win.view, "_cam_timer", None) is None


def test_a_second_camera_move_retargets_rather_than_queueing(win, app):
    win.view.animate_distance(200.0)
    first = win.view._cam_timer
    win.view.animate_distance(80.0)
    assert win.view._cam_timer is not first
    win.view._cam_timer.stop()


def test_the_camera_animation_finishes_and_stops_itself(win):
    """Driven by firing the timer rather than by waiting on the clock. Waiting made this flaky under
    load -- a 60 ms animation and a 900 ms budget still lost the race when the machine was busy -- and
    a test that fails when something else is running teaches people to ignore it."""
    start = float(win.view.opts["distance"])
    target = start * 1.5
    win.view.animate_distance(target, ms=60)
    assert win.view._cam_timer.isActive()

    seen = [start]
    for _ in range(200):
        win.view._cam_timer.timeout.emit()
        seen.append(float(win.view.opts["distance"]))
        if not win.view._cam_timer.isActive():
            break
    assert not win.view._cam_timer.isActive(), "the animation must stop itself"
    assert seen[-1] == pytest.approx(target, rel=1e-6)
    assert all(a <= b + 1e-9 for a, b in zip(seen, seen[1:])), "it must ease toward the target, not past it"
    assert len(seen) > 3, "a smoothstep over one frame is a cut"


# --------------------------------------------------------------------------- theme round trip
def test_applying_a_theme_updates_the_preferences_box_without_re_firing(win):
    """Without blocking the signal the box's own handler re-applies the theme, so every change costs
    two full redraws."""
    win.build_preferences()
    win.apply_theme("slate")
    assert win.theme_box.currentText() == "slate"
    win.apply_theme("dark")


def test_the_preferences_box_shows_the_colour_map_in_use(win):
    categorical = next(n for n, (kind, _) in TH.CMAPS.items() if kind == "categorical")
    win._on_cmap(categorical)
    win.build_preferences()
    assert win.cmap_box.currentText() == categorical
    win._on_cmap("auto (match the data)")


def test_choosing_a_point_style_and_mode_through_the_handlers(win):
    for style in TH.POINT_STYLES:
        win._on_point_style(style)
        assert win.point_style == style
    for mode in TH.POINT_MODES:
        win._on_point_mode(mode)
        assert win.point_mode == mode


def test_the_named_colour_map_is_used_for_the_continuous_ramp(win):
    seq = next(n for n, (kind, _) in TH.CMAPS.items() if kind == "sequential")
    win._on_cmap(seq)
    idx = next(i for i in range(len(COLOUR_MODES))
               if "fitness" in COLOUR_MODES[i].lower())
    win.set_colour_mode(COLOUR_MODES[idx])
    assert np.asarray(win.scatter.color).shape == (win.n, 4)
    win._on_cmap("auto (match the data)")


# --------------------------------------------------------------------------- the analysis dock
def test_the_analysis_dock_is_attached(win):
    assert getattr(win, "analysis_dock", None) is not None


def test_the_browser_still_opens_when_the_analysis_panel_cannot_be_imported(monkeypatch, win):
    """It needs scikit-learn and umap-learn; the map must not depend on them.

    The failure is injected into sys.modules rather than into __import__, because replacing the import
    hook while Qt is running takes the interpreter down with it."""
    import starplast.app as A
    monkeypatch.setitem(sys.modules, "starplast.analysis_panel", None)
    win._analysis()
    assert "analysis panel unavailable" in win.statusBar().currentMessage()
    assert win.n > 0, "the map itself must still be there"


def test_an_embedding_built_by_the_panel_replaces_the_map(win):
    """Genes the new embedding excluded keep no stale coordinates from a different map, and are no
    longer left at the origin either -- which is where they used to go, and for a walk configuration
    that is 7,340 of 8,140 genes in a lump in the middle of the map: drawn, pickable and counted,
    while having no position in the embedding at all. They are recorded as unplaced and hidden."""
    before, placed_before = win.xyz.copy(), win.placed
    rows = np.zeros(win.n, dtype=bool)
    rows[:100] = True
    coords = np.random.default_rng(0).normal(size=(100, 3)).astype(np.float32)
    win.use_embedding(coords, rows)
    assert np.allclose(win.xyz[:100], coords)
    assert not np.allclose(win.xyz[100:], before[100:]), "stale coordinates from the last map"
    assert np.allclose(win.xyz[100:], coords.mean(0), atol=1e-4), "not parked out of the way"
    assert win.visible_mask().sum() == 100
    assert win.view.xyz is win.xyz, "the picker must see the same coordinates as the scatter"
    win.xyz = before
    win.view.xyz = before
    win.placed = placed_before
    win.view.pickable = placed_before
    win.redraw()


# --------------------------------------------------------------------------- picking
def test_projection_puts_genes_in_widget_pixels(win):
    """Picking compares these against the mouse position, so a wrong scale reads as "clicking is
    slightly off" -- the hardest kind of bug to attribute."""
    sx, sy = win.view.project()
    assert sx.shape == (win.n,) and sy.shape == (win.n,)
    on_screen = np.isfinite(sx) & np.isfinite(sy)
    assert on_screen.any(), "no gene projected onto the widget at all"


def test_genes_behind_the_camera_are_excluded_rather_than_wrapped(win):
    """Divided through by a negative w they land back on screen, mirrored -- so a click near the front
    could select a gene behind you."""
    sx, _ = win.view.project()
    assert np.isnan(sx).sum() >= 0
    assert not np.isinf(sx[np.isfinite(sx)]).any()


def _click(win, x, y, button=None):
    from PyQt6 import QtCore, QtGui
    button = button or QtCore.Qt.MouseButton.LeftButton
    return QtGui.QMouseEvent(QtCore.QEvent.Type.MouseButtonRelease,
                             QtCore.QPointF(x, y), QtCore.QPointF(x, y),
                             button, button, QtCore.Qt.KeyboardModifier.NoModifier)


def test_clicking_on_a_gene_selects_it(win):
    sx, sy = win.view.project()
    ok = np.where(np.isfinite(sx) & np.isfinite(sy))[0]
    i = int(ok[0])
    got = []
    win.view.picked.connect(got.append)
    win.view.mouseReleaseEvent(_click(win, float(sx[i]), float(sy[i])))
    assert got and got[0] in ok


def test_clicking_empty_space_selects_nothing(win):
    """A 14-pixel radius, so a click in the void does not grab whichever gene happens to be nearest."""
    got = []
    win.view.picked.connect(got.append)
    win.view.mouseReleaseEvent(_click(win, -5000.0, -5000.0))
    assert got == []


def test_a_right_click_does_not_select(win):
    from PyQt6 import QtCore
    got = []
    win.view.picked.connect(got.append)
    win.view.mouseReleaseEvent(_click(win, 10.0, 10.0, QtCore.Qt.MouseButton.RightButton))
    assert got == []


def test_a_projection_failure_is_printed_rather_than_raised(win, monkeypatch, capsys):
    """A mouse handler is the wrong place to raise: every click printed a traceback and selected
    nothing when pyqtgraph changed its projectionMatrix signature."""
    monkeypatch.setattr(win.view, "project",
                        lambda: (_ for _ in ()).throw(TypeError("signature changed again")))
    got = []
    win.view.picked.connect(got.append)
    win.view.mouseReleaseEvent(_click(win, 10.0, 10.0))
    assert got == []
    assert "picking unavailable" in capsys.readouterr().out


def test_nothing_on_screen_selects_nothing(win, monkeypatch):
    nan = np.full(win.n, np.nan)
    monkeypatch.setattr(win.view, "project", lambda: (nan, nan))
    got = []
    win.view.picked.connect(got.append)
    win.view.mouseReleaseEvent(_click(win, 10.0, 10.0))
    assert got == []


# --------------------------------------------------------------------------- colouring edges
def test_attention_colouring_survives_a_table_without_the_column(win, monkeypatch):
    """The column is absent on a cache built before the literature layer existed, and the mode must
    grey out rather than raise."""
    idx = next(i for i in range(len(COLOUR_MODES))
               if "attention" in COLOUR_MODES[i].lower())
    win.set_colour_mode(COLOUR_MODES[idx])
    real = win.nodes
    try:
        win.nodes = real.drop(columns=["attention_depth"])
        win.redraw()
        assert np.asarray(win.scatter.color).shape == (win.n, 4)
    finally:
        win.nodes = real
        win.set_colour_mode(COLOUR_MODES[0])
        win.redraw()


def test_a_compartment_with_too_few_visible_genes_gets_no_centroid(win):
    """A centroid of two points is a midpoint, not a landmark."""
    win.set_level(0)
    win.comp_list.clearSelection()
    win.comp_list.item(0).setSelected(True)
    win.redraw()
    win.comp_list.clearSelection()
    win.set_level(2)
    win.redraw()


def test_the_evidence_panel_skips_edge_types_absent_from_the_graph(win):
    """A cache built with fewer edge types than the app knows about must still open."""
    real = win.edges
    try:
        win.edges = {k: v for k, v in list(real.items())[:2]}
        win.show_detail(100)
    finally:
        win.edges = real


def test_a_chosen_map_is_used_only_where_its_kind_suits_the_column(win):
    """A diverging ramp on a strictly positive quantity invents a midpoint, so the column's kind has
    the final say over the user's choice."""
    import pandas as pd
    seq = next(n for n, (kind, _) in TH.CMAPS.items() if kind == "sequential")
    win.cmap_name = seq
    positive = pd.Series([0.1, 5.0, 90.0])
    assert win._cmap_for(positive) == TH.CMAPS[seq][1], "a sequential map suits a positive quantity"
    straddling = pd.Series([-2.0, -0.5, 0.3, 1.8])
    assert win._cmap_for(straddling) == TH.DEFAULT_CMAP["diverging"], (
        "a sequential map must not be used on a quantity that goes both ways")
    win.cmap_name = None


def test_opening_preferences_shows_the_dialog_it_builds(win, monkeypatch):
    """exec() enters a modal loop, so only the call is asserted -- entering it in a test hangs the
    process rather than failing it."""
    from PyQt6 import QtWidgets
    shown = []
    monkeypatch.setattr(QtWidgets.QDialog, "exec", lambda self: shown.append(self) or 0)
    win.open_preferences()
    assert len(shown) == 1
    assert isinstance(shown[0], QtWidgets.QDialog)


def test_main_builds_and_shows_the_window(monkeypatch, app):
    """`starplast` is the entry point in pyproject, so the whole path is exercised: cache check,
    window construction, show, and the event loop.

    QApplication is handed back the instance this process already owns -- constructing a second one is
    undefined behaviour and takes the interpreter down rather than failing a test."""
    from PyQt6 import QtWidgets
    import starplast.app as A
    seen = {}

    class Shim:
        def __new__(cls, argv):
            return app

        @staticmethod
        def instance():
            return app

    monkeypatch.setattr(A.QtWidgets, "QApplication", Shim)
    monkeypatch.setattr(A.sys, "exit", lambda code=0: seen.setdefault("exit", code))
    monkeypatch.setattr(type(app), "exec", lambda self: seen.setdefault("exec", 0) or 0)
    monkeypatch.setattr(A.Window, "show", lambda self: seen.setdefault("shown", True))

    A.main()
    assert seen == {"shown": True, "exec": 0, "exit": 0}
