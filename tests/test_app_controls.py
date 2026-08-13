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

from PyQt6 import QtCore  # noqa: E402
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


# --------------------------------------------------------------------------- logging preferences
def test_logging_is_off_until_it_is_asked_for(win, tmp_path):
    """Writing files to someone's disk uninvited is how a tool loses trust. The console still gets
    warnings, because a warning nobody enabled a log to see is a silent failure with extra steps."""
    from starplast import logging_util
    win.apply_log_settings(enabled=False, directory=str(tmp_path))
    assert logging_util.log_file() == ""
    assert os.listdir(tmp_path) == []
    assert logging_util.state()["console_level"] in logging_util.LEVELS


def test_enabling_the_log_writes_a_file_and_says_where(win, tmp_path):
    from starplast import logging_util
    path = win.apply_log_settings(enabled=True, file_level="DEBUG", directory=str(tmp_path))
    assert path and os.path.exists(path)
    assert path in win.statusBar().currentMessage()
    logging_util.get_logger("test").info("a line")
    assert "a line" in open(path).read()
    win.apply_log_settings(enabled=False, directory=str(tmp_path))


def test_the_logging_choice_persists(win, tmp_path):
    """A preference that has to be set again every launch is a preference nobody sets."""
    win.apply_log_settings(enabled=True, console_level="DEBUG", directory=str(tmp_path))
    assert win.log_settings()["enabled"] is True
    assert win.log_settings()["console_level"] == "DEBUG"
    win.apply_log_settings(enabled=False, console_level="WARNING")
    assert win.log_settings()["enabled"] is False


def test_the_preferences_dialog_carries_the_logging_controls(win):
    win.build_preferences()
    for attr in ("log_box", "log_file_level", "log_console_level", "log_path"):
        assert hasattr(win, attr), attr
        tip = getattr(win, attr).toolTip()
        assert tip and len(tip.split()) >= 15, f"{attr} does not explain itself: {tip!r}"


def test_the_preference_controls_take_effect_when_used(win, tmp_path, monkeypatch):
    """A checkbox that stores a setting without applying it is the worst of both."""
    from starplast import logging_util
    monkeypatch.setattr(logging_util, "log_dir", lambda: str(tmp_path))
    win.build_preferences()
    win.log_console_level.setCurrentText("DEBUG")
    assert logging_util.state()["console_level"] == "DEBUG"
    win.log_box.setChecked(True)
    assert logging_util.log_file(), "ticking the box did not start a log"
    assert win.log_path.text() == logging_util.log_file()
    win.log_box.setChecked(False)
    assert logging_util.log_file() == "" and "no file" in win.log_path.text()


# --------------------------------------------------------------------------- navigate and select
def test_the_left_button_has_two_modes_and_says_which(win):
    """A silent mode change is a trap: the same drag rotates in one and gates in the other."""
    from starplast.app import INTERACTION_MODES
    for name in INTERACTION_MODES:
        win.set_interaction_mode(name)
        assert win.view.mode == name
        assert name in win.statusBar().currentMessage()
    win.set_interaction_mode("something else")
    assert win.view.mode == "navigate", "an unknown mode must fall back, not disable the mouse"


def test_a_constrained_orbit_moves_one_angle_and_leaves_the_other(win):
    """An unconstrained orbit never returns to the same view twice, which is exactly wrong for
    comparing two maps."""
    from PyQt6 import QtCore, QtGui

    def drag(dx, dy):
        start = QtCore.QPointF(300.0, 300.0)
        win.view.mousePos = start
        pos = QtCore.QPointF(300.0 + dx, 300.0 + dy)
        ev = QtGui.QMouseEvent(QtCore.QEvent.Type.MouseMove, pos, pos,
                               QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.MouseButton.LeftButton,
                               QtCore.Qt.KeyboardModifier.NoModifier)
        win.view.mouseMoveEvent(ev)

    win.set_interaction_mode("navigate")
    win.set_navigate_axis("z")
    before = dict(win.view.opts)
    drag(40, 40)
    assert win.view.opts["azimuth"] != before["azimuth"]
    assert win.view.opts["elevation"] == before["elevation"], "z rotation changed the elevation"

    win.set_navigate_axis("x")
    drag(40, 40)
    assert win.view.opts["azimuth"] == 0.0, "the x axis must pin the azimuth, or it is not an axis"
    win.set_navigate_axis("y")
    drag(0, 20)
    assert win.view.opts["azimuth"] == 90.0
    win.set_navigate_axis("free")


def test_free_orbit_is_still_the_default_behaviour(win):
    from PyQt6 import QtCore, QtGui
    win.set_navigate_axis("free")
    win.set_interaction_mode("navigate")
    before = dict(win.view.opts)
    pos = QtCore.QPointF(320.0, 320.0)
    win.view.mousePos = QtCore.QPointF(300.0, 300.0)
    win.view.mouseMoveEvent(QtGui.QMouseEvent(
        QtCore.QEvent.Type.MouseMove, pos, pos, QtCore.Qt.MouseButton.LeftButton,
        QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.KeyboardModifier.NoModifier))
    assert (win.view.opts["azimuth"] != before["azimuth"]
            or win.view.opts["elevation"] != before["elevation"])


def _on_screen(win):
    sx, sy = win.view.project()
    ok = np.flatnonzero(np.isfinite(sx) & np.isfinite(sy))
    assert len(ok), "no gene projected onto the widget"
    return sx, sy, ok


def test_a_lasso_takes_the_genes_inside_it_and_nothing_else(win):
    """The gate is a display claim -- it says "these genes are the ones you drew around" -- so it is
    checked against the projection rather than trusted."""
    win.set_interaction_mode("select")
    win.set_gate_shape("lasso (2D)")
    sx, sy, ok = _on_screen(win)
    cx, cy = float(sx[ok[0]]), float(sy[ok[0]])
    box = [(cx - 30, cy - 30), (cx + 30, cy - 30), (cx + 30, cy + 30), (cx - 30, cy + 30)]
    win.view.begin_gate(*box[0])
    for x, y in box[1:]:
        win.view.extend_gate(x, y)
    idx = win.view.finish_gate()
    assert len(idx), "a gate over a gene caught nothing"
    inside = (np.abs(sx[idx] - cx) <= 30 + 1) & (np.abs(sy[idx] - cy) <= 30 + 1)
    assert inside.all(), "a gene outside the lasso was gated"
    expected = np.flatnonzero((np.abs(sx - cx) <= 29) & (np.abs(sy - cy) <= 29))
    assert set(expected) <= set(idx), "a gene inside the lasso was missed"
    win.clear_gate()
    win.set_interaction_mode("navigate")


def test_a_lasso_of_two_points_is_not_a_polygon_and_gates_nothing(win):
    win.set_interaction_mode("select")
    win.set_gate_shape("lasso (2D)")
    sx, sy, ok = _on_screen(win)
    win.view.begin_gate(float(sx[ok[0]]), float(sy[ok[0]]))
    win.view.extend_gate(float(sx[ok[0]]) + 20, float(sy[ok[0]]))
    assert len(win.view.finish_gate()) == 0
    win.set_interaction_mode("navigate")


def test_a_brush_gates_a_ball_in_world_space_not_a_disc_on_screen(win):
    """The whole reason the 3D gate exists: deep in the cloud a lasso also catches the far side,
    which looks like a selection of one structure and is a selection of two."""
    win.set_interaction_mode("select")
    win.set_gate_shape("brush (3D)")
    sx, sy, ok = _on_screen(win)
    i = int(ok[0])
    cx, cy = float(sx[i]), float(sy[i])
    win.view.begin_gate(cx, cy)
    win.view.extend_gate(cx + 40, cy)
    idx = win.view.finish_gate()
    assert len(idx), "the brush caught nothing"
    anchor = win.xyz[i]
    d = np.linalg.norm(win.xyz[idx] - anchor, axis=1)
    assert d.max() <= np.linalg.norm(win.xyz - anchor, axis=1).max(), "gated beyond the cloud"
    # Every gated gene is nearer the anchor IN WORLD SPACE than the farthest one gated, which a
    # screen-space disc would not guarantee: the far side of the cloud projects into the same disc.
    outside = np.setdiff1d(np.arange(win.n), idx)
    assert np.linalg.norm(win.xyz[outside] - anchor, axis=1).min() >= d.max() - 1e-6
    win.clear_gate()
    win.set_interaction_mode("navigate")


def test_a_brush_that_is_barely_dragged_gates_nothing(win):
    """A click in select mode is not a one-gene gate; it is a click, and gating one gene by accident
    would be indistinguishable from a gate that failed."""
    win.set_interaction_mode("select")
    win.set_gate_shape("brush (3D)")
    sx, sy, ok = _on_screen(win)
    win.view.begin_gate(float(sx[ok[0]]), float(sy[ok[0]]))
    win.view.extend_gate(float(sx[ok[0]]) + 1, float(sy[ok[0]]))
    assert len(win.view.finish_gate()) == 0
    win.set_interaction_mode("navigate")


def test_finishing_without_starting_gates_nothing(win):
    assert len(win.view.finish_gate()) == 0


def test_a_gate_never_takes_a_gene_that_has_no_position(win):
    """The genes a walk configuration does not cover are not on the map; a gate over where they
    would have been must not collect them."""
    placed = np.zeros(win.n, bool)
    placed[:200] = True
    win.view.pickable = placed
    win.set_interaction_mode("select")
    win.set_gate_shape("lasso (2D)")
    sx, sy, _ = _on_screen(win)
    for x, y in [(-1e4, -1e4), (1e4, -1e4), (1e4, 1e4), (-1e4, 1e4)]:   # everything
        (win.view.begin_gate if x < 0 and y < 0 else win.view.extend_gate)(x, y)
    idx = win.view.finish_gate()
    assert placed[idx].all(), "a gene with no position was gated"
    win.view.pickable = None
    win.clear_gate()
    win.set_interaction_mode("navigate")


def test_a_gated_set_recedes_the_rest_of_the_map_without_recolouring_it(win):
    """A gate is a selection, not a claim about the data. Recolouring the gated genes would put a
    selection into the one channel that means measurement, inference or absence."""
    win.set_colour_mode("compartment")
    before = np.asarray(win.scatter.color).copy()
    win.on_gated(np.arange(50))
    after = np.asarray(win.scatter.color)
    assert np.allclose(after[:50, :3], before[:50, :3]), "the gated genes were recoloured"
    assert after[500, 3] < before[500, 3], "the rest of the map did not recede"
    assert np.asarray(win.scatter.size)[:50].max() > np.asarray(win.scatter.size)[500]
    win.clear_gate()


def test_a_gate_reports_what_is_in_it_and_leads_with_the_unlabelled(win):
    """The question a gate is drawn to answer is "what is this clump", and the genes with no label
    are the candidates it exists to produce."""
    win.on_gated(np.arange(300))
    html = win.detail.toHtml()
    assert "300" in html and win.category in html
    assert "candidates" in html
    win.clear_gate()
    assert "Click a gene" in win.detail.toHtml()


def test_an_empty_gate_says_it_is_empty_rather_than_looking_broken(win):
    win.on_gated(np.array([], dtype=int))
    assert "empty, not broken" in win.statusBar().currentMessage()
    win.clear_gate()


def test_a_gated_set_can_be_exported(win, tmp_path):
    """A gated set is the natural input to annotation, and it has to be able to leave the window."""
    import pandas as pd
    win.on_gated(np.arange(25))
    path = win.export_gated(str(tmp_path / "gated.csv"), columns=["compartment"])
    got = pd.read_csv(path)
    assert len(got) == 25 and {"gene_id", "compartment", "x", "y", "z"} <= set(got.columns)
    assert list(got.gene_id) == list(win.nodes.gene_id.iloc[:25])
    win.clear_gate()


def test_exporting_without_a_gate_says_what_to_do(win, tmp_path):
    win.gated = None
    assert win.export_gated(str(tmp_path / "none.csv")) is None
    assert "no gated selection" in win.statusBar().currentMessage()


def test_the_gate_is_reachable_from_the_menus(win):
    win.on_gated(np.arange(10))
    labels = [a.text() for a in win.build_context_menu().actions()]
    assert any("gated genes" in a for a in labels)
    assert any("Clear the gate" in a for a in labels)
    win.clear_gate()
    assert not any("Clear the gate" in a.text() for a in win.build_context_menu().actions())


def test_the_gate_being_drawn_is_painted_over_the_view(win):
    """A gate you cannot see while dragging it is a gate drawn by guesswork. Painted onto an image
    rather than rendered from the widget, which would bring the widget's background with it and
    make "drew a lasso" indistinguishable from "filled the rectangle"."""
    from PyQt6 import QtGui
    ov = win.view.overlay

    def painted():
        img = QtGui.QImage(80, 80, QtGui.QImage.Format.Format_ARGB32_Premultiplied)
        img.fill(0)
        p = QtGui.QPainter(img)
        ov.draw(p)
        p.end()
        return sum(QtGui.QColor.fromRgba(img.pixel(x, y)).alpha() > 0
                   for x in range(0, 80, 2) for y in range(0, 80, 2))

    ov.clear()
    assert painted() == 0, "something was drawn with no gate in progress"
    ov.show_lasso([(10, 10), (70, 10), (70, 70), (10, 70)])
    assert painted() > 0, "the lasso was not drawn"
    ov.show_brush(40, 40, 25)
    assert painted() > 0, "the brush was not drawn"
    ov.clear()
    assert painted() == 0, "clearing left the last gate on screen"


def test_the_overlay_follows_the_size_of_the_view(win):
    """Left at its original size it would clip the gate the moment the window is resized.

    The handler is driven directly: these windows are never shown, and Qt defers a resize event to
    the moment a widget becomes visible -- so `resize()` alone proves nothing here."""
    from PyQt6 import QtCore, QtGui
    win.view.resize(640, 480)
    win.view.resizeEvent(QtGui.QResizeEvent(QtCore.QSize(640, 480), QtCore.QSize(399, 950)))
    assert win.view.overlay.size() == win.view.size()


def test_the_mouse_drives_the_gate_end_to_end(win):
    """The handlers, not the helpers: connected to the wrong events, every gate test above would
    still pass while dragging on the map did nothing at all."""
    from PyQt6 import QtCore, QtGui

    def ev(kind, x, y):
        pos = QtCore.QPointF(x, y)
        return QtGui.QMouseEvent(kind, pos, pos, QtCore.Qt.MouseButton.LeftButton,
                                 QtCore.Qt.MouseButton.LeftButton,
                                 QtCore.Qt.KeyboardModifier.NoModifier)

    sx, sy, ok = _on_screen(win)
    cx, cy = float(sx[ok[0]]), float(sy[ok[0]])
    got = []
    win.view.gated.connect(got.append)
    win.set_interaction_mode("select")
    win.set_gate_shape("lasso (2D)")
    win.view.mousePressEvent(ev(QtCore.QEvent.Type.MouseButtonPress, cx - 30, cy - 30))
    for x, y in ((cx + 30, cy - 30), (cx + 30, cy + 30), (cx - 30, cy + 30)):
        win.view.mouseMoveEvent(ev(QtCore.QEvent.Type.MouseMove, x, y))
    win.view.mouseReleaseEvent(ev(QtCore.QEvent.Type.MouseButtonRelease, cx - 30, cy + 30))
    assert got and len(got[0]), "dragging in select mode gated nothing"
    win.clear_gate()
    win.set_interaction_mode("navigate")


def test_a_click_in_navigate_mode_still_selects_one_gene(win):
    """Select mode must not be the only way to use the map; the ordinary click has to survive it."""
    from PyQt6 import QtCore, QtGui
    win.set_interaction_mode("navigate")
    sx, sy, ok = _on_screen(win)
    i = int(ok[0])
    pos = QtCore.QPointF(float(sx[i]), float(sy[i]))
    got = []
    win.view.picked.connect(got.append)
    win.view.mouseReleaseEvent(QtGui.QMouseEvent(
        QtCore.QEvent.Type.MouseButtonRelease, pos, pos, QtCore.Qt.MouseButton.LeftButton,
        QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.KeyboardModifier.NoModifier))
    assert got


def test_gating_survives_a_projection_that_cannot_be_computed(win, monkeypatch, capsys):
    """A mouse handler is the wrong place to raise: this printed a traceback per click once already,
    through the picking path."""
    win.set_interaction_mode("select")
    win.set_gate_shape("lasso (2D)")
    win.view.begin_gate(10.0, 10.0)
    for pt in ((90.0, 10.0), (90.0, 90.0)):
        win.view.extend_gate(*pt)
    monkeypatch.setattr(win.view, "project", lambda: (_ for _ in ()).throw(RuntimeError("no GL")))
    assert len(win.view.finish_gate()) == 0
    assert "gating unavailable" in capsys.readouterr().out
    monkeypatch.undo()
    assert win.view.nearest(10.0, 10.0) is not None
    monkeypatch.setattr(win.view, "project", lambda: (_ for _ in ()).throw(RuntimeError("no GL")))
    assert win.view.nearest(10.0, 10.0) is None
    win.set_interaction_mode("navigate")


def test_the_brush_falls_back_when_too_little_is_on_screen(win):
    """The pixels-to-world scale is measured from the genes near the press. With almost nothing
    there it has to guess from the data radius rather than divide by an empty set."""
    scale = win.view._world_per_pixel(0, -5000.0, -5000.0, *win.view.project())
    assert scale > 0 and np.isfinite(scale)


def test_extending_a_gate_that_was_never_started_does_nothing(win):
    win.view._gate = None
    win.view.extend_gate(10.0, 10.0)
    assert win.view._gate is None


def test_the_overlay_paints_only_when_a_gate_is_in_progress(win):
    """The paint handler itself, driven directly: an overlay that painted on every frame would put
    the last gate back on screen after it had been cleared."""
    from PyQt6 import QtCore, QtGui
    ov = win.view.overlay
    ov.clear()
    ov.paintEvent(QtGui.QPaintEvent(QtCore.QRect(0, 0, 10, 10)))     # returns before painting
    ov.resize(60, 60)
    ov.show_brush(30, 30, 10)
    ov.paintEvent(QtGui.QPaintEvent(QtCore.QRect(0, 0, 60, 60)))
    ov.clear()


def test_the_brush_anchor_ignores_genes_with_no_position(win):
    """`nearest` is shared with picking, so the mask has to hold for both or a brush would anchor on
    a gene that is not in the displayed embedding."""
    placed = np.zeros(win.n, bool)
    placed[:5] = True
    win.view.pickable = placed
    sx, sy = win.view.project()
    ok = np.flatnonzero(np.isfinite(sx))
    i = win.view.nearest(float(sx[ok[0]]), float(sy[ok[0]]))
    assert i is None or placed[i]
    win.view.pickable = np.zeros(win.n, bool)          # nothing has a position at all
    assert win.view.nearest(10.0, 10.0) is None
    win.view.pickable = None


def test_nothing_is_returned_when_the_nearest_gene_is_too_far(win):
    sx, sy = win.view.project()
    assert win.view.nearest(-9999.0, -9999.0, within=5.0) is None


def test_a_gated_export_can_be_cancelled_at_either_dialog(win, monkeypatch, tmp_path):
    """Two dialogs stand between the menu item and the file, and cancelling either must write
    nothing rather than a file with default columns."""
    from PyQt6 import QtWidgets
    win.on_gated(np.arange(5))
    monkeypatch.setattr(win, "_ask_path", lambda *a, **k: "")
    assert win.export_gated() is None
    monkeypatch.setattr(win, "_ask_path", lambda *a, **k: str(tmp_path / "x.csv"))
    monkeypatch.setattr(QtWidgets.QDialog, "exec", lambda self: 0)
    assert win.export_gated() is None
    assert not (tmp_path / "x.csv").exists()
    win.clear_gate()


def test_pressing_in_navigate_mode_reaches_the_camera(win):
    """Select mode intercepts the press; navigate mode must still hand it to the view, or the map
    stops rotating the moment gating exists."""
    from PyQt6 import QtCore, QtGui
    win.set_interaction_mode("navigate")
    pos = QtCore.QPointF(120.0, 140.0)
    win.view.mousePressEvent(QtGui.QMouseEvent(
        QtCore.QEvent.Type.MouseButtonPress, pos, pos, QtCore.Qt.MouseButton.LeftButton,
        QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.KeyboardModifier.NoModifier))
    assert win.view.mousePos == pos, "the press never reached the camera"
    assert win.view._gate is None, "navigate mode started a gate"


def test_the_gated_export_uses_the_columns_that_were_ticked(win, monkeypatch, tmp_path):
    import pandas as pd
    from PyQt6 import QtWidgets
    win.on_gated(np.arange(8))
    monkeypatch.setattr(QtWidgets.QDialog, "exec", lambda self: 1)
    monkeypatch.setattr(win, "_ticked", staticmethod(lambda d: ["compartment", "n_publications"]))
    path = win.export_gated(str(tmp_path / "ticked.csv"))
    got = pd.read_csv(path)
    assert {"gene_id", "compartment", "n_publications", "x", "y", "z"} == set(got.columns)
    win.clear_gate()


# --------------------------------------------------------------------------- annotations on the map
def test_annotations_are_drawn_in_a_colour_of_their_own(win, tmp_path):
    """Measurement, inference and absence have one each. A proposal is a fourth thing, and reading
    as any of the three is the failure this application is built to prevent."""
    from starplast.annotations import ANNOTATION_COLOUR, Annotation, AnnotationStore
    store = AnnotationStore(str(tmp_path / "a.csv"))
    gene = str(win.nodes.gene_id.iloc[7])
    store.save(Annotation(gene_id=gene, proposed="dense granules", target="compartment", cluster=1,
                          precision=0.6, recall=0.4, date="2026-08-12"))
    win.annotations = store
    win.refresh_annotations()
    assert win.colour_mode == "annotations"
    c = np.asarray(win.scatter.color)
    assert np.allclose(c[7, :3], ANNOTATION_COLOUR, atol=1e-3)
    # Everything else is grey: not a category, not zero -- nobody has proposed anything for it.
    assert not np.allclose(c[8, :3], ANNOTATION_COLOUR)
    assert "proposals, not measurements" in win.statusBar().currentMessage()
    win.set_colour_mode("compartment")


def test_with_no_annotations_the_mode_draws_everything_as_unknown(win, tmp_path):
    from starplast.annotations import AnnotationStore
    win.annotations = AnnotationStore(str(tmp_path / "none.csv"))
    win.sel = None                     # a selected gene is drawn white, which is a third colour
    win.refresh_annotations()
    win.set_colour_mode("annotations")
    c = np.asarray(win.scatter.color)
    assert len({tuple(np.round(x, 3)) for x in c[:, :3]}) == 1, "something was coloured as annotated"
    win.set_colour_mode("compartment")


def test_an_unreadable_annotations_file_does_not_take_the_window_down(win, tmp_path, capsys):
    """These files are hand-edited and shared. A broken one must cost the colour, not the session."""
    class Broken:
        def mask(self, ids):
            raise ValueError("row 4 is not a row")

    win.annotations = Broken()
    win.refresh_annotations()
    assert "annotations unavailable" in capsys.readouterr().out


# --------------------------------------------------------------------------- color by
def test_the_panel_offers_columns_runs_and_binned_quantities_in_one_list(win):
    """One list rather than three controls: they answer the same question -- what should colour mean
    right now -- and having to know which of three places to look is the state this replaced."""
    from starplast.app import BIN_PREFIX, RUN_PREFIX
    win.keep_run(np.arange(win.n) % 4, name="test_run_a")
    sources = win.colour_sources()
    assert "compartment" in sources
    assert RUN_PREFIX + "test_run_a" in sources
    assert any(s.startswith(BIN_PREFIX) for s in sources)
    assert [win.category_box.itemText(i) for i in range(win.category_box.count())] == sources


def test_colouring_by_a_kept_run_uses_that_run(win):
    from starplast.app import RUN_PREFIX
    labels = np.where(np.arange(win.n) % 3 == 0, 0, 1)
    run = win.keep_run(labels, name="test_run_b")
    win.on_category_changed(RUN_PREFIX + run.name)
    vals = win.category_values()
    assert set(vals.unique()) == {"cluster 0", "cluster 1"}
    assert win.comp_list.count() == 2
    c = np.asarray(win.scatter.color)
    assert len({tuple(np.round(x, 3)) for x in c[:, :3]}) > 1, "one run drew one colour"


def test_two_runs_are_both_kept_and_can_be_switched_between(win):
    """The whole point: the second used to replace the first with no way back."""
    from starplast.app import RUN_PREFIX
    a = win.keep_run(np.zeros(win.n, int), name="test_two_a")
    b = win.keep_run(np.arange(win.n) % 5, name="test_two_b")
    assert {a.name, b.name} <= set(win.runs.names())
    win.on_category_changed(RUN_PREFIX + a.name)
    assert set(win.category_values().unique()) == {"cluster 0"}
    win.on_category_changed(RUN_PREFIX + b.name)
    assert len(set(win.category_values().unique())) == 5


def test_a_run_over_a_subsample_leaves_the_rest_absent_not_unclustered(win):
    """They were not put in no cluster; they were not in the map. Conflating the two would put
    thousands of genes into a category that means something else entirely."""
    from starplast.app import RUN_PREFIX
    placed = np.zeros(win.n, bool)
    placed[:300] = True
    win.placed = placed
    run = win.keep_run(np.arange(300) % 2, name="test_subsample")
    win.on_category_changed(RUN_PREFIX + run.name)
    vals = win.category_values()
    assert (vals[300:] == "").all() and set(vals[:300].unique()) == {"cluster 0", "cluster 1"}
    win.placed = None


def test_a_run_can_be_renamed_and_keeps_its_recipe(win, tmp_path):
    from starplast.app import RUN_PREFIX
    from starplast.runs import RunStore
    win.runs = RunStore(str(tmp_path))
    win.keep_run(np.zeros(win.n, int), recipe={"algorithm": "hdbscan", "min_cluster_size": 60},
                 name="test_rename_from")
    assert win.rename_run("test_rename_from", "gras_are_clean")
    assert win.runs.get("gras_are_clean").recipe["min_cluster_size"] == 60
    assert win.category_box.currentText() == RUN_PREFIX + "gras_are_clean"
    import json, os
    saved = json.load(open(str(tmp_path / "gras_are_clean.json")))
    assert saved["recipe"]["algorithm"] == "hdbscan"
    assert not os.path.exists(str(tmp_path / "test_rename_from.json"))


def test_a_name_already_in_use_is_refused(win, tmp_path):
    from starplast.runs import RunStore
    win.runs = RunStore(str(tmp_path))
    win.keep_run(np.zeros(win.n, int), name="taken")
    win.keep_run(np.zeros(win.n, int), name="other")
    assert win.rename_run("other", "taken") is False
    assert "name is taken" in win.statusBar().currentMessage()


def test_a_quantity_can_be_coloured_as_bins(win):
    """Binning makes a measurement behave like a category, which is what makes it comparable with a
    clustering -- the comparison this panel exists for."""
    from starplast.app import BIN_PREFIX
    win.bins_box.setValue(4)
    win.on_category_changed(BIN_PREFIX + "mean_plddt")
    vals = win.category_values()
    named = {v for v in vals.unique() if v}
    assert 1 < len(named) <= 4
    # The list carries the bins plus, where the column has gaps, one entry for absence -- which is
    # listed and sunk to the bottom rather than folded into the lowest bin.
    assert win.comp_list.count() == len(named) + int((vals == "").any())
    assert win.bins_box.isVisibleTo(win.bins_box.parentWidget())


def test_a_quantity_that_is_mostly_one_value_gets_fewer_bins_and_says_so(win):
    """`n_publications` is zero for most of this proteome, so its quartile edges are all zero.
    Splitting the tie by rank or by equal width would draw four colours over a column with one
    level -- a picture of a distinction that does not exist."""
    from starplast.app import BIN_PREFIX
    win.bins_box.setValue(4)
    win.on_category_changed(BIN_PREFIX + "n_publications")
    named = {v for v in win.category_values().unique() if v}
    assert len(named) < 4
    assert "share one value" in win.statusBar().currentMessage()


def test_changing_the_number_of_bins_recolours(win):
    from starplast.app import BIN_PREFIX
    win.on_category_changed(BIN_PREFIX + "mean_plddt")
    win.set_bins(3)
    three = {v for v in win.category_values().unique() if v}
    win.set_bins(8)
    eight = {v for v in win.category_values().unique() if v}
    assert len(eight) > len(three)


def test_genes_with_no_value_for_a_binned_quantity_are_absent_not_a_low_bin(win):
    from starplast.app import BIN_PREFIX
    col = "mean_plddt"
    win.on_category_changed(BIN_PREFIX + col)
    vals = win.category_values()
    missing = win.nodes[col].isna().to_numpy()
    if missing.any():
        assert (vals[missing] == "").all()
        c = np.asarray(win.scatter.color)
        import starplast.theme as TH
        assert np.allclose(c[missing][0, :3], TH.unknown_colour(win.theme)[:3], atol=1e-3)


def test_the_filter_and_the_fly_to_follow_the_chosen_source(win):
    """The list under the chooser filters what the chooser names, or the two controls describe
    different things while sitting on top of each other."""
    from starplast.app import RUN_PREFIX
    run = win.keep_run(np.where(np.arange(win.n) < 100, 0, 1), name="test_filter_run")
    win.on_category_changed(RUN_PREFIX + run.name)
    item = next(win.comp_list.item(i) for i in range(win.comp_list.count())
                if win.comp_list.item(i).data(QtCore.Qt.ItemDataRole.UserRole) == "cluster 0")
    item.setSelected(True)
    assert win.visible_mask().sum() == 100
    win.fly_to_compartment(item)
    win.comp_list.clearSelection()


def test_a_source_that_no_longer_exists_colours_nothing_rather_than_raising(win):
    win.on_category_changed("clustering: never_existed")
    assert (win.category_values() == "").all()
    win.on_category_changed("binned: not_a_column")
    assert (win.category_values() == "").all()
    win.on_category_changed("compartment")


# --------------------------------------------------------------------------- the cell diagram
def test_the_diagram_shows_for_a_localisation_category_and_hides_otherwise(win):
    """There is no sensible mapping from cell-cycle phase onto organelles, and colouring them by one
    would be a picture of a relationship that does not exist."""
    if win.diagram is None:
        pytest.skip("the artwork is not present in this checkout")
    win.on_category_changed("compartment")
    assert not win.diagram.isHidden()
    other = next(c for c in win.categories if c not in ("compartment", "compartment_best"))
    win.on_category_changed(other)
    assert win.diagram.isHidden()
    win.on_category_changed("compartment")


def test_the_diagram_and_the_map_are_coloured_from_the_same_dict(win):
    """Two palettes are two claims about what a colour means, and one of them will drift."""
    if win.diagram is None:
        pytest.skip("the artwork is not present in this checkout")
    from starplast.celldiagram import COMPARTMENT_SL
    win.on_category_changed("compartment")
    name = next(c for c in COMPARTMENT_SL if c in win.colour_of)
    win.select_compartment(name)
    fills = win.diagram.fills()
    assert list(fills) == [COMPARTMENT_SL[name]], "the selection is not the one thing coloured"
    assert np.allclose(fills[COMPARTMENT_SL[name]], win.colour_of[name], atol=1e-6)
    win.comp_list.clearSelection()
    win._refresh_diagram()
    assert win.diagram.fills() == {}, "with nothing selected the drawing is grey"


def test_selecting_in_the_list_colours_the_shape_it_shares(win):
    if win.diagram is None:
        pytest.skip("the artwork is not present in this checkout")
    win.on_category_changed("compartment")
    for name in ("rhoptries 1", "rhoptries 2"):
        item = next((win.comp_list.item(i) for i in range(win.comp_list.count())
                     if win.comp_list.item(i).data(QtCore.Qt.ItemDataRole.UserRole) == name), None)
        if item is None:
            pytest.skip("this cache has no rhoptry classes")
        win.comp_list.clearSelection()
        item.setSelected(True)
        assert np.allclose(win.diagram.fills()["SL0233"], win.colour_of[name], atol=1e-6)
        assert name in win.diagram_note.text()
    win.comp_list.clearSelection()


def test_clicking_the_diagram_selects_in_the_list(win):
    """Both directions, or the diagram is decoration."""
    if win.diagram is None:
        pytest.skip("the artwork is not present in this checkout")
    win.on_category_changed("compartment")
    win.select_compartment("apicoplast")
    chosen = [i.data(QtCore.Qt.ItemDataRole.UserRole) for i in win.comp_list.selectedItems()]
    assert chosen == ["apicoplast"]
    assert win.visible_mask().sum() < win.n, "selecting in the list must also filter the map"
    win.select_compartment("not a compartment")
    assert "not in this list" in win.statusBar().currentMessage()
    win.comp_list.clearSelection()


def test_the_compartments_with_no_organelle_are_named_beside_the_diagram(win):
    if win.diagram is None:
        pytest.skip("the artwork is not present in this checkout")
    win.on_category_changed("compartment")
    win.comp_list.clearSelection()
    note = win.diagram_note.text()
    assert "no organelle in this drawing" in note
    assert "unassigned" not in note, "absence is not a compartment that is missing a shape"


def test_changing_theme_recolours_the_diagram_with_the_map(win):
    if win.diagram is None:
        pytest.skip("the artwork is not present in this checkout")
    from starplast.celldiagram import COMPARTMENT_SL
    win.on_category_changed("compartment")
    name = next(c for c in COMPARTMENT_SL if c in win.colour_of)
    win.select_compartment(name)
    before = dict(win.diagram.fills())
    win.apply_theme("paper")
    win.select_compartment(name)
    after = dict(win.diagram.fills())
    assert before != after, "the diagram kept the old theme's colours"
    sl = COMPARTMENT_SL[name]
    assert np.allclose(after[sl], win.colour_of[name], atol=1e-6)
    win.apply_theme("dark")
    win.comp_list.clearSelection()
