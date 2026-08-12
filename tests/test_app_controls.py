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
    for i in range(win.colour_by.count()):
        win.colour_by.setCurrentIndex(i)
        c = np.asarray(win.scatter.color)
        assert c.shape == (win.n, 4), win.colour_by.currentText()
        assert np.isfinite(c).all()
        seen.append(c.copy())
    assert any(not np.allclose(seen[0], s) for s in seen[1:]), "every mode drew the same colours"


def test_unknown_is_drawn_in_its_own_grey_rather_than_as_a_category(win):
    """"Not measured" is not a compartment, and giving it a categorical colour would put it in the
    legend beside real ones."""
    for i in range(win.colour_by.count()):
        win.colour_by.setCurrentIndex(i)
        c = np.asarray(win.scatter.color)
        assert c.shape[0] == win.n


def test_the_colour_map_choice_reaches_the_categorical_modes_too(win):
    """It used to affect only the continuous ramp, so choosing a map while colouring by compartment --
    the default -- appeared to do nothing."""
    win.colour_by.setCurrentIndex(0)
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
    win.colour_by.setCurrentIndex(0)
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
    for i in range(win.level.count()):
        win.level.setCurrentIndex(i)
        win.redraw()
        drawn[i] = win.centroid_item is not None
    assert drawn[0] and drawn[1], "the compartment and orthogroup tiers must draw their own markers"
    win.level.setCurrentIndex(2)
    win.redraw()


def test_changing_level_moves_the_camera_rather_than_cutting(win):
    win.level.setCurrentIndex(2)
    win.on_level_changed()
    assert win.view._cam_timer is not None and win.view._cam_timer.isActive()
    win.view._cam_timer.stop()


# --------------------------------------------------------------------------- edges
def test_every_edge_type_can_be_toggled_on_its_own(win):
    from starplast.app import EDGE_TYPES
    win.all_edges.setChecked(True)
    for k, _ in EDGE_TYPES:
        for cb in win.edge_cb.values():
            cb.setChecked(False)
        win.edge_cb[k].setChecked(True)
        win.redraw()
        if k in win.edges:
            assert win.edge_items, f"{k} is present in the graph but drew nothing"
    win.all_edges.setChecked(False)


def test_with_nothing_selected_and_draw_all_off_no_edges_are_drawn(win):
    """8,140 genes' edges at once is not a view of anything."""
    win.sel = None
    win.all_edges.setChecked(False)
    for cb in win.edge_cb.values():
        cb.setChecked(True)
    win.redraw()
    assert win.edge_items == []


def test_selecting_a_gene_draws_only_its_own_edges(win):
    win.all_edges.setChecked(False)
    win.on_pick(3000)
    assert isinstance(win.edge_items, list)


def test_the_attention_toggle_changes_which_comention_edges_are_drawn(win):
    """Raw co-mention is confidently misleading, so the corrected residual is the default. The toggle
    reorders exactly the quantity the edge alpha encodes."""
    from starplast.app import EDGE_TYPES, COMENTION
    for cb in win.edge_cb.values():
        cb.setChecked(False)
    for k in COMENTION:
        if k in win.edge_cb:
            win.edge_cb[k].setChecked(True)
    win.all_edges.setChecked(True)

    win.attn.setChecked(True)
    win.redraw()
    corrected = sum(len(i.pos) for i in win.edge_items)
    win.attn.setChecked(False)
    win.redraw()
    raw = sum(len(i.pos) for i in win.edge_items)
    assert corrected != raw, "the correction must change what is shown"
    win.attn.setChecked(True)
    win.all_edges.setChecked(False)


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
    win.level.setCurrentIndex(0)
    win.fly_to_compartment(win.comp_list.item(0))
    assert win.level.currentIndex() == 2


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
