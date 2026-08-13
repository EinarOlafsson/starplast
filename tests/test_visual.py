#!/usr/bin/env python3
"""The visual layer, tested on values rather than on appearance.

Rendering bugs in this project have all passed array-level tests before: additive blending summed 8,140
points to white while every color array was correct. So these tests assert the things that *can* be
asserted numerically -- monotonicity, agreement between subsets, degenerate inputs -- and the appearance
claims they support are the ones that were also checked by rendering the app and looking at it.
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


# --------------------------------------------------------------------------- depth cueing
def test_depth_t_is_zero_at_the_nearest_point_and_one_at_the_farthest():
    xyz = np.array([[0, 0, 0], [0, 0, 5], [0, 0, 10]], float)
    t = TH.depth_t(xyz, (0, 0, 0))
    assert t[0] == pytest.approx(0.0)
    assert t[2] == pytest.approx(1.0)
    assert t[1] == pytest.approx(0.5)


def test_depth_cue_fades_and_shrinks_monotonically_with_distance():
    xyz = np.array([[0, 0, float(z)] for z in range(10)])
    a, s = TH.depth_cue(xyz, (0, 0, 0), np.full(10, 0.9), np.full(10, 6.0))
    assert np.all(np.diff(a) < 0), "alpha must decrease with distance"
    assert np.all(np.diff(s) < 0), "size must decrease with distance"
    assert a[-1] == pytest.approx(0.9 * TH.DEPTH_FADE)
    assert s[-1] == pytest.approx(6.0 * TH.DEPTH_SHRINK)


def test_equidistant_points_do_not_divide_by_zero():
    """A camera on the axis of a ring, or a single point: the range is zero and the cue is off."""
    ring = np.array([[np.cos(a), np.sin(a), 0.0] for a in np.linspace(0, 2 * np.pi, 12)])
    t = TH.depth_t(ring, (0, 0, 3.0))
    assert np.isfinite(t).all()
    assert t == pytest.approx(np.zeros(len(ring)))

    a, s = TH.depth_cue(np.zeros((1, 3)), (0, 0, 1.0), np.array([0.8]), np.array([5.0]))
    assert a[0] == pytest.approx(0.8) and s[0] == pytest.approx(5.0)


def test_a_subset_normalized_against_the_full_range_agrees_with_the_full_set():
    """Edges are cued against the whole cloud's range for exactly this reason.

    Normalize a subset against its own extent and an edge fades by a different amount than the point it
    touches, at the same place on screen -- which reads as flicker while orbiting.
    """
    rng_state = np.random.default_rng(0)
    xyz = rng_state.normal(size=(200, 3)) * 10
    cam = (30.0, 0.0, 0.0)
    d = np.linalg.norm(xyz - np.asarray(cam), axis=1)
    full_range = (float(d.min()), float(d.max()))

    full = TH.depth_t(xyz, cam, full_range)
    subset_idx = np.arange(0, 200, 7)
    shared = TH.depth_t(xyz[subset_idx], cam, full_range)
    own = TH.depth_t(xyz[subset_idx], cam)          # the wrong way, kept as the contrast

    assert shared == pytest.approx(full[subset_idx])
    assert not np.allclose(own, full[subset_idx]), "the bug this guards against must be reproducible"


# --------------------------------------------------------------------------- theme
@pytest.mark.parametrize("theme,light", [("dark", False), ("light", True),
                                         ("slate", False), ("paper", True)])
def test_ground_lightness_is_detected_per_theme(theme, light):
    """The halo brightens on dark grounds and darkens on light ones; this is what decides which."""
    assert TH.is_light(theme) is light
    assert TH.is_light(TH.palette_for(theme)) is light


def test_every_theme_defines_every_role():
    roles = set(TH.palette_for("dark"))
    for t in TH.THEMES:
        assert set(TH.palette_for(t)) == roles, f"{t} does not define the same roles as dark"


# --------------------------------------------------------------------------- edge ink budget
def _ink(n, target=None):
    from starplast.app import EDGE_INK_TARGET
    return float(np.clip(np.sqrt((target or EDGE_INK_TARGET) / max(n, 1)), 0.12, 1.0))


def test_a_small_neighbourhood_is_drawn_at_full_strength():
    """Selecting one gene shows a handful of edges. Those must not be dimmed."""
    assert _ink(1) == 1.0
    assert _ink(50) == 1.0


def test_ink_falls_as_edge_count_rises_so_dense_types_stay_transparent():
    from starplast.app import EDGE_INK_TARGET
    assert _ink(EDGE_INK_TARGET) == pytest.approx(1.0)
    assert _ink(20000) < 0.4, "20,000 edges at full alpha render as an opaque sheet"
    assert _ink(20000) > _ink(60000), "denser must be fainter"


def test_ink_has_a_floor_so_edges_never_vanish_completely():
    """Truncating to nothing would be a silent lie about the data; fading to a floor is not."""
    assert _ink(10 ** 9) == pytest.approx(0.12)


# --------------------------------------------------------------------------- the app's visual state
@pytest.fixture(scope="module")
def win():
    from PyQt6 import QtWidgets
    from starplast.app import Window
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    w = Window()
    yield w
    w.close()


def test_points_carry_a_depth_gradient_when_cueing_is_on(win):
    win.depth_cue = True
    win.redraw()
    alpha = np.asarray(win.scatter.color)[:, 3]
    assert alpha.max() - alpha.min() > 0.1, "the cloud should have a visible front and back"


def test_turning_depth_cueing_off_restores_flat_alpha(win):
    """It is a toggle because the fade changes apparent color with position, and comparing two
    points' colors exactly requires it off."""
    win.depth_cue = False
    win.redraw()
    alpha = np.asarray(win.scatter.color)[:, 3]
    win.depth_cue = True
    assert alpha.max() - alpha.min() < 1e-6


def test_the_halo_appears_on_selection_and_is_removed_on_clearing(win):
    win.on_pick(100)
    assert win.halo_item is not None
    assert len(win.halo_item.pos) == 3, "three concentric rings make the soft edge"
    win.sel = None
    win.redraw()
    assert win.halo_item is None


def test_the_horizon_grid_sits_below_the_cloud_not_through_it(win):
    """A grid at the origin cuts through the points, because the embedding is not centred on zero."""
    win.show_ground = True
    win.redraw()
    assert win.grid_item is not None
    z = float(win.grid_item.transform().matrix()[2][3])
    assert z < float(win.xyz[:, 2].min())


def test_the_grid_can_be_turned_off(win):
    win.show_ground = False
    win.redraw()
    assert win.grid_item is None
    win.show_ground = True


def test_changing_level_of_detail_retargets_the_camera(win):
    """The tiers are a change of scale, so the camera moves; it eases rather than cutting."""
    win.set_level(2)
    win.view.setCameraPosition(distance=win.view.data_radius())
    win.on_level_changed()
    assert win.view._cam_timer is not None and win.view._cam_timer.isActive()
    win.view._cam_timer.stop()
