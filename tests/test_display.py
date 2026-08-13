#!/usr/bin/env python3
"""The drifting background and the moving lights.

Both are scenery, and scenery in a data display has one hard rule: it may change how the map LOOKS
and must not change what it says. The lighting tests hold that line -- the flat colours survive
untouched under the shading, so nothing here can quietly darken a map into a different picture.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import ambient as A, lighting as L  # noqa: E402


# --------------------------------------------------------------------------- the blob field
def test_the_same_seed_gives_the_same_field():
    assert [b.x for b in A.field(7)] == [b.x for b in A.field(7)]
    assert [b.x for b in A.field(7)] != [b.x for b in A.field(8)]


def test_density_adds_blobs_rather_than_redrawing_them():
    """The pool is rolled at the top of the range and a prefix painted, so blob 3 keeps the numbers
    blob 3 always had. Otherwise nudging the density would rearrange the whole field."""
    few, many = A.field(0, 0.5), A.field(0, 2.0)
    assert len(few) < len(many)
    assert [b.x for b in few] == [b.x for b in many[:len(few)]]


def test_density_is_clamped_to_its_range():
    assert 1 <= len(A.field(0, -5.0)) <= len(A.field(0, 500.0))
    assert len(A.field(0, 500.0)) == len(A.field(0, A.DENSITY_RANGE[1]))


def test_position_is_a_function_of_the_clock():
    """Not a random walk: any frame can be drawn at any time, a dropped frame costs nothing, and the
    same t gives the same picture on any machine."""
    blobs = A.field(1)
    assert A.geometry(blobs, 12.5, 400, 300) == A.geometry(blobs, 12.5, 400, 300)
    assert A.geometry(blobs, 0.0, 400, 300) != A.geometry(blobs, 9.0, 400, 300)


def test_blobs_stay_a_sane_size_however_the_controls_are_set():
    blobs = A.field(2)
    for speed, size in ((0.0, 0.0), (99.0, 99.0), (1.0, 1.0)):
        for cx, cy, r, ci in A.geometry(blobs, 3.0, 400, 300, speed, size):
            assert r >= 1.0 and np.isfinite([cx, cy, r]).all()


def test_the_widget_paints_and_is_scenery(qapp):
    """Transparent to the mouse: a background that swallowed a click on the panel in front of it
    would be a bug nobody could describe."""
    from PyQt6 import QtCore
    w = A.AmbientWidget()
    w.resize(300, 200)
    w.set_time(4.0)
    assert not w.grab().isNull()
    assert w.testAttribute(QtCore.Qt.WidgetAttribute.WA_TransparentForMouseEvents)


def test_the_timer_runs_only_while_it_is_visible(qapp):
    """Sixty repaints a second of a picture nobody is looking at."""
    w = A.AmbientWidget()
    w.resize(100, 100)
    assert not w._timer.isActive()
    w.show()
    assert w._timer.isActive()
    w.hide()
    assert not w._timer.isActive()


def test_configuring_it_rebuilds_only_when_the_count_changed(qapp):
    w = A.AmbientWidget()
    before = list(w.blobs)
    w.configure(speed=2.0, size=1.5)
    assert w.blobs is before or [b.x for b in w.blobs] == [b.x for b in before]
    w.configure(density=2.5)
    assert len(w.blobs) > len(before)
    w.configure(colors=["#ff0000"])
    assert w.colors[0].name() == "#ff0000"
    # `size` under spaCR's name, stored where it cannot shadow QWidget.size().
    w.configure(size=2.0)
    assert w.blob_size == 2.0 and callable(w.size)


def test_the_settings_do_not_shadow_the_widget(qapp):
    """An attribute called `size` overrides QWidget.size(), and every caller asking this widget how
    big it is gets a float and 'float object is not callable' from somewhere that never mentioned
    size. Found exactly that way."""
    from PyQt6 import QtCore
    w = A.AmbientWidget()
    w.resize(120, 90)
    assert w.size() == QtCore.QSize(120, 90)
    assert isinstance(w.blob_size, float)


def test_a_tick_advances_the_clock(qapp):
    w = A.AmbientWidget()
    w.resize(80, 60)
    w.set_time(0.0)
    before = A.geometry(w.blobs, 0.0, 80, 60)
    w._tick()
    assert A.geometry(w.blobs, w._t, 80, 60) != before


# --------------------------------------------------------------------------- the lights
def test_lights_move_and_do_not_march_in_step():
    """Three lights orbiting in formation read as one light."""
    a = L.lights(0.0, 3)
    b = L.lights(6.0, 3)
    assert len(a) == 3
    assert not np.allclose(a[0]["pos"], b[0]["pos"])
    spread = np.linalg.norm(a[0]["pos"] - a[1]["pos"])
    assert spread > 1e-6


def test_the_number_of_lights_is_clamped():
    assert len(L.lights(0.0, 0)) == L.LIGHT_RANGE[0]
    assert len(L.lights(0.0, 99)) == L.LIGHT_RANGE[1]


def test_a_point_facing_a_light_is_brighter_than_one_facing_away():
    """The whole point of the feature: direction, on a cloud that has no surfaces of its own."""
    rng = np.random.default_rng(0)
    xyz = rng.normal(size=(500, 3)) * 20
    flat = np.full((500, 4), 0.5)
    lit = [{"pos": np.array([200.0, 0.0, 0.0]), "color": np.ones(3)}]
    out = L.shade(xyz, flat, lit)
    facing = xyz[:, 0] > np.percentile(xyz[:, 0], 90)
    away = xyz[:, 0] < np.percentile(xyz[:, 0], 10)
    assert out[facing, :3].mean() > out[away, :3].mean()


def test_nothing_is_ever_lit_to_invisibility():
    """An unlit half that went black would hide half the genes, and this is a data display before
    it is a rendering."""
    rng = np.random.default_rng(1)
    xyz = rng.normal(size=(300, 3)) * 20
    flat = np.full((300, 4), 0.6)
    out = L.shade(xyz, flat, L.lights(0.0, 3))
    assert out[:, :3].min() >= 0.6 * L.AMBIENT * 0.99


def test_transparency_is_never_touched():
    """Alpha carries the class filter and the depth cue. A lighting model that wrote it would hide
    genes the user chose to see."""
    rng = np.random.default_rng(2)
    xyz = rng.normal(size=(200, 3))
    flat = np.column_stack([np.full(200, 0.5)] * 3 + [rng.random(200)])
    out = L.shade(xyz, flat, L.lights(1.0))
    assert np.array_equal(out[:, 3], flat[:, 3])


def test_specular_adds_light_rather_than_replacing_it():
    rng = np.random.default_rng(3)
    xyz = rng.normal(size=(400, 3)) * 10
    flat = np.full((400, 4), 0.4)
    lit = L.lights(0.0, 2)
    assert L.shade(xyz, flat, lit, specular=True)[:, :3].sum() >= \
        L.shade(xyz, flat, lit, specular=False)[:, :3].sum()


def test_an_empty_map_shades_to_nothing_rather_than_raising():
    assert L.shade(np.zeros((0, 3)), np.zeros((0, 4)), L.lights(0.0)).shape == (0, 4)


# --------------------------------------------------------------------------- in the window
def test_the_display_tab_exists_with_both_controls(win):
    d = win.build_preferences()
    assert [win.pref_tabs.tabText(i) for i in range(win.pref_tabs.count())] == \
        ["Appearance", "Display"]
    assert win.ambient_box.currentText() in ("none", "blobs")
    assert win.light_box.currentText() in L.MODES
    d.close()


def test_turning_the_background_on_puts_it_behind_the_panel(win):
    assert win.set_ambient("blobs") == "blobs"
    assert win._ambient_widget is not None
    assert win._ambient_widget.parent() is win.left_panel
    win.set_ambient_option("density", 2.0)
    assert win._ambient["density"] == 2.0
    assert win.set_ambient("none") == "none"
    assert not win._ambient_widget.isVisible()


def test_an_unknown_background_is_refused_rather_than_drawn(win):
    assert win.set_ambient("aurora") == "none"


def test_lighting_runs_a_timer_only_while_it_is_on(win):
    assert win.set_lighting("lit") == "lit"
    assert win._light_timer.isActive()
    assert win.set_lighting("off") == "off"
    assert not win._light_timer.isActive()
    assert win.set_lighting("nonsense") == "off"


def test_shading_never_compounds(win):
    """Every frame shades the FLAT colours, not what is on screen. Shading an already-shaded array
    darkens it a little more each frame until the map goes black -- which looks like a slow fade and
    reads as the data changing."""
    win.set_lighting("lit")
    win.redraw()
    first = np.array(win._base_colors, copy=True)
    win._light_tick()
    win._light_tick()
    win._light_tick()
    assert np.array_equal(win._base_colors, first), "the flat colours were shaded in place"
    win.set_lighting("off")


def test_the_lights_orbit_outside_the_cloud(win):
    """Inside it, half the map is behind every light and the shading says nothing about shape."""
    reach = float(np.abs(win.xyz).max())
    assert win._light_radius() > reach


# --------------------------------------------------------------------------- menus and windows
def test_preferences_is_in_the_file_menu_and_about_is_in_help(win):
    from PyQt6 import QtWidgets
    labels = {m.title(): [a.text() for a in m.actions()]
              for m in win.menuBar().findChildren(QtWidgets.QMenu)}
    assert "Preferences…" in labels["&File"]
    assert "Preferences…" not in labels.get("&View", [])
    assert "About starplast" in labels["&Help"]


def test_preferences_is_a_window_of_its_own(win):
    """Modal, it froze the map behind it, so every setting had to be judged from memory of what the
    map looked like a moment ago."""
    from PyQt6 import QtCore
    d = win.open_preferences()
    assert not d.isModal()
    assert d.windowFlags() & QtCore.Qt.WindowType.Window
    assert win.open_preferences() is d, "a second call built a second dialog"
    d.close()


def test_about_says_what_is_running_rather_than_what_was_written(win):
    from starplast import __version__
    text = win.about_text()
    assert __version__ in text
    assert f"{len(win.nodes):,}" in text
    assert "hyperLOPIT" in text and "SwissBioPics" in text
    assert "never evidence on its own" in text


@pytest.fixture(scope="module")
def qapp():
    from PyQt6 import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture(scope="module")
def win(qapp):
    from starplast.app import Window
    w = Window()
    w.resize(1200, 800)
    yield w
    w.close()


def test_the_settings_are_remembered_across_a_restart(qapp, tmp_path, monkeypatch):
    """A background someone turned on, off again on the next launch, is a background they will not
    turn on twice."""
    from PyQt6 import QtCore
    from starplast import paths
    from starplast.app import Window
    monkeypatch.setenv(paths.ENV_STATE, str(tmp_path))
    s = QtCore.QSettings("starplast", "starplast")
    s.setValue("display/ambient", "blobs")
    s.setValue("display/ambient_density", 2.0)
    s.setValue("display/lighting", "lit")
    s.setValue("display/light_lights", 5)
    s.sync()
    try:
        w = Window()
        assert w._ambient_mode == "blobs" and w._ambient_widget is not None
        assert w._ambient["density"] == 2.0
        assert w._lighting["mode"] == "lit" and w._lighting["lights"] == 5
        assert w._light_timer.isActive(), "the lights were remembered but not started"
        w.set_lighting("off")
        w.close()
    finally:
        for key in ("display/ambient", "display/ambient_density", "display/lighting",
                    "display/light_lights"):
            s.remove(key)
        s.sync()


def test_the_background_follows_the_panel_when_it_is_resized(win):
    """It is a child of the panel, not laid out by it: without this it stays the size it was born."""
    from PyQt6 import QtCore, QtGui
    win.set_ambient("blobs")
    win.left_panel.resize(300, 500)
    win.eventFilter(win.left_panel, QtGui.QResizeEvent(QtCore.QSize(300, 500),
                                                       QtCore.QSize(200, 400)))
    assert win._ambient_widget.size() == win.left_panel.size()
    win.set_ambient("none")


def test_light_options_are_remembered(win):
    from PyQt6 import QtCore
    win.set_lighting_option("lights", 4)
    win.set_lighting_option("speed", 1.25)
    assert win._lighting["lights"] == 4 and win._lighting["speed"] == 1.25
    s = QtCore.QSettings("starplast", "starplast")
    assert int(s.value("display/light_lights")) == 4


def test_a_tick_with_the_lights_off_does_nothing(win):
    """The timer is stopped, but a queued tick can still arrive after it."""
    win.set_lighting("off")
    win._light_tick()


def test_about_opens(win, monkeypatch):
    from PyQt6 import QtWidgets
    shown = []
    monkeypatch.setattr(QtWidgets.QMessageBox, "about",
                        staticmethod(lambda *a: shown.append(a[-1])))
    text = win.about()
    assert shown and "starplast" in shown[0] and text == shown[0]
