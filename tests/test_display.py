#!/usr/bin/env python3
"""The drifting background and the moving lights.

Both are scenery, and scenery in a data display has one hard rule: it may change how the map LOOKS
and must not change what it says. The lighting tests hold that line -- the flat colours survive
untouched under the shading, so nothing here can quietly darken a map into a different picture.
"""
from __future__ import annotations

import os
import re
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import ambient as A, app as AW, lighting as L  # noqa: E402


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
    # A colored fill may lower one channel slightly, but never enough to hide the point.
    assert out[:, :3].min() >= 0.6 * L.AMBIENT * 0.90


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
        ["Appearance", "Display", "Window"]
    assert win.ambient_box.currentText() in ("none", "blobs")
    assert win.light_box.currentText() in L.MODES
    d.close()


# --------------------------------------------------------------------------- window size
def test_the_window_can_be_made_to_fit_a_1080p_screen(win):
    """The regression test for the bug that started this: the window would not fit a 1080p monitor.

    Not a size preference at all -- a MINIMUM. The Data tab's `feature blocks` group is one checkbox
    per slot block, so the tab's minimum height was 2,736 px; a QDockWidget passes its widget's
    minimum straight through, which made the window's own minimum 2,897. A window cannot be resized
    below its minimum, so on a 1920x1080 screen it opened with the bottom off the display and no
    amount of dragging brought it back. Nothing reported it, because a minimum size is not an error.

    The bound is the work area of a 1080p screen after a taskbar and a title bar. Asserting the
    minimum rather than the current size is the point: a size that merely happens to fit today would
    pass while the window was still unable to shrink.
    """
    assert win.minimumSizeHint().height() <= 1000, "the window cannot fit a 1080p screen"
    assert win.minimumSizeHint().width() <= 1900
    # And it must actually take the size, not snap back to a minimum.
    win.resize(1280, 720)
    assert (win.width(), win.height()) == (1280, 720)


def test_every_analysis_tab_scrolls_rather_than_forcing_the_window_taller(win):
    """The fix, stated where it can regress: a tab added later without a scroll area would put the
    window's minimum back where it was."""
    from PyQt6 import QtWidgets
    tabs = win.panel.findChild(QtWidgets.QTabWidget)
    # The count is asserted so the loop below cannot pass by iterating over nothing. It moves when a
    # tab is added -- 8 since 2026-08-18, when Questions was added.
    assert tabs.count() == 8
    for i in range(tabs.count()):
        page = tabs.widget(i)
        assert isinstance(page, QtWidgets.QScrollArea), f"tab {tabs.tabText(i)} does not scroll"
        assert page.widgetResizable(), f"tab {tabs.tabText(i)} would not reflow"


def test_the_default_is_the_screen_and_not_full_screen(win, monkeypatch):
    """What a fresh install does: fill the monitor it opens on, windowed."""
    win.settings().remove("window/size")
    win.settings().remove("window/fullscreen")
    cfg = win.window_settings()
    assert cfg == {"size": AW.MATCH_SCREEN, "fullscreen": False}
    assert not win.isFullScreen()


def test_a_size_larger_than_the_monitor_is_clamped_and_says_so(win):
    """Choosing 4K on a 1080p screen must give a window that fits. Silently ignoring the choice would
    read as a broken setting, so the clamp is reported."""
    note = win.apply_window_size(size="3840 × 2160", fullscreen=False)
    area = win.work_area()
    assert win.width() <= area.width() and win.height() <= area.height()
    if (area.width(), area.height()) < (3840, 2160):
        assert "clamped" in note


def test_match_screen_uses_the_work_area_not_the_raw_resolution(win):
    """The taskbar is the difference between fitting and a title bar off the top of the screen."""
    note = win.apply_window_size(size=AW.MATCH_SCREEN, fullscreen=False)
    area = win.work_area()
    assert (win.width(), win.height()) == (area.width(), area.height())
    assert f"{area.width()} × {area.height()}" in note


def test_full_screen_is_a_setting_that_survives_and_can_be_turned_off(win):
    win.apply_window_size(size="1280 × 720", fullscreen=True)
    assert win.isFullScreen()
    assert win.window_settings()["fullscreen"] is True
    win.apply_window_size(fullscreen=False)
    assert not win.isFullScreen()
    assert win.window_settings() == {"size": "1280 × 720", "fullscreen": False}


def test_a_stored_size_that_is_no_longer_offered_falls_back(win):
    """A hand-edited or stale config must not stop the window opening."""
    win.settings().setValue("window/size", "9999 x nonsense")
    assert win.window_settings()["size"] == AW.MATCH_SCREEN
    win.apply_window_size()
    assert win.width() > 0


def test_with_no_screen_at_all_the_window_still_gets_a_size(win, monkeypatch):
    """The defensive branch, exercised rather than assumed. `screen()` is None before a window is
    shown and `primaryScreen()` is None on a machine with no display attached -- and this runs during
    `__init__`, so returning None here would raise before the window exists at all."""
    from PyQt6 import QtGui
    monkeypatch.setattr(win, "screen", lambda: None)
    monkeypatch.setattr(QtGui.QGuiApplication, "primaryScreen", staticmethod(lambda: None))
    area = win.work_area()
    assert (area.width(), area.height()) == (1280, 720)
    assert win.apply_window_size(size=AW.MATCH_SCREEN, fullscreen=False) == "1280 × 720"


def test_the_controls_report_the_size_actually_used(win):
    d = win.build_preferences()
    assert win.window_size_box.currentText() in AW.WINDOW_SIZES
    note = win._on_window_change(size="1600 × 900")
    assert note in win.window_note.text()
    assert "this monitor allows up to" in win.window_note.text()
    d.close()


def test_turning_the_background_on_puts_it_behind_the_window(win):
    assert win.set_ambient("blobs") == "blobs"
    assert win._ambient_widget is not None
    assert win._ambient_widget.parent() is win
    win.set_ambient_option("density", 2.0)
    assert win._ambient["density"] == 2.0
    assert win.set_ambient("none") == "none"
    assert not win._ambient_widget.isVisible()


def test_an_unknown_background_is_refused_rather_than_drawn(win):
    assert win.set_ambient("aurora") == "none"


def test_lighting_runs_a_timer_only_while_it_is_on(win):
    assert win.set_lighting("soft") == "soft"
    assert win._light_timer.isActive()
    assert win.set_lighting("ray traced") == "ray traced"
    assert win._light_timer.isActive()
    assert win.set_lighting("off") == "off"
    assert not win._light_timer.isActive()
    assert win.set_lighting("nonsense") == "off"


def test_shading_never_compounds(win):
    """Every frame shades the FLAT colours, not what is on screen. Shading an already-shaded array
    darkens it a little more each frame until the map goes black -- which looks like a slow fade and
    reads as the data changing."""
    win.set_lighting("soft")
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
    assert "balanced feature recipe" in text and "SwissBioPics" in text
    assert "never evidence on its own" in text


@pytest.fixture(scope="module")
def qapp():
    from PyQt6 import QtWidgets
    return QtWidgets.QApplication.instance() or QtWidgets.QApplication([])


@pytest.fixture(scope="module")
def win(qapp):
    """The window, and the application font put back afterwards.

    The text-size tests move the font on the QApplication -- which is shared with every other test
    module in the run -- so leaving it scaled would resize widgets belonging to files that never
    asked for it, in whatever order the suite happens to run them.
    """
    from PyQt6 import QtGui
    from starplast.app import Window
    before = QtGui.QFont(qapp.font())
    w = Window()
    w.resize(1200, 800)
    yield w
    w.close()
    qapp.setFont(before)


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
    s.setValue("display/light_finish", "metallic")
    s.setValue("display/light_source", "mouse")
    s.sync()
    try:
        w = Window()
        assert w._ambient_mode == "blobs" and w._ambient_widget is not None
        assert w._ambient["density"] == 2.0
        # `target_marker` is not among the stored keys above, so it takes its default -- which is
        # now `halo` rather than `none`: the wandering and bouncing lights read as broken while it
        # was `none`, because the orb they are meant to show is drawn by the marker.
        assert w._lighting == {"mode": "soft", "source": "mouse flashlight",
                               "point_mode": "metallic 3D", "mood": "neutral",
                               "pointer_mode": "broad flashlight", "response": "smooth",
                               "target_marker": "halo"}
        assert w._light_timer.isActive(), "the lights were remembered but not started"
        w.set_lighting("off")
        w.close()
    finally:
        for key in ("display/ambient", "display/ambient_density", "display/lighting",
                    "display/light_finish", "display/light_source"):
            s.remove(key)
        s.sync()


def test_the_background_follows_the_window_when_it_is_resized(win):
    """It is a child of the window, not laid out by it: without this it stays the size it was born,
    and the blobs stop where the window used to end."""
    from PyQt6 import QtCore, QtGui
    win.set_ambient("blobs")
    win.resize(900, 640)
    win.eventFilter(win, QtGui.QResizeEvent(QtCore.QSize(900, 640), QtCore.QSize(1200, 800)))
    assert win._ambient_widget.size() == win.size()
    win.set_ambient("none")


def test_light_options_are_remembered(win):
    from PyQt6 import QtCore
    win.set_lighting_option("source", "selected gene")
    win.set_lighting_option("point_mode", "flat")
    win.set_lighting_option("mood", "warm")
    assert win._lighting["source"] == "selected gene"
    assert win._lighting["point_mode"] == "flat" and win._lighting["mood"] == "warm"
    s = QtCore.QSettings("starplast", "starplast")
    assert s.value("display/light_source") == "selected gene"
    assert s.value("display/light_point_mode") == "flat"
    assert s.value("display/light_mood") == "warm"


def test_preferences_offer_only_distinct_understandable_render_controls(win):
    from starplast import lighting as light
    dialog = win.build_preferences()
    try:
        choices = lambda box: [box.itemText(i) for i in range(box.count())]
        assert choices(win.light_box) == list(light.MODES)
        assert choices(win.light_source) == list(light.SOURCES)
        assert choices(win.light_mood) == list(light.LIGHT_MOODS)
        assert choices(win.point_render) == list(light.POINT_MODES)
        for removed in ("light_rays", "light_finish", "light_count", "light_speed"):
            assert not hasattr(win, removed), f"obsolete control {removed} returned"
        assert "GPU" in win.light_box.toolTip() and "volumetric" in win.light_box.toolTip()
        assert "Vulkan" in win.light_box.toolTip() and "path tracing" in win.light_box.toolTip()
    finally:
        dialog.close()


def test_a_tick_with_the_lights_off_does_nothing(win):
    """The timer is stopped, but a queued tick can still arrive after it."""
    win.set_lighting("off")
    win._light_tick()


def test_about_opens(win):
    """In a glass window beside the map now, not a modal box; the text shown is the text returned."""
    text = win.about()
    shown = win._about
    assert shown.isVisible() and shown.windowTitle() == "About starplast"
    assert "starplast" in shown.text and text == shown.text
    shown.close()


# --------------------------------------------------------------------------- the console, quiet
def test_the_artwork_no_longer_makes_qt_complain(qapp, capfd):
    """Three warnings per render, on every repaint of every organelle mask, is a console nobody can
    read -- and a real warning would be lost in it. UniProt nests hidden `<text>` descriptions and
    `<a>` links inside `<path>` elements, where the content model does not allow them."""
    from starplast import celldiagram as CD
    capfd.readouterr()
    d = CD.CellDiagram()
    d.resize(200, 280)
    d.set_palette({"rhoptries 1": (0.2, 0.6, 1.0)}, "rhoptries 1")
    d.masks()
    d.grab()
    err = capfd.readouterr().err
    assert "Could not add child element" not in err, err


def test_stripping_the_metadata_keeps_every_organelle(qapp):
    """It removes what is not drawn. An organelle lost with it would be a compartment that stops
    being clickable, which is the whole feature."""
    from starplast import celldiagram as CD
    raw = open(CD.icon_path(), encoding="utf8").read()
    stripped = CD.strip_metadata(raw)
    assert CD.groups_in(stripped) == CD.groups_in(raw)
    assert "<text" not in stripped and "<a " not in stripped
    d = CD.CellDiagram()
    d.resize(200, 280)
    assert len(d.masks()) == 14


def test_the_credit_survives_the_stripping(qapp):
    """The creator's name and the licence ARE that metadata. Stripping it and then looking for it is
    how a CC BY attribution silently disappears."""
    from starplast import celldiagram as CD
    d = CD.CellDiagram()
    assert d.credit["creator"] and "creativecommons.org" in d.credit["license"]
    assert d.credit["creator"] in d.toolTip()


# --------------------------------------------------------------------------- tooltips and text size
def test_a_tooltip_is_a_rectangle_of_prose(qapp):
    """A fixed table width does NOT constrain a tooltip -- measured at 1,588 px when it was tried.
    `white-space: pre` with the breaks already in is what holds it."""
    from PyQt6 import QtGui
    from starplast import theme as TH
    text = ("One long sentence about why a control exists and what choosing badly costs, which on "
            "one line would run off the side of any monitor and be unreadable because the start "
            "and the end are too far apart to hold in the head at once.")
    doc = QtGui.QTextDocument()
    doc.setHtml(TH.tip(text))
    assert doc.idealWidth() < 600, f"{doc.idealWidth():.0f}px is still a strip"
    lines = [x for x in doc.toPlainText().split("\n") if x.strip()]
    assert len(lines) > 2
    widths = [len(x.rstrip()) for x in lines]
    assert max(widths) - min(widths) < 70, "ragged: not a rectangle"


def test_the_tooltip_is_on_the_setting_not_only_the_field(qapp):
    """The label is the word people hover over; the box beside it is the thing they click."""
    from PyQt6 import QtWidgets
    from starplast import theme as TH
    w = QtWidgets.QWidget()
    form = QtWidgets.QFormLayout(w)
    box = QtWidgets.QComboBox()
    box.setToolTip("Why this setting exists and what it costs to choose it badly.")
    form.addRow("theme", box)
    assert TH.wrap_tooltips(w) == 1
    assert "white-space:pre" in box.toolTip()
    assert form.labelForField(box).toolTip() == box.toolTip()
    assert TH.wrap_tooltips(w) == 0, "wrapping twice would wrap the wrapper"


def test_every_control_in_preferences_explains_itself_in_a_block(win):
    from PyQt6 import QtWidgets
    d = win.build_preferences()
    tips = [q.toolTip() for q in d.findChildren(QtWidgets.QWidget) if q.toolTip()]
    assert tips
    assert all("white-space:pre" in t for t in tips), "a one-line tooltip survived"
    d.close()


def test_text_size_scales_the_whole_interface(win):
    """On the application's font, not a stylesheet: a stylesheet font-size does not change what a
    widget reports as its size hint, so the text grew and the boxes did not -- and labels were cut
    off at the old width."""
    from PyQt6 import QtWidgets
    from starplast.app import UI_SCALE_RANGE
    app = QtWidgets.QApplication.instance()
    base = win._base_font.pointSizeF()
    win.set_ui_scale(1.4)
    assert abs(app.font().pointSizeF() - base * 1.4) < 0.6
    win.set_ui_scale(1.0)
    assert abs(app.font().pointSizeF() - base) < 0.6
    assert win.set_ui_scale(99.0) == UI_SCALE_RANGE[1]
    assert win.set_ui_scale(0.01) == UI_SCALE_RANGE[0]
    win.set_ui_scale(1.0)


def test_the_stylesheet_grows_with_the_text_size(win):
    """The application font alone was not enough, and this is the bug that was reported as "text
    size does nothing".

    A stylesheet `font-size` BEATS the application font on every widget the sheet touches -- which
    here is nearly all of them. So the setting moved the app font, the sheet's hard-coded 13px
    immediately overrode it, and nothing on screen changed. The sheet has to be rebuilt at the new
    scale as well; the app font still matters, for the size hints the test above measures."""
    from PyQt6 import QtWidgets
    app = QtWidgets.QApplication.instance()
    sizes = lambda: [int(t) for t in re.findall(r"font-size:\s*(\d+)px", app.styleSheet())]
    win.set_ui_scale(1.0)
    small = sizes()
    assert small, "the stylesheet declares no font size at all"
    win.set_ui_scale(1.5)
    grown = sizes()
    win.set_ui_scale(1.0)
    assert len(grown) == len(small)
    assert all(g > s for g, s in zip(grown, small)), f"{small} -> {grown}"


def test_the_text_size_control_is_a_slider_that_says_where_it_is(win):
    """A number box asks for a figure nobody knows in advance. This is a thing you drag until it
    looks right -- and the read-out beside it is what lets you put it back."""
    from PyQt6 import QtWidgets
    from starplast.app import UI_SCALE_RANGE
    win.build_preferences()
    assert isinstance(win.zoom_box, QtWidgets.QSlider)
    assert (win.zoom_box.minimum(), win.zoom_box.maximum()) == (
        int(UI_SCALE_RANGE[0] * 100), int(UI_SCALE_RANGE[1] * 100))
    win.zoom_box.setValue(130)
    assert abs(win._ui_scale - 1.3) < 1e-6, "dragging the slider did not change the scale"
    assert "1.30" in win.zoom_label.text()
    win.zoom_box.setValue(100)


def test_scaling_never_compounds(win):
    """Applied to the desktop's own font every time: compounding 1.2 three times is 1.7, and the
    text creeps every time the dialog is opened."""
    from PyQt6 import QtWidgets
    app = QtWidgets.QApplication.instance()
    win.set_ui_scale(1.2)
    once = app.font().pointSizeF()
    win.set_ui_scale(1.2)
    win.set_ui_scale(1.2)
    assert abs(app.font().pointSizeF() - once) < 1e-6
    win.set_ui_scale(1.0)


def test_bigger_text_makes_a_wider_label_rather_than_a_clipped_one(win):
    """The whole point of moving the application font rather than a stylesheet: a widget's size hint
    grows with its text, so the layout gives it room instead of cutting it off.

    Measured on a label, because a QListWidget's hint is a fixed 256x192 whatever is in it -- which
    is a fact about QListWidget, not about whether the scaling worked."""
    from PyQt6 import QtWidgets
    win.set_ui_scale(1.0)
    QtWidgets.QApplication.processEvents()
    # Parentless and destroyed here rather than handed to Qt's deferred deletion: this file's
    # window is module-scoped, and a child queued for deletion outlives the test that made it.
    label = QtWidgets.QLabel("nucleus - non-chromatin  (461)")
    label.setFont(QtWidgets.QApplication.instance().font())
    small = label.sizeHint().width()
    win.set_ui_scale(1.5)
    QtWidgets.QApplication.processEvents()
    label.setFont(QtWidgets.QApplication.instance().font())
    grown = label.sizeHint().width()
    win.set_ui_scale(1.0)
    del label
    assert grown > small, "the text grew and its box did not"


def test_a_font_measured_in_pixels_scales_too(win, monkeypatch):
    """Some desktops hand out a font with no point size at all -- pointSizeF is -1 and the size is
    in pixels. Scaling by a negative number would make the text vanish."""
    from PyQt6 import QtGui, QtWidgets
    pixel_font = QtGui.QFont(win._base_font)
    pixel_font.setPixelSize(12)
    monkeypatch.setattr(win, "_base_font", pixel_font)
    win.set_ui_scale(1.5)
    assert QtWidgets.QApplication.instance().font().pixelSize() == 18
    monkeypatch.undo()
    win.set_ui_scale(1.0)


def test_the_caption_under_the_cell_has_no_card_of_its_own(win):
    """A scroll area paints its own background, and a black box under the caption covers whatever
    the theme is drawing behind the panel."""
    assert "transparent" in win.diagram_note_area.styleSheet()
    assert not win.diagram_note_area.viewport().autoFillBackground()
    assert "transparent" in win.diagram_note.styleSheet()


def test_fields_are_one_dark_grey_everywhere(qapp):
    from starplast import theme as TH
    for name in TH.THEMES:
        assert TH.FIELD_GREY in TH.stylesheet(name)


def test_a_widget_destroyed_mid_walk_is_skipped(qapp):
    """findChildren hands back wrappers around C++ objects, and a widget under construction can
    destroy one while this loop runs. PyQt raises RuntimeError for the ones it can catch; without
    the guard the walk stops and the rest of the window keeps its one-line tooltips."""
    from PyQt6 import QtWidgets
    from starplast import theme as TH

    class Vanishing(QtWidgets.QLabel):
        def toolTip(self):
            raise RuntimeError("wrapped C/C++ object of type QLabel has been deleted")

    parent = QtWidgets.QWidget()
    layout = QtWidgets.QFormLayout(parent)
    gone = Vanishing("x", parent)
    layout.addRow("gone", gone)
    ok = QtWidgets.QComboBox()
    ok.setToolTip("Why this setting exists.")
    layout.addRow("fine", ok)
    assert TH.wrap_tooltips(parent) == 1, "the walk stopped at the dead widget"
    assert "white-space:pre" in ok.toolTip()


# --------------------------------------------------------------------------- finishes and sources
def test_every_point_mode_changes_the_surface_rather_than_the_data(qapp):
    """A point mode is how a point answers the light. All have to look different, and none may
    move a point or touch its alpha -- the map's shape and its filter are not a rendering choice."""
    rng = np.random.default_rng(0)
    xyz = rng.normal(size=(400, 3)) * 15
    flat = np.column_stack([np.full(400, 0.5)] * 3 + [rng.random(400)])
    lit = L.lights(0.0, 2)
    seen = {}
    for name in L.POINT_MODES:
        out = L.shade(xyz, flat, lit, point_mode=name)
        assert np.array_equal(out[:, 3], flat[:, 3]), f"{name} touched alpha"
        seen[name] = round(float(out[:, :3].mean()), 4)
    assert len(set(seen.values())) == len(seen), f"two point modes render identically: {seen}"


def test_flat_has_no_highlight_and_glossy_has_a_tight_one(qapp):
    rng = np.random.default_rng(1)
    xyz = rng.normal(size=(600, 3)) * 15
    flat = np.full((600, 4), 0.5)
    lit = L.lights(0.0, 1)
    plain = L.shade(xyz, flat, lit, point_mode="flat")[:, :3]
    glossy = L.shade(xyz, flat, lit, point_mode="glossy 3D")[:, :3]
    assert glossy.max() > plain.max(), "glossy did not add a highlight"
    # Tight means the peak stands FURTHER ABOVE THE BODY, which is the thing that reads as gloss.
    # Counting points near the top of the range instead -- which this did -- measures the shape of
    # the range rather than the highlight, and it moves whenever the floor moves, so it broke as
    # soon as the finishes started separating on contrast as well.
    stands_out = lambda a: float(a.max() / max(np.median(a), 1e-9))
    assert stands_out(glossy) > stands_out(plain), "the highlight does not rise above the body"


def test_metallic_takes_its_highlight_from_the_point_not_the_light(qapp):
    """The difference between a copper bead and a white-glinting plastic one."""
    xyz = np.array([[0.0, 0.0, 10.0], [0.0, 0.0, -10.0]])
    red = np.array([[1.0, 0.0, 0.0, 1.0], [1.0, 0.0, 0.0, 1.0]])
    lit = [{"pos": np.array([0.0, 0.0, 60.0]), "color": np.array([1.0, 1.0, 1.0])}]
    metal = L.shade(xyz, red, lit, point_mode="metallic 3D")
    plastic = L.shade(xyz, red, lit, point_mode="glossy 3D")
    # A white light on a red point: the plastic highlight whitens the green channel, the metal's
    # stays red. This only works because specular is ADDED rather than multiplied into the colour --
    # multiplied, a pure red point has no green to raise and every finish looks identical on it.
    assert plastic[0, 1] > metal[0, 1] + 0.05, (plastic[0], metal[0])
    assert metal[0, 0] > metal[0, 1], "the metal highlight lost the point's own colour"


def test_light_moods_are_visibly_different_without_touching_alpha(qapp):
    rng = np.random.default_rng(22)
    xyz = rng.normal(size=(800, 3)) * 15
    colors = np.tile(np.array([[0.75, 0.22, 0.12, 0.8], [0.12, 0.35, 0.80, 0.6]]), (400, 1))
    rendered = {}
    for mood in L.LIGHT_MOODS:
        lit = L.fixed((0.4, 0.6, 1.0), 40.0, color=L.mood_color(mood))
        rendered[mood] = L.shade(xyz, colors, lit, point_mode="glossy 3D")
        assert np.array_equal(rendered[mood][:, 3], colors[:, 3])
    for a, b in (("neutral", "cool blue"), ("neutral", "warm"), ("cool blue", "warm")):
        steps = np.abs(rendered[a][:, :3] - rendered[b][:, :3]).max(axis=1) * 255
        assert np.median(steps) > 8, f"{a}/{b} mood difference is invisible"


def test_3d_points_are_large_enough_for_the_sphere_texture_to_survive_downsampling(win):
    old = win._lighting["point_mode"]
    try:
        win.set_lighting("off")
        win.set_lighting_option("point_mode", "flat")
        win.redraw()
        flat = np.asarray(win.scatter.size, dtype=float)
        win.set_lighting_option("point_mode", "glossy 3D")
        win.redraw()
        glossy = np.asarray(win.scatter.size, dtype=float)
        assert np.median(glossy) == pytest.approx(np.median(flat) * 1.30)
    finally:
        win.set_lighting_option("point_mode", old)


@pytest.mark.parametrize("source", L.SOURCES)
def test_every_source_produces_at_least_one_light(qapp, source):
    xyz = np.random.default_rng(2).normal(size=(50, 3))
    lit = L.light_at(xyz, source, 1.0, 3, 0.3, 30.0, pointer=(0.2, 0.4, 1.0), selected=3,
                     neighbours=[1, 2])
    assert lit and all("pos" in x and "color" in x for x in lit)


def test_a_source_with_nothing_to_follow_lights_from_the_front(qapp):
    """No pointer in the view yet, nothing selected. Guessing would put the light behind the map.

    The travelling lights are exempt, and not by oversight: they follow nothing, so there is nothing
    for them to be missing. They light from wherever they have reached inside the cloud, and they
    return as many lights as were asked for rather than one. Both have their own tests above.
    """
    xyz = np.random.default_rng(3).normal(size=(20, 3))
    travelling = {"wandering light", "bouncing light"}
    for source in L.SOURCES:
        lit = L.light_at(xyz, source, 0.0, 3, 0.3, 20.0)
        if source in travelling:
            assert len(lit) == 3, f"{source} ignored the light count"
            continue
        assert len(lit) == 1 and lit[0]["pos"][2] > 0


def test_the_orbiting_lights_stay_on_a_shell_outside_the_cloud(qapp):
    """Reported as "the orbit looks strange -- it lights things up and it moves, but it looks
    strange". The path was a Lissajous taken as a POSITION, which ranges over the whole cube it is
    inscribed in: the lights swung between 0.2 and 1.7 times the radius and dived through the cloud
    on about 4% of frames. The map pulsed, and every so often a light surfaced among the genes.

    A light that changes direction while keeping its distance is the thing that was wanted, and the
    guard has to be on the distance rather than on "it moved" -- the broken version moved too."""
    radius, seen = 90.0, []
    for t in np.arange(0.0, 60.0, 0.25):
        for light in L.lights(float(t), 3, 0.35, radius=radius):
            seen.append(float(np.linalg.norm(light["pos"])))
    seen = np.array(seen)
    assert np.allclose(seen, radius), f"distance ranged {seen.min():.1f}-{seen.max():.1f}"
    # And still moving: a fixed distance must not have frozen the direction as well.
    first = L.lights(0.0, 3, 0.35, radius=radius)
    later = L.lights(4.0, 3, 0.35, radius=radius)
    assert not np.allclose(first[0]["pos"], later[0]["pos"]), "the lights stopped orbiting"


def test_the_orbit_never_lands_on_the_origin(qapp):
    """All three sines cross zero together at t=0 for the first light, and a direction of zero
    length has no direction to normalise."""
    for light in L.lights(0.0, 6, 0.35, radius=50.0):
        assert np.isfinite(light["pos"]).all()
        assert abs(float(np.linalg.norm(light["pos"])) - 50.0) < 1e-9


def test_point_modes_move_the_median_point_by_something_a_person_can_see(qapp):
    """The reduced menu is justified by visible, whole-cloud differences, not tiny glints."""
    rng = np.random.default_rng(11)
    xyz = rng.normal(size=(2000, 3)) * 20
    flat = np.full((2000, 4), 0.5)
    lit = L.lights(0.0, 2)
    rendered = {name: L.shade(xyz, flat, lit, point_mode=name)[:, :3]
                for name in L.POINT_MODES}
    for a, b in (("flat", "glossy 3D"), ("flat", "metallic 3D"),
                 ("glossy 3D", "metallic 3D")):
        steps = np.abs(rendered[a] - rendered[b]).max(axis=1) * 255
        assert np.median(steps) > 8, f"{a}/{b} differ only {np.median(steps):.1f}/255"


def test_a_shiny_finish_is_darker_in_the_body_and_brighter_at_the_peaks(qapp):
    """The SHAPE of the difference, not just its size. Chalk is evenly lit all over; a bead is dark
    across most of itself with a few bright places. Both directions have to hold, or "glossy" is
    just "brighter", which is what a gamma slider is for."""
    rng = np.random.default_rng(12)
    xyz = rng.normal(size=(1500, 3)) * 20
    flat = np.full((1500, 4), 0.6)
    lit = L.lights(0.0, 2)
    plain = L.shade(xyz, flat, lit, point_mode="flat")[:, :3]
    glossy = L.shade(xyz, flat, lit, point_mode="glossy 3D")[:, :3]
    assert np.median(glossy) < np.median(plain), "the glossy body is no darker than flat"
    assert glossy.max() > plain.max(), "the glossy peaks are no brighter than flat"


def test_no_finish_can_take_a_gene_below_the_floor(qapp):
    """A metal reading properly as metal is nearly black away from its highlights, and a gene that
    is nearly black is a gene nobody can find. A data display before it is a rendering."""
    xyz = np.random.default_rng(13).normal(size=(300, 3)) * 20
    flat = np.full((300, 4), 1.0)
    for name in L.POINT_MODES:
        out = L.shade(xyz, flat, [], point_mode=name)[:, :3]  # nothing but the floor
        assert out.min() >= L.MIN_AMBIENT - 1e-9, f"{name} bottomed out at {out.min():.3f}"


def test_the_rim_lights_the_silhouette_and_not_the_dark_side(qapp):
    """What the rim is for, and what it must not do.

    For: the edge of the cloud, where the normal turns away from the viewer -- that outline is what
    makes a shape read as glossy, and on a scatter this sparse it is the only specular cue big
    enough to see. Must not: glow on the side no light reaches. A free-standing fresnel does exactly
    that, which on a data display is a point claiming a brightness nothing gave it, so the term is
    multiplied by the light that arrives. It keeps a 0.3 floor, deliberately -- a silhouette point
    is edge-on to the viewer, so under a light behind the viewer its diffuse is ~0 and a strict
    product would switch the rim off precisely where it is wanted."""
    eye = np.array([0.0, 0.0, 120.0])
    xyz = np.array([[0.0, 0.0, 12.0],      # facing the eye and the light
                    [12.0, 0.0, 0.0],      # on the silhouette, edge-on to the eye
                    [0.0, 0.0, -12.0]])    # facing away from both
    flat = np.full((3, 4), 0.5)
    # At the distance the map actually puts its lights -- a few times the radius of the cloud. Set
    # them thirty radii out instead and distance falloff swamps every finish, which measures the
    # falloff rather than the rim.
    behind = [{"pos": np.array([0.0, 0.0, 40.0]), "color": np.array([1.0, 1.0, 1.0])}]
    # Against the SAME finish with the rim switched off, rather than against matt. Matt is not a
    # rimless glossy -- it also sits on a raised floor, so a matt-vs-glossy difference is the floor
    # and the rim together and says nothing about either.
    lit_with = L.shade(xyz, flat, behind, point_mode="glossy 3D", eye=eye)[:, :3].max(axis=1)
    was = L.POINT_MODES["glossy 3D"]["rim"]
    try:
        L.POINT_MODES["glossy 3D"]["rim"] = 0.0
        lit_without = L.shade(xyz, flat, behind, point_mode="glossy 3D", eye=eye)[:, :3].max(axis=1)
    finally:
        L.POINT_MODES["glossy 3D"]["rim"] = was
    gained = lit_with - lit_without
    assert gained[1] > 0.02, f"the silhouette gained nothing ({gained[1]:.4f})"
    assert gained[2] < 0.005, f"the dark side lit itself ({gained[2]:.4f})"


def test_soft_light_lets_every_light_through_everything(win):
    """Soft mode deliberately avoids the ray-tracing cost and the cloud cannot shadow itself."""
    win._lighting["mode"] = "soft"
    lit = [{"pos": np.zeros(3), "color": np.ones(3)}]
    assert win.cast_rays(lit) is lit


def test_a_gene_behind_a_cluster_gets_less_light_than_one_in_front_of_it(win):
    """What "if nothing is between the mouse and the datapoint" means, measured on the real map:
    take the light's own line to every gene and compare the genes it can see with the genes it
    cannot."""
    win._lighting["mode"] = "ray traced"
    try:
        lamp = win.xyz.mean(axis=0) + np.array([0.0, 0.0, float(np.abs(win.xyz).max()) * 2.0])
        lit = win.cast_rays([{"pos": lamp, "color": np.ones(3)}])
        seen = lit[0]["shadow"][:, 0]
        assert seen.max() > 0.9, "nothing at all had a clear line to the light"
        assert seen.min() < 0.6, "nothing at all was blocked -- every gene saw the light"
        # And what is blocked is what has the map in front of it: genes far from the lamp along its
        # own axis are the ones behind everything else.
        along = win.xyz[:, 2]
        assert seen[along < np.percentile(along, 10)].mean() < \
            seen[along > np.percentile(along, 90)].mean()
    finally:
        win._lighting["mode"] = "soft"


def test_the_light_mode_never_decides_whether_an_emitter_is_drawn(win):
    """Ray tracing changes transport; it does not add confusing glowing ray ornaments.

    The assertion this replaces was "no emitter in either mode", which was only true because the
    marker defaulted to `none`. The property actually worth holding is that the MARKER decides and the
    mode never does -- otherwise turning on shadows would quietly add an ornament, which is the thing
    the original test was written to prevent.
    """
    try:
        win.set_lighting_option("target_marker", "none")
        for mode in ("soft", "ray traced", "deep ray traced"):
            win.set_lighting(mode)
            win._light_tick()
            assert win.emitter_item is None or not win.emitter_item.visible(), mode
        win.set_lighting_option("target_marker", "halo")
        drawn = []
        for mode in ("soft", "ray traced", "deep ray traced"):
            win.set_lighting(mode)
            win._light_tick()
            drawn.append(win.emitter_item is not None and win.emitter_item.visible())
        assert len(set(drawn)) == 1, f"the mode changed whether the marker was drawn: {drawn}"
    finally:
        win.set_lighting("off")


def test_the_occupancy_grid_is_built_once_for_a_cloud_that_has_not_moved(win):
    """Rebuilding this every frame costs more than the shading it informs."""
    first = win.occupancy()
    assert win.occupancy() is first
    was = win.xyz
    try:
        win.xyz = np.array(win.xyz, copy=True)
        assert win.occupancy() is not first, "a new cloud reused the old map's shape"
    finally:
        win.xyz = was


def test_an_unknown_light_transport_is_refused_rather_than_mislabelled(win):
    assert win.set_lighting("path traced") == "off"
    assert win.set_lighting("ray traced") == "ray traced"
    win.set_lighting("off")


def test_mouse_light_is_a_camera_flashlight_not_the_gene_under_the_pointer(win, monkeypatch):
    """Overlapping depth planes must not change the illumination source discontinuously."""
    from PyQt6 import QtCore, QtGui
    win.set_lighting("soft")
    win._lighting["source"] = "mouse flashlight"
    try:
        monkeypatch.setattr(win.view, "under_pointer",
                            lambda: pytest.fail("flashlight queried the nearest gene"))
        at = QtCore.QPointF(win.view.width() * 0.7, win.view.height() * 0.35)
        win.view.mouseMoveEvent(QtGui.QMouseEvent(
            QtCore.QEvent.Type.MouseMove, at, at, QtCore.Qt.MouseButton.NoButton,
            QtCore.Qt.MouseButton.NoButton, QtCore.Qt.KeyboardModifier.NoModifier))
        lit = win.frame_lights()
        assert len(lit) == 1 and lit[0].get("spot") and not lit[0].get("local")
        assert np.allclose(lit[0]["pos"], win._eye()), "flashlight did not originate at camera"
        assert np.isfinite(lit[0]["target"]).all()
    finally:
        win.set_lighting("off")


def test_a_torch_lights_the_same_pool_whatever_the_zoom(win):
    """The pool is a fraction of the map's own radius, not of the camera's distance, so leaning in
    does not turn a light into a spotlight or a floodlight."""
    basis = win.view.camera_basis()
    anchor = win.xyz[10]
    close = L.torch(anchor, basis, win._light_radius())
    far_basis = (np.asarray(basis[0]) * 4.0,) + tuple(basis[1:])
    stepped_back = L.torch(anchor, far_basis, win._light_radius())
    assert np.linalg.norm(close[0]["pos"] - anchor) == pytest.approx(
        np.linalg.norm(stepped_back[0]["pos"] - anchor), rel=1e-6)


def test_a_light_inside_the_cloud_lifts_what_is_near_it_whichever_way_it_faces(qapp):
    """The reason a torch works at all here. These normals are radial -- a gene on the far side has
    one pointing away from the viewer -- so a light between the viewer and that gene lands on its
    back and lights nothing. Measured before this term existed, pointing into a far cluster left it
    DIMMER than the near face of the map."""
    xyz = np.random.default_rng(21).normal(size=(800, 3)) * 20
    flat = np.full((800, 4), 0.5)
    target = xyz[0]
    facing_away = {"pos": target - np.array([0.0, 0.0, 6.0]), "color": np.ones(3), "local": True}
    distant = dict(facing_away, local=False)
    near = np.linalg.norm(xyz - target, axis=1) < 12.0
    with_glow = L.shade(xyz, flat, [facing_away])[:, :3].max(axis=1)
    without = L.shade(xyz, flat, [distant])[:, :3].max(axis=1)
    assert with_glow[near].mean() > without[near].mean() * 1.3
    # And it stays local: the far side of the map is not lifted with it.
    far = np.linalg.norm(xyz - target, axis=1) > 40.0
    assert with_glow[far].mean() < with_glow[near].mean() * 0.8


def test_nothing_under_the_pointer_is_not_an_anchor(win):
    """Empty space beside the map. The light falls back to hanging off the pointer's direction
    rather than anchoring on whichever gene happened to be least far away."""
    win.view.pointer_px = (-5000.0, -5000.0)
    assert win.view.under_pointer() is None
    win.view.pointer_px = None
    assert win.view.under_pointer() is None


def test_the_sprite_is_rebuilt_when_the_point_mode_changes_and_not_otherwise(win):
    """It runs on the light timer, so "has anything actually moved" is the whole of its cost."""
    win.set_lighting("off")
    win._lighting["point_mode"] = "glossy 3D"
    assert win._refresh_sprite() is True
    assert win._refresh_sprite() is False, "rebuilt a sprite nothing had changed"
    assert win._sprite_state[0] == "glossy 3D"
    win._lighting["point_mode"] = "flat"
    assert win._refresh_sprite() is True


def test_the_ball_is_lit_from_above_left_when_there_is_no_camera_yet(win, monkeypatch):
    """During startup, and for the whole of an offscreen render, there is no camera to convert a
    world light into a screen direction. Above and to the left is where every reader assumes light
    comes from, and it reads as convex -- a ball lit from below reads as a hollow."""
    from starplast import sprite as SP

    def no_camera():
        raise RuntimeError("no GL context")
    monkeypatch.setattr(win.view, "camera_basis", no_camera)
    win._sprite_state = None
    assert win._refresh_sprite() is True
    assert win._sprite_state[1] == SP.DEFAULT_LIGHT


def test_the_ball_is_lit_from_where_the_light_is(win):
    """The CPU fallback highlight agrees with the scene light on old OpenGL contexts."""
    from starplast import sprite as SP
    old_sel, old_source, failed = win.sel, win._lighting["source"], win.scatter._gpu_failed
    win.set_lighting("soft")
    win._lighting["source"] = "selected gene"
    win.sel = 3
    try:
        win.scatter._gpu_failed = True
        win._sprite_state = None
        win._refresh_sprite()
        finish, where = win._sprite_state
        basis = win.view.camera_basis()
        lit = win.frame_lights()
        expected = SP.to_screen(np.asarray(lit[0]["pos"]) - win.xyz.mean(axis=0), basis)
        assert where == pytest.approx(expected)
        assert np.linalg.norm(where) > 0, f"the selected-gene light had no direction: {where}"
    finally:
        win.scatter._gpu_failed = failed
        win.sel, win._lighting["source"] = old_sel, old_source
        win.set_lighting("off")


def test_mouse_ray_direction_eases_continuously_without_moving_its_origin(win):
    """Smoothing applies to cursor coordinates; the light origin remains the camera."""
    old_source, old_mode = win._lighting["source"], win._lighting["mode"]
    try:
        win._lighting["source"] = "mouse flashlight"
        win._lighting["mode"] = "soft"
        win._lighting["response"] = "smooth"
        win._smoothed_pointer = None
        win.view.pointer = (-1.0, 0.0, 1.0)
        first = win.frame_lights()[0]
        win.view.pointer = (1.0, 0.0, 1.0)
        second = win.frame_lights()[0]
        assert np.allclose(first["pos"], second["pos"])
        right = win.view.camera_basis()[1]
        first_x, second_x = float(first["direction"] @ right), float(second["direction"] @ right)
        assert first_x < second_x < 0.0, "beam jumped to cursor"
    finally:
        win._lighting["source"], win._lighting["mode"] = old_source, old_mode
        win._smoothed_pointer = None


def test_the_pointer_is_followed_without_a_button_held(win):
    """Qt delivers a move event only while a button is down unless the widget asks for tracking. The
    default light source follows the pointer, and the drag that used to be the only way to move it
    is also what orbits the camera -- so the light looked stuck to the cloud, and the rest of the
    time it did not move at all. This is the whole of "I can't see the mouse light"."""
    from PyQt6 import QtCore, QtGui
    assert win.view.hasMouseTracking(), "hover events are not being delivered at all"
    seen = []
    for x in (0.2, 0.8):
        at = QtCore.QPointF(win.view.width() * x, win.view.height() * 0.5)
        win.view.mouseMoveEvent(QtGui.QMouseEvent(
            QtCore.QEvent.Type.MouseMove, at, at, QtCore.Qt.MouseButton.NoButton,
            QtCore.Qt.MouseButton.NoButton, QtCore.Qt.KeyboardModifier.NoModifier))
        seen.append(win.view.pointer)
    assert seen[0] is not None and seen[0][0] < seen[1][0], f"the pointer did not move: {seen}"


def test_a_drag_off_the_edge_does_not_throw_the_light_off_the_map(win):
    """A drag that leaves the widget keeps delivering moves, with coordinates outside it. Unclamped,
    the pointer read several widths out and took its light with it."""
    from PyQt6 import QtCore, QtGui
    # A left button held is the whole point -- that is what a drag is -- and it means pyqtgraph
    # orbits the camera by the same event. The window is shared with every other test in this file,
    # and the corner lights are placed in the camera's frame, so leaving the view spun round quietly
    # breaks a test three functions further down. Put it back.
    # The options dict rather than setCameraParams: that setter refuses to be given rotation and
    # elevation together, and cameraParams() hands back both.
    was = dict(win.view.opts)
    try:
        at = QtCore.QPointF(win.view.width() * 12.0, -400.0)
        win.view.mouseMoveEvent(QtGui.QMouseEvent(
            QtCore.QEvent.Type.MouseMove, at, at, QtCore.Qt.MouseButton.LeftButton,
            QtCore.Qt.MouseButton.LeftButton, QtCore.Qt.KeyboardModifier.NoModifier))
        x, y, _ = win.view.pointer
        assert -1.0 <= x <= 1.0 and -1.0 <= y <= 1.0, f"pointer left the widget: {(x, y)}"
    finally:
        win.view.opts.update(was)
        win.view.update()


def test_a_corner_light_keeps_its_corner_when_the_map_is_turned(qapp):
    """"Top left" is a statement about the picture, not about the data. Placed in world coordinates
    it drifted to the far side of the cloud as soon as anything was orbited, which is why the corner
    presets read as "something changes but the light doesn't move"."""
    centre = np.zeros(3)
    facing = (np.array([0.0, 0.0, 90.0]), np.array([1.0, 0.0, 0.0]),
              np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, -1.0]))
    turned = (np.array([90.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]),
              np.array([0.0, 0.0, 1.0]), np.array([-1.0, 0.0, 0.0]))
    a = L.on_screen(L.CORNERS["top left"], facing, centre, 40.0)[0]["pos"]
    b = L.on_screen(L.CORNERS["top left"], turned, centre, 40.0)[0]["pos"]
    assert not np.allclose(a, b), "the light stayed put while the camera moved"
    # Still top left of the PICTURE in both: left of the camera's right axis, above its up axis.
    for pos, (_, right, up, _) in ((a, facing), (b, turned)):
        assert float(pos @ right) < 0 and float(pos @ up) > 0


def test_the_pointer_light_follows_the_pointer(qapp):
    xyz = np.random.default_rng(5).normal(size=(40, 3)) * 10
    basis = (np.array([0.0, 0.0, 90.0]), np.array([1.0, 0.0, 0.0]),
             np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, -1.0]))
    left = L.light_at(xyz, "mouse", 0.0, 1, 0.3, 30.0, pointer=(-0.9, 0.5, 0.0), basis=basis)
    right = L.light_at(xyz, "mouse", 0.0, 1, 0.3, 30.0, pointer=(0.9, 0.5, 0.0), basis=basis)
    assert left[0]["direction"][0] < right[0]["direction"][0], "the light ignored the pointer"
    assert np.allclose(left[0]["pos"], right[0]["pos"]), "flashlight origin left the camera"


def test_the_highlight_lands_where_the_viewer_actually_is(qapp):
    """Computed against a fixed +Z, the specular term stays on the side of the cloud that faced the
    reader before they orbited -- so from anywhere else, every finish looks matt."""
    xyz = np.random.default_rng(6).normal(size=(500, 3)) * 15
    flat = np.full((500, 4), 0.5)
    lit = [{"pos": np.array([80.0, 0.0, 0.0]), "color": np.array([1.0, 1.0, 1.0])}]
    # Satin, and measured at a percentile rather than at the maximum: glossy's gain saturates, so
    # both views peg a point at pure white and comparing the maxima compares two clipped numbers.
    front = L.shade(xyz, flat, lit, finish="satin", eye=np.array([0.0, 0.0, 200.0]))[:, :3]
    side = L.shade(xyz, flat, lit, finish="satin", eye=np.array([200.0, 0.0, 0.0]))[:, :3]
    assert not np.allclose(front, side), "the highlight ignored the camera"
    assert np.percentile(side, 99) > np.percentile(front, 99), "looking down the light, no glint"


def test_lighting_a_gene_lights_the_gene(qapp):
    xyz = np.random.default_rng(4).normal(size=(30, 3)) * 5
    lit = L.at_points(xyz, [7], 20.0)
    assert np.allclose(lit[0]["pos"], xyz[7])


def test_the_pointer_is_read_in_the_view_s_own_frame(win):
    """The map rotates. A light fixed in world space swings away from the pointer the moment
    anything moves, which reads as the light being broken."""
    from PyQt6 import QtCore, QtGui
    win.view.resize(200, 100)
    ev = QtGui.QMouseEvent(QtCore.QEvent.Type.MouseMove, QtCore.QPointF(150, 25),
                           QtCore.QPointF(150, 25), QtCore.Qt.MouseButton.NoButton,
                           QtCore.Qt.MouseButton.NoButton, QtCore.Qt.KeyboardModifier.NoModifier)
    win.view.mouseMoveEvent(ev)
    x, y, z = win.view.pointer
    assert x > 0 and y > 0 and z > 0, "top right of the widget should be up and to the right"


def test_the_edges_that_light_a_neighbour_are_the_ones_being_drawn(win):
    """A light on a relationship the reader cannot see answers a question they did not ask."""
    with_edges = win.edge_neighbours(10)
    was = dict(win.edge_on)
    try:
        for k in win.edge_on:
            win.edge_on[k] = False
        assert win.edge_neighbours(10) == [], "a switched-off edge type still lit a neighbour"
    finally:
        win.edge_on.update(was)
    assert win.edge_neighbours(10) == with_edges


def test_lighting_still_draws_when_there_is_no_camera_to_ask(win, monkeypatch):
    """Every screen-relative source needs the camera's own axes, and there is a window of startup --
    and the whole of an offscreen render -- where the renderer has no camera to give. Falling over
    there would take the map with it, so the fallback is a light from the front."""
    def no_camera():
        raise RuntimeError("no GL context")
    monkeypatch.setattr(win.view, "camera_basis", no_camera)
    win._lighting["source"] = "mouse flashlight"
    assert win._eye() is None
    lit = win.frame_lights()
    assert len(lit) == 1 and lit[0]["pos"][2] > 0, "no camera left the map unlit"
    win._light_tick()                                  # and the frame still draws


def test_the_grid_is_lit_by_the_same_lights_as_the_points(win):
    """A lit cloud over an unlit grid reads as two pictures: the horizon is what the eye uses to
    judge where the light is coming from."""
    win.show_ground = True
    win.redraw()
    assert win.grid_item is not None
    before = win.grid_item.color().getRgb()
    win.set_lighting("soft")
    win._lighting["source"] = "mouse flashlight"
    win._lighting["mood"] = "neutral"
    win._light_tick()
    after = win.grid_item.color().getRgb()
    assert after[:3] != before[:3], "the grid ignored the light"
    # Alpha moves too, and it did not always. Holding it at the theme's value was defensible in
    # principle -- how solid the grid is IS a theme decision -- but the grid is drawn so faint that
    # a colour-only change was below the threshold of visible, which is how it was reported: "I
    # can't see the grid being lit". Light moves both now, and the band keeps it inside the range
    # the theme set rather than letting a bright light turn the horizon into a wall.
    assert after[3] != before[3], "the grid's alpha ignored the light"
    assert 15 <= after[3] <= 210, "the lit grid left the theme's range"

    # The floor shares the mood as well as the brightness; otherwise warm genes would float above
    # a neutral horizon and read as a separate composited layer.
    win._lighting["mood"] = "warm"
    win._light_tick()
    warm = win.grid_item.color().getRgb()
    assert warm[:3] != after[:3], "the grid ignored the light mood"
    win.set_lighting("off")


# --------------------------------------------------------------------------- containers
def test_the_background_sits_behind_every_container(win):
    """Behind one panel it was scenery for one corner of the screen."""
    win.set_ambient("blobs")
    assert win._ambient_widget.parent() is win
    assert win._ambient_widget.size() == win.size()
    win.set_ambient("none")


def test_panel_opacity_lets_the_background_through_but_not_the_fields(win):
    from starplast import theme as TH
    assert win.set_container_opacity(0.6) == 0.6
    sheet = TH.stylesheet(win.theme, 0.6)
    assert "rgba(" in sheet, "the containers are still opaque"
    field_line = [line for line in sheet.splitlines() if "QLineEdit" in line and "background" in line]
    assert not any("rgba(" in line for line in field_line), \
        "a translucent field puts moving colour behind text somebody is reading"
    assert win.set_container_opacity(0.0) == 0.35, "clamped: below this the list is unreadable"
    win.set_container_opacity(1.0)


def test_a_source_or_point_mode_that_does_not_exist_is_refused(win):
    """A settings file written by a later version, or a typo in one written by hand."""
    win.set_lighting_option("source", "selected gene")
    assert win.set_lighting_option("source", "from behind the sofa") == "selected gene"
    win.set_lighting_option("point_mode", "glossy 3D")
    assert win.set_lighting_option("point_mode", "velvet") == "glossy 3D"


def test_changing_a_light_setting_shows_at_once_while_lit(win):
    """Turning a knob and seeing nothing until the next frame reads as the knob doing nothing."""
    win.set_lighting("soft")
    win.redraw()
    before = np.array(win.scatter.color if hasattr(win.scatter, "color") else [], copy=True)
    win.set_lighting_option("mood", "cool blue")
    assert win._lighting["mood"] == "cool blue"
    after = np.array(win.scatter.color if hasattr(win.scatter, "color") else [], copy=True)
    if win.scatter.gpu_material_enabled(win._lighting["point_mode"]):
        assert np.allclose(win.scatter._scene["mood"], L.mood_color("cool blue"))
        assert np.array_equal(after, before), "the GPU material rewrote its data-color albedo"
    else:
        assert not np.array_equal(after, before), "the visible points did not update immediately"
    win.set_lighting("off")


def test_the_grid_is_left_alone_when_there_is_none(win):
    """`show_ground` off, or a table with no genes: the lighting must not assume a grid exists."""
    was = win.show_ground
    win.show_ground = False
    win.redraw()
    assert win.grid_item is None
    win._light_ground(win.frame_lights())          # must not raise
    win.show_ground = was
    win.redraw()


def test_a_pointer_mode_with_no_cone_falls_back_to_a_plain_pointer_light(qapp):
    """One of the pointer modes is not a flashlight at all -- it is the old point light that simply
    follows the cursor. It shares the branch with the cone modes and must not try to build one."""
    xyz = np.random.default_rng(11).normal(size=(40, 3)) * 10
    basis = (np.array([0.0, 0.0, 90.0]), np.array([1.0, 0.0, 0.0]),
             np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, -1.0]))
    plain = [name for name, cone in L.POINTER_MODES.items() if cone is None]
    assert plain, "no coneless pointer mode exists to exercise this path"
    lit = L.light_at(xyz, "mouse flashlight", 0.0, 1, 0.3, 30.0,
                     pointer=(-0.8, 0.4, 0.0), basis=basis, pointer_mode=plain[0])
    assert len(lit) == 1 and not lit[0].get("spot"), "a coneless mode produced a spotlight"
    # And with no pointer at all it lights from the front rather than dividing by nothing.
    ahead = L.light_at(xyz, "mouse flashlight", 0.0, 1, 0.3, 30.0,
                       pointer=None, basis=basis, pointer_mode=plain[0])
    assert len(ahead) == 1 and np.isfinite(ahead[0]["pos"]).all()


# --------------------------------------------------------------------------- lighting settings
def test_every_lighting_setting_refuses_a_value_it_does_not_know(win):
    """Each of these reaches the map from a settings file that a newer version may have written, so
    an unknown value is the ordinary case rather than a programming error. Refusing it keeps the
    previous value; falling through would set the map to a mode that does not exist."""
    for key in ("pointer_mode", "response", "target_marker", "mood"):
        was = win._lighting[key]
        assert win.set_lighting_option(key, "a mode from the future") == str(was)
        assert win._lighting[key] == was, f"{key} accepted a value it does not know"
    assert win.set_lighting_option("not_a_setting", "anything") == ""


def test_the_target_marker_draws_and_stops_drawing(win):
    """The marker says where the light is aimed. It is scenery -- it encodes nothing about the data
    -- so the only claims to hold are that each shape produces points, that beacon adds a lifted
    twin, and that switching it off hides the item rather than leaving it on screen."""
    win.set_lighting("lit")
    try:
        lit = win.frame_lights()
        counts = {}
        for marker in ("halo", "beacon", "pulse"):
            win._lighting["target_marker"] = marker
            win._draw_emitter(lit)
            assert win.emitter_item is not None and win.emitter_item.visible(), marker
            counts[marker] = len(win.emitter_item.pos)
        assert counts["beacon"] == 2 * counts["halo"], "the beacon lost its lifted twin"
        win._lighting["target_marker"] = "none"
        win._draw_emitter(lit)
        assert not win.emitter_item.visible(), "switching the marker off left it on screen"
        win._draw_emitter([])
        assert not win.emitter_item.visible()
    finally:
        win._lighting["target_marker"] = "none"
        win.set_lighting("off")


def test_the_cpu_shading_path_runs_when_the_gpu_material_is_not_in_use(win, monkeypatch):
    """Two ways a gene gets its colour: the GPU material shades per fragment and wants the flat
    colours untouched, or the CPU shades per point. Both must draw, because the second is what runs
    on any driver that refuses the shader -- and it is the one every headless test sees."""
    win.set_lighting("lit")
    try:
        monkeypatch.setattr(win.scatter, "gpu_material_enabled", lambda mode: False)
        win.redraw()
        cpu = np.array(win.scatter.color, dtype=float)
        monkeypatch.setattr(win.scatter, "gpu_material_enabled", lambda mode: True)
        win.redraw()
        gpu_side = np.array(win.scatter.color, dtype=float)
        assert cpu.shape == gpu_side.shape
        # The GPU path hands the flat colours through untouched; the CPU path has shaded them.
        assert not np.allclose(cpu[:, :3], gpu_side[:, :3]), "the two paths produced one picture"
    finally:
        win.set_lighting("off")


# --------------------------------------------------------------------------- lights that travel
def test_a_wandering_light_stays_inside_the_cloud_and_keeps_moving(qapp):
    """The point of this source is that it travels THROUGH the map rather than around it. A light
    that drifted outside the hull would be an ordinary exterior light that happened to move."""
    xyz = np.random.default_rng(3).normal(size=(600, 3)) * np.array([20.0, 8.0, 35.0])
    centre = xyz.mean(axis=0)
    extent = np.abs(xyz - centre).max(axis=0)
    seen = []
    for t in np.arange(0.0, 400.0, 7.0):
        lit = L.wandering(xyz, float(t), radius=60.0)
        assert len(lit) == 1 and lit[0]["local"], "an interior light that is not local lights nothing"
        where = lit[0]["pos"]
        assert (np.abs(where - centre) <= extent).all(), f"left the cloud at t={t}: {where}"
        seen.append(where)
    seen = np.array(seen)
    # It has to actually go somewhere: a light that jitters around one spot is not wandering.
    assert np.linalg.norm(seen.max(axis=0) - seen.min(axis=0)) > np.linalg.norm(extent) * 0.5


def test_the_wander_is_reproducible_and_never_repeats(qapp):
    """Deterministic because everything here is a function of the clock -- a dropped frame costs
    nothing and two machines agree. Non-repeating because the frequencies share no common multiple,
    so it does not settle into a loop that leaves half the map permanently dark."""
    xyz = np.random.default_rng(4).normal(size=(200, 3)) * 15
    assert np.allclose(L.wandering(xyz, 12.5, 40.0)[0]["pos"], L.wandering(xyz, 12.5, 40.0)[0]["pos"])
    early = np.array([L.wandering(xyz, float(t), 40.0)[0]["pos"] for t in range(0, 60)])
    late = np.array([L.wandering(xyz, float(t), 40.0)[0]["pos"] for t in range(600, 660)])
    assert not np.allclose(early, late, atol=1.0), "the path closed into a loop"


def test_a_bouncing_light_turns_around_at_the_walls(qapp):
    """A ball in a box, in closed form. The property that matters is the turn: the path must reverse
    at the wall rather than pass through it or wrap around to the other side."""
    xyz = np.random.default_rng(5).normal(size=(400, 3)) * 12
    centre = xyz.mean(axis=0)
    extent = np.abs(xyz - centre).max(axis=0)
    path = np.array([L.bouncing(xyz, float(t) * 0.25, 40.0)[0]["pos"] for t in range(1200)])
    assert (np.abs(path - centre) <= extent + 1e-9).all(), "the ball left the box"
    x = path[:, 0]
    turns = np.sum(np.diff(np.sign(np.diff(x))) != 0)
    assert turns > 4, f"only {turns} turns in 1,200 steps -- it is not bouncing"
    # And it reaches the walls rather than staying safely in the middle. The walls are at
    # INTERIOR of the extent, deliberately inside the hull -- a light sitting exactly on the hull
    # reads as an ordinary outside light that stopped moving.
    assert np.abs(x - centre[0]).max() > extent[0] * L.INTERIOR * 0.95


def test_the_bounce_is_closed_form_rather_than_simulated(qapp):
    """Asking for t=900 directly must give the same answer as arriving there frame by frame would.
    That is what having no state buys: a search that jumps the clock lands where the light really
    is, and two machines never diverge."""
    xyz = np.random.default_rng(6).normal(size=(150, 3)) * 10
    assert np.allclose(L.bouncing(xyz, 900.0, 30.0)[0]["pos"],
                       L.bouncing(xyz, 900.0, 30.0)[0]["pos"])
    # A triangle wave of period 4 over [-1, 1]: centre, wall, centre, other wall, centre.
    assert L._bounce_axis(0.0) == pytest.approx(0.0)
    assert L._bounce_axis(1.0) == pytest.approx(-1.0)
    assert L._bounce_axis(3.0) == pytest.approx(1.0)
    assert L._bounce_axis(4.0) == pytest.approx(0.0), "the period is not 4"
    assert all(-1.0 <= L._bounce_axis(v) <= 1.0 for v in np.arange(-20, 20, 0.37))


def test_both_travelling_lights_light_what_they_pass(qapp):
    """The claim a reader makes about them: genes near the light are brighter than genes far from
    it. Without the local glow an interior light lands on the back of everything and lights almost
    nothing, which was measured once already with the pointer torch."""
    xyz = np.random.default_rng(7).normal(size=(800, 3)) * 18
    flat = np.full((800, 4), 0.55)
    for source in (L.wandering, L.bouncing):
        lit = source(xyz, 31.0, 45.0)
        shaded = L.shade(xyz, flat, lit)[:, :3].max(axis=1)
        near = np.linalg.norm(xyz - lit[0]["pos"], axis=1)
        close, far = near < np.percentile(near, 10), near > np.percentile(near, 60)
        assert shaded[close].mean() > shaded[far].mean() * 1.15, source.__name__


def test_several_travelling_lights_do_not_move_as_one(qapp):
    xyz = np.random.default_rng(8).normal(size=(300, 3)) * 14
    for source in (L.wandering, L.bouncing):
        lit = source(xyz, 19.0, 40.0, n=3)
        assert len(lit) == 3
        spread = [np.linalg.norm(a["pos"] - b["pos"]) for a in lit for b in lit if a is not b]
        assert min(spread) > 1e-6, f"{source.__name__} lights are stacked on each other"


def test_a_travelling_light_with_no_map_to_travel_through_lights_from_the_front(qapp):
    for source in (L.wandering, L.bouncing):
        lit = source(np.zeros((0, 3)), 5.0, 30.0)
        assert len(lit) == 1 and lit[0]["pos"][2] > 0


def test_the_travelling_lights_are_selectable_and_need_no_camera(qapp):
    """Both are facts about the data's own extent rather than about where the reader stands, so they
    light the same genes whether the map has been orbited or not -- and they work offscreen, where
    every other source falls back to a light from the front."""
    xyz = np.random.default_rng(9).normal(size=(300, 3)) * 16
    assert "wandering light" in L.SOURCES and "bouncing light" in L.SOURCES
    for source in ("wandering light", "bouncing light"):
        with_camera = L.light_at(xyz, source, 21.0, 2, 0.35, 40.0,
                                 basis=(np.array([0.0, 0.0, 90.0]), np.array([1.0, 0.0, 0.0]),
                                        np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, -1.0])))
        without = L.light_at(xyz, source, 21.0, 2, 0.35, 40.0, basis=None)
        assert len(with_camera) == 2 and len(without) == 2, source
        assert np.allclose(with_camera[0]["pos"], without[0]["pos"]), \
            f"{source} moved when the camera did"
        assert all(light.get("local") for light in with_camera), source


# --------------------------------------------------------------------------- the requested defaults
def test_the_five_defaults_are_what_was_asked_for(qapp, tmp_path, monkeypatch):
    """A fresh install, with nothing stored. These were chosen by the person who uses the program;
    the test exists so a later refactor cannot quietly move them back."""
    from PyQt6 import QtCore
    from starplast import paths
    from starplast.app import Window, DEFAULT_POINT_SIZE
    monkeypatch.setenv(paths.ENV_STATE, str(tmp_path))
    QtCore.QSettings("starplast", "starplast").clear()
    w = Window()
    try:
        assert w._lighting["mode"] == "deep ray traced"
        assert w._lighting["source"] == "selected gene and its edges"
        assert w._lighting["target_marker"] == "halo"
        assert w._lighting["point_mode"] == "glossy 3D"
        assert w.point_size == DEFAULT_POINT_SIZE == 7.0
        assert dict(POINT_SIZES_BY_LABEL())["Medium (7 px)"] == DEFAULT_POINT_SIZE
    finally:
        w.close()


def POINT_SIZES_BY_LABEL():
    from starplast.app import POINT_SIZES
    return POINT_SIZES


# --------------------------------------------------------------------------- finishes worth telling apart
def test_every_finish_is_visibly_different_from_every_other():
    """The complaint that removed three finishes: `glossy 3D` and `pearl 3D` looked the same, and so
    did two other pairs. Measured rather than judged -- shade the whole cloud under one light and take
    the mean absolute RGB difference between each pair. `glossy 3D` against `pearl 3D` was 0.0204,
    which is the number this floor is set against."""
    import itertools
    import numpy as np
    rng = np.random.default_rng(0)
    # A calibration sphere samples all pseudo-normal directions. A particular
    # UMAP can collapse most directions and is not a stable material fixture.
    xyz = rng.normal(size=(4096, 3))
    xyz = xyz / np.linalg.norm(xyz, axis=1, keepdims=True) * 50
    radius = 50.
    colors = rng.random((len(xyz), 3)) * 0.7 + 0.2
    lit = L.fixed((0.4, 0.6, 0.7), radius)
    shaded = {n: L.shade(xyz, colors, lit, True, point_mode=n)[:, :3] for n in L.POINT_MODES}
    worst = min((float(np.abs(shaded[a] - shaded[b]).mean()), a, b)
                for a, b in itertools.combinations(shaded, 2))
    assert worst[0] >= 0.05, f"{worst[1]} and {worst[2]} differ by only {worst[0]:.4f}"


def test_the_retired_finishes_still_resolve_to_something_kept():
    """A stored setting or a saved recipe naming a removed finish must keep drawing, not fall back to
    the default and silently change a figure."""
    assert L.normalize_point_mode("pearl 3D") == "glossy 3D"
    assert L.normalize_point_mode("brushed metal 3D") == "metallic 3D"
    assert L.normalize_point_mode("silver 3D") == "metallic 3D"
    assert L.normalize_point_mode("nonsense") == L.DEFAULT_POINT_MODE
    for name in ("pearl 3D", "brushed metal 3D", "silver 3D"):
        assert name not in L.POINT_MODES


# --------------------------------------------------------------------------- the travelling orb
def test_a_travelling_light_actually_travels(win):
    """The report was that the wandering and bouncing lights "do not work". Two separate causes, and
    this is the first: at the old speed the slowest term of `_drift` had a period of 184 seconds, so
    the orb moved about 4% of the map's radius per second and the brightest gene did not change for
    seconds at a time. One unit of `t` is one second."""
    import numpy as np
    xyz = win.xyz
    radius = float(np.abs(xyz - xyz.mean(axis=0)).max())
    for source in ("wandering light", "bouncing light"):
        steps, brightest = [], []
        previous = None
        colors = np.tile(np.array([[0.4, 0.6, 0.9]]), (len(xyz), 1))
        for t in range(7):
            lit = L.light_at(xyz, source, float(t), 1, L.DEFAULT_SPEED, radius)
            pos = np.asarray(lit[0]["pos"], dtype=float)
            if previous is not None:
                steps.append(float(np.linalg.norm(pos - previous)))
            previous = pos
            brightest.append(int(L.shade(xyz, colors, lit, True,
                                         point_mode="glossy 3D")[:, :3].sum(axis=1).argmax()))
        assert np.mean(steps) > 0.10 * radius, f"{source} crawls: {np.mean(steps):.2f} of {radius}"
        assert len(set(brightest)) >= 4, f"{source} lights the same genes all along: {brightest}"


def test_a_travelling_light_is_brighter_and_tighter_than_a_held_one():
    """The second cause: at gain 1.0 the pool measured 1.2 to 2.0 times the median gene, which reads as
    no light at all against 8,000 other bright points. Raising the gain ALONE made it worse -- a pool
    that wide lifts the median too -- so the travelling lights get their own narrower width."""
    import numpy as np
    # Symmetric probes around the source isolate the local pool from unrelated
    # directional highlights and RGB clipping. A map-wide maximum can be a
    # distant specular highlight, so its ratio to the median did not measure
    # the travelling light's pool at all after the layout changed.
    xyz = np.array([[x,0.,0.] for x in (-50.,-25.,-10.,10.,25.,50.)])
    radius = 50.
    colors = np.full((len(xyz),3),0.1)
    lit = L.wandering(xyz, 3.0, radius)
    lit[0]['pos'] = np.zeros(3)
    assert lit[0]["gain"] == L.TRAVEL_GAIN > 1.0
    assert lit[0]["width"] == L.TRAVEL_WIDTH < L.LOCAL_WIDTH
    tight = L.shade(xyz, colors, lit, False)[:, :3].sum(axis=1)
    plain = L.shade(xyz, colors, [{k: v for k, v in lit[0].items() if k not in ("gain", "width")}],
                    False)[:, :3].sum(axis=1)
    assert tight[2] > plain[2] * 1.3
    assert tight[2] / tight[0] > plain[2] / plain[0]


def test_the_finishes_are_distinct_on_the_renderer_the_user_actually_sees(win):
    """The CPU-fallback measurement was the wrong path, and this test exists because of that mistake.

    The eight finishes were culled to five on `lighting.shade` numbers -- the CPU fallback -- which put
    the closest surviving pair at 0.070 and looked safely clear. Rendered through the PRODUCTION GPU
    shader and measured over the pixels that actually carry genes, `glass 3D` and `glossy 3D` differed
    on only **49%** of them while every other pair differed on 89-100%. Roughness alone does not
    separate two dielectrics at seven pixels across: the highlight is a couple of pixels either way,
    and what separates them has to be the body.

    Asserted here on the shader's own parameters rather than by rendering, because rendering needs a
    display and this suite must not. The rendered check lives in `scripts/benchmark_lighting_lab.py`,
    whose recorded run is in `results/`.
    """
    from starplast import sprite as SP
    glass, glossy = SP.SPRITES["glass 3D"], SP.SPRITES["glossy 3D"]
    # The body, which is what carries at this size.
    assert glass["ambient"] < glossy["ambient"], "glass is no darker in the body than glossy"
    assert glass["diffuse"] < glossy["diffuse"] * 0.6, "glass scatters too much to read as glass"
    assert glass["rim"] > glossy["rim"] * 3, "glass has no silhouette to be bright at"
    source = SP._SPHERE_FRAGMENT
    assert "albedo * 0.08" in source, "the glass body is no longer darkened in the shader"
    assert "envFresnel * 0.75" in source, "the glass rim is no longer strengthened in the shader"
