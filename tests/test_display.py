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
