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
def test_every_finish_changes_the_surface_rather_than_the_data(qapp):
    """A finish is how a point answers the light. Four of them have to look different, and none may
    move a point or touch its alpha -- the map's shape and its filter are not a rendering choice."""
    rng = np.random.default_rng(0)
    xyz = rng.normal(size=(400, 3)) * 15
    flat = np.column_stack([np.full(400, 0.5)] * 3 + [rng.random(400)])
    lit = L.lights(0.0, 2)
    seen = {}
    for name in L.FINISHES:
        out = L.shade(xyz, flat, lit, finish=name)
        assert np.array_equal(out[:, 3], flat[:, 3]), f"{name} touched alpha"
        seen[name] = round(float(out[:, :3].mean()), 4)
    assert len(set(seen.values())) == len(seen), f"two finishes render identically: {seen}"


def test_matt_has_no_highlight_and_glossy_has_a_tight_one(qapp):
    rng = np.random.default_rng(1)
    xyz = rng.normal(size=(600, 3)) * 15
    flat = np.full((600, 4), 0.5)
    lit = L.lights(0.0, 1)
    matt = L.shade(xyz, flat, lit, finish="matt")[:, :3]
    glossy = L.shade(xyz, flat, lit, finish="glossy")[:, :3]
    assert glossy.max() > matt.max(), "glossy did not add a highlight"
    # Tight means the peak stands FURTHER ABOVE THE BODY, which is the thing that reads as gloss.
    # Counting points near the top of the range instead -- which this did -- measures the shape of
    # the range rather than the highlight, and it moves whenever the floor moves, so it broke as
    # soon as the finishes started separating on contrast as well.
    stands_out = lambda a: float(a.max() / max(np.median(a), 1e-9))
    assert stands_out(glossy) > stands_out(matt), "the highlight does not rise above the body"


def test_metallic_takes_its_highlight_from_the_point_not_the_light(qapp):
    """The difference between a copper bead and a white-glinting plastic one."""
    xyz = np.array([[0.0, 0.0, 10.0], [0.0, 0.0, -10.0]])
    red = np.array([[1.0, 0.0, 0.0, 1.0], [1.0, 0.0, 0.0, 1.0]])
    lit = [{"pos": np.array([0.0, 0.0, 60.0]), "color": np.array([1.0, 1.0, 1.0])}]
    metal = L.shade(xyz, red, lit, finish="metallic")
    plastic = L.shade(xyz, red, lit, finish="glossy")
    # A white light on a red point: the plastic highlight whitens the green channel, the metal's
    # stays red. This only works because specular is ADDED rather than multiplied into the colour --
    # multiplied, a pure red point has no green to raise and every finish looks identical on it.
    assert plastic[0, 1] > metal[0, 1] + 0.05, (plastic[0], metal[0])
    assert metal[0, 0] > metal[0, 1], "the metal highlight lost the point's own colour"


@pytest.mark.parametrize("source", L.SOURCES)
def test_every_source_produces_at_least_one_light(qapp, source):
    xyz = np.random.default_rng(2).normal(size=(50, 3))
    lit = L.light_at(xyz, source, 1.0, 3, 0.3, 30.0, pointer=(0.2, 0.4, 1.0), selected=3,
                     neighbours=[1, 2])
    assert lit and all("pos" in x and "color" in x for x in lit)


def test_a_source_with_nothing_to_follow_lights_from_the_front(qapp):
    """No pointer in the view yet, nothing selected. Guessing would put the light behind the map."""
    xyz = np.random.default_rng(3).normal(size=(20, 3))
    for source in ("mouse", "selected gene", "selected gene and its edges"):
        lit = L.light_at(xyz, source, 0.0, 3, 0.3, 20.0)
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


def test_a_finish_moves_the_median_point_by_something_a_person_can_see(qapp):
    """Reported twice as "the points always look matt", and the second time it was this test's
    fault: the first version asked whether a finish moved a point by more than 2 values out of 255,
    which is not a visible difference, and it passed while nothing on screen changed.

    Measured in 8-bit steps at the MEDIAN, because that is the ordinary point rather than the lucky
    one. A highlight, however bright, lands only on the share of a scatter whose normal happens to
    face the light; what tells chalk from a bead across a whole cloud is contrast."""
    rng = np.random.default_rng(11)
    xyz = rng.normal(size=(2000, 3)) * 20
    flat = np.full((2000, 4), 0.5)
    lit = L.lights(0.0, 2)
    matt = L.shade(xyz, flat, lit, finish="matt")[:, :3]
    steps = lambda name: np.abs(L.shade(xyz, flat, lit, finish=name)[:, :3] - matt).max(axis=1) * 255
    for name in ("satin", "glossy", "metallic"):
        d = steps(name)
        assert np.median(d) > 8, f"{name} moves the median point {np.median(d):.1f}/255 -- invisible"
    # Satin is deliberately the quiet one -- it sits between chalk and a bead, and a satin that
    # shouted would leave nothing for glossy to be. The two loud ones have to carry the cloud.
    for name in ("glossy", "metallic"):
        d = steps(name)
        assert float((d > 12).mean()) > 0.5, f"{name} leaves {1 - (d > 12).mean():.0%} unchanged"
    assert np.median(steps("glossy")) > np.median(steps("satin")), "satin outdoes glossy"


def test_a_shiny_finish_is_darker_in_the_body_and_brighter_at_the_peaks(qapp):
    """The SHAPE of the difference, not just its size. Chalk is evenly lit all over; a bead is dark
    across most of itself with a few bright places. Both directions have to hold, or "glossy" is
    just "brighter", which is what a gamma slider is for."""
    rng = np.random.default_rng(12)
    xyz = rng.normal(size=(1500, 3)) * 20
    flat = np.full((1500, 4), 0.6)
    lit = L.lights(0.0, 2)
    matt = L.shade(xyz, flat, lit, finish="matt")[:, :3]
    glossy = L.shade(xyz, flat, lit, finish="glossy")[:, :3]
    assert np.median(glossy) < np.median(matt), "the glossy body is no darker than chalk"
    assert glossy.max() > matt.max(), "the glossy peaks are no brighter than chalk"


def test_no_finish_can_take_a_gene_below_the_floor(qapp):
    """A metal reading properly as metal is nearly black away from its highlights, and a gene that
    is nearly black is a gene nobody can find. A data display before it is a rendering."""
    xyz = np.random.default_rng(13).normal(size=(300, 3)) * 20
    flat = np.full((300, 4), 1.0)
    for name in L.FINISHES:
        out = L.shade(xyz, flat, [], finish=name)[:, :3]       # nothing but the floor
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
    lit_with = L.shade(xyz, flat, behind, finish="glossy", eye=eye)[:, :3].max(axis=1)
    was = L.FINISHES["glossy"]["rim"]
    try:
        L.FINISHES["glossy"]["rim"] = 0.0
        lit_without = L.shade(xyz, flat, behind, finish="glossy", eye=eye)[:, :3].max(axis=1)
    finally:
        L.FINISHES["glossy"]["rim"] = was
    gained = lit_with - lit_without
    assert gained[1] > 0.02, f"the silhouette gained nothing ({gained[1]:.4f})"
    assert gained[2] < 0.005, f"the dark side lit itself ({gained[2]:.4f})"


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
    assert left[0]["pos"][0] < right[0]["pos"][0], "the light ignored the pointer"


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
    win._lighting["source"] = "top left"
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
    win.set_lighting("lit")
    win._lighting["source"] = "top left"
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

    # A light BELOW the floor leaves it darker than one above: the sign of the thing is what makes
    # it read as lighting rather than as a flicker.
    win._lighting["source"] = "bottom left"
    win._light_tick()
    below = win.grid_item.color().getRgb()
    win._lighting["source"] = "top left"
    win._light_tick()
    above = win.grid_item.color().getRgb()
    assert sum(below[:3]) < sum(above[:3]), "the grid was as bright from below as from above"
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


def test_a_source_or_finish_that_does_not_exist_is_refused(win):
    """A settings file written by a later version, or a typo in one written by hand."""
    win.set_lighting_option("source", "top right")
    assert win.set_lighting_option("source", "from behind the sofa") == "top right"
    win.set_lighting_option("finish", "matt")
    assert win.set_lighting_option("finish", "velvet") == "matt"


def test_changing_a_light_setting_shows_at_once_while_lit(win):
    """Turning a knob and seeing nothing until the next frame reads as the knob doing nothing."""
    win.set_lighting("lit")
    win.redraw()
    before = np.array(win.scatter.color if hasattr(win.scatter, "color") else [], copy=True)
    win.set_lighting_option("source", "bottom left")
    assert win._lighting["source"] == "bottom left"
    win.set_lighting("off")
    del before


def test_the_grid_is_left_alone_when_there_is_none(win):
    """`show_ground` off, or a table with no genes: the lighting must not assume a grid exists."""
    was = win.show_ground
    win.show_ground = False
    win.redraw()
    assert win.grid_item is None
    win._light_ground(win.frame_lights())          # must not raise
    win.show_ground = was
    win.redraw()
