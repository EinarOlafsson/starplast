"""Translucent black, rounded windows: the colours, what is dressed, and the opaque fallback.

Each claim the theme makes is checked on the object it is made about: a menu that is said to be
translucent carries the attribute that makes it so, a dialog said to be rounded has no square frame
and paints nothing in its corners, and a display without a compositor gets opaque near-black and a
rounded mask rather than a translucent window it cannot show.
"""
import os
import re

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")

from PyQt6 import QtCore, QtTest, QtWidgets  # noqa: E402

from starplast import glass  # noqa: E402
from starplast import theme as TH  # noqa: E402

W = QtCore.Qt.WidgetAttribute


@pytest.fixture(scope="module")
def qapp():
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    glass.install(app)
    app.setProperty("starplastTheme", "dark")
    app.setStyleSheet(TH.stylesheet("dark"))
    yield app


# --------------------------------------------------------------------------- colours
@pytest.mark.parametrize("theme,ink", [("dark", 0), ("slate", 0), ("light", 255), ("paper", 255)])
def test_glass_is_black_on_a_dark_theme_and_white_on_a_light_one(theme, ink):
    for role in ("menu", "tooltip", "popup", "dialog"):
        colour = TH.glass(theme, role)
        assert colour.startswith(f"rgba({ink}, {ink}, {ink}, ")
        alpha = float(colour.rsplit(",", 1)[1].rstrip(")"))
        assert 0.8 <= alpha < 1.0, "translucent, but not so thin the text competes with the map"


def test_without_a_compositor_the_same_glass_is_opaque_and_near_the_theme_ground():
    dark = TH.glass("dark", "menu", translucent=False)
    assert re.fullmatch(r"#[0-9a-f]{6}", dark)
    assert max(int(dark[i:i + 2], 16) for i in (1, 3, 5)) < 0x10, "near-black, not grey"
    light = TH.glass("light", "menu", translucent=False)
    assert min(int(light[i:i + 2], 16) for i in (1, 3, 5)) > 0xf0


def _rule(sheet, selector):
    m = re.search(re.escape(selector) + r"\s*\{([^}]*)\}", sheet)
    assert m, f"no rule for {selector}"
    return m.group(1)


def test_menus_tooltips_and_lists_are_rounded_glass_in_the_stylesheet():
    sheet = TH.stylesheet("dark")
    for selector in ("QMenu", "QToolTip", "QComboBox QAbstractItemView",
                     "QFrame#HelpSearchResults"):
        rule = _rule(sheet, selector)
        assert "rgba(0, 0, 0, " in rule, f"{selector} is not black glass"
        assert "border-radius" in rule, f"{selector} has square corners"
    item = _rule(sheet, "QMenu::item")
    assert "border-radius" in item, "menu items are rounded pills, as in spaCR"
    bar = _rule(sheet, "QMenuBar::item:selected, QMenuBar::item:pressed")
    assert TH.palette_for("dark")["accent"] in bar, "the pointed-at word lights, spaCR's convention"
    assert "transparent" not in bar


def test_an_opaque_display_gets_no_translucent_popups():
    sheet = TH.stylesheet("dark", translucent=False)
    for selector in ("QMenu", "QToolTip", "QComboBox QAbstractItemView"):
        assert "rgba(" not in _rule(sheet, selector)


def test_fields_stay_readable_on_every_theme():
    """The field grey is dark in every theme, so its text must be light in every theme."""
    def luminance(h):
        lin = [c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
               for c in TH.rgbf(h)[:3]]
        return 0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2]
    for theme in TH.THEMES:
        rule = _rule(TH.stylesheet(theme), "QComboBox, QSpinBox, QDoubleSpinBox, QLineEdit, "
                                           "QPlainTextEdit, QTextEdit, QAbstractSpinBox")
        text = re.search(r"color:\s*(#[0-9a-fA-F]{6})", rule).group(1)
        ratio = (luminance(text) + 0.05) / (luminance(TH.FIELD_GREY) + 0.05)
        assert ratio > 7, f"{theme}: field text contrast {ratio:.1f}:1"  # WCAG AAA


# --------------------------------------------------------------------------- the display
def test_the_override_decides_whether_windows_are_translucent(monkeypatch):
    monkeypatch.setenv(glass.ENV_OVERRIDE, "0")
    assert glass.compositing_available() is False
    monkeypatch.setenv(glass.ENV_OVERRIDE, "1")
    assert glass.compositing_available() is True


def test_an_x11_display_that_cannot_be_asked_is_not_assumed_to_composite(monkeypatch):
    monkeypatch.delenv("DISPLAY", raising=False)
    assert glass._x11_compositor() is None
    monkeypatch.setenv("DISPLAY", ":987654")          # nothing listens there
    assert glass._x11_compositor() is None


# --------------------------------------------------------------------------- popups
def test_a_menu_is_translucent_frameless_and_shadowless_when_polished(qapp):
    menu = QtWidgets.QMenu()
    menu.addAction("one")
    menu.ensurePolished()
    assert glass.is_dressed(menu)
    assert menu.testAttribute(W.WA_TranslucentBackground)
    flags = menu.windowFlags()
    assert flags & QtCore.Qt.WindowType.FramelessWindowHint
    assert flags & QtCore.Qt.WindowType.NoDropShadowWindowHint


def test_a_shown_menu_paints_nothing_in_its_corners(qapp):
    menu = QtWidgets.QMenu()
    for text in ("one", "two", "three"):
        menu.addAction(text)
    menu.popup(QtCore.QPoint(40, 40))
    qapp.processEvents()
    img = menu.grab().toImage()
    assert img.pixelColor(0, 0).alpha() == 0, "a square corner survived"
    middle = img.pixelColor(img.width() // 2, img.height() // 2)
    assert 150 < middle.alpha() < 255, "the body is glass: mostly opaque, not solid"
    menu.hide()


def test_a_tooltip_is_translucent(qapp):
    host = QtWidgets.QLabel("host")
    host.show()
    host.activateWindow()                       # Qt shows tooltips for the active window only
    QtTest.QTest.qWaitForWindowActive(host, 2000)
    QtWidgets.QToolTip.showText(host.mapToGlobal(QtCore.QPoint(2, 2)), "a tooltip", host)
    for _ in range(5):
        qapp.processEvents()
    tips = [w for w in qapp.allWidgets() if w.metaObject().className() == "QTipLabel"]
    assert tips and all(t.testAttribute(W.WA_TranslucentBackground) for t in tips)
    QtWidgets.QToolTip.hideText()
    host.close()


def test_a_drop_down_list_has_no_square_frame_behind_it(qapp):
    host = QtWidgets.QWidget()
    lay = QtWidgets.QVBoxLayout(host)
    combo = QtWidgets.QComboBox()
    combo.addItems(["dark", "light", "slate", "paper"])
    lay.addWidget(combo)
    host.show()
    combo.showPopup()
    qapp.processEvents()
    pops = [w for w in qapp.topLevelWidgets() if w.isVisible()
            and w.metaObject().className() == "QComboBoxPrivateContainer"]
    assert pops and pops[0].testAttribute(W.WA_TranslucentBackground)
    assert pops[0].grab().toImage().pixelColor(0, 0).alpha() == 0
    combo.hidePopup()
    host.close()


def test_without_a_compositor_a_popup_is_rounded_by_a_mask_instead(qapp, monkeypatch):
    monkeypatch.setenv(glass.ENV_OVERRIDE, "0")
    menu = QtWidgets.QMenu()
    menu.addAction("one")
    menu.ensurePolished()
    assert not menu.testAttribute(W.WA_TranslucentBackground)
    menu.resize(120, 60)
    menu.show()
    qapp.processEvents()
    mask = menu.mask()
    assert not mask.isEmpty()
    assert not mask.contains(QtCore.QPoint(0, 0)), "the corner is cut"
    assert mask.contains(QtCore.QPoint(60, 30))
    menu.hide()


def test_a_popup_that_already_has_a_window_is_left_alone(qapp):
    menu = QtWidgets.QMenu()
    menu.setProperty(glass.NO_GLASS, True)
    menu.ensurePolished()
    assert not glass.is_dressed(menu)
    other = QtWidgets.QMenu()
    other.winId()                                  # its native window exists now
    assert glass.dress_popup(other) is False


# --------------------------------------------------------------------------- dressed windows
def _dialog():
    d = QtWidgets.QDialog()
    d.setWindowTitle("A dressed window")
    lay = QtWidgets.QVBoxLayout(d)
    lay.addWidget(QtWidgets.QLabel("Some words, " * 6))
    box = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Close)
    box.rejected.connect(d.reject)
    lay.addWidget(box)
    return d


def test_a_dressed_dialog_is_a_rounded_translucent_card(qapp):
    d = _dialog()
    before = d.layout().contentsMargins()
    assert glass.dress(d) is True
    assert glass.is_dressed(d) and d.testAttribute(W.WA_TranslucentBackground)
    assert d.windowFlags() & QtCore.Qt.WindowType.FramelessWindowHint
    after = d.layout().contentsMargins()
    assert after.left() > before.left() and after.top() > before.top(), "room for rim and title"
    d.resize(420, 220)
    d.show()
    qapp.processEvents()
    img = d.grab().toImage()
    assert img.pixelColor(0, 0).alpha() == 0, "no square corner round the rounded card"
    body = img.pixelColor(img.width() // 2, img.height() - glass.RIM_ROOM // 2 - 1)
    assert 200 < body.alpha() < 255 and body.red() < 20, "black glass, not solid grey"
    d.close()


def test_a_dressed_dialog_keeps_what_its_title_bar_did(qapp):
    d = _dialog()
    glass.dress(d)
    d.resize(420, 220)
    d.show()
    qapp.processEvents()
    chrome = d._glass_chrome
    close = d.findChild(QtWidgets.QToolButton, glass.CLOSE_NAME)
    assert close is chrome.close and close.isVisible()
    assert d.rect().contains(close.geometry()), "the close mark is inside the card"
    assert "white-space:pre" in close.toolTip(), "the close mark explains itself in a block"
    start = d.pos()
    QtTest.QTest.mousePress(d, QtCore.Qt.MouseButton.LeftButton, pos=QtCore.QPoint(200, 12))
    QtTest.QTest.mouseMove(d, QtCore.QPoint(230, 42))
    QtTest.QTest.mouseRelease(d, QtCore.Qt.MouseButton.LeftButton, pos=QtCore.QPoint(230, 42))
    assert d.pos() != start, "dragging the background moves the window"
    assert chrome.edges_at(QtCore.QPoint(2, d.height() // 2)) == QtCore.Qt.Edge.LeftEdge
    assert chrome.edges_at(QtCore.QPoint(d.width() - 1, d.height() - 1)) == \
        (QtCore.Qt.Edge.RightEdge | QtCore.Qt.Edge.BottomEdge)
    close.click()
    assert not d.isVisible()


def test_dressing_is_idempotent_and_refuses_a_window_already_shown(qapp):
    d = _dialog()
    assert glass.dress(d) and not glass.dress(d)
    assert len(d.findChildren(glass.GlassCard)) == 1
    shown = _dialog()
    shown.show()
    qapp.processEvents()
    assert glass.dress(shown) is False, "rebuilding a visible window's native surface is refused"
    shown.close()
    kept = _dialog()
    kept.setProperty(glass.NO_GLASS, True)
    assert glass.dress(kept) is False


def test_without_a_compositor_a_dressed_dialog_is_opaque_and_masked(qapp, monkeypatch):
    monkeypatch.setenv(glass.ENV_OVERRIDE, "0")
    d = _dialog()
    glass.dress(d)
    assert not d.testAttribute(W.WA_TranslucentBackground)
    d.resize(300, 160)
    d.show()
    qapp.processEvents()
    assert not d.mask().isEmpty() and not d.mask().contains(QtCore.QPoint(0, 0))
    card = d.findChild(glass.GlassCard)
    assert card.body_color().alpha() == 255
    d.close()


# --------------------------------------------------------------------------- on the window
@pytest.fixture(scope="module")
def win(qapp):
    from starplast import app as A
    w = A.Window()
    yield w
    w.console.remove()


def test_every_window_starplast_builds_is_dressed(win, tmp_path):
    import pandas as pd
    table = tmp_path / "t.csv"
    pd.DataFrame({"gene": list(win.nodes.gene_id.iloc[:4]), "x": range(4)}).to_csv(table,
                                                                                  index=False)
    windows = {
        "preferences": win.preferences_dialog(),
        "scoring": win.explain_scoring(),
        "columns": win.choose_export_columns(),
        "import": win.build_import_dialog(str(table)),
        "comparison": win.build_comparison({"note": ""}),
    }
    for name, w in windows.items():
        assert glass.is_dressed(w), name
    explained = win.explain_map()
    assert glass.is_dressed(explained) and explained.isVisible()
    explained.close()
    tree = win.open_slot_tree()
    assert glass.is_dressed(tree)
    win.slot_tree_act.setChecked(False)
    win._open_workflows(0)
    assert glass.is_dressed(win.workflows_dialog)
    win.workflows_dialog.close()


def test_preferences_is_still_a_window_of_its_own_when_dressed(win):
    d = win.open_preferences()
    assert d.windowFlags() & QtCore.Qt.WindowType.Window
    assert d.windowFlags() & QtCore.Qt.WindowType.FramelessWindowHint
    assert not d.isModal()
    d.close()


def test_the_window_applies_the_theme_the_painted_parts_read(win, qapp):
    win.apply_theme("light")
    assert qapp.property("starplastTheme") == "light"
    assert "rgba(255, 255, 255, " in _rule(qapp.styleSheet(), "QMenu")
    win.apply_theme("dark")
    assert "rgba(0, 0, 0, " in _rule(qapp.styleSheet(), "QMenu")
