"""Translucent black, rounded-corner windows: spaCR's glass on Starplast's palette.

Menus, tooltips, drop-down lists and Starplast's own dialogs are drawn as rounded panes of
translucent black (translucent white on a light theme) rather than solid grey rectangles. This is the
look spaCR already has, and the one the user asked for: "use the opacity, rounded corner windows
whenever you can, black opacity is better than solid gray".

WHAT CAN BE MADE TRANSLUCENT, AND WHEN. A top-level window gets an alpha channel only if it asks for
one BEFORE its native window exists; asking afterwards is silently ignored. Measured on this Qt:

* a menu, a tooltip and a combo box's drop-down list are polished before their native window is
  created, so :class:`GlassStyle` can dress every one of them, whoever built it, from the
  application style's polish hook;
* a dialog is polished only after its native window exists, so a dialog cannot be dressed from the
  style. :func:`dress` is called on each window Starplast builds, before its first show. A
  `QMessageBox` raised through Qt's static functions keeps its native frame and takes the black
  palette from the stylesheet instead of being rebuilt behind the user's back.

WITHOUT A COMPOSITOR, TRANSLUCENCY IS NOT AVAILABLE, and a translucent window there shows garbage or
black squares where its corners should be. :func:`compositing_available` asks the display: Wayland,
Windows and macOS always composite; on X11 the answer is whether anything owns the
``_NET_WM_CM_S<screen>`` selection. When nothing composites, the same windows are drawn OPAQUE
near-black and their corners are cut with a window mask, which rounds them on any display. The
``STARPLAST_TRANSLUCENT`` environment variable (``0`` or ``1``) overrides the check.

A dressed window has no title bar -- the rounded card is the window -- so :func:`dress` gives it back
what the title bar did: the title painted at the top, a close mark, dragging by the background, and
resizing from the edges. Escape still closes a dialog.
"""
from __future__ import annotations

import os
import sys

from PyQt6 import QtCore, QtGui, QtWidgets

#: Environment override for :func:`compositing_available`: ``1`` forces translucency, ``0`` forbids it.
ENV_OVERRIDE = "STARPLAST_TRANSLUCENT"

#: Dynamic property carried by a dressed window. The stylesheet keys its transparent background and
#: its see-through containers off it, so a window that is not dressed keeps its ordinary painting.
GLASS = "glass"

#: A window that must keep its own painting carries this property and is never dressed.
NO_GLASS = "starplastNoGlass"

#: Corner radius of a dressed window, in pixels. The card and the mask share it, so the cut edge and
#: the painted edge are one line.
CARD_RADIUS = 14

#: Corner radius of menus, tooltips and drop-down lists. Mirrors `theme.RADIUS["md"]`.
POPUP_RADIUS = 8

#: Extra margin added round a dressed window's contents. Without it the rim is painted underneath the
#: outermost controls and the edges that resize the window are covered by them.
RIM_ROOM = 10

#: Height of the band a dressed window paints its title in, and where the close mark sits.
TITLE_BAND = 30

#: How close to an edge, in pixels, a press resizes rather than moves.
EDGE = 7

#: Object names, so the stylesheet and the tests can find the parts.
CARD_NAME = "GlassCard"
CLOSE_NAME = "GlassClose"

#: Popup classes the filter dresses, by Qt class name. `QComboBoxPrivateContainer` is the frame that
#: holds a combo box's drop-down list; it is private, so it is recognised by name.
_POPUP_CLASSES = ("QMenu", "QTipLabel", "QComboBoxPrivateContainer")

_CACHE: dict = {}


def _x11_compositor() -> bool | None:
    """Whether an X11 compositing manager is running, or None when the display cannot be asked.

    The EWMH answer: a compositor owns the ``_NET_WM_CM_S<n>`` selection for the screen it manages.
    Read through libX11 directly because Qt 6 no longer exposes the X connection, and read on a
    connection of our own that is closed again at once.
    """
    display_name = os.environ.get("DISPLAY")
    if not display_name:
        return None
    try:
        import ctypes
        import ctypes.util
        name = ctypes.util.find_library("X11") or "libX11.so.6"
        x11 = ctypes.cdll.LoadLibrary(name)
        x11.XOpenDisplay.restype = ctypes.c_void_p
        x11.XOpenDisplay.argtypes = [ctypes.c_char_p]
        x11.XDefaultScreen.argtypes = [ctypes.c_void_p]
        x11.XInternAtom.restype = ctypes.c_ulong
        x11.XInternAtom.argtypes = [ctypes.c_void_p, ctypes.c_char_p, ctypes.c_int]
        x11.XGetSelectionOwner.restype = ctypes.c_ulong
        x11.XGetSelectionOwner.argtypes = [ctypes.c_void_p, ctypes.c_ulong]
        x11.XCloseDisplay.argtypes = [ctypes.c_void_p]
        dpy = x11.XOpenDisplay(display_name.encode())
        if not dpy:
            return None
        try:
            screen = x11.XDefaultScreen(dpy)
            atom = x11.XInternAtom(dpy, f"_NET_WM_CM_S{screen}".encode(), 0)
            return bool(atom and x11.XGetSelectionOwner(dpy, atom))
        finally:
            x11.XCloseDisplay(dpy)
    except Exception:
        return None


def compositing_available(refresh: bool = False) -> bool:
    """Whether top-level windows can be translucent on this display.

    Decided once per process and cached: the answer only changes if the user starts or stops a
    compositor while Starplast is running, and a window already built keeps what it was built with.

    :param refresh: ask the display again instead of using the cached answer.
    """
    forced = os.environ.get(ENV_OVERRIDE, "").strip()
    if forced in ("0", "1"):
        return forced == "1"
    if not refresh and "compositing" in _CACHE:
        return _CACHE["compositing"]
    platform = QtGui.QGuiApplication.platformName() if QtGui.QGuiApplication.instance() else ""
    platform = (platform or os.environ.get("QT_QPA_PLATFORM", "")).lower()
    if platform.startswith("wayland") or platform in ("windows", "cocoa", "offscreen"):
        # Offscreen renders into images that keep their alpha channel, which is what the suite and
        # the screenshot builders grab; it is treated as compositing so they exercise the real look.
        answer = True
    elif platform == "xcb" or (not platform and sys.platform.startswith("linux")):
        answer = bool(_x11_compositor())
    else:
        answer = False
    _CACHE["compositing"] = answer
    return answer


def rounded_region(rect: QtCore.QRect, radius: float) -> QtGui.QRegion:
    """A rounded rectangle as a region, for cutting a window's corners without a compositor.

    Built at four times the size and scaled down, so the polygon Qt fills has enough points on each
    arc that the corner reads as a curve rather than as three steps.
    """
    step = 4.0
    path = QtGui.QPainterPath()
    path.addRoundedRect(QtCore.QRectF(0, 0, rect.width() * step, rect.height() * step),
                        radius * step, radius * step)
    polygon = QtGui.QTransform().scale(1 / step, 1 / step).map(path.toFillPolygon())
    return QtGui.QRegion(polygon.toPolygon())


def _palette() -> tuple:
    """The palette and whether it is dark, for the theme the application is showing now."""
    from . import theme as TH
    app = QtWidgets.QApplication.instance()
    name = str(app.property("starplastTheme") or "dark") if app is not None else "dark"
    return TH.palette_for(name), TH.is_dark(name), name


class _MaskKeeper(QtCore.QObject):
    """Re-cuts a popup's rounded mask every time it is resized, for the opaque fallback."""

    def __init__(self, widget: QtWidgets.QWidget, radius: float):
        super().__init__(widget)
        self._radius = radius
        widget.installEventFilter(self)

    def eventFilter(self, obj, event):
        """Cut the corners again at the new size; never consumes the event."""
        if event.type() in (QtCore.QEvent.Type.Resize, QtCore.QEvent.Type.Show):
            rect = obj.rect()
            if rect.width() > 0 and rect.height() > 0:
                obj.setMask(rounded_region(rect, self._radius))
        return False


class GlassStyle(QtWidgets.QProxyStyle):
    """The application style, unchanged, plus one thing: it dresses every popup as Qt polishes it.

    One hook for the whole program rather than a line at every `QMenu(...)`: menus are built all over
    it -- the menu bar, the map's right-click menu, the category lists, the jobs panel -- and by Qt
    itself for a combo box or a tooltip, and a look applied by hand is a look missing from the next.

    A STYLE, NOT AN APPLICATION EVENT FILTER. A filter sees every event of every object, including
    the ones a widget receives while it is being destroyed, and PyQt must wrap each of those objects
    to hand it to Python -- which, for an object half torn down, crashed the process. It did, in the
    suite, on every module that closes a window. `polish` is only ever called on a live widget, at
    the moment it is about to need its style, which is exactly the moment wanted here.

    Only popups whose native window does not exist yet are touched, which is every one of them the
    first time it is polished; one that already has a window is left as it is rather than rebuilt.
    """

    def polish(self, target):
        """Dress a popup being polished, then let the wrapped style polish it as it always did."""
        if isinstance(target, QtWidgets.QWidget):
            try:
                if target.isWindow() and target.metaObject().className() in _POPUP_CLASSES:
                    dress_popup(target)
            except RuntimeError:
                pass
        return super().polish(target)


def dress_popup(widget: QtWidgets.QWidget) -> bool:
    """Make one popup translucent and frameless, or give it a rounded mask. True if dressed.

    Must run before the popup's native window exists; afterwards the alpha channel cannot be added,
    so a created popup is left alone and False is returned.
    """
    if widget.property(GLASS) or widget.property(NO_GLASS):
        return False
    if widget.testAttribute(QtCore.Qt.WidgetAttribute.WA_WState_Created):
        return False
    widget.setProperty(GLASS, True)
    if compositing_available():
        widget.setAttribute(QtCore.Qt.WidgetAttribute.WA_TranslucentBackground, True)
        if widget.metaObject().className() == "QComboBoxPrivateContainer":
            # The frame round a drop-down list paints itself as a square menu panel, whatever the
            # application sheet says about it -- measured, only a rule on the frame itself clears
            # it. The rounded glass is the list inside, styled by `QComboBox QAbstractItemView`.
            # Opaque displays keep the frame's fill: there is nothing behind it to show through.
            widget.setStyleSheet(
                "QComboBoxPrivateContainer { background: transparent; border: none; }")
        if isinstance(widget, QtWidgets.QMenu):
            # A menu is already frameless; the shadow hint stops a window manager drawing a square
            # shadow round the rounded pane, which reads as the grey rectangle this replaces.
            widget.setWindowFlags(widget.windowFlags()
                                  | QtCore.Qt.WindowType.FramelessWindowHint
                                  | QtCore.Qt.WindowType.NoDropShadowWindowHint)
    else:
        _MaskKeeper(widget, POPUP_RADIUS)
    return True


def install(app: QtWidgets.QApplication | None = None) -> GlassStyle | None:
    """Put the popup-dressing style under the application, once. Returns it, or None without an app.

    The wrapped style is the platform's default -- what the application would have used anyway --
    and the application palette is put back afterwards, because `setStyle` otherwise replaces it
    with the style's standard palette and every palette-drawn colour would move.
    """
    app = app or QtWidgets.QApplication.instance()
    if app is None:
        return None
    existing = getattr(app, "_starplast_glass_style", None)
    if existing is not None:
        return existing
    palette = QtGui.QPalette(app.palette())
    style = GlassStyle()
    app.setStyle(style)
    app.setPalette(palette)
    app._starplast_glass_style = style
    return style


# --------------------------------------------------------------------------- dressed windows
class GlassCard(QtWidgets.QWidget):
    """The rounded body of a dressed window: translucent fill, a hairline rim, and the title.

    A child laid BEHIND the window's own contents and kept at its size, holding no layout, so it
    cannot disturb one. It ignores the mouse; presses on it reach the window, which is what makes the
    background draggable.
    """

    def __init__(self, window: QtWidgets.QWidget, radius: int = CARD_RADIUS, title: bool = True):
        super().__init__(window)
        self.setObjectName(CARD_NAME)
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_TransparentForMouseEvents, True)
        self.radius = radius
        self.show_title = title

    def body_color(self) -> QtGui.QColor:
        """Black glass on a dark theme, white glass on a light one; opaque without a compositor."""
        from . import theme as TH
        _p, _dark, name = _palette()
        # From the float tuple, not the stylesheet string: QColor does not parse `rgba(...)`, and an
        # unparsed colour paints as opaque black -- the one failure that looks almost right here.
        return QtGui.QColor.fromRgbF(*TH.glass_rgba(name, "dialog", compositing_available()))

    def paintEvent(self, _event):
        """Fill the rounded body, stroke the rim, and draw the window title in the top band."""
        from . import theme as TH
        p, _dark, name = _palette()
        painter = QtGui.QPainter(self)
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        rect = QtCore.QRectF(self.rect())
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        painter.setBrush(self.body_color())
        painter.drawRoundedRect(rect, self.radius, self.radius)
        painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
        painter.setPen(QtGui.QPen(QtGui.QColor(TH.rim(name)), 1.0))
        painter.drawRoundedRect(rect.adjusted(0.5, 0.5, -0.5, -0.5), self.radius, self.radius)
        title = self.window().windowTitle() if self.show_title else ""
        if title:
            font = QtGui.QFont(self.font())
            font.setBold(True)
            painter.setFont(font)
            painter.setPen(QtGui.QColor(p["fg_muted"]))
            band = QtCore.QRectF(RIM_ROOM + 6, 0, rect.width() - 2 * RIM_ROOM - 48, TITLE_BAND)
            painter.drawText(band, int(QtCore.Qt.AlignmentFlag.AlignVCenter
                                       | QtCore.Qt.AlignmentFlag.AlignLeft),
                             QtGui.QFontMetrics(font).elidedText(
                                 title, QtCore.Qt.TextElideMode.ElideRight, int(band.width())))
        painter.end()


class _Chrome(QtCore.QObject):
    """What the title bar used to do for a dressed window: keep the card sized, move, resize.

    An event filter on the window rather than overridden handlers, because the window is somebody
    else's class -- a `QDialog`, a `QMainWindow`, a dialog built in another module -- and dressing it
    must not require subclassing it.
    """

    def __init__(self, window: QtWidgets.QWidget, card: GlassCard, close: QtWidgets.QToolButton):
        super().__init__(window)
        self.window, self.card, self.close = window, card, close
        self._drag = None          # (global press point, window top-left) during a manual move
        self._resize = None        # (edges, global press point, geometry) during a manual resize
        window.setMouseTracking(True)
        window.installEventFilter(self)
        self.fit()

    def fit(self):
        """Keep the card covering the window, the close mark in its corner, and the mask cut."""
        w = self.window
        rect = w.rect()
        self.card.setGeometry(rect)
        self.card.lower()
        if self.close is not None:
            size = self.close.sizeHint()
            self.close.setGeometry(rect.width() - size.width() - RIM_ROOM,
                                   max((TITLE_BAND - size.height()) // 2, 2),
                                   size.width(), size.height())
            self.close.raise_()
        if not compositing_available() and rect.width() > 0 and rect.height() > 0:
            w.setMask(rounded_region(rect, self.card.radius))

    def edges_at(self, pos: QtCore.QPoint) -> QtCore.Qt.Edge:
        """Which edges a point in window coordinates is close enough to resize."""
        w, h = self.window.width(), self.window.height()
        edges = QtCore.Qt.Edge(0)
        if pos.x() <= EDGE:
            edges |= QtCore.Qt.Edge.LeftEdge
        elif pos.x() >= w - EDGE:
            edges |= QtCore.Qt.Edge.RightEdge
        if pos.y() <= EDGE:
            edges |= QtCore.Qt.Edge.TopEdge
        elif pos.y() >= h - EDGE:
            edges |= QtCore.Qt.Edge.BottomEdge
        return edges

    @staticmethod
    def cursor_for(edges) -> QtCore.Qt.CursorShape:
        """The resize cursor for a set of edges, or the ordinary arrow."""
        E, C = QtCore.Qt.Edge, QtCore.Qt.CursorShape
        diagonal = {(E.LeftEdge | E.TopEdge), (E.RightEdge | E.BottomEdge)}
        anti = {(E.RightEdge | E.TopEdge), (E.LeftEdge | E.BottomEdge)}
        if edges in diagonal:
            return C.SizeFDiagCursor
        if edges in anti:
            return C.SizeBDiagCursor
        if edges & (E.LeftEdge | E.RightEdge):
            return C.SizeHorCursor
        if edges & (E.TopEdge | E.BottomEdge):
            return C.SizeVerCursor
        return C.ArrowCursor

    def eventFilter(self, obj, event):
        """Resize and show refit the card; presses on the background move or resize the window."""
        if obj is not self.window:
            return False
        t = event.type()
        T = QtCore.QEvent.Type
        if t in (T.Resize, T.Show):
            self.fit()
        elif t == T.WindowTitleChange:
            self.card.update()
        elif t == T.MouseButtonPress and event.button() == QtCore.Qt.MouseButton.LeftButton:
            return self._press(event)
        elif t == T.MouseMove:
            return self._move(event)
        elif t == T.MouseButtonRelease:
            self._drag = self._resize = None
        return False

    def _press(self, event) -> bool:
        pos = event.position().toPoint()
        gpos = event.globalPosition().toPoint()
        edges = self.edges_at(pos)
        handle = self.window.windowHandle()
        resizable = self.window.minimumSize() != self.window.maximumSize()
        if edges and resizable:
            if handle is not None and handle.startSystemResize(edges):
                return True
            self._resize = (edges, gpos, QtCore.QRect(self.window.geometry()))
            return True
        # The compositor's own move where it offers one: it snaps and crosses screens properly.
        if handle is not None and handle.startSystemMove():
            return True
        self._drag = (gpos, self.window.pos())
        return True

    def _move(self, event) -> bool:
        gpos = event.globalPosition().toPoint()
        if self._drag is not None:
            start, origin = self._drag
            self.window.move(origin + (gpos - start))
            return True
        if self._resize is not None:
            edges, start, geo = self._resize
            d = gpos - start
            g = QtCore.QRect(geo)
            E = QtCore.Qt.Edge
            if edges & E.LeftEdge:
                g.setLeft(g.left() + d.x())
            if edges & E.RightEdge:
                g.setRight(g.right() + d.x())
            if edges & E.TopEdge:
                g.setTop(g.top() + d.y())
            if edges & E.BottomEdge:
                g.setBottom(g.bottom() + d.y())
            if g.width() >= self.window.minimumWidth() and g.height() >= self.window.minimumHeight():
                self.window.setGeometry(g)
            return True
        edges = self.edges_at(event.position().toPoint())
        self.window.setCursor(self.cursor_for(edges))
        return False


def _widen_margins(window: QtWidgets.QWidget, title: bool) -> None:
    """Push the contents in far enough for the rim, the resize band and the title band."""
    top = TITLE_BAND if title else RIM_ROOM
    layout = window.layout()
    if isinstance(window, QtWidgets.QMainWindow) or layout is None:
        m = window.contentsMargins()
        window.setContentsMargins(m.left() + RIM_ROOM, m.top() + top,
                                  m.right() + RIM_ROOM, m.bottom() + RIM_ROOM)
        return
    m = layout.contentsMargins()
    layout.setContentsMargins(m.left() + RIM_ROOM, max(m.top(), 4) + top - 4,
                              m.right() + RIM_ROOM, m.bottom() + RIM_ROOM)


def dress(window: QtWidgets.QWidget, title: bool = True, closable: bool = True) -> bool:
    """Make a window Starplast built into a rounded glass card. True if it was dressed.

    Call it after the window's layout exists and before its first show. A window already shown (its
    native window exists) is left with its frame: the alpha channel cannot be added any more, and
    rebuilding the native window under a visible dialog is how a dialog vanishes on some window
    managers. Idempotent, and a window carrying the ``starplastNoGlass`` property is never touched.

    :param window: a top-level `QDialog` or `QMainWindow`.
    :param title: paint the window title in the top band (dropping the title bar would lose it).
    :param closable: add a close mark in the top-right corner.
    """
    if window.property(GLASS) or window.property(NO_GLASS):
        return False
    if window.testAttribute(QtCore.Qt.WidgetAttribute.WA_WState_Created):
        return False
    translucent = compositing_available()
    if translucent:
        window.setAttribute(QtCore.Qt.WidgetAttribute.WA_TranslucentBackground, True)
    window.setWindowFlags(window.windowFlags() | QtCore.Qt.WindowType.FramelessWindowHint)
    window.setProperty(GLASS, True)
    _widen_margins(window, title)
    card = GlassCard(window, title=title)
    card.lower()
    close = None
    if closable:
        from . import theme as TH
        close = QtWidgets.QToolButton(window)
        close.setObjectName(CLOSE_NAME)
        close.setText("✕")
        close.setAutoRaise(True)
        close.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        close.setToolTip(TH.tip("Close this window. Escape does the same, and the window can be "
                                "moved by dragging its background and resized from its edges."))
        close.clicked.connect(window.close)
    window._glass_chrome = _Chrome(window, card, close)
    # Re-read the stylesheet now that the window carries the property its rules are keyed on.
    style = window.style()
    style.unpolish(window)
    style.polish(window)
    return True


def is_dressed(widget: QtWidgets.QWidget) -> bool:
    """Whether a window or popup has been given the glass treatment."""
    return bool(widget.property(GLASS))


def readable_links(widget: QtWidgets.QWidget) -> None:
    """Colour a widget's links with the theme's accent.

    Qt's default link colour is a dark blue chosen for a white page; on black glass it is close to
    invisible, and a link nobody can see is a reference nobody follows.
    """
    p, _dark, _name = _palette()
    pal = widget.palette()
    for group in (QtGui.QPalette.ColorGroup.Active, QtGui.QPalette.ColorGroup.Inactive):
        pal.setColor(group, QtGui.QPalette.ColorRole.Link, QtGui.QColor(p["accent"]))
        pal.setColor(group, QtGui.QPalette.ColorRole.LinkVisited, QtGui.QColor(p["accent_hi"]))
    widget.setPalette(pal)


def message(parent: QtWidgets.QWidget | None, title: str, text: str,
            rich: bool = False) -> QtWidgets.QDialog:
    """A glass window carrying a passage of explanation, shown beside the map and returned.

    Not modal: an explanation is read while looking at the thing it explains, and a modal box froze
    the map behind it -- the reason Preferences stopped being modal as well.
    """
    d = QtWidgets.QDialog(parent)
    d.setWindowTitle(title)
    d.setObjectName("GlassMessage")
    lay = QtWidgets.QVBoxLayout(d)
    body = QtWidgets.QLabel(text)
    body.setObjectName("GlassMessageText")
    body.setWordWrap(True)
    body.setTextFormat(QtCore.Qt.TextFormat.RichText if rich else QtCore.Qt.TextFormat.PlainText)
    body.setOpenExternalLinks(True)
    body.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextBrowserInteraction)
    readable_links(body)
    body.setMinimumWidth(460)
    body.setMaximumWidth(620)
    lay.addWidget(body)
    buttons = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Close)
    buttons.rejected.connect(d.reject)
    lay.addWidget(buttons)
    d.text = text
    dress(d)
    d.show()
    return d
