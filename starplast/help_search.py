"""The search field beside the Help menu, and the exact place each result takes you.

spaCR has a box directly to the right of its Help menu that searches everything the program can do;
this is Starplast's, built the same way so the two tools work alike. Type, and a list drops down of
menu commands, panels and tabs, strategies, settings, information slots, datasets and pages of the
guide (`starplast.help_index` decides what is in it and how it is ranked). Down and Return, or a
click, go to the thing itself:

``action``    the menu command is carried out; a submenu is opened where it lives in the menu bar
``panel``     the dock is shown and raised, on the named tab; Preferences, a guided workflow or the
              slot tree is opened on the right page
``setting``   the panel or Preferences tab holding the control is opened, the control scrolled into
              view, focused and outlined for a moment; a display choice is applied
``strategy``  the Strategies tab is raised with that strategy selected and its Guide showing
``slot``      the slot tree opens on the slot's organism with the slot selected
``dataset``   a card with what the dataset provides, its coverage, citation and source link
``doc``       the page, from this checkout when there is one, scrolled to the heading
``tutorial``  the tutorial page in the browser, at the section

A CHILD OF THE WINDOW, NOT A POPUP WINDOW, for spaCR's reason: a `Qt.Popup` grabs the keyboard the
moment it appears, so the next letter typed would go to the list instead of the box. As a child, the
caret stays in the box, and the list is a glass pane drawn inside the window -- translucent on any
display, compositor or not. Ctrl+Shift+H puts the caret in the box, spaCR's key for the same field.

The index is built the first time something is typed, not when the window opens: reading every
registry costs a fraction of a second nobody should pay for a window they only wanted to look at.
"""
from __future__ import annotations

import os

from PyQt6 import QtCore, QtGui, QtWidgets

from . import help_index as HI
from . import theme as TH

#: objectNames, so the stylesheet can reach the parts and tests can find them.
FIELD_NAME = "HelpSearchField"
POPUP_NAME = "HelpSearchResults"
LIST_NAME = "HelpSearchResultList"
NOTE_NAME = "HelpSearchNote"

#: spaCR's key for the same field. Ctrl+F already means "find a gene" in the find panel.
SHORTCUT = "Ctrl+Shift+H"

#: Pause after a keystroke before searching, so a word typed quickly is searched once.
DEBOUNCE_MS = 120

#: How many rows the list shows.
RESULT_LIMIT = 24

#: Role the entry is stored under on its list item.
ENTRY_ROLE = int(QtCore.Qt.ItemDataRole.UserRole) + 1

#: How long a control found by the search stays outlined.
HIGHLIGHT_MS = 1600


# --------------------------------------------------------------------------- openers
def _status(window, message: str) -> str:
    try:
        window.statusBar().showMessage(message, 5000)
    except Exception:
        pass
    return message


def _find_action(window, path: str):
    """The menu-bar action at a " ▸ "-separated path, and the chain of menus leading to it."""
    parts = path.split(" ▸ ")
    actions, chain = window.menuBar().actions(), []
    found = None
    for depth, part in enumerate(parts):
        found = next((a for a in actions if not a.isSeparator()
                      and HI.clean_label(a.text()) == part), None)
        if found is None:
            return None, chain
        chain.append(found)
        if depth < len(parts) - 1:
            if found.menu() is None:
                return None, chain
            actions = found.menu().actions()
    return found, chain


def open_action(window, e) -> str:
    """Carry out a menu command, or open a submenu where it sits in the menu bar."""
    data = e.data
    action, chain = _find_action(window, data.get("path", ""))
    if action is None:
        return _status(window, f"{e.title} is no longer in the menu bar.")
    if action.menu() is not None:
        bar = window.menuBar()
        bar.setActiveAction(chain[0])
        menu = chain[0].menu()
        for step in chain[1:]:
            menu.setActiveAction(step)
            menu = step.menu() or menu
        return _status(window, f"{data['path']}")
    if not action.isEnabled():
        return _status(window, f"{data['path']} is not available right now.")
    if action.property("helpSearchOpens") and action.isCheckable():
        action.setChecked(True)            # a "show this" toggle: finding it must not close it
    else:
        action.trigger()
    state = (" -- now on" if action.isChecked() else " -- now off") if action.isCheckable() \
        and not (action.actionGroup() and action.actionGroup().isExclusive()) else ""
    return _status(window, f"{data['path']}{state}")


def _dock(window, name: str):
    return next((d for d in window.findChildren(QtWidgets.QDockWidget)
                 if d.windowTitle() == name), None)


def _raise_dock(window, name: str):
    d = _dock(window, name)
    if d is not None:
        d.show()
        d.raise_()
    return d


def _select_tab(tabs, label: str) -> int:
    for i in range(tabs.count()):
        if tabs.tabText(i) == label:
            tabs.setCurrentIndex(i)
            return i
    return -1


def _preference_tabs(window):
    """Open Preferences and return ITS tab widget.

    Read off the open dialog rather than `window.pref_tabs`, which names the tabs of whichever
    Preferences was built last -- not necessarily the one on screen.
    """
    return window.open_preferences().findChild(QtWidgets.QTabWidget)


def open_panel(window, e) -> str:
    """Show a dock on the named tab, or open Preferences, a workflow or the slot tree."""
    data = e.data
    if "preferences" in data:
        _select_tab(_preference_tabs(window), data["preferences"])
        return _status(window, f"Preferences ▸ {data['preferences']}")
    if "workflow" in data:
        window._open_workflows(int(data["workflow"]))
        return _status(window, f"Guided workflows ▸ {e.title}")
    if "slot_tree" in data:
        window.slot_tree_act.setChecked(True)
        tree = getattr(window, "_slot_tree", None) or window.open_slot_tree()
        tree.raise_()
        return _status(window, "Slot tree")
    d = _raise_dock(window, data.get("dock", ""))
    if d is None:
        return _status(window, f"The {e.title} panel is not available.")
    tab = data.get("tab")
    if tab:
        owner = window.strategy_panel if data.get("dock") == "strategies" else window.panel
        _select_tab(owner.tabs, tab)
    return _status(window, f"{e.subtitle} ▸ {e.title}" if tab else f"{e.title} panel")


def _form_field(root, label: str):
    """The field of the form row whose label reads `label`, anywhere under `root`."""
    for form in root.findChildren(QtWidgets.QFormLayout):
        for row in range(form.rowCount()):
            li = form.itemAt(row, QtWidgets.QFormLayout.ItemRole.LabelRole)
            fi = form.itemAt(row, QtWidgets.QFormLayout.ItemRole.FieldRole)
            lw = li.widget() if li is not None else None
            if isinstance(lw, QtWidgets.QLabel) and fi is not None and fi.widget() is not None:
                if HI.clean_label(HI.plain(lw.text())).rstrip(":") == label:
                    return fi.widget()
    return None


def highlight(widget: QtWidgets.QWidget, ms: int = HIGHLIGHT_MS) -> None:
    """Outline a control for a moment, so the eye lands on what the search found."""
    def mark(on):
        try:
            widget.setProperty("helpSearchHit", bool(on))
            widget.style().unpolish(widget)
            widget.style().polish(widget)
            widget.update()
        except RuntimeError:
            pass
    mark(True)
    QtCore.QTimer.singleShot(ms, lambda: mark(False))


def reveal(widget: QtWidgets.QWidget) -> None:
    """Scroll every scroll area holding `widget` until it is in view, focus it, and outline it."""
    parent = widget.parentWidget()
    while parent is not None:
        if isinstance(parent, QtWidgets.QScrollArea):
            parent.ensureWidgetVisible(widget, 40, 60)
        parent = parent.parentWidget()
    widget.setFocus(QtCore.Qt.FocusReason.OtherFocusReason)
    highlight(widget)


def open_setting(window, e) -> str:
    """Land on the control that holds a setting, or apply a display choice."""
    data = e.data
    if "display" in data:
        for label, _options, _current, apply in window.display_choices():
            if label == data["display"]:
                apply(data["option"])
                return _status(window, f"{label}: {data['option']}")
        return _status(window, f"{data['display']} is not a display setting any more.")
    if "preferences" in data:
        tabs = _preference_tabs(window)
        _select_tab(tabs, data["preferences"])
        root = tabs.currentWidget()
        where = f"Preferences ▸ {data['preferences']}"
    else:
        _raise_dock(window, data.get("dock", "analysis"))
        _select_tab(window.panel.tabs, data.get("tab", ""))
        root = window.panel.tabs.currentWidget()
        where = f"Analysis ▸ {data.get('tab', '')}"
    field = _form_field(root, data.get("label", e.title)) if root is not None else None
    if field is None:
        return _status(window, f"{e.title} is not on {where} any more.")
    reveal(field)
    return _status(window, f"{e.title} in {where}")


def open_strategy(window, e) -> str:
    """Raise the Strategies tab with the strategy selected and its Guide showing."""
    panel = getattr(window, "strategy_panel", None)
    if panel is None:
        return _status(window, "The Strategies tab is not available.")
    key = e.data["key"]
    item = panel.items.get(key)
    if item is not None and item.isHidden():
        panel.filter.clear()              # the panel's own filter hid it; the user asked for it
    _raise_dock(window, "strategies")
    panel.select(key)
    panel.tabs.setCurrentIndex(0)
    if item is not None:
        panel.tree.scrollToItem(item)
    return _status(window, f"Strategy {e.title}")


def open_slot(window, e) -> str:
    """Open the slot tree on the slot's organism, with the slot selected and scrolled to."""
    data = e.data
    window.slot_tree_act.setChecked(True)
    tree_window = getattr(window, "_slot_tree", None) or window.open_slot_tree()
    view = tree_window.view
    if view.organism.currentText() != data["organism"]:
        view.organism.setCurrentText(data["organism"])
    if view.filter.text():
        view.filter.clear()
    root = view.tree.invisibleRootItem()
    stack = [root]
    while stack:
        node = stack.pop()
        for i in range(node.childCount()):
            child = node.child(i)
            if child.data(0, QtCore.Qt.ItemDataRole.UserRole) == data["name"]:
                parent = child.parent()
                while parent is not None:
                    parent.setExpanded(True)
                    parent = parent.parent()
                view.tree.setCurrentItem(child)
                view.tree.scrollToItem(child, QtWidgets.QAbstractItemView.ScrollHint.PositionAtCenter)
                tree_window.raise_()
                return _status(window, f"Slot {data['organism']}_{data['name']}")
            stack.append(child)
    tree_window.raise_()
    return _status(window, f"{e.title} is not in the slot tree.")


def dataset_card(window, key: str) -> QtWidgets.QDialog:
    """A glass card saying what one dataset is: provides, coverage, citation, columns, source."""
    from html import escape
    from . import datasets as D
    from . import glass
    d = D.get(key)
    rows = [("Provides", d.provides), ("Type", d.kind), ("Level", d.level.replace("_", " ")),
            ("Coverage", d.coverage), ("Reference", d.citation), ("Accession", d.accession),
            ("PMID", d.pmid), ("Note", d.note)]
    body = "".join(f"<p><b>{escape(k)}</b><br>{escape(str(v))}</p>" for k, v in rows if v)
    if d.columns:
        present = [c for c in d.columns if c in getattr(window, "nodes", {}).columns] \
            if hasattr(window, "nodes") else []
        body += (f"<p><b>Columns</b> ({len(d.columns)}, {len(present)} in the table on the map)"
                 f"<br>{escape(', '.join(d.columns[:24]))}{' …' if len(d.columns) > 24 else ''}</p>")
    if d.url:
        body += f"<p><a href='{escape(d.url)}'>Source data</a> &nbsp; "
    body += (f"<a href='{HI.SITE}datasets.html#dataset-{escape(d.key)}'>Dataset catalogue</a></p>")
    card = glass.message(window, d.name, body, rich=True)
    card.setObjectName("DatasetCard")
    return card


def open_dataset(window, e) -> str:
    """Show the dataset's card."""
    window._dataset_card = dataset_card(window, e.data["key"])
    return _status(window, f"Dataset {e.title}")


class DocViewer(QtWidgets.QDialog):
    """A page of the guide, read from this checkout and scrolled to one heading.

    Local rather than the browser: it works offline, it is the version that matches the code
    running, and it can land on the heading itself. The published page is one button away.
    """

    def __init__(self, path: str, heading: str, url: str, parent=None):
        super().__init__(parent)
        self.setObjectName("DocViewer")
        self.path, self.heading, self.url = path, heading, url
        self.setWindowTitle(f"{os.path.basename(path)} — {heading}")
        lay = QtWidgets.QVBoxLayout(self)
        self.view = QtWidgets.QTextBrowser()
        self.view.setOpenExternalLinks(True)
        from .glass import readable_links
        readable_links(self.view)
        with open(path, encoding="utf8") as fh:
            self.view.setMarkdown(fh.read())
        lay.addWidget(self.view, 1)
        row = QtWidgets.QHBoxLayout()
        web = QtWidgets.QPushButton("Open the published page")
        web.setToolTip(TH.tip("Open this section of the guide on the documentation site, in the "
                              "browser. The copy shown here is the one in this checkout."))
        web.clicked.connect(lambda: QtGui.QDesktopServices.openUrl(QtCore.QUrl(self.url)))
        row.addStretch(1)
        row.addWidget(web)
        lay.addLayout(row)
        self.resize(820, 680)

    def heading_block(self):
        """The text block holding the heading, or an invalid block when it is not on the page."""
        doc = self.view.document()
        block = doc.begin()
        while block.isValid():
            if block.blockFormat().headingLevel() > 0 and \
                    " ".join(block.text().split()) == " ".join(self.heading.split()):
                return block
            block = block.next()
        return doc.findBlockByNumber(-1)

    def scroll_to_heading(self) -> bool:
        """Put the heading at the top of the view. True when the heading was found."""
        block = self.heading_block()
        if not block.isValid():
            return False
        cursor = QtGui.QTextCursor(block)
        self.view.setTextCursor(cursor)
        top = self.view.document().documentLayout().blockBoundingRect(block).top()
        self.view.verticalScrollBar().setValue(int(top))
        return True


def open_doc(window, e) -> str:
    """The guide page at the heading: from this checkout when present, else the published page."""
    from . import glass
    data = e.data
    path = data.get("path", "")
    if path and os.path.exists(path):
        viewer = DocViewer(path, data.get("heading", e.title), data.get("url", ""), window)
        glass.dress(viewer)
        viewer.show()
        QtCore.QTimer.singleShot(0, viewer.scroll_to_heading)
        viewer.scroll_to_heading()
        window._doc_viewer = viewer
        return _status(window, f"{e.subtitle} ▸ {e.title}")
    QtGui.QDesktopServices.openUrl(QtCore.QUrl(data.get("url", HI.SITE)))
    return _status(window, f"Opened {e.title} in the browser")


def open_tutorial(window, e) -> str:
    """The tutorial page in the browser, at its section."""
    data = e.data
    url = QtCore.QUrl.fromLocalFile(data["path"])
    if data.get("anchor"):
        url.setFragment(data["anchor"])
    QtGui.QDesktopServices.openUrl(url)
    return _status(window, f"Opened {e.title}")


#: What choosing each kind of result does. A new kind needs a provider in `help_index` and an
#: opener here, and nothing else changes.
OPENERS = {"action": open_action, "panel": open_panel, "setting": open_setting,
           "strategy": open_strategy, "slot": open_slot, "dataset": open_dataset,
           "doc": open_doc, "tutorial": open_tutorial}


def open_entry(window, e) -> str:
    """Take the user to `e`; returns what the status bar says about it."""
    opener = OPENERS.get(e.kind)
    if opener is None:
        return ""
    try:
        return opener(window, e)
    except Exception as exc:                     # a broken result must not take the window down
        return _status(window, f"Could not open {e.title}: {exc}")


# --------------------------------------------------------------------------- the field
class _RowDelegate(QtWidgets.QStyledItemDelegate):
    """Paints a result as its title, where it lives in a quieter colour, and its kind at the right."""

    def sizeHint(self, option, index):
        """One line of body text plus padding."""
        h = QtGui.QFontMetrics(option.font).height()
        return QtCore.QSize(option.rect.width(), h + 12)

    def paint(self, painter, option, index):
        """Draw the pill for the current row, then the three pieces of text."""
        e = index.data(ENTRY_ROLE)
        app = QtWidgets.QApplication.instance()
        p = TH.palette_for(str(app.property("starplastTheme") or "dark")) if app else TH.DARK
        painter.save()
        painter.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        r = QtCore.QRectF(option.rect).adjusted(1, 1, -1, -1)
        if option.state & QtWidgets.QStyle.StateFlag.State_Selected:
            painter.setPen(QtCore.Qt.PenStyle.NoPen)
            painter.setBrush(QtGui.QColor(p["accent_soft"]))
            painter.drawRoundedRect(r, 5, 5)
        font = QtGui.QFont(option.font)
        small = QtGui.QFont(font)
        if font.pointSizeF() > 0:
            small.setPointSizeF(max(font.pointSizeF() * 0.88, 6.0))
        else:
            small.setPixelSize(max(int(font.pixelSize() * 0.88), 6))
        # Float metrics, rounded UP, with a pixel to spare: an integer advance is rounded to the
        # nearest pixel while eliding compares the exact width, so a title measured at its own
        # rounded width was elided -- "Clusters" came out as "Clu…ers".
        fm, sfm = QtGui.QFontMetricsF(font), QtGui.QFontMetricsF(small)
        kind = HI.KIND_LABEL.get(e.kind, e.kind) if e is not None else ""
        kind_w = sfm.horizontalAdvance(kind) + 6
        x, right = r.left() + 9, r.right() - 9
        title = e.title if e is not None else index.data()
        where = e.subtitle if e is not None else ""
        avail = right - x - kind_w - 12
        title_w = min(fm.horizontalAdvance(title) + 2, avail * 0.62 if where else avail)
        painter.setFont(font)
        painter.setPen(QtGui.QColor(p["fg"]))
        painter.drawText(QtCore.QRectF(x, r.top(), title_w, r.height()),
                         int(QtCore.Qt.AlignmentFlag.AlignVCenter),
                         fm.elidedText(title, QtCore.Qt.TextElideMode.ElideMiddle, title_w))
        if where:
            painter.setFont(small)
            painter.setPen(QtGui.QColor(p["fg_muted"]))
            wx = x + title_w + 12
            painter.drawText(QtCore.QRectF(wx, r.top(), max(right - kind_w - 8 - wx, 0), r.height()),
                             int(QtCore.Qt.AlignmentFlag.AlignVCenter),
                             sfm.elidedText(where, QtCore.Qt.TextElideMode.ElideRight,
                                            max(right - kind_w - 8 - wx, 0.0)))
        painter.setFont(small)
        painter.setPen(QtGui.QColor(p["accent"] if option.state
                                    & QtWidgets.QStyle.StateFlag.State_Selected else p["fg_dim"]))
        painter.drawText(QtCore.QRectF(right - kind_w, r.top(), kind_w, r.height()),
                         int(QtCore.Qt.AlignmentFlag.AlignVCenter
                             | QtCore.Qt.AlignmentFlag.AlignRight), kind)
        painter.restore()


class HelpSearchField(QtWidgets.QLineEdit):
    """The box beside the Help menu: it owns the result list, the index and the keyboard.

    Up, Down, Page Up and Page Down steer the list while the caret stays in the box; Return opens the
    highlighted row; Escape closes the list, and a second Escape gives the window its focus back.
    """

    #: Emitted after a result was opened, with the entry.
    opened = QtCore.pyqtSignal(object)

    def __init__(self, window: QtWidgets.QMainWindow, parent=None):
        super().__init__(parent)
        self.setObjectName(FIELD_NAME)
        self.setClearButtonEnabled(True)
        self.setPlaceholderText("Search Starplast…")
        self.setAccessibleName("Search Starplast")
        self.setToolTip(TH.tip(
            "Search every menu command, panel and tab, strategy, setting, information slot, dataset "
            f"and page of the guide. {SHORTCUT} puts the cursor here from anywhere in the window. "
            "Up and Down choose a result and Return goes to it: the command is carried out, the "
            "panel raised, the strategy selected, the setting scrolled to and outlined."))
        self._window = window
        self._index = None
        self._results = []

        self._popup = QtWidgets.QFrame(window)
        self._popup.setObjectName(POPUP_NAME)
        self._popup.setAttribute(QtCore.Qt.WidgetAttribute.WA_StyledBackground, True)
        self._popup.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        self._popup.hide()
        column = QtWidgets.QVBoxLayout(self._popup)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(0)
        self._note = QtWidgets.QLabel(self._popup)
        self._note.setObjectName(NOTE_NAME)
        self._note.setWordWrap(True)
        self._note.hide()
        column.addWidget(self._note)
        self._list = QtWidgets.QListWidget(self._popup)
        self._list.setObjectName(LIST_NAME)
        self._list.setFocusPolicy(QtCore.Qt.FocusPolicy.NoFocus)
        self._list.setItemDelegate(_RowDelegate(self._list))
        self._list.setVerticalScrollMode(QtWidgets.QAbstractItemView.ScrollMode.ScrollPerPixel)
        self._list.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._list.setMouseTracking(True)
        self._list.itemClicked.connect(self._activate)
        column.addWidget(self._list, 1)

        self._debounce = QtCore.QTimer(self)
        self._debounce.setSingleShot(True)
        self._debounce.setInterval(DEBOUNCE_MS)
        self._debounce.timeout.connect(self.refresh)
        self.textEdited.connect(lambda _t: self._debounce.start())
        self.returnPressed.connect(self.activate_current)

    # ---------------------------------------------------------------- state
    def index(self) -> list:
        """The index, built on first use."""
        if self._index is None:
            self._index = HI.build_index(self._window)
        return self._index

    def set_index(self, entries) -> None:
        """Use `entries` instead of building an index, and search again."""
        self._index = list(entries)
        if self.text().strip():
            self.refresh()

    def invalidate(self) -> None:
        """Forget the index, so the next search reads the window again (after menus change)."""
        self._index = None

    def results(self) -> list:
        """What the list is showing, best first."""
        return list(self._results)

    def popup(self) -> QtWidgets.QFrame:
        """The result pane, so a caller can ask whether it is up."""
        return self._popup

    def note(self) -> str:
        """The line above the list, shown when nothing matched."""
        return self._note.text()

    def type_and_search(self, text: str) -> list:
        """Put `text` in the box and search at once, skipping the debounce. The seam tests drive."""
        self.setText(text)
        self._debounce.stop()
        self.refresh()
        return self.results()

    # ---------------------------------------------------------------- searching
    def refresh(self) -> None:
        """Search for whatever is in the box and show the list."""
        query = self.text().strip()
        if not query:
            self._results = []
            self._list.clear()
            self.hide_popup()
            return
        self._results = HI.search(self.index(), query, limit=RESULT_LIMIT)
        self._list.clear()
        for e in self._results:
            item = QtWidgets.QListWidgetItem(f"{e.title}    {e.subtitle}")
            item.setData(ENTRY_ROLE, e)
            if e.description:
                item.setToolTip(TH.tip(e.description[:600]))
            self._list.addItem(item)
        if self._results:
            self._list.setCurrentRow(0)
            self._note.hide()
            self._note.setText("")
        else:
            self._note.setText(f"Nothing called “{query}”. Try a shorter word, a strategy "
                               "number, a column name or a method.")
            self._note.show()
        self._list.setVisible(bool(self._results))
        self.place_popup()

    def place_popup(self) -> None:
        """Size the pane to its rows and hang it under the box, inside the window."""
        window = self._popup.parentWidget()
        if window is None:
            return
        rows = min(max(self._list.count(), 1), 12)
        row_h = self._list.sizeHintForRow(0) if self._list.count() else 0
        height = (rows * max(row_h, 18) + 12) if self._list.count() else 0
        if self._note.isVisible():
            height += self._note.sizeHint().height()
        width = min(max(self.width() * 2 + 120, 560), max(window.width() - 16, 200))
        corner = self.mapTo(window, QtCore.QPoint(0, self.height() + 4))
        left = max(8, min(corner.x(), window.width() - width - 8))
        self._popup.setGeometry(int(left), int(corner.y()), int(width), int(height))
        self._popup.show()
        self._popup.raise_()

    def hide_popup(self) -> None:
        """Put the list away."""
        self._popup.hide()

    # ---------------------------------------------------------------- choosing
    def activate_current(self) -> None:
        """Open whatever row is highlighted."""
        item = self._list.currentItem()
        if item is not None and self._popup.isVisible():
            self._activate(item)

    def _activate(self, item) -> None:
        e = item.data(ENTRY_ROLE) if item is not None else None
        if e is None:
            return
        self.hide_popup()
        open_entry(self._window, e)
        self.opened.emit(e)

    def keyPressEvent(self, event: QtGui.QKeyEvent) -> None:
        """Steer the list from the box, and let Escape give the window back."""
        key = event.key()
        K = QtCore.Qt.Key
        if key in (K.Key_Down, K.Key_Up, K.Key_PageDown, K.Key_PageUp) and self._list.count():
            if not self._popup.isVisible():
                self.place_popup()
            QtWidgets.QApplication.sendEvent(self._list, event)
            return
        if key == K.Key_Escape:
            if self._popup.isVisible():
                self.hide_popup()
            else:
                self.clear()
                self._window.setFocus()
            event.accept()
            return
        super().keyPressEvent(event)

    def focusOutEvent(self, event) -> None:
        """Close the list when the caret leaves the box; the list never takes focus itself."""
        QtCore.QTimer.singleShot(150, self._hide_if_unfocused)
        super().focusOutEvent(event)

    def _hide_if_unfocused(self) -> None:
        if not self.hasFocus():
            self.hide_popup()


class _FieldPlacer(QtCore.QObject):
    """Keeps the field immediately to the right of the last menu -- which is Help.

    spaCR tried the menu bar's corner widget first and moved away from it: the corner is the far
    right edge of the window, away from the menus the field belongs with. A widget in a menu-bar
    action is driven as a menu item and never takes focus. So the field is an ordinary child of the
    bar, moved to sit after Help whenever the bar is laid out again.
    """

    #: Space between Help and the field.
    GAP = 10

    def __init__(self, bar: QtWidgets.QMenuBar, field: HelpSearchField):
        super().__init__(bar)
        self.bar, self.field = bar, field
        bar.installEventFilter(self)

    def place(self) -> None:
        """Move the field to just after the last menu, as tall as the bar allows."""
        bar, field = self.bar, self.field
        try:
            edges = [bar.actionGeometry(a).right() for a in bar.actions()
                     if not a.isSeparator() and a.isVisible()]
            left = (max(edges) if edges else 0) + self.GAP
            height = max(min(field.sizeHint().height(), bar.height() - 6), 16)
            top = max((bar.height() - height) // 2, 0)
            room = bar.width() - left - 12
            width = max(min(field.maximumWidth(), room), 0)
            field.setGeometry(left, top, width, height)
            field.setVisible(width >= field.minimumWidth())
            field.raise_()
        except RuntimeError:
            pass

    def eventFilter(self, obj, event):
        """Re-place the field whenever the bar changes shape or its menus change."""
        T = QtCore.QEvent.Type
        if obj is self.bar and event.type() in (T.Resize, T.Show, T.LayoutRequest, T.ActionAdded,
                                                T.ActionRemoved, T.ActionChanged, T.FontChange,
                                                T.StyleChange):
            QtCore.QTimer.singleShot(0, self.place)
        return False


def field_of(window) -> HelpSearchField | None:
    """The field installed on `window`, if there is one."""
    field = getattr(window, "_help_search", None)
    return field if isinstance(field, HelpSearchField) else None


def install(window: QtWidgets.QMainWindow) -> HelpSearchField:
    """Put the search field directly to the right of the Help menu. Idempotent."""
    existing = field_of(window)
    if existing is not None:
        return existing
    bar = window.menuBar()
    field = HelpSearchField(window, bar)
    field.setMinimumWidth(150)
    field.setMaximumWidth(300)
    window._help_search = field
    window._help_search_placer = _FieldPlacer(bar, field)
    window._help_search_placer.place()
    field.show()
    return field


def focus_field(window) -> bool:
    """Put the caret in the search box and select what is in it. True when there was a box."""
    field = field_of(window)
    if field is None:
        return False
    if not field.isVisible():
        field.show()
    field.setFocus(QtCore.Qt.FocusReason.ShortcutFocusReason)
    field.selectAll()
    return True
