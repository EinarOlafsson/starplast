"""The search beside the Help menu: what it can find, where it sits, and where each result lands.

The index half is tested without a window, because it is meant to work without one; the field half is
tested on a real window, because "the exact place" is only meaningful there -- a strategy result that
reports success while the Strategies tab shows a different strategy is the failure this is for.
"""
import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")

from PyQt6 import QtCore, QtGui, QtTest, QtWidgets  # noqa: E402

from starplast import help_index as HI  # noqa: E402


# --------------------------------------------------------------------------- the index, headless
def test_every_strategy_is_in_the_index_with_its_key():
    from starplast import strategies as S
    rows = HI.strategy_entries()
    assert {e.data["key"] for e in rows} == {s.key for s in S.catalog()}
    one = next(e for e in rows if e.data["key"] == S.catalog()[0].key)
    s = S.catalog()[0]
    assert one.title == f"{s.number:02d} · {s.name}"
    assert s.family in one.subtitle
    assert s.method in one.description and s.key in one.description


def test_every_slot_and_every_dataset_is_in_the_index():
    from starplast import datasets as D
    from starplast import slots as SL
    slots = HI.slot_entries()
    assert len(slots) == len(SL.all_slots())
    assert {(e.data["organism"], e.data["name"]) for e in slots} == \
        {(s.organism, s.name) for s in SL.all_slots()}
    assert {e.data["key"] for e in HI.dataset_entries()} == {d.key for d in D.REGISTRY}


def test_the_guide_is_indexed_by_heading_and_code_blocks_are_not_headings():
    rows = HI.doc_entries()
    pages = {os.path.basename(e.data.get("path", "")) for e in rows}
    assert "guide.md" in pages and "strategies.md" in pages
    guide = [e for e in rows if e.data.get("path", "").endswith("guide.md")]
    first = open(os.path.join(HI.DOCS_DIR, "guide.md"), encoding="utf8").readline().lstrip("# ").strip()
    assert any(e.title == HI._heading_text(first) for e in guide)
    assert not any(e.title.startswith(("import ", "python ", "pip ")) for e in rows)


def test_without_a_docs_directory_the_published_pages_stand_in(tmp_path):
    rows = HI.doc_entries(str(tmp_path / "nowhere"))
    assert rows and all(e.data["url"].startswith(HI.SITE) for e in rows)
    assert any(e.title == "User guide" for e in rows)


def test_tutorials_are_indexed_with_their_sections():
    rows = HI.tutorial_entries()
    assert any(e.title == "The complete guide to Starplast" for e in rows)
    assert any(e.data.get("anchor") for e in rows), "the complete guide's sections have anchors"


def test_a_name_typed_in_full_beats_a_word_in_a_description():
    rows = [HI.entry("action", "Spin", "View", "rotate"),
            HI.entry("setting", "spin speed", "Preferences", "how fast"),
            HI.entry("doc", "Rotation", "Guide", "press spin to rotate")]
    got = HI.search(rows, "spin")
    assert [e.title for e in got] == ["Spin", "spin speed", "Rotation"]


def test_a_second_word_narrows_the_list():
    rows = HI.static_entries()
    one, two = HI.search(rows, "hdbscan", per_kind=None, limit=999), \
        HI.search(rows, "hdbscan gene", per_kind=None, limit=999)
    assert two and len(two) < len(one)
    assert set(two) <= set(one)


def test_a_method_finds_every_strategy_that_uses_it():
    from starplast import strategies as S
    got = HI.search(HI.static_entries(), "hdbscan", per_kind=None, limit=999)
    found = {e.data["key"] for e in got if e.kind == "strategy"}
    assert found == {s.key for s in S.catalog() if "HDBSCAN" in s.method}


def test_one_kind_cannot_take_every_seat():
    got = HI.search(HI.static_entries(), "transcription", limit=24)
    kinds = [e.kind for e in got]
    assert max(kinds.count(k) for k in set(kinds)) <= HI.PER_KIND_LIMIT
    assert len(set(kinds)) > 1


def test_a_misspelled_name_is_still_found():
    rows = [HI.entry("setting", "cell_diameter", "Analysis"), HI.entry("setting", "save_map", "")]
    assert HI.search(rows, "clldiam")[0].title == "cell_diameter"


def test_tooltip_markup_is_not_what_gets_matched():
    assert HI.plain('<div style="white-space:pre">Spin&#160;it&#160;&#160;</div>') == "Spin it"
    assert HI.clean_label("&Export image…") == "Export image"


# --------------------------------------------------------------------------- on a window
@pytest.fixture(scope="module")
def win():
    from starplast import app as A
    qapp = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])  # noqa: F841 -- kept
    w = A.Window()
    w.resize(1400, 900)
    w.show()
    for _ in range(5):
        QtWidgets.QApplication.processEvents()
    yield w
    w.console.remove()
    w.close()


@pytest.fixture
def field(win):
    f = win._help_search
    yield f
    f.hide_popup()
    f.clear()


def _choose(field, entry):
    """Put `entry` in the list and open it the way a user does: highlight it, press Return."""
    field.set_index([entry])
    field.type_and_search(entry.title.split()[0] if entry.title.split() else entry.title)
    assert field.results() and field.results()[0] == entry
    QtTest.QTest.keyClick(field, QtCore.Qt.Key.Key_Return)
    QtWidgets.QApplication.processEvents()


def test_the_field_sits_in_the_menu_bar_directly_right_of_help(win, field):
    bar = win.menuBar()
    assert field.parentWidget() is bar
    menus = [a for a in bar.actions() if a.menu() is not None]
    assert "Help" in menus[-1].text(), "Help is the last menu, so the field follows it"
    help_right = bar.actionGeometry(menus[-1]).right()
    assert field.isVisible()
    assert 0 < field.geometry().left() - help_right <= 16, "beside Help, not at the far edge"
    assert 0 <= field.geometry().top() and field.geometry().bottom() <= bar.height()
    assert field.objectName() == "HelpSearchField"
    assert len(field.toolTip().split()) >= 15, "the field explains itself like every control"


def test_the_shortcut_focuses_the_field(win, field):
    from starplast import help_search as HS
    assert win.search_act.shortcut() == QtGui.QKeySequence(HS.SHORTCUT) == \
        QtGui.QKeySequence("Ctrl+Shift+H")
    win.view.setFocus()
    win.activateWindow()
    QtTest.QTest.qWaitForWindowActive(win, 2000)
    QtTest.QTest.keyClick(win, QtCore.Qt.Key.Key_H,
                          QtCore.Qt.KeyboardModifier.ControlModifier
                          | QtCore.Qt.KeyboardModifier.ShiftModifier)
    QtWidgets.QApplication.processEvents()
    assert field.hasFocus()


def test_the_help_menu_entry_focuses_the_field_too(win, field):
    win.view.setFocus()
    win.search_act.trigger()
    assert QtWidgets.QApplication.focusWidget() is field or field.hasFocus()


def test_the_window_index_holds_every_menu_command(win, field):
    field.invalidate()
    index = field.index()
    titles = {(e.kind, e.title) for e in index}
    docks = {d.toggleViewAction() for d in win.findChildren(QtWidgets.QDockWidget)}

    def walk(menu):
        for act in menu.actions():
            if act.isSeparator() or act in docks:
                continue
            assert ("action", HI.clean_label(act.text())) in titles, act.text()
            if act.menu():
                walk(act.menu())
    for top in win.menuBar().actions():
        walk(top.menu())


def test_the_window_index_holds_the_panels_tabs_and_settings(win, field):
    field.invalidate()
    index = field.index()
    panels = {e.title for e in index if e.kind == "panel"}
    assert {"Strategies", "Analysis", "Evidence", "Console", "Jobs", "Slot tree"} <= panels
    for i in range(win.panel.tabs.count()):
        assert win.panel.tabs.tabText(i).split("·")[-1].strip() in panels
    settings = {e.title for e in index if e.kind == "setting"}
    assert "spin speed" in settings, "a Preferences control"
    assert "Theme: paper" in settings, "a display choice from the right-click menu"
    assert {"strategy", "slot", "dataset", "doc"} <= {e.kind for e in index}


def test_choosing_a_strategy_selects_it_in_the_strategies_tab(win, field):
    from starplast import strategies as S
    target = S.catalog()[3]
    entry = next(e for e in HI.strategy_entries() if e.data["key"] == target.key)
    win.strategy_panel.select(S.catalog()[0].key)
    win.right_dock.raise_()
    _choose(field, entry)
    assert win.strategy_panel.current.key == target.key
    assert win.strategy_panel.tree.currentItem() is win.strategy_panel.items[target.key]
    assert win.strategy_panel.tabs.currentIndex() == 0, "lands on its Guide"
    assert win.strategies_dock.isVisible()
    assert target.key in win.strategy_panel.guide.toHtml() or \
        target.title.split()[0] in win.strategy_panel.guide.toPlainText()
    assert not field.popup().isVisible()


def test_a_strategy_hidden_by_the_panel_filter_is_still_reached(win, field):
    from starplast import strategies as S
    target = S.catalog()[5]
    win.strategy_panel.filter.setText("no strategy is called this")
    assert win.strategy_panel.items[target.key].isHidden()
    _choose(field, next(e for e in HI.strategy_entries() if e.data["key"] == target.key))
    assert win.strategy_panel.current.key == target.key
    assert not win.strategy_panel.items[target.key].isHidden()


def test_choosing_a_command_carries_it_out(win, field):
    entry = next(e for e in HI.action_entries(win) if e.title == "clusters"
                 and e.subtitle == "View ▸ Color by")
    win.set_color_mode("compartment")
    _choose(field, entry)
    assert win.color_mode == "clusters"
    win.set_color_mode("compartment")


def test_choosing_a_panel_raises_it_on_the_named_tab(win, field):
    entry = next(e for e in HI.panel_entries(win) if e.subtitle == "Analysis ▸ tab"
                 and e.title == "Clusters")
    win.right_dock.raise_()
    _choose(field, entry)
    assert win.panel.tabs.tabText(win.panel.tabs.currentIndex()).endswith("Clusters")
    assert win.analysis_dock.isVisible()


def test_choosing_a_setting_lands_on_the_control_and_outlines_it(win, field):
    entry = next(e for e in HI.setting_entries(win) if e.title == "spin speed")
    _choose(field, entry)
    from starplast import help_search as HS
    prefs = win.preferences_dialog()
    assert prefs.isVisible()
    tabs = prefs.findChild(QtWidgets.QTabWidget)
    assert tabs.tabText(tabs.currentIndex()) == entry.data["preferences"]
    control = HS._form_field(tabs.currentWidget(), "spin speed")
    assert control is not None and control.property("helpSearchHit") is True
    prefs.hide()


def test_choosing_a_display_choice_applies_it(win, field):
    entry = next(e for e in HI.setting_entries(win) if e.title == "Theme: paper")
    _choose(field, entry)
    assert win.theme == "paper"
    win.apply_theme("dark")


def test_choosing_a_slot_opens_the_slot_tree_on_it(win, field):
    from starplast import slots as SL
    slot = next(s for s in SL.all_slots("Pf"))
    entry = next(e for e in HI.slot_entries() if e.data == {"organism": "Pf", "name": slot.name})
    _choose(field, entry)
    tree = win._slot_tree
    assert tree.isVisible() and win.slot_tree_act.isChecked()
    assert tree.view.organism.currentText() == "Pf"
    assert tree.view.tree.currentItem().data(0, QtCore.Qt.ItemDataRole.UserRole) == slot.name
    win.slot_tree_act.setChecked(False)


def test_choosing_a_dataset_shows_what_it_is(win, field):
    from starplast import datasets as D
    d = D.REGISTRY[0]
    _choose(field, next(e for e in HI.dataset_entries() if e.data["key"] == d.key))
    card = win._dataset_card
    assert card.isVisible() and card.windowTitle() == d.name
    assert d.provides.split()[0] in card.text
    card.close()


def test_choosing_a_guide_heading_opens_the_page_at_that_heading(win, field):
    entry = next(e for e in HI.doc_entries() if e.data.get("path", "").endswith("workflows.md")
                 and e.data.get("level") == "2")
    _choose(field, entry)
    viewer = win._doc_viewer
    assert viewer.isVisible()
    assert viewer.scroll_to_heading()
    assert viewer.view.textCursor().block().text().strip() == entry.title
    viewer.close()


def test_choosing_a_tutorial_opens_it_at_its_section(win, field, monkeypatch):
    opened = []
    monkeypatch.setattr(QtGui.QDesktopServices, "openUrl", staticmethod(opened.append))
    entry = next(e for e in HI.tutorial_entries() if e.data.get("anchor"))
    _choose(field, entry)
    assert opened and opened[0].isLocalFile() and opened[0].fragment() == entry.data["anchor"]


def test_the_keyboard_steers_the_list_and_escape_puts_it_away(win, field):
    field.invalidate()
    got = field.type_and_search("export")
    assert len(got) > 2 and field.popup().isVisible()
    # The caret stays in the box while the list is steered: the list is a child of the window that
    # never takes focus, not a popup window that would grab the keyboard from the box.
    assert not field.popup().isWindow()
    assert field._list.focusPolicy() == QtCore.Qt.FocusPolicy.NoFocus
    first = field._list.currentRow()
    QtTest.QTest.keyClick(field, QtCore.Qt.Key.Key_Down)
    assert field._list.currentRow() == first + 1
    QtTest.QTest.keyClick(field, QtCore.Qt.Key.Key_Escape)
    assert not field.popup().isVisible()


def test_a_query_that_finds_nothing_says_so(win, field):
    field.type_and_search("zzqqxx")
    assert not field.results()
    assert "zzqqxx" in field.note()


def test_every_menu_command_explains_itself(win):
    """The rule every control follows here, extended to the menu bar: a tooltip that says more than
    the label, and a status tip for the status bar while the command is pointed at."""
    bare = []

    def walk(menu):
        for act in menu.actions():
            if act.isSeparator():
                continue
            if act.menu() is not None:
                walk(act.menu())
                continue
            plain = HI.plain(act.toolTip())
            if not plain or plain == HI.clean_label(act.text()) or not act.statusTip():
                bare.append(act.text())
    for top in win.menuBar().actions():
        walk(top.menu())
    assert not bare, f"menu commands with no explanation: {bare}"


def test_keyboard_shortcuts_are_listed_from_the_menus(win):
    rows = dict(win.shortcut_rows())
    assert rows[QtGui.QKeySequence("Ctrl+Shift+H").toString(
        QtGui.QKeySequence.SequenceFormat.NativeText)] == "Help ▸ Search Starplast"
    dialog = win.show_shortcuts()
    assert "Search Starplast" in dialog.text
    dialog.close()
