"""The simplified app preserves gene inspection and exposes help for its settings."""
import numpy as np
import pytest
from PyQt6 import QtWidgets
from starplast import app as A


@pytest.fixture(scope="module")
def window():
    qapp = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    win = A.Window(species="Toxoplasma gondii")
    yield win
    win.close()
    qapp.processEvents()


def test_summary_modes_are_absent_and_gene_inspection_still_works(window):
    menus = window.menuBar().findChildren(QtWidgets.QMenu) + window.build_context_menu().findChildren(QtWidgets.QMenu)
    assert not any(m.title() == "Level of detail" for m in menus)
    assert len(window.scatter.pos) == len(window.nodes)
    window.on_pick(0)
    assert window.nodes.gene_id.iloc[0] in window.detail.toPlainText()
    assert "orthogroup" in window.edges
    assert "View: individual genes" in window.describe_state()


def test_display_menu_settings_have_help(window):
    menu = window.build_context_menu()
    display = next(a.menu() for a in menu.actions() if a.text() == "Display")
    for action in display.actions():
        if action.menu():
            assert action.menu().toolTipsVisible()
            for choice in action.menu().actions():
                assert choice.toolTip()
                assert choice.statusTip()
        elif action.isCheckable():
            assert action.toolTip()


def test_preferences_cover_previously_unexplained_numeric_settings(window):
    dialog = window.build_preferences()
    for name in ("spin_speed", "ambient_speed", "ambient_size", "ambient_density"):
        assert len(getattr(window, name).toolTip()) > 40
    dialog.close()


def test_application_uses_a_renderable_icon(window):
    assert not window.windowIcon().isNull()
    assert not window.windowIcon().pixmap(32, 32).isNull()


def test_flying_to_a_category_keeps_individual_gene_coordinates(window):
    before = np.array(window.scatter.pos, copy=True)
    window.fly_to_compartment(window.comp_list.item(0))
    np.testing.assert_array_equal(window.scatter.pos, before)


@pytest.mark.parametrize("species,database", [("Toxoplasma gondii", "toxodb.org/toxo"),
                                               ("Plasmodium falciparum", "plasmodb.org/plasmo")])
def test_gene_evidence_handles_each_species_schema(window, species, database):
    win = A.Window(species=species)
    try:
        for index in (0, win.n // 2, win.n - 1):
            win.on_pick(index)
            gid = win.nodes.gene_id.iloc[index]
            assert gid in win.detail.toPlainText()
            assert f"{database}/app/record/gene/{gid}" in win.detail.toHtml()
    finally:
        win.close()
