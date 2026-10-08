"""Gene entry never guesses an alias collision or borrows a changed-table scorecard."""
from types import SimpleNamespace

import pandas as pd
from PyQt6 import QtCore, QtWidgets

from starplast import app as A, organisms as O
from starplast.gene_evidence import build_resolver
from starplast.identity import GeneIndex


def window(qapp, *, subset=False):
    nodes = pd.DataFrame({'gene_id': ['TGME49_100001', 'TGME49_100002'],
                          'symbol': ['SAME', 'SAME'], 'product': ['kinase one', 'kinase two']})
    index = GeneIndex(canonical=set(nodes.gene_id), lookup={gene.lower(): (gene, 'accession') for gene in nodes.gene_id},
                      ambiguous={'same': set(nodes.gene_id)})
    if subset:
        nodes = nodes.iloc[:1].copy()
    resolver = build_resolver(nodes, O.TOXOPLASMA, index=index, index_source='authored-index')
    messages = []
    w = SimpleNamespace(nodes=nodes, species=O.get(O.TOXOPLASMA).species, sel=99,
        search=QtWidgets.QLineEdit(), detail=QtWidgets.QTextBrowser(),
        status=SimpleNamespace(showMessage=messages.append), redraw=lambda: None,
        _gene_resolver=lambda: resolver)
    w._select_search_row = lambda row: setattr(w, 'sel', row)
    return w, messages


def test_alias_collision_requires_choice_even_with_one_available_candidate(qapp):
    w, _ = window(qapp, subset=True)
    w.search.setText('SAME')
    A.Window.do_search(w)
    assert w.sel is None
    assert 'Choose a gene' in w.detail.toPlainText()
    assert 'TGME49_100002 — absent from the current table' in w.detail.toPlainText()
    assert 'authored-index' in w.detail.toPlainText()
    A.Window._detail_link(w, QtCore.QUrl('starplast://select/0'))
    assert w.sel == 0


def test_product_matches_require_choice_and_invalid_row_cannot_select(qapp):
    w, _ = window(qapp)
    w.search.setText('kinase')
    A.Window.do_search(w)
    assert w.sel is None and '2 matches' in w.detail.toPlainText()
    A.Window._detail_link(w, QtCore.QUrl('starplast://select/200'))
    assert w.sel is None
    A.Window._detail_link(w, QtCore.QUrl('starplast://select/1'))
    assert w.sel == 1


def test_unavailable_alias_target_does_not_fall_back_to_product(qapp):
    w, messages = window(qapp, subset=True)
    w.search.setText('TGME49_100002')
    A.Window.do_search(w)
    assert w.sel is None and 'outside the current gene table' in messages[-1]
    assert 'TGME49_100002' in w.detail.toPlainText()


def test_changed_table_cannot_reuse_installed_class_scorecard(qapp, monkeypatch):
    w, messages = window(qapp)
    installed = w.nodes.copy()
    w.nodes.loc[0, 'product'] = 'changed value'
    monkeypatch.setattr(A.pd, 'read_parquet', lambda _: installed)
    called = []
    w._detail_link = called.append
    dialog = SimpleNamespace(close=lambda: called.append('closed'))
    A.Window._gene_label_route(w, dialog, 0, 'symbol', 'SAME')
    assert not called and 'validation required' in messages[-1]


def test_original_table_routes_label_and_quoted_class(qapp, monkeypatch):
    w, _ = window(qapp)
    monkeypatch.setattr(A.pd, 'read_parquet', lambda _: w.nodes.copy())
    called = []
    w._detail_link = lambda url: called.append(url.toString())
    opened = []
    w.open_label_browser = lambda column, label: opened.append((column, label))
    dialog = SimpleNamespace(close=lambda: None)
    A.Window._gene_label_route(w, dialog, 1, 'label/name', 'a/b')
    assert called == ['starplast://class/label%2Fname/a%2Fb']
    assert w._detail_row == 1
    assert opened == [('label/name', 'a/b')]


def test_label_source_filters_organism_and_clears_ambiguous_selection(qapp):
    filters = {key: QtWidgets.QComboBox() for key in ('organism', 'unit')}
    filters['organism'].addItem('All', None)
    filters['organism'].addItem(O.TOXOPLASMA, O.TOXOPLASMA)
    filters['organism'].addItem(O.FALCIPARUM, O.FALCIPARUM)
    filters['unit'].addItem('All', None)
    refreshed = []
    browser = SimpleNamespace(search=QtWidgets.QLineEdit(), filters=filters,
        shown_rows=[{'question': 'original'}, {'question': 'refused'}],
        source_table=QtWidgets.QTableWidget(2, 1), _show_selected=lambda: refreshed.append(True))
    browser.source_table.setCurrentCell(0, 0)
    w = SimpleNamespace(species=O.get(O.TOXOPLASMA).species, open_dataset_browser=lambda: browser)
    A.Window._label_source(w, 'authored-source')
    assert filters['organism'].currentData() == O.TOXOPLASMA
    assert browser.search.text() == 'authored-source'
    assert browser.source_table.currentRow() == -1 and refreshed


def test_reopening_gene_evidence_releases_the_previous_window(qapp, monkeypatch):
    from PyQt6 import sip
    from starplast import dataset_space_inputs, gene_evidence_browser

    class Panel(QtWidgets.QWidget):
        source_requested = QtCore.pyqtSignal(object)
        label_requested = QtCore.pyqtSignal(str)
        class_requested = QtCore.pyqtSignal(str, str)

        def __init__(self, nodes, organism, gene, space, parent):
            super().__init__(parent)

    monkeypatch.setattr(dataset_space_inputs, 'from_sources', lambda *args, **kwargs: object())
    monkeypatch.setattr(gene_evidence_browser, 'GeneEvidenceBrowser', Panel)
    w = QtWidgets.QWidget()
    w.nodes = pd.DataFrame({'gene_id': ['TGME49_100001']})
    w.species, w.sel, w.imports = O.get(O.TOXOPLASMA).species, 0, []
    w._gene_source = w._gene_label_route = lambda *args: None
    A.Window.open_gene_evidence(w)
    previous = w._gene_evidence_dialog
    A.Window.open_gene_evidence(w)
    QtCore.QCoreApplication.sendPostedEvents(None, QtCore.QEvent.Type.DeferredDelete)
    assert sip.isdeleted(previous) and not sip.isdeleted(w._gene_evidence_dialog)
    w.close()
