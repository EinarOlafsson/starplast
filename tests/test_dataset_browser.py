"""Tiny dataset UI fixtures preserve scoped filters, values, links and paging."""
import json

import numpy as np
import pandas as pd
import pytest
from PyQt6 import QtCore, QtGui, QtWidgets

from starplast import datasets as D, inventory as I, organisms as O
from starplast.dataset_browser import DatasetBrowser
from starplast.dataset_space import DatasetSpace


@pytest.fixture(scope='session')
def app(qapp):
    # One QApplication per process: pytest-qt also owns the later slot-view tests.
    yield qapp


def fixture(count=3, *, no_origins=False):
    ids = [f'TGME49_{100000+i}' for i in range(count)]
    measured = [0., np.nextafter(1., 2.), None] + [float(i) for i in range(3, count)]
    genes = pd.DataFrame({'gene_id': ids, 'measurement': measured[:count],
        'flag': pd.array(([False, True, None] + [None] * count)[:count], dtype='boolean'),
        'imported_measurement': ([7., None, 0.] + [None] * count)[:count]})
    host = pd.DataFrame({'host_id': ['fixture_host1', 'fixture_host2'], 'host_measurement': [0., None]})
    sources = [D.Dataset(key, name, level, 'software_fixture', 'Authored test only',
        organism=organism, columns=(column,), url='https://example.org/' + key)
        for key, name, level, organism, column in [
            ('numeric_fixture', 'Numeric fixture', 'transcription', O.TOXOPLASMA, 'measurement'),
            ('flag_fixture', 'Flag fixture', 'reference', O.TOXOPLASMA, 'flag'),
            ('import_fixture', 'Imported fixture', 'reference', O.TOXOPLASMA, 'imported_measurement'),
            ('host_fixture', 'Host fixture', 'translation', O.HUMAN, 'host_measurement')]]
    tables = {(O.TOXOPLASMA, 'gene'): genes, (O.HUMAN, 'protein'): host}
    refusal = {'source_id': 'numeric_fixture', 'organism': O.TOXOPLASMA, 'unit': 'gene',
               'question': 'different synthetic quantity', 'reason': 'Wrong synthetic quantity; scoped refusal'}
    inventory = I.build_inventory(tables, sources=sources, catalog=[], refusals=[refusal])
    inventory['contexts'] = inventory.source_id.map(lambda key: ['host context' if key == 'host_fixture' else 'fixture context'])
    origins = {(key, organism, unit, ''): {'kind': kind, 'evidence_grade': grade,
        'lineage': 'Explicit synthetic fixture declaration; no authenticated biological source'}
        for key, organism, unit, kind, grade in [
            ('numeric_fixture', O.TOXOPLASMA, 'gene', 'measured', 'direct_experiment'),
            ('flag_fixture', O.TOXOPLASMA, 'gene', 'transferred', 'orthology_transfer'),
            ('import_fixture', O.TOXOPLASMA, 'gene', 'user_imported', 'unresolved'),
            ('host_fixture', O.HUMAN, 'protein', 'predicted', 'prediction')]}
    return DatasetSpace(inventory, tables=tables, sources=sources, catalog=[],
                        origins=None if no_origins else origins), genes, host


def select(browser, source_id, question=''):
    index = next(i for i, row in enumerate(browser.shown_rows)
                 if row['source_id'] == source_id and row['question'] == question)
    browser.source_table.setCurrentCell(index, 0)
    browser._show_selected()


def set_filter(browser, key, value):
    index = browser.filters[key].findData(value)
    assert index >= 0
    browser.filters[key].setCurrentIndex(index)


def test_all_facets_and_text_reproduce_inventory_counts(app):
    space, _, _ = fixture()
    browser = DatasetBrowser(space, organism=O.TOXOPLASMA)
    try:
        for widget in [*browser.filters.values(), browser.search, browser.page_previous, browser.page_next]:
            assert widget.toolTip().strip()
        assert browser.source_table.rowCount() == len(space.filter_rows(organism=O.TOXOPLASMA))
        for key, value in [('organism', O.HUMAN), ('unit', 'protein'),
            ('family', 'translation / software_fixture'), ('context', 'host context'), ('status', 'installed')]:
            set_filter(browser, key, value)
            current = {name: combo.currentData() for name, combo in browser.filters.items()}
            assert browser.shown_rows == space.filter_rows(**current)
            assert browser.source_table.rowCount() == len(browser.shown_rows) == 1
        browser.search.setText('not a recorded source')
        assert browser.source_table.rowCount() == 0 and browser.card.view is None
        assert browser.entity_table.model().rowCount() == 0
        browser.search.setText('HOST FIXTURE')
        assert browser.source_table.rowCount() == 1
    finally:
        browser.close()


def test_source_selection_keeps_exact_values_counts_zero_false_and_missingness(app):
    space, genes, _ = fixture()
    browser = DatasetBrowser(space, organism=O.TOXOPLASMA)
    try:
        select(browser, 'numeric_fixture')
        pd.testing.assert_frame_equal(browser.entity_frame.reset_index(drop=True),
            genes[['gene_id', 'measurement']], check_exact=True, check_flags=False)
        model = browser.entity_table.model()
        assert model.data(model.index(0, 1)) == '0.0'
        assert model.data(model.index(1, 1)) == str(np.nextafter(1., 2.))
        assert model.data(model.index(2, 1)) == 'Unavailable'
        assert dict(browser.card.view.counts)['stored_any_rows'] == 2
        assert dict(browser.card.view.counts)['missing_unknown_rows'] == 1
        assert browser.card.view.kind == 'evidence' and not browser.card.view.metrics
        assert 'accuracy' in browser.summary.text()
        select(browser, 'flag_fixture')
        assert model is not browser.entity_table.model()
        model = browser.entity_table.model()
        assert model.data(model.index(0, 1)) == 'False'
        assert dict(browser.card.view.counts)['false_cells'] == 1
        assert browser.entity_frame.flag.tolist() == genes.flag.tolist()
    finally:
        browser.close()


def test_scoped_refusal_cannot_borrow_installed_records_or_old_card(app):
    space, _, _ = fixture()
    browser = DatasetBrowser(space, organism=O.TOXOPLASMA)
    try:
        select(browser, 'numeric_fixture')
        old_snapshot = json.loads(browser.card.export_json())['snapshot']
        select(browser, 'numeric_fixture', 'different synthetic quantity')
        assert browser.entity_frame.empty and browser.entity_frame.attrs['status'] == 'unavailable'
        assert 'Wrong synthetic quantity' in browser.value_note.text()
        assert browser.card.view.status == 'rejected'
        assert json.loads(browser.card.export_json())['snapshot'] != old_snapshot
        assert browser.page_note.text() == '0 shown; stored record total unavailable'
        assert not browser.page_previous.isEnabled() and not browser.page_next.isEnabled()
        opened = []
        browser.entity_requested.connect(lambda *values: opened.append(values))
        browser._open_entity(QtCore.QModelIndex())
        assert not opened
    finally:
        browser.close()


def test_origin_labels_source_link_navigation_and_qualified_entity_signal(app, monkeypatch):
    space, _, _ = fixture()
    browser = DatasetBrowser(space)
    try:
        assert {row['origin'] for row in browser.shown_rows} >= {'measured', 'transferred', 'predicted', 'user_imported'}
        opened_urls = []
        monkeypatch.setattr(QtGui.QDesktopServices, 'openUrl', lambda url: opened_urls.append(url.toString()) or True)
        select(browser, 'host_fixture')
        assert 'Double-click a gene identifier' in browser.value_note.text()
        assert 'other units remain in this table' in browser.value_note.text()
        assert browser.card.navigate('https://example.org/host_fixture')
        assert opened_urls == ['https://example.org/host_fixture']
        assert not browser.card.navigate('https://example.org/numeric_fixture')
        assert browser.card.navigate('scorecard:detail/source')
        requested = []
        browser.entity_requested.connect(lambda *values: requested.append(values))
        browser._open_entity(browser.entity_table.model().index(0, 1))
        assert requested == [(O.HUMAN, 'protein', 'fixture_host1')]
        browser.entity_frame.attrs['organism'] = O.TOXOPLASMA
        browser._open_entity(browser.entity_table.model().index(0, 0))
        assert len(requested) == 1
    finally:
        browser.close()


def test_pagination_makes_every_stored_record_reachable_without_changing_counts(app):
    space, genes, _ = fixture(401)
    browser = DatasetBrowser(space, organism=O.TOXOPLASMA)
    try:
        select(browser, 'numeric_fixture')
        selected = []
        while True:
            assert len(browser.entity_frame) <= 200
            selected.extend(browser.entity_frame.gene_id.tolist())
            assert browser.entity_frame.attrs['total_rows'] == 401
            assert dict(browser.card.view.counts)['table_rows'] == 401
            if not browser.page_next.isEnabled():
                break
            browser.page_next.click()
        assert selected == genes.gene_id.tolist()
        assert browser.page_note.text() == '401–401 of 401 stored records'
        browser.page_previous.click()
        assert browser.page_start == 200 and len(browser.entity_frame) == 200
        select(browser, 'flag_fixture')
        assert browser.page_start == 0
    finally:
        browser.close()


def test_refresh_preserves_available_filters_and_never_borrows_an_unknown_organism(app):
    space, _, _ = fixture()
    browser = DatasetBrowser(space, organism=O.TOXOPLASMA)
    try:
        select(browser, 'import_fixture')
        browser.set_space(space)
        assert browser.selected_row['source_id'] == 'import_fixture'
        browser.set_space(space, organism='unavailable fixture organism')
        assert browser.filters['organism'].currentData() == 'unavailable fixture organism'
        assert not browser.shown_rows and browser.card.view is None
        assert browser.entity_frame.empty
    finally:
        browser.close()


def test_readable_widths_and_exact_tooltips_preserve_unverified_metadata(app):
    space, _, _ = fixture(no_origins=True)
    browser = DatasetBrowser(space, organism=O.TOXOPLASMA)
    try:
        select(browser, 'numeric_fixture')
        row = browser.source_table.currentRow()
        origin = browser.source_table.item(row, 3)
        assert origin.text() == 'Unverified registry'
        assert origin.toolTip() == 'registry_asserted_unverified'
        assert browser.selected_row['origin'] == 'registry_asserted_unverified'
        assert 'registry_asserted_unverified' in browser.card.export_json()
        assert browser.source_table.columnWidth(0) == 220
        assert browser.source_table.horizontalHeader().stretchLastSection()
        assert browser.source_table.horizontalHeaderItem(7).toolTip() == 'Stored / table rows'
        assert browser.entity_table.columnWidth(0) == 170
        model = browser.entity_table.model()
        assert model.data(model.index(1, 1), QtCore.Qt.ItemDataRole.ToolTipRole) == str(np.nextafter(1., 2.))
        assert model.headerData(0, QtCore.Qt.Orientation.Horizontal, QtCore.Qt.ItemDataRole.ToolTipRole) == 'gene_id'
    finally:
        browser.close()


def test_clearing_ambiguous_source_selection_clears_old_evidence(app):
    space, _, _ = fixture()
    browser = DatasetBrowser(space, organism=O.TOXOPLASMA)
    try:
        select(browser, 'numeric_fixture')
        assert browser.card.view is not None and not browser.entity_frame.empty
        browser.source_table.setCurrentCell(-1, -1)
        browser._show_selected()
        assert browser.selected_row is None and browser.card.view is None
        assert browser.entity_frame.empty and browser.page_start == 0
        assert not browser.page_previous.isEnabled() and not browser.page_next.isEnabled()
    finally:
        browser.close()
