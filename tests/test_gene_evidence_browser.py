"""Exercise condensed evidence, exact values and scoped navigation with tiny records."""
from copy import deepcopy
import json

import numpy as np
import pandas as pd
import pytest
from PyQt6 import QtWidgets

from starplast import gene_evidence_browser as G
from starplast.scorecard_view import build_scorecard_view


@pytest.fixture(scope='session')
def app():
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield application


class Space:
    rows = []

    def source_card(self, row):
        return build_scorecard_view({'source': {'name': row['source_id']},
            'scope': {'organism': row['organism'], 'unit': row['unit'], 'target': row['question']},
            'counts': {'table_rows': 4, 'stored_any_rows': 3}, 'evidence': {'inventory': row}},
            kind='evidence')


def fixture(monkeypatch):
    nodes = pd.DataFrame({'gene_id': ['a', 'b', 'c', 'd'], 'product': ['word'] * 4,
        'label': ['one', 'one', 'two', None], 'quantity': [np.nextafter(1., 2.), 0., None, 3.],
        'flag': [False, True, None, False], 'missing': [None] * 4,
        **{f'measurement_{i}': [float(i)] * 4 for i in range(15)}})
    sources = [{'source_id': key, 'organism': 'Tg', 'unit': 'gene', 'question': 'recorded question'}
               for key in ('source_one', 'source_two')]
    rows = [{'column': column, 'value': deepcopy(nodes.iloc[0][column]), 'question': 'recorded question',
             'contexts': ['recorded context'], 'quantity_unit': 'supplied unit',
             'origins': ['measured'], 'gaps': [],
             'sources': deepcopy(sources) if column in {'label', 'quantity'} else []}
            for column in nodes if column != 'gene_id']
    monkeypatch.setattr(G, 'evidence_rows', lambda *args: deepcopy(rows))
    return G.GeneEvidenceBrowser(nodes, 'Tg', 'a', Space()), nodes, rows


def select(browser, column):
    index = next(i for i, row in enumerate(browser.shown_rows) if row['column'] == column)
    browser.table.setCurrentCell(index, 0)
    browser._selected()
    return index


def test_condensed_and_expanded_views_retain_every_value_and_missingness(app, monkeypatch):
    browser, nodes, rows = fixture(monkeypatch)
    try:
        assert len(browser.shown_rows) == 12
        assert not any(row['column'] == 'missing' for row in browser.shown_rows)
        browser.show_all.setChecked(True)
        assert {row['column'] for row in browser.shown_rows} == set(nodes) - {'gene_id'}
        assert len(browser.shown_rows) == len(rows)
        index = select(browser, 'quantity')
        assert browser.table.item(index, 2).text() == str(np.nextafter(1., 2.))
        assert browser.table.item(index, 2).toolTip() == repr(browser.selected_row['value'])
        index = select(browser, 'flag')
        assert browser.table.item(index, 2).text() == 'False'
        index = select(browser, 'measurement_0')
        assert browser.table.item(index, 2).text() == '0.0'
        index = select(browser, 'missing')
        assert browser.table.item(index, 2).text() == 'Unknown / missing'
        assert browser.selected_row['value'] is None
        assert browser.table.item(index, 2).toolTip() == 'None'
        assert browser.table.item(index, 3).text() == "['recorded context']"
    finally:
        browser.close()


def test_sources_are_separately_scoped_and_exact_exports_survive_navigation(app, monkeypatch):
    browser, _, _ = fixture(monkeypatch)
    try:
        browser.show_all.setChecked(True)
        select(browser, 'quantity')
        assert browser.source_choice.count() == 2
        assert browser.card.view.kind == 'evidence' and not browser.card.view.metrics
        requested = []
        browser.source_requested.connect(requested.append)
        browser.source_choice.setCurrentIndex(1)
        browser.source_open.click()
        assert requested == [{'source_id': 'source_two', 'organism': 'Tg', 'unit': 'gene', 'question': 'recorded question'}]
        snapshot = json.loads(browser.card.export_json())['snapshot']['card']
        assert snapshot['evidence']['inventory'] == requested[0]
        assert snapshot['counts']['table_rows'] == 4
        requested[0]['organism'] = 'other'
        assert browser.source_choice.currentData()['organism'] == 'Tg'
        select(browser, 'product')
        assert browser.card.view is None and not browser.source_open.isEnabled()
        assert not browser.source_choice.isEnabled()
    finally:
        browser.close()


def test_label_and_class_routes_exclude_quantities_free_text_and_missing_classes(app, monkeypatch):
    browser, nodes, rows = fixture(monkeypatch)
    try:
        browser.show_all.setChecked(True)
        labels, classes = [], []
        browser.label_requested.connect(labels.append)
        browser.class_requested.connect(lambda *values: classes.append(values))
        select(browser, 'label')
        assert browser.label_open.isEnabled() and browser.class_open.isEnabled()
        browser.label_open.click()
        browser.class_open.click()
        assert labels == ['label'] and classes == [('label', 'one')]
        for column in ('product', 'quantity', 'missing', 'flag'):
            select(browser, column)
            assert not browser.label_open.isEnabled() and not browser.class_open.isEnabled()
        for row in rows:
            if row['column'] == 'label':
                row['value'] = None
        monkeypatch.setattr(G, 'evidence_rows', lambda *args: deepcopy(rows))
        missing = G.GeneEvidenceBrowser(nodes, 'Tg', 'd', Space())
        missing.show_all.setChecked(True)
        select(missing, 'label')
        assert missing.label_open.isEnabled() and not missing.class_open.isEnabled()
        missing.close()
        assert browser.table.columnWidth(0) == 240
        assert browser.table.columnWidth(1) == 220
        assert browser.table.horizontalHeaderItem(2).toolTip() == 'Original value'
        assert all(widget.toolTip() for widget in (browser.show_all, browser.source_choice,
            browser.source_open, browser.label_open, browser.class_open))
    finally:
        browser.close()


def test_real_backend_preserves_columns_and_unavailable_target_gap(app):
    nodes = pd.DataFrame({'gene_id': ['TGME49_100000', 'TGME49_100001'],
                          'product': ['fixture only', None], 'quantity': [0., None],
                          'flag': [False, None], 'labels': [['one', 'two'], []]})
    browser = G.GeneEvidenceBrowser(nodes, 'Tg', 'TGME49_100000', Space())
    unavailable = G.GeneEvidenceBrowser(nodes, 'Tg', 'TGME49_999999', Space())
    try:
        browser.show_all.setChecked(True)
        assert {row['column'] for row in browser.rows} == set(nodes) - {'gene_id'}
        select(browser, 'labels')
        assert browser.selected_row['value'] == ['one', 'two']
        assert 'Source attribution unavailable' in str(browser.selected_row['gaps'])
        assert browser.selected_row['question'] == 'Unassigned'
        assert not unavailable.rows and 'unavailable' in unavailable.summary.text()
        assert unavailable.card.view is None
        assert not any(button.isEnabled() for button in (unavailable.source_open, unavailable.label_open, unavailable.class_open))
        with pytest.raises(ValueError, match='registered organism'):
            G.GeneEvidenceBrowser(nodes, 'unregistered organism', 'TGME49_999999', Space())
    finally:
        browser.close()
        unavailable.close()
