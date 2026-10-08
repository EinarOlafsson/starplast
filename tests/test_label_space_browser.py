"""Check native class keys, complete paging and separate recorded evaluation scopes."""
import json

import numpy as np
import pandas as pd
import pytest
from PyQt6 import QtCore, QtWidgets

from starplast import scorecard as SC
from starplast.label_space_browser import LabelSpaceBrowser
from starplast.scorecard_view import build_scorecard_view


@pytest.fixture(scope='session')
def app():
    application = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])
    yield application


class Model:
    organism = 'Tg'

    def __init__(self):
        self.labels = pd.DataFrame([{'target': target, 'title': title, 'family': family,
            'genes': 402, 'annotated_genes': 401, 'unannotated_genes': 1,
            'source_ids': sources, 'truth_interpretation': 'Source annotation; missing membership is unknown'}
            for target, title, family, sources in [('location', 'Location', 'Localization', 'unresolved'),
                ('function_flag', 'Function flag', 'Function', 'fixture_one; fixture_two')]])
        self.cards = [build_scorecard_view({'task': SC.T_LABEL,
            'scope': {'organism': 'Tg', 'target': 'function_flag', 'strategy': 'fixture_strategy',
                      'settings': str(seed), 'seed': seed, 'partition': 'fixture_holdout'},
            'counts': {'eligible': 2, 'answered': 2, 'correct': 1},
            'accuracy': value, 'precision_of_calls': value, 'coverage': 1.},
            status='recorded', details={'rows': {'seed': seed}})
            for seed, value in ((1, np.nextafter(.5, 1.)), (2, .75))]

    def classes(self, target):
        return pd.DataFrame([{'value': value, 'native_value': native, 'description': 'Recorded class',
            'annotated_genes': size, 'class_fraction': size / 402, 'unknown_genes': 1}
            for value, native, size in [('False', False, 1), ('True', True, 400)]])

    def members(self, target, value=None):
        frame = pd.DataFrame({'organism': ['Tg'] * 401, 'gene_id': [f'fixture_gene_{i}' for i in range(401)],
            'target': [target] * 401, 'value': ['False'] + ['True'] * 400,
            'description': ['Recorded class'] * 401})
        return frame if value is None else frame[frame.value.eq(value)].reset_index(drop=True)

    def evidence_card(self, target, value=None):
        return build_scorecard_view({'scope': {'organism': 'Tg', 'target': target, 'unit': 'gene'},
            'counts': {'table_rows': 402, 'stored_any_rows': 401, 'missing_unknown_rows': 1},
            'evidence': {'native_value': False if value == 'False' else value,
                         'hierarchy': 'Native classes only; hierarchy not supplied'}}, kind='evidence')

    def mechanisms(self, target, value=None):
        return [{'strategy': f'fixture_strategy_{i}', 'title': f'Fixture strategy {i}',
            'applicable': True, 'status': 'recorded' if i == 0 else 'unavailable',
            'metrics': None, 'card': None, 'gaps': [] if i == 0 else ['No recorded evaluation'],
            'evaluations': [{'scope': {'organism': 'Tg', 'target': target, 'value': value,
                'settings': str(seed), 'seed': seed, 'mode': 'held_out', 'partition': 'fixture_holdout'},
                'card': card, 'class_metrics': {'precision': np.nextafter(.5, 1.)},
                'confusion': [{'truth': 'False', 'prediction': 'True', 'count': 1}],
                'gaps': []} for seed, card in enumerate(self.cards, 1)] if i == 0 else []}
            for i in range(39)]


def test_function_first_filters_and_explicit_selection_keep_native_false(app):
    browser = LabelSpaceBrowser(Model())
    try:
        assert browser.target_choice.currentData() == 'function_flag'
        browser.search.setText('location')
        assert browser.target_choice.count() == 1 and browser.target == 'location'
        browser.select('function_flag', 'False')
        assert not browser.search.text() and browser.value == 'False'
        assert browser.class_frame.iloc[0].native_value == False
        assert browser.member_frame.gene_id.tolist() == ['fixture_gene_0']
        snapshot = json.loads(browser.evidence_card.export_json())['snapshot']['card']
        assert snapshot['evidence']['native_value'] is False
        assert dict(browser.evidence_card.view.counts)['missing_unknown_rows'] == 1
        assert 'unknown' in browser.summary.text() and 'may overlap' in browser.summary.text()
        with pytest.raises(ValueError, match='Native class'):
            browser.select('function_flag', 'new guessed class')
        browser.search.setText('nothing matches')
        assert browser.target is None and browser.evidence_card.view is None
        assert browser.member_frame.empty and browser.strategy_card.view is None
    finally:
        browser.close()


def test_all_members_are_reachable_and_organism_qualified_routes_are_scoped(app):
    browser = LabelSpaceBrowser(Model())
    try:
        requested, sources = [], []
        browser.gene_requested.connect(lambda *values: requested.append(values))
        browser.source_requested.connect(sources.append)
        seen = []
        while True:
            assert len(browser.page_frame) <= 200
            seen.extend(browser.page_frame.gene_id)
            if not browser.page_next.isEnabled():
                break
            browser.page_next.click()
        assert seen == browser.member_frame.gene_id.tolist()
        assert browser.page_note.text() == '401–401 of 401 membership records'
        browser._open_gene(0, 0)
        assert requested == [('Tg', 'fixture_gene_400')]
        browser.page_frame = browser.page_frame.copy()
        browser.page_frame['organism'] = 'Hs'
        browser._open_gene(0, 0)
        assert len(requested) == 1
        browser.source_choice.setCurrentIndex(1)
        browser.source_open.click()
        assert sources == ['fixture_two']
        browser.select('location')
        assert not browser.source_open.isEnabled()
        assert browser.member_table.columnWidth(browser.page_frame.columns.get_loc('gene_id')) == 170
    finally:
        browser.close()


def test_all_strategies_and_separate_exact_cohort_exports_remain_visible(app):
    model = Model()
    browser = LabelSpaceBrowser(model)
    try:
        assert browser.strategy_choice.count() == 39
        assert browser.evaluation_choice.count() == 2
        assert browser.strategy_card.view == model.cards[0]
        details = json.loads(browser.evaluation_details.toPlainText())
        assert details['Class precision, recall and F1']['precision'] == np.nextafter(.5, 1.)
        assert details['Confusion records'] == [{'truth': 'False', 'prediction': 'True', 'count': 1}]
        assert 'seed' in browser.evaluation_choice.itemData(0, QtCore.Qt.ItemDataRole.ToolTipRole)
        browser.evaluation_choice.setCurrentIndex(1)
        assert browser.strategy_card.view == model.cards[1]
        snapshot = json.loads(browser.strategy_card.export_json())['snapshot']['card']
        assert snapshot['scope']['seed'] == 2 and snapshot['accuracy'] == .75
        browser.evaluation_choice.setCurrentIndex(0)
        assert json.loads(browser.strategy_card.export_json())['snapshot']['card']['accuracy'] == np.nextafter(.5, 1.)
        browser.strategy_choice.setCurrentIndex(38)
        assert browser.strategy_card.view is None
        assert 'metrics unavailable' in browser.strategy_note.text()
        assert 'Class metrics and confusions unavailable' in browser.evaluation_details.toPlainText()
        assert not browser.evaluation_choice.isEnabled()
        for control in (browser.family, browser.target_choice, browser.class_choice, browser.source_choice,
                browser.source_open, browser.strategy_choice, browser.evaluation_choice,
                browser.page_previous, browser.page_next):
            assert control.toolTip()
    finally:
        browser.close()


def test_real_backend_keeps_native_overlaps_false_unknowns_and_no_frozen_accuracy(app):
    from starplast.label_space import LabelSpace

    nodes = pd.DataFrame({'gene_id': [f'TGME49_{100000+i}' for i in range(4)],
        'ec_number': ['1.1.1.1; 2.2.2.2', '1.1.1.1', None, '2.2.2.2'],
        'is_member': pd.array([False, True, None, True], dtype='boolean')})
    browser = LabelSpaceBrowser(LabelSpace(nodes, 'Tg'))
    try:
        browser.select('ec_number')
        assert len(browser.member_frame) == 4
        assert browser.member_frame.gene_id.nunique() == 3
        assert set(browser.class_frame.value) == {'1.1.1.1', '2.2.2.2'}
        assert '1 unknown' in browser.summary.text()
        evidence = json.loads(browser.evidence_card.export_json())['snapshot']['card']['evidence']
        assert evidence['hierarchy']['status'] == 'unavailable'
        browser.select('is_member', 'False')
        assert browser.member_frame.gene_id.tolist() == ['TGME49_100000']
        assert browser.class_frame[browser.class_frame.value.eq('False')].iloc[0].native_value == False
        assert browser.strategy_choice.count() == 39
        assert browser.strategy_card.view is None and 'metrics unavailable' in browser.strategy_note.text()
    finally:
        browser.close()
