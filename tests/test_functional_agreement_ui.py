"""Paired desktop cards retain both calls and refuse stale source/report scopes."""
import json
from pathlib import Path

import pandas as pd
import pytest
from PyQt6 import sip

from starplast import claims as C, functional_agreement as A, functional_results as F
from starplast import strategies as S, track_record as T
from starplast.discoveries_panel import DiscoveriesPanel


@pytest.fixture
def panel(monkeypatch):
    monkeypatch.setattr(C, 'shipped', lambda *args: pd.DataFrame())
    monkeypatch.setattr(C, 'recipes', lambda *args: pd.DataFrame())
    monkeypatch.setattr(T, 'shipped', lambda *args: pd.DataFrame())
    panel = DiscoveriesPanel('Pf')
    yield panel
    dialog = getattr(panel, 'functional_outcome_dialog', None)
    if dialog is not None and not sip.isdeleted(dialog): dialog.close()
    panel.close()


def select(panel, strategy, address):
    index = next(i for i in range(panel.functional_choice.count())
                 if panel.functional_choice.itemData(i).metadata['strategy'] == strategy)
    panel.functional_choice.setCurrentIndex(index)
    index = panel.functional_view.findData(address)
    assert index >= 0
    panel.functional_view.setCurrentIndex(index)


def test_pair_from_both_strategies_retains_exact_export_rows_and_gene_outcomes(panel):
    source = Path(A.__file__).with_name('data') / 'functional_agreement.json'
    report = json.loads(source.read_text())
    for strategy in ('feature_knn', 'random_forest'):
        select(panel, strategy, 'comparison:paired')
        assert panel.functional_export.isEnabled()
        assert panel.functional_rows.rowCount() == 152
        assert panel.functional_shown_rows == report['rows']
        assert json.loads(panel.functional_card.export_json())['snapshot'] == report
        assert '33 / 94 (35.1%)' in panel.functional_card.toPlainText()
        assert panel.functional_rows.item(0, 5).text() == 'Same correct profile'
        for index, row in enumerate(report['rows']):
            assert panel.functional_rows.item(index, 3).text() == (row['left_prediction'] or 'abstained')
            assert panel.functional_rows.item(index, 4).text() == (row['right_prediction'] or 'abstained')
        chosen = []
        panel.gene_chosen.connect(chosen.append)
        panel._functional_gene_clicked(0, 0)
        assert chosen[-1] == report['rows'][0]['entity']
        panel._functional_outcome(0, 0)
        browser = panel.functional_outcome_dialog.findChild(type(panel.functional_card))
        assert browser.view.kind == 'outcome' and not browser.view.metrics
        assert json.loads(browser.export_json())['snapshot']['outcome'] == report['rows'][0]
        panel.functional_outcome_dialog.close()
        select(panel, strategy, 'strategy:')
        assert panel.functional_comparison_report is None
        assert panel.functional_rows.horizontalHeaderItem(3).text() == 'recovered profile'
        assert panel.functional_rows.horizontalHeaderItem(5).text() == 'method support (uncalibrated)'
        assert panel.functional_card.view.kind == 'performance'


def test_unavailable_report_clears_comparison_and_original_method_still_opens(panel, monkeypatch):
    select(panel, 'feature_knn', 'comparison:paired')
    def unavailable(*args, **kwargs): raise ValueError('Paired report checksum changed')
    monkeypatch.setattr(A, 'load_report', unavailable)
    panel._functional_view_selected()
    assert panel.functional_card.view is None and not panel.functional_export.isEnabled()
    assert not panel.functional_shown_rows and panel.functional_rows.rowCount() == 0
    assert panel.functional_comparison_report is None
    assert 'checksum changed' in panel.functional_card.toPlainText()
    select(panel, 'feature_knn', 'strategy:')
    assert panel.functional_export.isEnabled() and panel.functional_rows.rowCount() == 152


def test_changed_live_context_clears_a_previously_displayed_comparison(panel):
    select(panel, 'feature_knn', 'comparison:paired')
    nodes = panel.context.nodes.copy()
    nodes.loc[nodes.index[0], 'gene_id'] = 'altered_source'
    panel.context = S.Context(nodes, graph={}, organism='Pf')
    panel._functional_view_selected()
    assert panel.functional_card.view is None and not panel.functional_export.isEnabled()
    assert panel.functional_rows.rowCount() == 0 and not panel.functional_shown_rows
    assert 'values or entities differ' in panel.functional_card.toPlainText()
    select(panel, 'feature_knn', 'strategy:')
    assert panel.functional_card.view is None and not panel.functional_export.isEnabled()
    assert panel.functional_rows.rowCount() == 0


def test_report_checksum_and_cross_organism_are_refused(panel, tmp_path):
    pair = A.find_pair(panel.functional_choice.currentData(), panel.functional_benchmarks)
    path = tmp_path / 'changed.json'; path.write_text('{}')
    with pytest.raises(ValueError, match='checksum'):
        A.load_report(*pair, panel.context, path=path)
    wrong = S.Context(panel.context.nodes, graph={}, organism='Tg')
    with pytest.raises(ValueError, match='organism'):
        A.load_report(*pair, wrong)
    with pytest.raises(ValueError, match='outside'):
        A.build_gene_scorecard(A.load_report(*pair, panel.context), 'not_in_test')
    assert A.find_pair(pair[0], [pair[0]]) is None
    assert A.find_pair(pair[0], [*pair, pair[0]]) is None
