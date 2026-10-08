"""Synthetic dataset browsing checks protect source scopes and original values."""
from dataclasses import FrozenInstanceError, replace
import json

import numpy as np
import pandas as pd
import pytest

from starplast import datasets as D, organisms as O, slots as S
from starplast.dataset_space import DatasetSpace
from starplast.inventory import build_inventory
from starplast.provenance import MeasurementTrace


class Graph(dict):
    @property
    def files(self):
        return list(self)


def source(key='example', columns=('value',), **kwargs):
    return D.Dataset(key, 'Synthetic evidence', 'DNA', 'synthetic', 'Stored values',
                     organism=O.TOXOPLASMA, columns=columns, **kwargs)


def table():
    return pd.DataFrame({'gene_id': ['TGME49_100001', 'TGME49_100002', 'TGME49_100003'],
                         'value': [0., 2., np.nan], 'detected': pd.Series([False, True, pd.NA], dtype='boolean')})


def space(data=None, sources=None, **kwargs):
    return DatasetSpace(None, tables={(O.TOXOPLASMA, 'gene'): table() if data is None else data},
                        sources=[source()] if sources is None else sources, catalog=[], **kwargs)


def test_original_values_counts_pagination_and_nulls_are_preserved():
    browser = space(sources=[source(columns=('value', 'detected'))])
    row = browser.rows[0]
    assert row['table_rows'] == 3 and row['stored_any_rows'] == 2
    assert row['false_cells'] == 1 and row['missing_unknown_rows'] == 1
    result = browser.entities(row, start=0, limit=1)
    assert result.value.iloc[0] == 0 and not result.detected.iloc[0]
    assert result.attrs['entity_column'] == 'gene_id'
    assert result.attrs['total_rows'] == 3 and result.attrs['displayed_rows'] == 1
    last = browser.entities(row, start=2, limit=1)
    assert pd.isna(last.value.iloc[0]) and pd.isna(last.detected.iloc[0])
    assert browser.entities(row, start=3, limit=1).empty
    assert browser.entities(row).gene_id.tolist() == table().gene_id.tolist()


def test_rows_table_and_cards_are_copied_and_evidence_has_no_accuracy():
    supplied = table()
    browser = space(supplied)
    row = browser.rows[0]
    row['columns'].append('invented')
    supplied.loc[0, 'value'] = 999
    assert browser.entities(browser.rows[0]).value.iloc[0] == 0
    card = browser.source_card(row)
    assert card.kind == 'evidence' and card.task is None and card.metrics == ()
    assert dict(card.counts)['stored_any_rows'] == 2
    assert card.source.grade == 'unresolved'
    with pytest.raises(FrozenInstanceError):
        card.status = 'changed'
    assert 'invented' not in card.snapshot_json


def test_object_cells_retain_explicit_nulls_and_nested_values_without_aliases():
    supplied = table()
    supplied['value'] = pd.Series([['original'], None, False], dtype=object)
    browser = space(supplied)
    supplied.value.iloc[0].append('changed')
    result = browser.entities(browser.rows[0])
    assert result.value.iloc[0] == ['original'] and result.value.iloc[1] is None
    assert result.value.iloc[2] is False and result.value.dtype == object
    result.value.iloc[0].append('also changed')
    assert browser.entities(browser.rows[0]).value.iloc[0] == ['original']


def test_filters_and_facets_reproduce_scoped_inventory():
    slot = S.Slot(O.TOXOPLASMA, 'assay', 'DNA', 'culture', 'gene', ('value',), 'one',
                  context_path=('culture', 'stage'), evidence_path=('measurements', 'synthetic'))
    sources = [source(), source('other', ('absent',))]
    inventory = build_inventory({(O.TOXOPLASMA, 'gene'): table()}, sources=sources, catalog=[slot])
    browser = DatasetSpace(inventory, tables={(O.TOXOPLASMA, 'gene'): table()}, sources=sources, catalog=[slot])
    assert len(browser.rows) == len(inventory) == 2
    assert len(browser.filter_rows(organism=O.TOXOPLASMA, unit='gene', status='installed')) == 1
    assert len(browser.filter_rows(context='culture / stage', family='measurements / synthetic')) == 1
    assert browser.filter_rows(text='OTHER')[0]['source_id'] == 'other'
    assert browser.filter_rows(organism=O.FALCIPARUM) == []
    assert browser.facets()['status'] == ['installed', 'unavailable']


def test_question_refusal_never_borrows_same_source_installed_values():
    refusal = {'source_id': 'example', 'organism': O.TOXOPLASMA, 'unit': 'gene',
               'question': 'wrong quantity', 'reason': 'Does not measure this question'}
    inventory = build_inventory({(O.TOXOPLASMA, 'gene'): table()}, sources=[source()], catalog=[], refusals=[refusal])
    browser = DatasetSpace(inventory, tables={(O.TOXOPLASMA, 'gene'): table()}, sources=[source()], catalog=[])
    installed, rejected = browser.rows
    assert len(browser.entities(installed)) == 3
    result = browser.entities(rejected)
    assert result.empty and result.attrs['gap'] == refusal['reason']
    assert browser.source_card(rejected).status == 'rejected'
    with pytest.raises(ValueError, match='outside'):
        browser.entities(dict(installed, organism=O.FALCIPARUM))


def test_unavailable_partial_and_missing_identifiers_remain_explicit():
    browser = space(sources=[source(columns=('value', 'absent'))])
    assert browser.rows[0]['status'] == 'partial'
    assert browser.entities(browser.rows[0]).columns.tolist() == ['gene_id', 'value']
    assert 'absent' in browser.entities(browser.rows[0]).attrs['gap']
    browser = space(pd.DataFrame({'value': [0.]}))
    unavailable = browser.entities(browser.rows[0])
    assert unavailable.empty and 'identifiers unavailable' in unavailable.attrs['gap']
    browser = DatasetSpace(None, sources=[source(columns=())], catalog=[])
    assert browser.rows[0]['unit'] == 'source_records'
    assert browser.entities(browser.rows[0]).attrs['status'] == 'unavailable'


def test_host_proteins_and_metabolites_do_not_project_to_genes():
    host_source = source(path='starplast/data/' + O.HOST_TABLES[O.HUMAN])
    proteins = pd.DataFrame({'host_id': ['P12345'], 'value': [0.]})
    browser = DatasetSpace(None, tables={(O.HUMAN, 'protein'): proteins}, sources=[host_source], catalog=[])
    assert browser.rows[0]['organism'] == O.HUMAN
    assert browser.entities(browser.rows[0]).attrs['entity_column'] == 'host_id'
    metabolite = S.Slot(O.TOXOPLASMA, 'metabolites', 'DNA', 'culture', 'metabolite', ('value',), 'one')
    browser = DatasetSpace(None, tables={(O.TOXOPLASMA, 'metabolite'): pd.DataFrame({'metabolite': ['M1'], 'value': [0.]})},
                           sources=[source()], catalog=[metabolite])
    assert browser.entities(browser.rows[0]).attrs['entity_column'] == 'metabolite'


@pytest.mark.parametrize(('grade', 'letter'), [('direct_experiment', 'A'), ('orthology_transfer', 'B'),
                                              ('prediction', 'C'), ('derived_quantity', 'C'), ('curation', None)])
def test_supplied_grades_are_distinct_from_installed_status(grade, letter):
    key = ('example', O.TOXOPLASMA, 'gene', '')
    browser = space(origins={key: {'evidence_grade': grade}})
    assert browser.rows[0]['slot_grades'] == ([letter] if letter else [])
    assert browser.source_card(browser.rows[0]).source.grade == grade


def test_import_origin_is_scoped_and_table_origin_cannot_grade_all_sources():
    native, imported = source(), source('session:1', ('detected',))
    origin = {('session:1', O.TOXOPLASMA, 'gene', ''): {'kind': 'user_imported', 'authorship': 'User supplied'}}
    browser = space(sources=[native, imported], origins=origin)
    assert browser.rows[0]['origin'] == 'registry_asserted_unverified'
    assert browser.rows[1]['origin'] == 'user_imported'
    assert browser.rows[1]['evidence_grade'] == 'unresolved'
    browser = space(origins={(O.TOXOPLASMA, 'gene'): {'sha256': 'a' * 64, 'evidence_grade': 'direct_experiment'}})
    assert browser.source_card(browser.rows[0]).source.sha256 == 'a' * 64
    assert browser.rows[0]['slot_grades'] == []


def test_typed_trace_retains_units_mappings_and_source_lineage_without_admission():
    trace = MeasurementTrace('example', O.TOXOPLASMA, O.TOXOPLASMA, 'gene', 'value',
                             evidence_grade='prediction', quantity_unit='arbitrary units', gaps=('mapping_unresolved',))
    browser = space(origins={('example', O.TOXOPLASMA, 'gene', ''): {'traces': (trace,)}})
    card = browser.source_card(browser.rows[0])
    saved = json.loads(card.snapshot_json)
    assert saved['card']['evidence']['quantity_units'] == {'value': 'arbitrary units'}
    assert 'mapping_unresolved' in saved['details']['gaps']
    assert card.source.grade == 'prediction'
    with pytest.raises(ValueError, match='another'):
        space(origins={('example', O.TOXOPLASMA, 'gene', ''): {'traces': (replace(trace, source_id='other'),)}})


def test_registry_declared_derivation_is_unverified_and_https_links_are_typed():
    browser = space(sources=[source(url='https://example.org/source', derived_from=('input',))])
    row = browser.rows[0]
    assert row['declared_derivation'] == ['input'] and row['slot_grades'] == []
    assert browser.source_card(row).links[0].url == 'https://example.org/source'
    browser = space(sources=[source(url='javascript:bad')])
    assert browser.source_card(browser.rows[0]).links == ()
    assert 'not an allowed HTTPS' in browser.source_card(browser.rows[0]).snapshot_json


def test_pair_paging_preserves_original_endpoints_weights_and_full_counts():
    graph = Graph(layer__a=np.array([0, 1, 2]), layer__b=np.array([1, 2, 0]), layer__w=np.array([0., .5, 1.]))
    browser = DatasetSpace(None, graphs={O.TOXOPLASMA: graph}, sources=[source(columns=('edge:layer',))], catalog=[])
    row = browser.rows[0]
    frame = browser.entities(row, start=1, limit=1)
    assert frame.layer__a.tolist() == [1] and frame.layer__w.tolist() == [.5]
    assert frame.attrs['entity_column'] is None and frame.attrs['total_rows'] == 3
    assert frame.attrs['endpoint_mapping'] == 'unresolved'
    assert browser.entities(row, start=3, limit=1).columns.tolist() == frame.columns.tolist()
    assert row['pair_records'] == 3 and row['assayed_pair_denominator'] is None


def test_pair_endpoint_names_require_matching_explicit_gene_order():
    graph = Graph(layer__a=np.array([0]), layer__b=np.array([1]), gene_ids=table().gene_id.to_numpy())
    browser = space(sources=[source(columns=('edge:layer',))], graphs={O.TOXOPLASMA: graph})
    frame = browser.entities(browser.rows[0], limit=1)
    assert frame.endpoint_a_gene_id.iloc[0] == 'TGME49_100001'
    assert frame.attrs['endpoint_mapping'] == 'verified_gene_order'
    graph['gene_ids'] = graph['gene_ids'][::-1]
    assert 'endpoint_a_gene_id' not in browser.entities(browser.rows[0])


def test_pair_page_crosses_layers_without_truncating_full_population():
    graph = Graph(first__a=np.array([0, 1]), first__b=np.array([1, 0]),
                  second__a=np.array([2, 3]), second__b=np.array([3, 2]))
    browser = DatasetSpace(None, graphs={O.TOXOPLASMA: graph},
                           sources=[source(columns=('edge:first', 'edge:second'))], catalog=[])
    result = browser.entities(browser.rows[0], start=1, limit=2)
    assert result.source_column.tolist() == ['edge:first', 'edge:second']
    assert result.attrs['total_rows'] == 4 and len(result) == 2


def test_empty_graph_is_installed_zero_records_and_bad_weights_are_refused():
    graph = Graph(layer__a=np.array([], dtype=int), layer__b=np.array([], dtype=int))
    browser = DatasetSpace(None, graphs={O.TOXOPLASMA: graph}, sources=[source(columns=('edge:layer',))], catalog=[])
    frame = browser.entities(browser.rows[0], limit=1)
    assert frame.empty and frame.attrs['total_rows'] == 0 and frame.attrs['status'] == 'installed'
    graph = Graph(layer__a=np.array([0, 1]), layer__b=np.array([1, 0]), layer__w=np.array([1.]))
    browser = DatasetSpace(None, graphs={O.TOXOPLASMA: graph}, sources=[source(columns=('edge:layer',))], catalog=[])
    with pytest.raises(ValueError, match='complete stored'):
        browser.entities(browser.rows[0], limit=1)


def test_bridges_require_source_attribution_and_page_only_selected_paper():
    bridges = pd.DataFrame({'bridge': ['host'] * 4, 'source': ['example', 'other', 'example', 'example'],
                            'host_id': ['P1', 'P2', 'P3', 'P4'], 'value': [0., 5., 1., 2.]})
    browser = DatasetSpace(None, bridges={O.TOXOPLASMA: bridges}, sources=[source(columns=('bridge:host',))], catalog=[])
    frame = browser.entities(browser.rows[0], start=1, limit=1)
    assert frame.host_id.tolist() == ['P3'] and frame.attrs['total_rows'] == 3
    browser = DatasetSpace(None, bridges={O.TOXOPLASMA: bridges.drop(columns='source')},
                           sources=[source(columns=('bridge:host',))], catalog=[])
    assert browser.rows[0]['status'] == 'unattributed'
    assert browser.entities(browser.rows[0]).attrs['status'] == 'unavailable'


@pytest.mark.parametrize(('start', 'limit'), [(-1, None), (True, 1), (0, 0), (0, -1), (0, True)])
def test_invalid_page_parameters_are_refused(start, limit):
    browser = space()
    with pytest.raises(ValueError, match='Page'):
        browser.entities(browser.rows[0], start=start, limit=limit)


def test_invalid_table_addresses_scopes_and_grades_are_refused():
    with pytest.raises(ValueError, match='explicit'):
        DatasetSpace(None, tables={'Tg': table()}, sources=[])
    with pytest.raises(ValueError, match='explicit organism'):
        space(table().assign(gene_id=['PF3D7_0100001', 'PF3D7_0100002', 'PF3D7_0100003']))
    with pytest.raises(ValueError, match='grade'):
        space(origins={('example', O.TOXOPLASMA, 'gene', ''): {'evidence_grade': 'A'}})
    report = build_inventory(sources=[source()], catalog=[])
    with pytest.raises(ValueError, match='duplicate'):
        DatasetSpace(pd.concat([report, report]), sources=[source()], catalog=[])
    assert DatasetSpace(None, sources=[], catalog=[]).rows == []
