"""Synthetic gene evidence contracts preserve scoped identities and original cells."""
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

from starplast import datasets as D, organisms as O, slots as S
from starplast.dataset_space import DatasetSpace
from starplast.gene_evidence import build_resolver, evidence_rows
from starplast.identity import GeneIndex, norm
from starplast.provenance import MeasurementTrace

G1, G2, G3 = 'TGME49_100001', 'TGME49_100002', 'TGME49_100003'


def nodes():
    return pd.DataFrame({'gene_id': [G1, G2], 'symbol': ['SHARED', 'SHARED'],
                         'gene_name': ['Short', 'Other'], 'value': [0., 3.],
                         'flag': pd.Series([False, True], dtype=object),
                         'nested': pd.Series([['PF00001', 'PF00002'], []], dtype=object),
                         'missing': pd.Series([None, None], dtype=object)})


def source(key='example', columns=('value', 'flag')):
    return D.Dataset(key, 'Synthetic source', 'DNA', 'synthetic', 'Original stored values',
                     organism=O.TOXOPLASMA, columns=columns)


def space(table, sources=None, origins=None):
    return DatasetSpace(None, tables={(O.TOXOPLASMA, 'gene'): table},
                        sources=[source()] if sources is None else sources, catalog=[], origins=origins)


def test_exact_canonical_and_explicit_aliases_preserve_ambiguities_sources():
    resolver = build_resolver(nodes(), O.TOXOPLASMA)
    assert resolver.resolve(G1).entity.identifier == G1
    assert resolver.resolve('short').entity.identifier == G1
    assert resolver.resolve('SHARED').status == 'ambiguous'
    assert {choice.entity.identifier for choice in resolver.resolve('SHARED').choices} == {G1, G2}
    assert resolver.resolve('SHARED').entity is None
    assert resolver.resolve('short-extra').status == 'unresolved'
    assert resolver.resolve('TgShort').status == 'unresolved'
    assert resolver.resolve('TGGT1_100001').status == 'unresolved'
    assert all(choice.source.startswith('table_alias_mapping_sha256:') for choice in resolver.resolve('short').choices)


def test_explicit_symbol_collections_are_not_delimiter_or_prose_guesses():
    table = nodes()
    table['symbol'] = pd.Series([['First', 'Second'], 'Whole; cell'], dtype=object)
    resolver = build_resolver(table, O.TOXOPLASMA)
    assert resolver.resolve('First').entity.identifier == G1
    assert resolver.resolve('Second').entity.identifier == G1
    assert resolver.resolve('Whole; cell').entity.identifier == G2
    assert resolver.resolve('Whole').status == 'unresolved'


def test_mapping_identity_changes_with_aliases_and_order_without_mutating_input():
    table = nodes()
    original = deepcopy(table)
    first = build_resolver(table, O.TOXOPLASMA).resolve('Short').choices[0].source
    table.loc[0, 'gene_name'] = 'Renamed'
    second = build_resolver(table, O.TOXOPLASMA).resolve('Renamed').choices[0].source
    assert first != second
    third = build_resolver(original.iloc[::-1], O.TOXOPLASMA).resolve('Short').choices[0].source
    assert first != third
    assert original.gene_name.iloc[0] == 'Short'


def test_same_alias_in_another_organism_never_changes_current_scope():
    table = pd.DataFrame({'gene_id': ['PF3D7_0100001'], 'symbol': ['Short']})
    assert build_resolver(table, O.FALCIPARUM).resolve('Short').entity.organism == O.FALCIPARUM
    assert build_resolver(nodes(), O.TOXOPLASMA).resolve('PF3D7_0100001').status == 'unresolved'
    with pytest.raises(ValueError):
        build_resolver(table, O.TOXOPLASMA)


def test_supplied_index_preserves_collisions_with_unavailable_same_organism_targets():
    index = GeneIndex(canonical={G1, G3}, lookup={norm('Historical'): (G3, 'accession_prev')},
                      ambiguous={norm('IndexCollision'): {G1, G3}})
    resolver = build_resolver(nodes(), O.TOXOPLASMA, index=index, index_source='index_sha256:verified')
    collision = resolver.resolve('IndexCollision')
    assert collision.status == 'ambiguous'
    assert {choice.entity.identifier for choice in collision.choices} == {G1, G3}
    assert all(choice.source == 'index_sha256:verified' for choice in collision.choices)
    assert resolver.resolve('Historical').entity.identifier == G3
    assert resolver.resolve(G3).entity.identifier == G3
    assert resolver.resolve(G2).entity.identifier == G2
    with pytest.raises(ValueError, match='exact canonical'):
        evidence_rows(nodes(), O.TOXOPLASMA, G3, space(nodes()))


def test_index_requires_source_and_canonical_scoped_references():
    with pytest.raises(ValueError, match='index_source'):
        build_resolver(nodes(), O.TOXOPLASMA, index=GeneIndex())
    foreign = GeneIndex(canonical={'PF3D7_0100001'})
    with pytest.raises(ValueError):
        build_resolver(nodes(), O.TOXOPLASMA, index=foreign, index_source='source')
    malformed = GeneIndex(canonical={G1}, lookup={'BAD': (G2, 'symbol')})
    with pytest.raises(ValueError, match='membership'):
        build_resolver(nodes(), O.TOXOPLASMA, index=malformed, index_source='source')


def test_every_column_original_zero_false_nested_and_missing_are_reachable(monkeypatch):
    monkeypatch.setattr(S, 'all_slots', lambda organism: ())
    table = nodes()
    rows = evidence_rows(table, O.TOXOPLASMA, G1, space(table))
    assert [row['column'] for row in rows] == list(table.columns[1:])
    by_column = {row['column']: row for row in rows}
    assert by_column['value']['value'] == 0
    assert by_column['flag']['value'] is False
    assert by_column['missing']['value'] is None
    assert by_column['nested']['value'] == ['PF00001', 'PF00002']
    by_column['nested']['value'].append('changed')
    assert table.nested.iloc[0] == ['PF00001', 'PF00002']
    assert all(row['question'] == 'Unassigned' for row in rows)
    assert by_column['missing']['quantity_unit'] == 'unresolved'
    assert 'Source attribution unavailable for stored column' in by_column['missing']['gaps']


def test_all_overlapping_questions_and_contexts_survive(monkeypatch):
    one = S.Slot(O.TOXOPLASMA, 'Question one', 'DNA', 'Culture', 'gene', ('value',), 'one')
    two = S.Slot(O.TOXOPLASMA, 'Question two', 'DNA', 'Stage', 'gene', ('value',), 'one')
    foreign = S.Slot(O.FALCIPARUM, 'Foreign question', 'DNA', 'Wrong', 'gene', ('value',), 'one')
    pair = S.Slot(O.TOXOPLASMA, 'Pair question', 'DNA', 'Wrong', 'pair', ('value',), 'one')
    monkeypatch.setattr(S, 'all_slots', lambda organism: (one, two, foreign, pair))
    row = next(row for row in evidence_rows(nodes(), O.TOXOPLASMA, G1, space(nodes())) if row['column'] == 'value')
    assert row['questions'] == ['Question one', 'Question two']
    assert row['question'] == 'Question one / Question two'
    assert row['contexts'] == ['Culture', 'Stage']


def test_original_scoped_sources_units_and_gaps_come_from_evidence_cards(monkeypatch):
    monkeypatch.setattr(S, 'all_slots', lambda organism: ())
    trace = MeasurementTrace('example', O.TOXOPLASMA, O.TOXOPLASMA, 'gene', 'value',
                             evidence_grade='prediction', quantity_unit='stored units', gaps=('mapping_unresolved',))
    browser = space(nodes(), origins={('example', O.TOXOPLASMA, 'gene', ''): {'traces': (trace,)}})
    row = next(row for row in evidence_rows(nodes(), O.TOXOPLASMA, G1, browser) if row['column'] == 'value')
    assert row['sources'] == browser.rows
    assert row['quantity_unit'] == 'stored units'
    assert 'mapping_unresolved' in row['gaps']
    row['sources'][0]['columns'].append('changed')
    assert 'changed' not in browser.rows[0]['columns']


def test_conflicting_source_units_never_choose_convert_or_invent_one(monkeypatch):
    monkeypatch.setattr(S, 'all_slots', lambda organism: ())
    sources = [source(), source('second')]
    origins = {('example', O.TOXOPLASMA, 'gene', ''): {'quantity_unit': 'unit A'},
               ('second', O.TOXOPLASMA, 'gene', ''): {'quantity_unit': 'unit B', 'kind': 'user_imported'}}
    row = next(row for row in evidence_rows(nodes(), O.TOXOPLASMA, G1, space(nodes(), sources, origins)) if row['column'] == 'value')
    assert row['quantity_units'] == ['unit A', 'unit B']
    assert row['quantity_unit'] == 'conflicting: unit A / unit B'
    assert row['origins'] == ['registry_asserted_unverified', 'user_imported']
    assert any('conflict' in gap for gap in row['gaps'])


@pytest.mark.parametrize('gene_id', ['SHARED', 'Short', G3, None])
def test_evidence_requires_one_existing_exact_canonical_gene(gene_id):
    with pytest.raises(ValueError, match='exact canonical'):
        evidence_rows(nodes(), O.TOXOPLASMA, gene_id, space(nodes()))


@pytest.mark.parametrize('change', ['missing', 'duplicate', 'foreign', 'column'])
def test_invalid_canonical_tables_are_refused(change):
    table = nodes()
    if change == 'missing':
        table = table.drop(columns='gene_id')
    elif change == 'duplicate':
        table.loc[1, 'gene_id'] = G1
    elif change == 'foreign':
        table.loc[1, 'gene_id'] = 'PF3D7_0100001'
    else:
        table = pd.concat([table, table[['symbol']]], axis=1)
    with pytest.raises(ValueError):
        build_resolver(table, O.TOXOPLASMA)
