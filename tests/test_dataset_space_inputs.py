"""Local browser assembly preserves explicit source addresses and import scope."""
import numpy as np
import pandas as pd
import pytest

from starplast import datasets as D, organisms as O, source_refusals as R, slots as S
from starplast.dataset_space_inputs import from_sources
from starplast.scorecard_view import export_scorecard


def test_import_records_remain_separate_and_processed_cache_is_not_raw(monkeypatch):
    source = D.Dataset('example', 'Example', 'gene', 'synthetic', 'Values',
        organism=O.TOXOPLASMA, columns=('value',), path='starplast/data/nodes.parquet')
    monkeypatch.setattr(D, 'registry', lambda: [source])
    monkeypatch.setattr(D, 'local_path', lambda key: 'located_processed_cache')
    monkeypatch.setattr(R, 'records', lambda: [])
    nodes = pd.DataFrame({'gene_id': ['TGME49_100001'], 'value': [0.], 'imported_value': [False]})
    model = from_sources({O.TOXOPLASMA: {'nodes': nodes}},
        imports=[{'columns': ['imported_value'], 'quantification': 'mean', 'genes': 1}],
        imported_organism=O.TOXOPLASMA)
    original, imported = model.rows
    assert original['origin'] == 'registry_asserted_unverified' and original['raw_available'] is False
    assert original['processed_available'] is True
    assert imported['origin'] == 'user_imported' and imported['stored_any_rows'] == 1
    assert not model.entities(imported).imported_value.iloc[0]
    assert model.source_card(imported).metrics == ()
    assert 'original source-file identity' in model.source_card(imported).source.lineage.lower()
    pd.testing.assert_frame_equal(nodes, pd.DataFrame({'gene_id': ['TGME49_100001'], 'value': [0.], 'imported_value': [False]}), check_exact=True)


def test_host_addresses_and_conflicting_tables_are_checked(monkeypatch):
    source = D.Dataset('host', 'Host', 'gene', 'synthetic', 'Values',
        organism=O.HUMAN, columns=('value',))
    monkeypatch.setattr(D, 'registry', lambda: [source])
    monkeypatch.setattr(R, 'records', lambda: [])
    host = pd.DataFrame({'host_id': ['P00001'], 'value': [2.]})
    model = from_sources({O.TOXOPLASMA: {'tables': {'host_gene': {O.HUMAN: host}}}})
    assert model.rows[0]['organism'] == O.HUMAN and model.rows[0]['unit'] == 'protein'
    pd.testing.assert_frame_equal(model.entities(model.rows[0]), host, check_exact=True)
    with pytest.raises(ValueError, match='Conflicting host'):
        from_sources({O.TOXOPLASMA: {'tables': {'host_gene': {O.HUMAN: host}}},
            O.FALCIPARUM: {'tables': {'host_gene': {O.HUMAN: host.assign(value=3.)}}}})


def test_graph_mapping_retains_exact_index_binding_and_publication_route(monkeypatch):
    source = D.Dataset('pair', 'Pairs', 'protein', 'synthetic', 'Links',
        organism=O.TOXOPLASMA, columns=(S.EDGE_PREFIX + 'example',), pmid='33053376', url='http://example.invalid/source')
    monkeypatch.setattr(D, 'registry', lambda: [source])
    monkeypatch.setattr(R, 'records', lambda: [])
    nodes = pd.DataFrame({'gene_id': ['TGME49_100001', 'TGME49_100002']})
    graph = {'gene_ids': nodes.gene_id.to_numpy(), 'example__a': np.array([0]),
        'example__b': np.array([1]), 'example__w': np.array([0.])}
    model = from_sources({O.TOXOPLASMA: {'nodes': nodes, 'graph': graph}})
    row = model.rows[0]
    assert row['pair_records'] == 1
    assert model.entities(row).iloc[0]['example__w'] == 0.
    card = model.source_card(row)
    assert any(link.url == 'https://pubmed.ncbi.nlm.nih.gov/33053376/' for link in card.links)
    assert 'not an allowed HTTPS' in export_scorecard(card)


def test_scoped_refusals_match_the_canonical_generator_dictionary():
    assert len(R.REFUSED_CANDIDATES) == len(R.records()) == 16
    assert all(row['question'] and row['reason'] for row in R.records())
    assert len(R.records([])) == 0


def test_located_processed_source_remains_visible_without_an_entity_mapping(monkeypatch):
    source = D.Dataset('embedding', 'Embedding', 'protein', 'prediction', 'Encoded values',
        organism=O.TOXOPLASMA, path='starplast/data/embedding.parquet')
    monkeypatch.setattr(D, 'registry', lambda: [source])
    monkeypatch.setattr(D, 'local_path', lambda key: 'located_processed_cache')
    monkeypatch.setattr(R, 'records', lambda: [])
    model = from_sources({})
    row = model.rows[0]
    assert row['status'] == 'source_only' and row['processed_available'] is True
    assert row['raw_available'] is False and row['unit'] == 'source_records'
    assert model.entities(row).attrs['total_rows'] is None
    assert model.source_card(row).metrics == ()
