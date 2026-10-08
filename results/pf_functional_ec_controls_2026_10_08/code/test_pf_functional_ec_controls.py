"""Synthetic complete EC-profile controls preserve source and native group semantics."""
from copy import deepcopy

import numpy as np
import pandas as pd
import pytest

from scripts import freeze_pf_functional_ec_controls as C
from starplast import functional_ontology as F


def fixture():
    entries = F.parse_enzyme('''ID   1.1.1.1
DE   First activity.
//
ID   2.7.1.1
DE   Kinase activity.
//
ID   3.6.3.1
DE   Transferred entry: 7.6.2.1.
//
ID   7.6.2.1
DE   Transport activity.
//
ID   4.1.1.1
DE   Transferred entry: 1.1.1.1 and 2.7.1.1.
//
ID   5.1.1.1
DE   Deleted entry.
//
''')
    values = ['1.1.1.1;2.7.1.1' if i % 3 else '1.1.1.1' for i in range(96)]
    values[-6:] = [None, None, '4.1.1.1', '5.1.1.1', '6.1.1.1', 'malformed']
    nodes = pd.DataFrame({'gene_id': [f'fixture_{i}' for i in range(96)], 'ec_number': values,
        'ec_number_orthology': ['3.6.3.1'] * 96, 'orthogroup': [f'group_{i // 2}' for i in range(96)]})
    nodes.at[0, 'orthogroup'] = None
    nodes.at[1, 'orthogroup'] = ''
    return nodes, entries


def test_complete_direct_profiles_unknowns_and_transfer_census_remain_separate():
    nodes, entries = fixture()
    original = nodes.copy(deep=True)
    prepared = C.prepare_profiles(nodes, entries)
    assert prepared['direct'].profile.iloc[1] == '["1","2"]'
    assert prepared['direct'].status.value_counts().to_dict() == {'complete': 90, 'unresolved': 4, 'unannotated': 2}
    assert prepared['orthology'].profile.eq('["7"]').all()
    assert prepared['direct'].profile.iloc[-6:].isna().all()
    assert len(prepared['split'].assignments) == 90
    assert prepared['groups'].gene_id.tolist() == nodes.gene_id.tolist()
    assert prepared['groups'].group.iloc[:2].tolist() == ['__0', '__1']
    assert prepared['groups'].fallback.sum() == 2
    assert prepared['groups'].role.iloc[-6:].isna().all()
    assert prepared['exclusions'].fit_entities == prepared['split'].entities('train')
    assert 'ec_number_orthology' in prepared['exclusions'].columns
    assert 'orthogroup' in prepared['exclusions'].columns
    pd.testing.assert_frame_equal(nodes, original, check_exact=True)
    roles = prepared['groups'].dropna(subset='role').groupby('group').role.nunique()
    assert roles.eq(1).all()


def test_control_calls_replay_complete_train_profiles_exactly():
    nodes, entries = fixture()
    prepared = C.prepare_profiles(nodes, entries)
    rows, scores, cards = C.control_records(prepared)
    labels = prepared['labels']
    counts = labels.value_counts().sort_index()
    majority = min(counts.index[counts == counts.max()])
    assert majority == '["1","2"]'
    assert rows.entity.tolist() == list(prepared['split'].entities('test'))
    assert rows.majority.eq(majority).all()
    assert rows.prevalence_call.tolist() == np.random.default_rng(C.SEED).choice(
        counts.index.to_numpy(), len(rows), p=(counts / counts.sum()).to_numpy()).tolist()
    for name, bundle in cards.items():
        card = bundle['card']
        assert card['scope']['strategy'] == name
        assert card['extra']['classifier_fitting_performed'] is False
        assert card['extra']['biological_accuracy'] is None
        assert card['extra']['unknown_states'] == {'unresolved': 4, 'unannotated': 2}
        assert card['counts']['eligible'] == len(rows) == card['counts']['answered']
        assert card['counts']['abstained'] == 0
        assert sum(item['count'] for item in card['extra']['confusion']) == len(rows)
        assert len(bundle['major_class_cards']) == 7
        assert all(item['biological_precision'] is None for item in bundle['major_class_cards'])
        assert scores[name].columns.tolist() == ['["1","2"]', '["1"]']


def test_unsupported_test_profiles_retained_and_never_choose_control_parameters():
    nodes, entries = fixture()
    first = C.prepare_profiles(nodes, entries)
    rows, scores, _ = C.control_records(first)
    test_gene = first['split'].entities('test')[0]
    changed_nodes = nodes.copy()
    changed_nodes.loc[changed_nodes.gene_id.eq(test_gene), 'ec_number'] = '3.6.3.1'
    changed = C.prepare_profiles(changed_nodes, entries)
    second, second_scores, cards = C.control_records(changed)
    assert changed['split'] == first['split']
    pd.testing.assert_series_equal(changed['labels'], first['labels'], check_exact=True)
    pd.testing.assert_frame_equal(rows[['entity', 'majority', 'prevalence_call']],
        second[['entity', 'majority', 'prevalence_call']], check_exact=True)
    assert int((~second.training_supported).sum()) == 1
    assert second.truth[second.entity.eq(test_gene)].iloc[0] == '["7"]'
    for name in scores:
        pd.testing.assert_frame_equal(scores[name], second_scores[name][scores[name].columns], check_exact=True)
        assert second_scores[name]['["7"]'].eq(0).all()
        class_card = next(item for item in cards[name]['profile_class_cards'] if item['class'] == '["7"]')
        assert class_card['test_truth_genes'] == 1 and class_card['recall'] == 0


def test_orthology_changes_cannot_modify_direct_split_or_controls():
    nodes, entries = fixture()
    first = C.prepare_profiles(nodes, entries)
    nodes.ec_number_orthology = None
    second = C.prepare_profiles(nodes, entries)
    assert first['split'] == second['split']
    pd.testing.assert_frame_equal(first['direct'], second['direct'], check_exact=True)
    first_rows, first_scores, first_cards = C.control_records(first)
    second_rows, second_scores, second_cards = C.control_records(second)
    pd.testing.assert_frame_equal(first_rows, second_rows, check_exact=True)
    assert first_cards == second_cards
    for name in first_scores:
        pd.testing.assert_frame_equal(first_scores[name], second_scores[name], check_exact=True)


def test_native_group_fallback_collision_is_explicitly_refused():
    nodes, entries = fixture()
    nodes.at[2, 'orthogroup'] = '__0'
    with pytest.raises(ValueError, match='fallback collides'):
        C.prepare_profiles(nodes, entries)


def test_unknown_source_and_array_cell_replay_preserves_original_meanings(tmp_path):
    nodes, entries = fixture()
    prepared = C.prepare_profiles(nodes, entries)
    path = tmp_path / 'profiles.parquet'
    prepared['direct'].to_parquet(path, index=False)
    assert C._records(prepared['direct']) == C._records(pd.read_parquet(path))
    assert C._plain({'source': None, 'value': False, 'zero': 0, 'list': np.array(['1', '2'])}) == {
        'source': None, 'value': False, 'zero': 0, 'list': ['1', '2']}


def test_control_training_guard_refuses_hidden_annotation_inputs():
    nodes, entries = fixture()
    prepared = C.prepare_profiles(nodes, entries)
    contaminated = dict(prepared)
    contaminated['labels'] = prepared['direct'].set_index('gene_id').profile.dropna()
    with pytest.raises(ValueError, match='Forbidden'):
        C.control_records(contaminated)


@pytest.mark.parametrize('change', ['duplicate_gene', 'missing_source', 'missing_gene_id'])
def test_original_source_address_failures_are_refused(change):
    nodes, entries = fixture()
    if change == 'duplicate_gene':
        nodes.at[1, 'gene_id'] = nodes.gene_id.iloc[0]
    else:
        nodes = nodes.drop(columns='ec_number' if change == 'missing_source' else 'gene_id')
    with pytest.raises(ValueError, match='original gene/source'):
        C.prepare_profiles(nodes, entries)
