"""Synthetic exact native replay and source closure for a fixed complete-profile candidate."""
from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from scripts import freeze_functional_profile_knn as K
from starplast import datasets as D, functional_domain_profiles as F, functional_exclusions as E
from starplast.functional_profile_targets import freeze_target
from starplast.splits import Assignment, SplitManifest


def fixture(monkeypatch):
    genes = tuple(f'fixture_{i}' for i in range(64))
    values = ['PF00001' if i % 2 == 0 else 'PF00001;PF00002' for i in range(64)]
    values[44:53] = ['PF00003;PF00004'] * 9
    values[-2:] = [None, 'PF00005;bad']
    roles = ['train'] * 36 + ['tune'] * 4 + ['calibration'] * 4 + ['test'] * 18
    profiles = pd.DataFrame([{'organism': 'Tg', 'gene_id': gene, 'source_target': 'pfam_id',
        **F.parse_source_profile(value, 'pfam_id'), 'source_ids': ['synthetic-pfam-source'],
        'source_release': 'fixture-only', 'evidence_grade': 'synthetic_control',
        'nomenclature_releases': ['Pfam:synthetic'], 'negative_semantics': F.NEGATIVE_SEMANTICS,
        'biological_admission': False} for gene, value in zip(genes, values)])
    split = SplitManifest('Tg', 'synthetic-knn-fixture', 'entity',
        tuple(Assignment(gene, gene, role) for gene, role in zip(genes, roles)), 7)
    contract = freeze_target(profiles, organism='Tg', source_target='pfam_id', universe=genes,
        protected_groups=genes, split=split, target='fixture_complete_profile')
    nodes = pd.DataFrame({'gene_id': genes, 'pfam_id': values,
        'raw_x': [float(i * 17 % 101) / 100 for i in range(64)],
        'raw_y': [float(i * i % 113) / 100 for i in range(64)],
        'opaque_domain_alias': [float(i % 2) for i in range(64)],
        'pfam_annotation_digest': [float(i % 2) for i in range(64)],
        'citations_total': [float(i % 4) for i in range(64)],
        'unregistered_map': [float(i) for i in range(64)],
        'unresolved_registered_embedding': [float(i % 3) for i in range(64)]})
    original = D.provenance

    def provenance(column, organism=None):
        if column == 'unresolved_registered_embedding':
            return SimpleNamespace(key='synthetic-opaque', derived_from=('unresolved_operator',))
        if column in {'raw_x', 'raw_y', 'opaque_domain_alias', 'citations_total'}:
            return SimpleNamespace(key='synthetic-numeric-source', derived_from=())
        return original(column, organism)

    monkeypatch.setattr(D, 'provenance', provenance)
    lineage = (E.DerivedInput('column', 'opaque_domain_alias', (('column', 'pfam_id'),)),)
    features, excluded, metadata = K.prepare_inputs(nodes, contract, derived_inputs=lineage)
    return nodes, contract, features, excluded, metadata


def test_source_closed_registered_numeric_inputs_and_train_selection(monkeypatch):
    nodes, contract, features, excluded, metadata = fixture(monkeypatch)
    assert features.columns.tolist() == ['raw_x', 'raw_y']
    assert tuple(features.index) == tuple(a.entity for a in contract.split.assignments)
    assert excluded.fit_entities == contract.split.entities('train')
    assert metadata['selection_fit_entities'] == list(contract.split.entities('train'))
    for column in ('pfam_id', 'opaque_domain_alias', 'citations_total', 'unregistered_map', 'unresolved_registered_embedding'):
        assert column in excluded.columns
    assert metadata['sources']['raw_x']['evidence_grade'] == 'unresolved'
    changed = nodes.copy()
    changed.loc[36:, ['raw_x', 'raw_y']] = 10000.0
    second, _, second_metadata = K.prepare_inputs(changed, contract,
        derived_inputs=(E.DerivedInput('column', 'opaque_domain_alias', (('column', 'pfam_id'),)),))
    assert second.columns.tolist() == features.columns.tolist()
    assert second_metadata == metadata


def test_exact_native_ranks_votes_classes_and_unsupported_denominators(monkeypatch):
    _, contract, features, excluded, _ = fixture(monkeypatch)
    result = K.fit_recovery(features, contract, excluded)
    assert result['rows'].entity.tolist() == list(contract.split.entities('test'))
    assert len(result['rows']) == 18 and int((~result['rows'].training_supported).sum()) == 9
    assert result['batch'].model_state['k'] == 15 and result['batch'].model_state['min_share'] == 0.3
    assert result['batch'].model_state['native_classes'] == ['["PF00001","PF00002"]', '["PF00001"]']
    assert result['class_scores']['["PF00003","PF00004"]'].eq(0).all()
    assert result['batch'].class_scores.columns.tolist() != result['class_scores'].columns.tolist()
    assert result['card']['counts']['eligible'] == 18
    assert result['card']['extra']['unknown_states'] == {'malformed': 1, 'unannotated': 1}
    assert result['card']['extra']['biological_accuracy'] is None
    assert result['card']['extra']['deployment_applicability'] == 'unknown'
    assert all(len(record['neighbors']) == 15 for record in result['raw_votes'])
    assert all(not (set(record['neighbors']) & set(contract.split.entities('test'))) for record in result['raw_votes'])
    assert result['rows'].calibrated_confidence.isna().all()
    for comparison in result['baseline_comparisons'].values():
        assert comparison['same_entity_order'] and comparison['same_truth'] and comparison['same_class_order']
    np.testing.assert_array_equal(result['batch'].model_state['training_vectors'],
        (features.loc[list(contract.split.entities('train'))].rank(pct=True) - 0.5).to_numpy())


def test_changed_test_truth_cannot_change_model_votes_or_native_scores(monkeypatch):
    _, contract, features, excluded, _ = fixture(monkeypatch)
    first = K.fit_recovery(features, contract, excluded)
    changed_rows = list(contract.rows)
    changed_rows[44] = replace(changed_rows[44], profile='["PF99999"]', recorded_identifiers=('PF99999',))
    changed = replace(contract, rows=tuple(changed_rows))
    second = K.fit_recovery(features, changed, excluded)
    assert first['batch'].model_state == second['batch'].model_state
    assert first['raw_votes'] == second['raw_votes']
    pd.testing.assert_frame_equal(first['batch'].rows, second['batch'].rows, check_exact=True)
    pd.testing.assert_frame_equal(first['batch'].class_scores, second['batch'].class_scores, check_exact=True)


def test_frozen_training_distributions_ignore_heldout_values(monkeypatch):
    _, contract, features, excluded, _ = fixture(monkeypatch)
    first = K.fit_recovery(features, contract, excluded)
    changed = features.copy()
    changed.loc[list(contract.split.entities('test'))] = -10000.0
    second = K.fit_recovery(changed, contract, excluded)
    for key in ('training_distributions', 'training_vectors', 'kept_columns', 'training_labels', 'fit_entities'):
        assert first['batch'].model_state[key] == second['batch'].model_state[key]


def test_native_replay_is_exact_and_rejects_changed_rank_state(monkeypatch):
    _, contract, features, excluded, _ = fixture(monkeypatch)
    first = K.fit_recovery(features, contract, excluded)
    replay = K.fit_recovery(features, contract, excluded)
    assert first['raw_votes'] == replay['raw_votes'] and first['card'] == replay['card']
    changed = dict(first['batch'].model_state)
    changed['training_vectors'] = np.zeros_like(changed['training_vectors']).tolist()
    with pytest.raises(AssertionError):
        K.native_replay(features, contract.training_labels(), contract.split, changed)
    with pytest.raises(ValueError, match='exact split order'):
        K.native_replay(features.iloc[::-1], contract.training_labels(), contract.split, first['batch'].model_state)


def test_no_features_retains_abstentions_and_unknown_biology(monkeypatch):
    _, contract, features, excluded, _ = fixture(monkeypatch)
    result = K.fit_recovery(features.iloc[:, :0], contract, excluded)
    assert len(result['rows']) == 18 and result['rows'].abstained.all()
    assert result['class_scores'].isna().all().all()
    assert result['card']['counts']['answered'] == 0 and result['card']['counts']['abstained'] == 18
    assert all(record['raw_call'] is None for record in result['raw_votes'])
    assert result['model_available'] is False
    assert result['card']['extra']['model_status'] == 'unavailable_no_supported_native_vote'
    payload = K._json(K._score_payload(result['class_scores']))
    assert payload['classes'] == result['class_scores'].columns.tolist()
    assert payload['scores'] == [[None] * len(payload['classes']) for _ in range(18)]


def test_registered_descendants_of_late_vetoed_inputs_are_closed_with_reasons(monkeypatch):
    nodes, contract, _, _, _ = fixture(monkeypatch)
    nodes['registered_child'] = nodes.raw_x * 2
    nodes['registered_grandchild'] = nodes.raw_x * 3
    nodes['missing_parent_child'] = nodes.raw_x * 4
    nodes['missing_parent_grandchild'] = nodes.raw_x * 5
    original = D.provenance
    dependencies = {'registered_child': ('unregistered_map',), 'registered_grandchild': ('registered_child',),
        'missing_parent_child': ('unresolved_registered_embedding',), 'missing_parent_grandchild': ('missing_parent_child',)}

    def provenance(column, organism=None):
        if column in dependencies:
            assert organism == 'Tg'
            return SimpleNamespace(key='synthetic-registered-child', derived_from=dependencies[column])
        return original(column, organism)

    monkeypatch.setattr(D, 'provenance', provenance)
    features, excluded, metadata = K.prepare_inputs(nodes, contract,
        derived_inputs=(E.DerivedInput('column', 'opaque_domain_alias', (('column', 'pfam_id'),)),))
    assert features.columns.tolist() == ['raw_x', 'raw_y']
    for column, parents in dependencies.items():
        assert column in excluded.columns
        assert metadata['dependency_closure_withheld'][column]['vetoed_parents'] == list(parents)


def test_late_registry_closure_is_organism_scoped(monkeypatch):
    nodes, contract, _, _, _ = fixture(monkeypatch)
    original = D.derived_dependents
    addresses = []

    def dependents(columns, organism=None):
        addresses.append(organism)
        return original(columns, organism)

    monkeypatch.setattr(D, 'derived_dependents', dependents)
    K.prepare_inputs(nodes, contract)
    # Existing shared guard may conservatively call its global closure; the
    # final late-input closure explicitly repeats it for this organism.
    assert addresses[-1] == 'Tg'


@pytest.mark.parametrize('change', ['order', 'missing_source', 'existing_target', 'foreign_organism'])
def test_source_address_or_original_gene_order_mismatch_refused(monkeypatch, change):
    nodes, contract, _, _, _ = fixture(monkeypatch)
    if change == 'order':
        nodes = nodes.iloc[::-1]
    elif change == 'missing_source':
        nodes = nodes.drop(columns='pfam_id')
    elif change == 'existing_target':
        nodes[contract.target] = 'invented'
    else:
        contract = replace(contract, organism='Pf', split=replace(contract.split, organism='Pf'))
    with pytest.raises(ValueError):
        K.prepare_inputs(nodes, contract)
