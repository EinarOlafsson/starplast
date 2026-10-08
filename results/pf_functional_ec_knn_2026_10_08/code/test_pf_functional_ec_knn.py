"""Exact native Pf EC-profile replay, full source closure and artifact scope checks."""
from copy import deepcopy
from dataclasses import asdict, replace
import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from scripts import freeze_pf_functional_ec_controls as C, freeze_pf_functional_ec_knn as K
from starplast import datasets as D, functional_exclusions as E, functional_ontology as F


def fixture(monkeypatch):
    entries = F.parse_enzyme('ID   1.1.1.1\nDE   First.\n//\nID   2.7.1.1\nDE   Kinase.\n//\nID   7.6.2.1\nDE   Transport.\n//\n')
    genes = [f'PF3D7_fixture_{i}' for i in range(96)]
    nodes = pd.DataFrame({'gene_id': genes, 'orthogroup': genes,
        'ec_number': ['1.1.1.1;2.7.1.1' if i % 3 else '1.1.1.1' for i in range(96)],
        'ec_number_orthology': ['7.6.2.1'] * 96,
        'raw_x': [float(i * 17 % 101) / 100 for i in range(96)],
        'raw_y': [float(i * i % 113) / 100 for i in range(96)],
        'opaque_domain_alias': [float(i % 3) for i in range(96)], 'papers_total': [float(i % 7) for i in range(96)],
        'unregistered_map': [float(i % 13) for i in range(96)], 'registered_child': [float(i % 11) for i in range(96)],
        'registered_grandchild': [float(i % 5) for i in range(96)], 'unknown_parent_embedding': [float(i % 4) for i in range(96)]})
    nodes.loc[94, 'ec_number'] = None
    nodes.loc[95, 'ec_number'] = 'malformed'
    first = C.prepare_profiles(nodes, entries)
    nodes.loc[nodes.gene_id.eq(first['split'].entities('test')[0]), 'ec_number'] = '7.6.2.1'
    prepared = C.prepare_profiles(nodes, entries)
    control_rows, _, controls = C.control_records(prepared)
    packet = {'profiles': prepared['direct'], 'groups': prepared['groups'], 'split': prepared['split'],
        'control_rows': control_rows, 'controls': controls, 'gaps': prepared['gaps'],
        'target_identity': C.P._identity(C._records(prepared['direct'])), 'cohort_identity': prepared['cohort_identity'],
        'manifest_hash': 'b' * 64}
    original = D.provenance
    dependencies = {'registered_child': ('unregistered_map',), 'registered_grandchild': ('registered_child',),
                    'unknown_parent_embedding': ('absent_unreviewed_source',)}

    def provenance(column, organism=None):
        assert organism == 'Pf'
        if column in dependencies or column in {'raw_x', 'raw_y', 'opaque_domain_alias', 'papers_total'}:
            return SimpleNamespace(key='synthetic-feature-source', derived_from=dependencies.get(column, ()))
        return original(column, organism)

    monkeypatch.setattr(D, 'provenance', provenance)
    lineage = (E.DerivedInput('column', 'opaque_domain_alias', (('column', 'ec_number_orthology'),)),)
    features, labels, excluded, provenance = K.source_closed_inputs(nodes, packet, derived_inputs=lineage)
    return nodes, packet, features, labels, excluded, provenance


def test_full_source_closure_withholds_direct_transfer_attention_and_late_descendants(monkeypatch):
    _, packet, features, labels, excluded, provenance = fixture(monkeypatch)
    assert features.columns.tolist() == ['raw_x', 'raw_y']
    assert tuple(labels.index) == packet['split'].entities('train')
    assert excluded.fit_entities == tuple(labels.index)
    assert tuple(features.index) == tuple(a.entity for a in packet['split'].assignments)
    for column in ('ec_number', 'ec_number_orthology', 'orthogroup', 'opaque_domain_alias', 'papers_total',
                   'unregistered_map', 'registered_child', 'registered_grandchild', 'unknown_parent_embedding'):
        assert column in excluded.columns
    assert provenance['dependency_closure_withheld']['registered_grandchild']['vetoed_parents'] == ['registered_child']
    assert provenance['sources']['raw_x']['evidence_grade'] == 'unresolved'


def test_native_predictions_scores_votes_and_atomic_unsupported_profiles_are_exact(monkeypatch):
    _, packet, features, labels, excluded, _ = fixture(monkeypatch)
    result = K.fit_candidate(features, labels, packet, excluded)
    assert result['rows'].entity.tolist() == packet['control_rows'].entity.tolist()
    assert result['rows'].truth.tolist() == packet['control_rows'].truth.tolist()
    assert int((~result['rows'].training_supported).sum()) == 1
    assert result['class_scores']['["7"]'].eq(0).all()
    assert '["7"]' not in result['batch'].class_scores.columns
    assert result['batch'].model_state['k'] == 15 and result['batch'].model_state['min_share'] == 0.3
    assert result['batch'].model_state['fit_entities'] == list(packet['split'].entities('train'))
    assert result['card']['extra']['unknown_states'] == {'unannotated': 1, 'unresolved': 1}
    assert result['card']['extra']['biological_accuracy'] is None
    assert sum(item['count'] for item in result['card']['extra']['confusion']) == len(result['rows'])
    assert len(result['major_class_cards']) == 7
    for name, comparison in result['baseline_comparisons'].items():
        assert comparison['source_scope'] == packet['controls'][name]['card']['scope']
        assert comparison['source_scope']['strategy'] == name
        assert comparison['same_entity_order_truth_groups_and_support']
    again = K.fit_candidate(features, labels, packet, excluded)
    assert result['raw_votes'] == again['raw_votes'] and result['card'] == again['card']
    assert result['batch'].model_state == again['batch'].model_state
    pd.testing.assert_frame_equal(result['batch'].class_scores, again['batch'].class_scores, check_exact=True)


def test_heldout_features_cannot_choose_columns_or_training_transforms(monkeypatch):
    nodes, packet, features, labels, excluded, _ = fixture(monkeypatch)
    first = K.fit_candidate(features, labels, packet, excluded)
    nodes.loc[nodes.gene_id.isin(packet['split'].entities('test')), ['raw_x', 'raw_y']] = -10000.0
    second_features, second_labels, second_excluded, _ = K.source_closed_inputs(nodes, packet,
        derived_inputs=(E.DerivedInput('column', 'opaque_domain_alias', (('column', 'ec_number_orthology'),)),))
    assert features.columns.tolist() == second_features.columns.tolist()
    second = K.fit_candidate(second_features, second_labels, packet, second_excluded)
    for key in ('training_vectors', 'training_distributions', 'training_labels', 'kept_columns'):
        assert first['batch'].model_state[key] == second['batch'].model_state[key]


def test_changed_hidden_truth_never_changes_native_model_or_calls(monkeypatch):
    _, packet, features, labels, excluded, _ = fixture(monkeypatch)
    first = K.fit_candidate(features, labels, packet, excluded)
    changed = deepcopy(packet)
    gene = packet['split'].entities('test')[0]
    changed['profiles'].loc[changed['profiles'].gene_id.eq(gene), 'profile'] = '["6"]'
    changed['control_rows'].loc[changed['control_rows'].entity.eq(gene), 'truth'] = '["6"]'
    changed['target_identity'] = 'c' * 64
    second = K.fit_candidate(features, labels, changed, excluded)
    assert first['batch'].model_state == second['batch'].model_state and first['raw_votes'] == second['raw_votes']
    pd.testing.assert_frame_equal(first['batch'].rows, second['batch'].rows, check_exact=True)


def test_frozen_control_truth_or_scope_mismatch_is_refused(monkeypatch):
    _, packet, features, labels, excluded, _ = fixture(monkeypatch)
    changed = deepcopy(packet)
    changed['control_rows'].loc[0, 'truth'] = '["6"]'
    with pytest.raises(AssertionError, match='control scope differs'):
        K.fit_candidate(features, labels, changed, excluded)
    changed = deepcopy(packet)
    changed['controls']['training_majority']['card']['scope']['protocol'] = 'unrelated'
    with pytest.raises(ValueError, match='control card scope'):
        K.fit_candidate(features, labels, changed, excluded)


def test_no_feature_candidate_keeps_abstentions_and_null_reporting_scores(monkeypatch):
    _, packet, features, labels, excluded, _ = fixture(monkeypatch)
    result = K.fit_candidate(features.iloc[:, :0], labels, packet, excluded)
    assert result['model_available'] is False
    assert result['rows'].abstained.all() and result['class_scores'].isna().all().all()
    assert result['card']['counts']['abstained'] == len(result['rows'])
    assert all(item['raw_call'] is None for item in result['raw_votes'])


@pytest.mark.parametrize('available', [True, False])
def test_shared_exact_artifact_factory_roundtrips_scope_and_original_controls(monkeypatch, tmp_path, available):
    from starplast import artifacts as A
    _, packet, _, _, _, _ = fixture(monkeypatch)
    dependencies = tuple(A.Dependency(kind, kind, 'a' * 64) for kind in ('code', 'table', 'truth', 'exclusions', 'model')) + (
        A.Dependency('split', 'prepared', packet['split'].identity),)
    spec = K.artifact_spec(packet, dependencies, model_available=available)
    scope = json.loads(spec.evaluation_scope_json)
    assert scope['unit'] == 'gene' and scope['eligible_population'] == len(packet['split'].entities('test'))
    assert scope['target_identity'] == packet['target_identity'] and scope['prepared_packet_manifest_sha256'] == packet['manifest_hash']
    assert scope['biological_admission'] is False and spec.confidence_kind == 'method_support'
    assert spec.status == ('ok' if available else 'unavailable')
    assert spec.fit_entities == packet['split'].entities('train')
    scores = pd.DataFrame(1.0 if available else np.nan, index=spec.entity_order, columns=['["1","2"]'])
    payload = C._plain({'baseline_cards.json': packet['controls'], 'class_scores.json': K.N._score_payload(scores)})
    identity = A.write_artifact(tmp_path / 'held_out', spec, payload, split=packet['split'])
    restored = A.read_artifact(tmp_path / 'held_out', expected=spec, split=packet['split'])
    assert restored.payloads == payload and restored.identity == identity
    assert restored.payloads['baseline_cards.json'] == packet['controls']
    with pytest.raises(ValueError, match='complete ordered outer-test cohort'):
        A.write_artifact(tmp_path / 'reordered', replace(spec, entity_order=spec.entity_order[::-1]), payload, split=packet['split'])


@pytest.mark.parametrize('change', ['order', 'foreign_split', 'transfer_as_truth', 'missing_source'])
def test_original_source_universe_and_target_guards(monkeypatch, change):
    nodes, packet, _, _, _, _ = fixture(monkeypatch)
    if change == 'order':
        nodes = nodes.iloc[::-1]
    elif change == 'foreign_split':
        packet['split'] = replace(packet['split'], organism='Tg')
    elif change == 'transfer_as_truth':
        packet['profiles'].source_target = 'ec_number_orthology'
    else:
        nodes = nodes.drop(columns='ec_number')
    with pytest.raises(ValueError, match='unchanged Pf direct'):
        K.source_closed_inputs(nodes, packet)


def test_preparation_manifest_pin_rejects_arbitrary_packets(tmp_path):
    (tmp_path / 'SHA256SUMS.txt').write_text('untrusted\n')
    with pytest.raises(ValueError, match='pinned source'):
        K.load_packet(tmp_path)
