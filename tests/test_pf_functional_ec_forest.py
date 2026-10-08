"""Exact Pf native forest scope, unsupported profiles and shared artifact checks."""
from copy import deepcopy
from dataclasses import replace
import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from scripts import freeze_pf_functional_ec_controls as C, freeze_pf_functional_ec_knn as K
from scripts import freeze_pf_functional_ec_forest as RF
from starplast import datasets as D, functional_exclusions as E, functional_ontology as F


def fixture(monkeypatch):
    entries = F.parse_enzyme('ID   1.1.1.1\nDE   First.\n//\nID   2.7.1.1\nDE   Kinase.\n//\nID   7.6.2.1\nDE   Transport.\n//\n')
    genes = [f'PF3D7_fixture_{i}' for i in range(96)]
    nodes = pd.DataFrame({'gene_id': genes, 'orthogroup': genes,
        'ec_number': ['1.1.1.1;2.7.1.1' if i % 3 else '1.1.1.1' for i in range(96)],
        'ec_number_orthology': ['7.6.2.1'] * 96,
        'raw_x': [float(i * 17 % 101) / 100 for i in range(96)],
        'raw_y': [float(i * i % 113) / 100 for i in range(96)],
        'unregistered_map': [float(i % 13) for i in range(96)],
        'registered_child': [float(i % 11) for i in range(96)], 'papers_total': [float(i % 7) for i in range(96)]})
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

    def provenance(column, organism=None):
        assert organism == 'Pf'
        if column in {'raw_x', 'raw_y', 'registered_child', 'papers_total'}:
            return SimpleNamespace(key='synthetic-feature-source', derived_from=('unregistered_map',) if column == 'registered_child' else ())
        return original(column, organism)

    monkeypatch.setattr(D, 'provenance', provenance)
    monkeypatch.setitem(RF.SETTINGS, 'trees', 7)
    features, labels, excluded, provenance = K.source_closed_inputs(nodes, packet)
    return nodes, packet, features, labels, excluded, provenance


def test_original_source_closure_and_all_unsupported_profiles_are_retained(monkeypatch):
    _, packet, features, labels, excluded, provenance = fixture(monkeypatch)
    assert features.columns.tolist() == ['raw_x', 'raw_y']
    assert {'ec_number', 'ec_number_orthology', 'unregistered_map', 'registered_child', 'papers_total'} <= set(excluded.columns)
    assert provenance['policy'] == E.POLICY_VERSION
    result = RF.fit_candidate(features, labels, packet, excluded)
    assert result['rows'].entity.tolist() == packet['control_rows'].entity.tolist()
    assert int((~result['rows'].training_supported).sum()) == 1
    assert result['class_scores']['["7"]'].eq(0).all()
    assert '["7"]' not in result['batch'].class_scores.columns
    assert result['card']['scope']['strategy'] == 'random_forest'
    assert result['card']['extra']['unknown_states'] == {'unannotated': 1, 'unresolved': 1}
    assert result['card']['extra']['biological_accuracy'] is None
    assert result['model_available'] and len(result['major_class_cards']) == 7
    assert sum(item['count'] for item in result['card']['extra']['confusion']) == len(result['rows'])
    assert result['batch'].model_state['execution_n_jobs'] == 1
    for name, comparison in result['baseline_comparisons'].items():
        assert comparison['source_scope'] == packet['controls'][name]['card']['scope']
        assert comparison['source_scope']['strategy'] == name
        assert comparison['same_entity_order_truth_groups_and_support']
    repeated = RF.fit_candidate(features, labels, packet, excluded)
    assert repeated['batch'].model_state == result['batch'].model_state
    assert repeated['card'] == result['card']
    pd.testing.assert_frame_equal(repeated['class_scores'], result['class_scores'], check_exact=True)


def test_heldout_truth_does_not_enter_native_training_or_calls(monkeypatch):
    _, packet, features, labels, excluded, _ = fixture(monkeypatch)
    first = RF.fit_candidate(features, labels, packet, excluded)
    changed = deepcopy(packet)
    gene = packet['split'].entities('test')[0]
    changed['profiles'].loc[changed['profiles'].gene_id.eq(gene), 'profile'] = '["6"]'
    changed['control_rows'].loc[changed['control_rows'].entity.eq(gene), 'truth'] = '["6"]'
    changed['target_identity'] = 'c' * 64
    second = RF.fit_candidate(features, labels, changed, excluded)
    assert first['batch'].model_state == second['batch'].model_state
    pd.testing.assert_frame_equal(first['batch'].rows, second['batch'].rows, check_exact=True)
    pd.testing.assert_frame_equal(first['batch'].class_scores, second['batch'].class_scores, check_exact=True)


def test_control_cohort_or_protocol_mismatch_is_refused(monkeypatch):
    _, packet, features, labels, excluded, _ = fixture(monkeypatch)
    changed = deepcopy(packet)
    changed['control_rows'].loc[0, 'truth'] = '["6"]'
    with pytest.raises(AssertionError, match='control scope differs'):
        RF.fit_candidate(features, labels, changed, excluded)
    changed = deepcopy(packet)
    changed['controls']['training_majority']['card']['scope']['protocol'] = 'unrelated'
    with pytest.raises(ValueError, match='control card scope'):
        RF.fit_candidate(features, labels, changed, excluded)


def test_unavailable_features_retain_null_scores_and_every_abstention(monkeypatch):
    _, packet, features, labels, excluded, _ = fixture(monkeypatch)
    result = RF.fit_candidate(features.iloc[:, :0], labels, packet, excluded)
    assert not result['model_available'] and result['rows'].abstained.all()
    assert result['class_scores'].isna().all().all()
    assert result['card']['counts']['abstained'] == len(result['rows'])
    assert result['card']['extra']['model_status'] == 'unavailable_native_forest'
    assert all(value is None for row in K.N._score_payload(result['class_scores'])['scores'] for value in row)


@pytest.mark.parametrize('available', [True, False])
def test_exact_shared_artifact_factory_restores_native_strategy_and_controls(monkeypatch, tmp_path, available):
    from starplast import artifacts as A, capabilities as Cap
    _, packet, _, _, _, _ = fixture(monkeypatch)
    dependencies = tuple(A.Dependency(kind, kind, 'a' * 64) for kind in ('code', 'table', 'truth', 'exclusions', 'model')) + (
        A.Dependency('split', 'prepared', packet['split'].identity),)
    spec = RF.artifact_spec(packet, dependencies, model_available=available)
    scope = json.loads(spec.evaluation_scope_json)
    assert spec.strategy == 'random_forest' and spec.output == Cap.get('random_forest').outputs[0]
    assert scope['unit'] == 'gene' and scope['eligible_population'] == len(packet['split'].entities('test'))
    assert scope['target_identity'] == packet['target_identity'] and scope['prepared_packet_manifest_sha256'] == packet['manifest_hash']
    assert scope['native_n_jobs'] == 4 and scope['execution_n_jobs'] == 1
    assert not scope['biological_admission'] and spec.confidence_kind == 'method_support'
    assert spec.status == ('ok' if available else 'unavailable')
    assert spec.fit_entities == packet['split'].entities('train')
    scores = pd.DataFrame(1.0 if available else np.nan, index=spec.entity_order, columns=['["1","2"]'])
    payload = C._plain({'baseline_cards.json': packet['controls'], 'class_scores.json': K.N._score_payload(scores)})
    identity = A.write_artifact(tmp_path / 'held_out', spec, payload, split=packet['split'])
    restored = A.read_artifact(tmp_path / 'held_out', expected=spec, split=packet['split'])
    assert restored.payloads == payload and restored.identity == identity
    with pytest.raises(ValueError, match='complete ordered outer-test cohort'):
        A.write_artifact(tmp_path / 'reordered', replace(spec, entity_order=spec.entity_order[::-1]), payload, split=packet['split'])
