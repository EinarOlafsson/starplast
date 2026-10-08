"""Native EC controls retain independent scopes and unknown biological accuracy."""
from copy import deepcopy
import json

import pytest

from starplast import functional_results as F
from tests.test_functional_results import _document, _repair_source_hashes, _write


def _native_document():
    document=_document();entry=document['benchmarks'][0]
    metadata=entry['metadata'];payload=entry['payloads'];spec=metadata['source_manifest']['spec']
    metadata.update(control_format='native_ec_controls_v1',organism='Pf',target='ec_direct_complete_major_profile')
    scope=payload['card.json']['scope'];scope.update(organism='Pf',target=metadata['target'])
    for card in payload['profile_class_cards.json']:card['scope']=deepcopy(scope)
    spec['target']=metadata['target']
    query=json.loads(spec['query_json']);query.update(organism='Pf',target=metadata['target']);spec['query_json']=F._canonical(query)
    evaluation=json.loads(spec['evaluation_scope_json'])
    evaluation.update(target_identity='1'*64,cohort_identity='2'*64);spec['evaluation_scope_json']=F._canonical(evaluation)
    controls={};comparisons={}
    for old,name in (('majority','training_majority'),('prevalence_call','training_prevalence')):
        card=payload['baseline_cards.json'][old];card['scope']=deepcopy(scope);card['scope']['strategy']=name
        card['extra']={'control_name':name,'baseline_parameters_estimated_from_training_only':True,
            'classifier_fitting_performed':False,'test_support_used_for_selection':False,
            'biological_admission':False,'biological_accuracy':None}
        controls[name]={'card':card}
        comparisons[name]={'source_scope':deepcopy(card['scope']),'counts':deepcopy(card['counts']),
            'metrics':deepcopy(card['metrics']),'same_entity_order_truth_groups_and_support':True}
    payload['baseline_cards.json']=controls
    for row in payload['rows.json']:row['training_supported']=True
    metadata['source_summary']=deepcopy(metadata['summary'])
    metadata['source_summary'].update(target_identity='1'*64,cohort_identity='2'*64,
        split_identity=scope['protocol'],unsupported_test_genes=0,baseline_comparisons=comparisons)
    _repair(document)
    return document


def _repair(document):
    _repair_source_hashes(document)
    metadata=document['benchmarks'][0]['metadata']
    metadata['source_summary']['artifact_identity']=metadata['summary']['artifact_identity']


def test_native_control_scopes_and_payload_names_remain_original(tmp_path):
    document=_native_document();path,digest=_write(tmp_path,document)
    benchmark,=F.load(path,expected_sha256=digest)
    assert set(benchmark.baseline_cards)=={'training_majority','training_prevalence'}
    assert benchmark.payloads==document['benchmarks'][0]['payloads']
    assert all(card['scope']['strategy']==name for name,card in benchmark.baseline_cards.items())


@pytest.mark.parametrize('mutation', ['format','source','protocol','truth','classifier','selection','biology',
    'unsupported','comparison','cohort'])
def test_native_control_refusals_survive_repaired_envelope(tmp_path,mutation):
    document=_native_document();entry=document['benchmarks'][0];metadata=entry['metadata']
    card=entry['payloads']['baseline_cards.json']['training_majority']['card']
    if mutation=='format':metadata['control_format']='unknown'
    elif mutation=='source':metadata['source_targets']=['ec_number_orthology']
    elif mutation=='protocol':card['scope']['protocol']='0'*64
    elif mutation=='truth':card['scope']['truth_grade']='measured'
    elif mutation=='classifier':card['extra']['classifier_fitting_performed']=True
    elif mutation=='selection':card['extra']['test_support_used_for_selection']=True
    elif mutation=='biology':card['extra']['biological_admission']=True
    elif mutation=='unsupported':entry['payloads']['rows.json'][0]['training_supported']=False
    elif mutation=='comparison':metadata['source_summary']['baseline_comparisons']['training_majority']['metrics']={}
    elif mutation=='cohort':metadata['source_summary']['cohort_identity']='3'*64
    _repair(document);path,digest=_write(tmp_path,document)
    with pytest.raises(ValueError):F.load(path,expected_sha256=digest)
