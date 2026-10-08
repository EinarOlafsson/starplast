"""Namespace-aware readers retain exact reference profiles and reject invented biological scope."""
from copy import deepcopy
import hashlib
import json
from pathlib import Path

import pytest

from starplast import functional_results as F


def _document():
    shipped=json.loads((Path(F.__file__).with_name('data')/'functional_results.json').read_text())
    entry=deepcopy(next(entry for entry in shipped['benchmarks'] if entry['metadata'].get('profile_namespace','ec_major')=='ec_major'))
    metadata=entry['metadata'];spec=metadata['source_manifest']['spec']
    scope=deepcopy(entry['payloads']['card.json']['scope'])
    metadata.update(profile_namespace='pfam',source_targets=['pfam_id'],target='pfam_id_complete_profile',benchmark_id='software-pfam-fixture')
    scope.update(target=metadata['target'],benchmark_id=metadata['benchmark_id'])
    spec.update(target=metadata['target'],benchmark_id=metadata['benchmark_id'])
    entities=spec['entity_order'][:3];spec['entity_order']=entities
    split={'organism':metadata['organism'],'benchmark_id':metadata['benchmark_id'],'group_kind':'homology',
        'seed':scope['seed'],'feature_access':'inductive','protocol_version':1,
        'assignments':[{'entity':entity,'group':entity,'role':'train'} for entity in spec['fit_entities']]
            +[{'entity':'software-tune','group':'software-tune','role':'tune'},
              {'entity':'software-cal','group':'software-cal','role':'calibration'}]
            +[{'entity':entity,'group':entity,'role':'test'} for entity in entities]}
    split_identity=F._hash(split);scope['protocol']=split_identity
    spec['dependencies']=[d for d in spec['dependencies'] if d['kind']!='split']+[{'kind':'split','name':'fixture','sha256':split_identity}]
    truths=['["PF00001"]','["PF00002"]','["PF00001","PF00002"]']
    predictions=['["PF00001"]','["PF00001"]',None]
    rows=[{'entity':entity,'truth':truth,'prediction':pred,'abstained':pred is None,
        'calibrated_confidence':None,'training_supported':i!=2} for i,(entity,truth,pred) in enumerate(zip(entities,truths,predictions))]
    card=deepcopy(entry['payloads']['card.json']);card['scope']=scope
    card['counts']={'eligible':3,'answered':2,'abstained':1,'correct':1,'wrong':1}
    profile_cards=[]
    for term in sorted(truths):
        actual=[row['truth']==term for row in rows];called=[row['prediction']==term for row in rows]
        tp=sum(t and p for t,p in zip(actual,called));fp=sum(not t and p for t,p in zip(actual,called));fn=sum(t and not p for t,p in zip(actual,called))
        precision=tp/(tp+fp) if tp+fp else 0.;recall=tp/(tp+fn)
        local=deepcopy(card);local['class']=term;local['counts']['eligible']=sum(actual)
        local['class_metrics']={'precision':precision,'recall':recall,'f1':2*precision*recall/(precision+recall) if precision+recall else 0.,
            'true_positive':tp,'false_positive':fp,'false_negative':fn,'evaluation_prevalence':sum(actual)/3}
        profile_cards.append(local)
    evaluation=json.loads(spec['evaluation_scope_json']);evaluation.update(target_identity='a'*64,cohort_identity='b'*64)
    spec['evaluation_scope_json']=F._canonical(evaluation)
    controls={}
    for name in ('training_majority','training_prevalence'):
        control=deepcopy(card)
        control['scope']={key:scope[key] for key in ('organism','target','seed','task','truth_grade','unit','negative_semantics')}
        control['scope'].update(strategy=name,partition='test',benchmark_id='software-control-fixture')
        control['extra']={'baseline_parameters_estimated_from_training_only':True,'classifier_or_feature_fitting_performed':False,
            'test_support_used_for_selection':False,'biological_accuracy':None,'split_identity':scope['protocol'],
            'target_identity':'a'*64,'cohort_identity':'b'*64}
        controls[name]={'card':control,'view_snapshot':{}}
    payloads={'rows.json':rows,'card.json':card,'class_cards.json':profile_cards,'baseline_cards.json':controls,
        'major_class_cards.json':F.member_cards(rows,'pfam')}
    entry['payloads']=payloads
    original={'train_genes':len(spec['fit_entities']),'test_genes':3,'eligible_genes':len(spec['fit_entities'])+5,
        'unsupported_test_genes':1,'artifact_identity':'pending','counts':card['counts'],'metrics':card['metrics'],
        'baseline_comparisons':{name:{'metrics':control['card']['metrics']} for name,control in controls.items()}}
    mapping={'train':'train_genes','test':'test_genes','eligible_profiles':'eligible_genes'}
    metadata.update(source_summary=original,summary={**original,**{k:original[v] for k,v in mapping.items()},'tune':1,'calibration':1},summary_mapping=mapping,
        split_manifest=split,source_split_receipt={'path':'synthetic-split.json','sha256':'c'*64,'bytes':123,'identity':split_identity},
        summary_role_mapping={role:role for role in ('train','tune','calibration','test')})
    document={'schema_version':1,'benchmarks':[entry]}
    _repair(document)
    return document


def _repair(document):
    entry=document['benchmarks'][0];metadata=entry['metadata'];manifest=metadata['source_manifest']
    for name in ('rows.json','card.json','class_cards.json','baseline_cards.json'):
        data=(F._canonical(entry['payloads'][name])+'\n').encode()
        manifest['files'][name]={'sha256':hashlib.sha256(data).hexdigest(),'bytes':len(data)}
    manifest['key']=F._hash(manifest['spec']);manifest['identity']=F._hash({k:v for k,v in manifest.items() if k!='identity'})
    metadata['source_artifact_identity']=manifest['identity']
    metadata['summary']['artifact_identity']=metadata['source_summary']['artifact_identity']=manifest['identity']
    metadata['membership_derivation']={'recipe':'recorded_profile_member_cards_v1','namespace':'pfam',
        'source_rows_sha256':manifest['files']['rows.json']['sha256'],'payload_sha256':F._hash(entry['payloads']['major_class_cards.json'])}


def _load(tmp_path,document):
    path=tmp_path/'fixture.json';path.write_text(json.dumps(document,allow_nan=False))
    return F.load(path,expected_sha256=hashlib.sha256(path.read_bytes()).hexdigest())


def test_complete_pfam_namespace_keeps_overlaps_support_and_distinct_control_scopes(tmp_path):
    document=_document();benchmark=_load(tmp_path,document)[0]
    assert benchmark.namespace=='pfam' and benchmark.source_matches('pfam_id') and not benchmark.source_matches('ec_number')
    assert len(benchmark.rows)==3 and len(benchmark.profile_class_cards)==3
    assert list(benchmark.baseline_cards)==['training_majority','training_prevalence']
    assert all(card['scope']!=benchmark.card['scope'] for card in benchmark.baseline_cards.values())
    cards={card['class']:card for card in benchmark.major_class_cards}
    assert cards['PF00001']['true_positive']==1 and cards['PF00001']['false_positive']==1 and cards['PF00001']['false_negative']==1
    assert cards['PF00002']['recall']==0 and cards['PF00002']['precision'] is None
    assert 'not verified biological absence' in cards['PF00002']['complement_semantics']
    assert all(card['biological_precision'] is None for card in cards.values())
    assert sum(row['training_supported'] is False for row in benchmark.rows)==1
    assert benchmark.payloads['class_cards.json']==document['benchmarks'][0]['payloads']['class_cards.json']


@pytest.mark.parametrize('mutation',['wrong_source','wrong_namespace','missing_domain','biological_member','member_rate',
    'missing_profile','profile_precision','support_count','missing_control','control_scope','control_fit','control_selection',
    'role_count','split_order','split_identity'])
def test_self_consistent_pfam_envelopes_do_not_override_recorded_semantics(tmp_path,mutation):
    document=_document();entry=document['benchmarks'][0];metadata=entry['metadata'];payloads=entry['payloads']
    if mutation=='wrong_source':metadata['source_targets']=['ec_number']
    elif mutation=='wrong_namespace':metadata['profile_namespace']='ec_major'
    elif mutation=='missing_domain':payloads['major_class_cards.json'].pop()
    elif mutation=='biological_member':payloads['major_class_cards.json'][0]['biological_precision']=.9
    elif mutation=='member_rate':payloads['major_class_cards.json'][0]['precision']=.99
    elif mutation=='missing_profile':payloads['class_cards.json'].pop()
    elif mutation=='profile_precision':payloads['class_cards.json'][0]['class_metrics']['precision']=.99
    elif mutation=='support_count':payloads['rows.json'][0]['training_supported']=False
    elif mutation=='missing_control':payloads['baseline_cards.json'].pop('training_majority')
    elif mutation=='control_scope':payloads['baseline_cards.json']['training_majority']['card']['extra']['split_identity']='0'*64
    elif mutation=='control_fit':payloads['baseline_cards.json']['training_majority']['card']['extra']['classifier_or_feature_fitting_performed']=True
    elif mutation=='control_selection':payloads['baseline_cards.json']['training_majority']['card']['extra']['test_support_used_for_selection']=True
    elif mutation=='role_count':metadata['summary']['calibration']=2
    elif mutation=='split_order':metadata['split_manifest']['assignments'].reverse()
    elif mutation=='split_identity':metadata['source_split_receipt']['identity']='0'*64
    _repair(document)
    with pytest.raises((ValueError,KeyError)):_load(tmp_path,document)


def test_profile_namespace_validation_rejects_versions_duplicates_and_ec_pfam_mix():
    assert F.profile_classes('["PF00001","PF00002"]','pfam')==('PF00001','PF00002')
    assert F.profile_classes(None,'pfam') is None
    for value in ('[]','["PF00001.2"]','["PF00001","PF00001"]','["PF00002","PF00001"]','["2"]'):
        with pytest.raises(ValueError):F.profile_classes(value,'pfam')
    with pytest.raises(ValueError):F.profile_classes('["PF00001"]')
    with pytest.raises(ValueError):F.profile_classes(None,'unreviewed')


def test_domain_titles_use_only_verified_current_nomenclature(monkeypatch):
    from starplast import functional_domains as D
    class Lookup:
        def lookup(self,_value):return {'status':'current_metadata','name':'Reference family name'}
    monkeypatch.setattr(D,'shipped',lambda:Lookup())
    assert F.class_title('PF00001','pfam')=='PF00001 — Reference family name (current nomenclature)'
    assert F.class_title('["PF00001"]','pfam')=='["PF00001"]'
    assert F.class_title('2')=='2 — Transferases'
