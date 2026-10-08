"""Synthetic train-fit operators preserve native ranks, source identities and DAG roles."""
from dataclasses import replace
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from starplast import artifacts as A,capabilities as C,organisms as O,precompute_features as F,precompute_jobs as P
from starplast.query import Query
from starplast.splits import Assignment,ExclusionManifest,SplitManifest

IDS=tuple(f'TGME49_{100001+i}' for i in range(10))
SPLIT=SplitManifest(O.TOXOPLASMA,'numeric_software_fixture','entity',tuple(
    Assignment(entity,entity,role) for entity,role in zip(IDS,['train']*5+['tune','calibration']+['test']*3)),3)
EXCLUSIONS=ExclusionManifest(SPLIT.benchmark_id,('target',),('target','sibling_source'),(),IDS[:5],SPLIT.identity)
LIMITS=F.FeatureLimits(20,10,200,100000)
BUDGET=P.Budget(10,30,2000000,2000000,400*1024*1024)


@pytest.fixture(autouse=True)
def _small_operator_rss(monkeypatch):
    """Test small-job semantics independently of an earlier Qt process's resident memory.

    The entire operator suite also passes in a real external 400 MB scope. A
    combined application suite retains Qt/context memory before these tiny jobs.
    Production current-RSS admission and external caps remain unchanged.
    """
    monkeypatch.setattr(P,'_rss',lambda:128*1024*1024)


def _features():
    return pd.DataFrame({'x':[1.,1.,3.,None,7.,123.,123.,0.,2.,8.],
        'y':[None,2.,2.,4.,8.,123.,123.,None,3.,100.],
        'missing_train':[None]*5+[1.,2.,3.,4.,5.]},index=IDS)


def _builder(tmp_path,features=None,*,split=SPLIT,exclusions=EXCLUSIONS,limits=LIMITS,name='source.json'):
    features=_features() if features is None else features
    path=tmp_path/name;path.write_text(A.canonical_object(F.snapshot_payload(features,organism=split.organism))+'\n')
    snapshot=P.Snapshot(A.Dependency('table','numeric_inputs',hashlib.sha256(path.read_bytes()).hexdigest()),str(path))
    return F.NumericFeatureBuilder(features,snapshot,split=split,exclusions=exclusions,limits=limits)


def _spec(builder):
    cap=C.get('feature_knn');available=builder.availability()
    return A.ArtifactSpec(Query(O.TOXOPLASMA,'label',target='target').to_json(),'feature_knn',cap.benchmark_tasks[0],
        cap.outputs[0],'target',tuple(assignment.entity for assignment in builder.split.assignments),'reusable_base',
        'numeric-software-fixture:'+builder.split.identity,
        tuple(snapshot.dependency for snapshot in builder.snapshots)+(A.Dependency('split','nested',builder.split.identity),
         A.Dependency('exclusions','source_closure',F.exclusions_identity(builder.exclusions))),
        A.canonical_object(F.settings(builder.limits)),A.canonical_object({'unit':'gene','eligible_population':len(IDS),
        'truth_grade':'synthetic_control','negative_semantics':'Software fixture only','context':{},'limitations':['No biological validation']}),
        3,'software-fixture-v1',fit_entities=builder.split.entities('train'),fit_role='train',
        benchmark_id=builder.split.benchmark_id,status=available['status'],gaps=tuple(available['gaps']))


def _context(builder,*,spec=None,checkpoint=lambda:None):
    spec=_spec(builder) if spec is None else spec
    return P.Context(P.Job('features',spec,split=builder.split),spec,{},checkpoint)


def test_native_training_ranks_and_heldout_ecdf_are_exact_and_missingness_is_explicit(tmp_path):
    builder=_builder(tmp_path)
    payload=builder(_context(builder))['operator.json']
    assert payload['columns']==['x','y'] and payload['withheld_all_missing_training_columns']==['missing_train']
    native=(_features().loc[list(IDS[:5]),['x','y']].rank(pct=True)-.5).fillna(0).to_numpy()
    np.testing.assert_array_equal(np.array(payload['matrix'])[:5],native)
    np.testing.assert_array_equal(np.array(payload['matrix'])[-3:],np.array([[-.5,0.],[0.,0.],[.5,.5]]))
    assert payload['observed'][-3:]==[[True,False],[True,True],[True,True]]
    assert payload['entity_order']==list(IDS) and payload['fit_entities']==list(IDS[:5])
    assert payload['inference_outputs'] is False and payload['biological_accuracy'] is None and payload['calibrated_confidence'] is None
    assert json.loads(json.dumps(payload,allow_nan=False))==payload


def test_heldout_feature_sentinels_cannot_change_fitted_column_selection_or_training_state(tmp_path):
    first=_builder(tmp_path,name='first.json');a=first(_context(first))['operator.json']
    changed=_features();changed.loc[list(IDS[5:]),:]=10**12
    second=_builder(tmp_path,changed,name='second.json');b=second(_context(second))['operator.json']
    for field in ('columns','withheld_all_missing_training_columns','training_distributions','fit_entities'):
        assert a[field]==b[field]
    np.testing.assert_array_equal(np.array(a['matrix'])[:5],np.array(b['matrix'])[:5])
    assert a['source_sha256']!=b['source_sha256']


@pytest.mark.parametrize('case',['target','sibling','attention','publication_attention','source_content','source_bytes',
    'source_dependency','exclusion_dependency','recipe','order','fit_role'])
def test_input_source_closure_identity_and_fit_role_mismatches_are_refused(tmp_path,case):
    features=_features()
    if case=='target':features['target']=1.
    elif case=='sibling':features['sibling_source']=1.
    elif case=='attention':features['attention_count']=1.
    elif case=='publication_attention':features['n_publications']=1.
    builder=_builder(tmp_path,features)
    if case in {'target','sibling'}:
        with pytest.raises(ValueError,match='contamination'):builder.availability()
        return
    spec=_spec(builder)
    if case=='source_content':builder.features.iloc[0,0]=np.nextafter(1.,2.)
    elif case=='source_bytes':
        path=Path(builder.snapshot.path)
        path.write_bytes(path.read_bytes()+b' ')
    elif case=='source_dependency':spec=replace(spec,dependencies=tuple(replace(d,sha256='a'*64) if d.kind=='table' else d for d in spec.dependencies))
    elif case=='exclusion_dependency':spec=replace(spec,dependencies=tuple(replace(d,sha256='a'*64) if d.kind=='exclusions' else d for d in spec.dependencies))
    elif case=='recipe':spec=replace(spec,settings_json=A.canonical_object({'phase':'other'}))
    elif case=='order':builder.features=builder.features.iloc[::-1]
    elif case=='fit_role':spec=replace(spec,fit_entities=SPLIT.entities('calibration'),fit_role='calibration')
    with pytest.raises(ValueError):builder(_context(builder,spec=spec))


def test_all_missing_training_publishes_unavailable_and_blocks_dependents(tmp_path):
    builder=_builder(tmp_path,_features()[['missing_train']])
    spec=_spec(builder);assert spec.status=='unavailable'
    payload=builder(_context(builder))['operator.json']
    assert payload['status']=='unavailable' and payload['columns']==[]
    assert payload['entity_order']==list(IDS) and payload['matrix']==[[]]*len(IDS)
    child=replace(spec,status='ok',gaps=(),dependencies=spec.dependencies+(A.Dependency('artifact','operator','0'*64),))
    plan=P.Plan((P.Job('features',spec,split=SPLIT),P.Job('child',child,(P.Upstream('features','artifact','operator'),),SPLIT)),builder.snapshots)
    result=P.run(plan,{'features':builder},tmp_path/'run',BUDGET)
    assert dict(result.states)=={'features':'unavailable','child':'blocked'}
    assert result.artifacts['features'].payloads['operator.json']==payload
    incorrect=replace(spec,status='ok',gaps=())
    with pytest.raises(ValueError,match='availability'):builder(_context(builder,spec=incorrect))


def test_typed_graph_resumes_exact_operator_and_child_binding_without_refitting(tmp_path):
    builder=_builder(tmp_path);spec=_spec(builder)
    child=replace(spec,dependencies=spec.dependencies+(A.Dependency('artifact','operator','0'*64),))
    plan=P.Plan((P.Job('features',spec,split=SPLIT),P.Job('child',child,(P.Upstream('features','artifact','operator'),),SPLIT)),builder.snapshots)
    def consume(context):
        context.checkpoint()
        parent=context.parents['features'];parent.verify_contents()
        return {'reference.json':{'upstream_identity':parent.identity,'columns':parent.payloads['operator.json']['columns']}}
    result=P.run(plan,{'features':builder,'child':consume},tmp_path/'run',BUDGET)
    assert set(result.states.values())=={'succeeded'}
    assert result.artifacts['child'].payloads['reference.json']['upstream_identity']==result.artifacts['features'].identity
    replay=P.run(plan,{},tmp_path/'run',BUDGET)
    assert {key:value.identity for key,value in replay.artifacts.items()}=={key:value.identity for key,value in result.artifacts.items()}


def test_cooperative_checkpoint_stops_before_publishing_and_shape_byte_caps_are_checked(tmp_path):
    builder=_builder(tmp_path)
    def stop():raise P.Cancelled('Fixture stop')
    with pytest.raises(P.Cancelled):builder(_context(builder,checkpoint=stop))
    with pytest.raises(P.BudgetExceeded):_builder(tmp_path,limits=replace(LIMITS,max_cells=1),name='small.json')
    bounded=_builder(tmp_path,limits=replace(LIMITS,max_input_bytes=1),name='bytes.json')
    with pytest.raises(P.BudgetExceeded):bounded(_context(bounded))


def test_inf_and_non_numeric_data_and_signed_zero_snapshot_mismatches_are_not_silently_changed(tmp_path):
    features=_features();features.iloc[0,0]=np.inf
    with pytest.raises(ValueError,match='Infinite'):_builder(tmp_path,features)
    features=_features();features['text']='1.0'
    with pytest.raises(ValueError,match='numeric dtypes'):_builder(tmp_path,features)
    features=_features();features.iloc[0,0]=-0.
    builder=_builder(tmp_path,features);spec=_spec(builder);builder.features.iloc[0,0]=0.
    with pytest.raises(ValueError,match='exactly match'):builder(_context(builder,spec=spec))


@pytest.mark.parametrize('column',['n_publications_total','n_fulltext_total','publication_count','fulltext_count','citation_count'])
def test_literature_attention_aliases_are_not_features(tmp_path,column):
    features=_features();features[column]=1.
    builder=_builder(tmp_path,features)
    with pytest.raises(ValueError,match='attention'):builder(_context(builder))


def test_complex_input_is_refused_before_lossy_float_cast(tmp_path):
    features=_features();features['x']=features['x'].astype(complex)+1j
    with pytest.raises(ValueError,match='Complex'):F.snapshot_payload(features,organism=O.TOXOPLASMA)
    reference=_builder(tmp_path)
    with pytest.raises(ValueError,match='Complex'):
        F.NumericFeatureBuilder(features,reference.snapshot,split=SPLIT,exclusions=EXCLUSIONS,limits=LIMITS)


def test_support_guard_code_is_an_explicit_snapshot_and_dependency(tmp_path):
    builder=_builder(tmp_path);spec=_spec(builder)
    assert {s.dependency.name for s in builder.guard_code}=={
        'numeric_guard_strategies','numeric_guard_splits','numeric_guard_artifacts','numeric_guard_precompute_jobs'}
    changed=replace(spec,dependencies=tuple(d for d in spec.dependencies if d.name!='numeric_guard_strategies'))
    with pytest.raises(ValueError,match='source/code'):builder(_context(builder,spec=changed))


def test_nullable_numeric_input_preserves_observed_flags_and_exact_float_operator(tmp_path):
    features=_features().astype('Float64')
    builder=_builder(tmp_path,features)
    payload=builder(_context(builder))['operator.json']
    assert payload['observed'][0]==[True,False]
    native=(_features().iloc[:5][['x','y']].rank(pct=True)-.5).fillna(0).to_numpy()
    np.testing.assert_array_equal(np.array(payload['matrix'])[:5],native)
