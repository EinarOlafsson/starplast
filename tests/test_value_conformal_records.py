"""Known-truth interval capacity, native parity and frozen calibration boundaries."""
from dataclasses import replace
import json

import numpy as np
import pandas as pd
import pytest

from starplast import artifacts as A, capabilities as C, organisms as O, record_scorecards as R, scorecard as SC, strategies as S, strategy_learning as N, strategy_catalog as L, value_conformal_records as F
from starplast.query import Query
from starplast.splits import Assignment, SplitManifest

IDS = tuple(f'TGME49_{i:06d}' for i in range(100001, 100036))
ROLES = ['train']*24 + ['tune']*3 + ['calibration']*3 + ['test']*5
SPLIT = SplitManifest(O.TOXOPLASMA, 'interval-fixture', 'entity', tuple(Assignment(e,e,r) for e,r in zip(IDS,ROLES)),17)


def model(tmp_path):
    deps = tuple(A.Dependency(kind,name,digest) for kind,name,digest in (
        ('code','fixture','a'*64),('table','fixture','a'*64),('exclusions','fixture','a'*64),('split','nested',SPLIT.identity)))
    scope = A.canonical_object({'unit':'gene','eligible_population':24,'truth_grade':'synthetic_control',
        'negative_semantics':'observed fixture values only','context':{},'limitations':['not biological']})
    spec = A.ArtifactSpec(Query(O.TOXOPLASMA,'trait',target='fixture').to_json(),'trait_regression',SC.T_VALUES,
        next(o for o in C.get('trait_regression').outputs if o.kind=='numeric_estimates'),'fixture',
        SPLIT.entities('train'),'fitted_model','train:'+SPLIT.identity,deps,A.canonical_object({'model':'ridge','alpha':1}),
        scope,17,'fixture',fit_entities=SPLIT.entities('train'),fit_role='train',benchmark_id=SPLIT.benchmark_id)
    A.write_artifact(tmp_path/'base',spec,{'model_state.json':{'fit_entities':list(SPLIT.entities('train'))}},split=SPLIT)
    return A.read_artifact(tmp_path/'base',expected=spec,split=SPLIT)


def predictions():
    matrix = np.column_stack([np.linspace(-.5,.5,35),np.sin(np.arange(35.))])
    values = pd.Series(np.cos(np.arange(35.)))
    visible = values.copy(); visible.iloc[24:] = np.nan
    native = L._fit_predict(matrix,visible,'ridge')
    return matrix,values,pd.Series(native[27:30],index=SPLIT.entities('calibration')),pd.Series(native[30:],index=SPLIT.entities('test'))


@pytest.mark.parametrize('alpha',[.1,.5])
def test_exact_native_predictions_quantile_and_intervals_on_frozen_roles(tmp_path,monkeypatch,alpha):
    matrix,values,cal,test = predictions()
    truth = pd.Series(values.iloc[27:30].to_numpy(),index=cal.index)
    batch = F.conformal_values(cal,test,truth,split=SPLIT,base_model=model(tmp_path),alpha=alpha)
    monkeypatch.setattr(N,'_group_split',lambda *args:(np.arange(24),np.arange(27,30)))
    context = S.Context(pd.DataFrame({'gene_id':IDS}),graph={},organism=O.TOXOPLASMA)
    visible = values.copy(); visible.iloc[24:27] = np.nan; visible.iloc[30:] = np.nan
    native,lower,upper,half = N._intervals(context,matrix,visible,'ridge',alpha)
    np.testing.assert_array_equal(batch.rows.prediction,native[30:])
    status = 'finite' if np.isfinite(half) else 'unbounded'
    assert batch.rows.interval_status.eq(status).all()
    assert batch.model_state['half_width']=={'status':status,'value':half if np.isfinite(half) else None}
    if status=='finite':
        np.testing.assert_array_equal(batch.rows.lower,lower[30:])
        np.testing.assert_array_equal(batch.rows.upper,upper[30:])
    else:
        assert batch.rows[['lower','upper']].isna().all().all()
    json.dumps(batch.model_state,allow_nan=False)


def test_missing_calibration_prediction_cannot_silently_reduce_population(tmp_path):
    _,values,cal,test = predictions()
    cal.iloc[0] = np.nan; test.iloc[-1] = np.nan
    batch = F.conformal_values(cal,test,pd.Series(values.iloc[27:30].to_numpy(),index=cal.index),
        split=SPLIT,base_model=model(tmp_path),alpha=.5)
    assert batch.rows.entity.tolist()==list(SPLIT.entities('test'))
    assert batch.rows.interval_status.eq('unavailable').all()
    assert batch.rows.abstained.tolist()==[False,False,False,False,True]
    assert batch.model_state['missing_calibration_predictions']==1
    assert len(batch.model_state['absolute_residuals'])==3
    assert batch.model_state['half_width']=={'status':'unavailable','value':None}
    assert batch.rows.calibrated_confidence.isna().all()


def test_model_and_cohort_guards_refuse_reordering_contamination_and_mutation(tmp_path):
    _,values,cal,test = predictions()
    truth = pd.Series(values.iloc[27:30].to_numpy(),index=cal.index)
    base = model(tmp_path)
    for args in ((cal.iloc[::-1],test,truth),(cal,test.iloc[::-1],truth),(cal,test,truth.iloc[::-1])):
        with pytest.raises(ValueError,match='complete ordered'):
            F.conformal_values(*args,split=SPLIT,base_model=base)
    with pytest.raises(ValueError,match='fitted-model artifact'):
        F.conformal_values(cal,test,truth,split=SPLIT,base_model=replace(base,spec=replace(base.spec,role='reusable_base')))
    base.payloads['model_state.json']['fit_entities'] = list(test.index)
    with pytest.raises(ValueError,match='changed after loading'):
        F.conformal_values(cal,test,truth,split=SPLIT,base_model=base)


@pytest.mark.parametrize('alpha',[True,0,1,np.nan,np.inf])
def test_invalid_alpha_is_refused(tmp_path,alpha):
    _,values,cal,test = predictions()
    with pytest.raises(ValueError,match='alpha'):
        F.conformal_values(cal,test,pd.Series(values.iloc[27:30].to_numpy(),index=cal.index),split=SPLIT,base_model=model(tmp_path),alpha=alpha)


def test_unknown_truth_infinite_inputs_and_residual_overflow_are_refused(tmp_path):
    _,values,cal,test = predictions()
    truth = pd.Series(values.iloc[27:30].to_numpy(),index=cal.index)
    base = model(tmp_path)
    unknown = truth.copy(); unknown.iloc[0] = np.nan
    with pytest.raises(ValueError,match='observed finite'):
        F.conformal_values(cal,test,unknown,split=SPLIT,base_model=base)
    with pytest.raises(ValueError,match='Infinite'):
        F.conformal_values(cal,test*0+np.inf,truth,split=SPLIT,base_model=base)
    with pytest.raises(ValueError,match='residuals overflowed'):
        F.conformal_values(cal*0+1e308,test,truth*0-1e308,split=SPLIT,base_model=base)


def test_known_truth_coverage_width_and_baselines_use_declared_denominators():
    rows = pd.DataFrame({'entity':['a','b','c','d'],'truth':[1.,3.,4.,100.],
        'prediction':[1.,2.,4.,np.nan],'lower':[0.,1.,np.nan,np.nan],'upper':[2.,2.5,np.nan,np.nan],
        'interval_status':['finite','finite','unbounded','unavailable'],'baseline_prediction':[2.]*4})
    metrics = SC.value_intervals(rows.prediction,rows.truth,rows.lower,rows.upper,rows.interval_status)
    assert metrics['interval_coverage']==.5 and metrics['interval_coverage_answered']==pytest.approx(2/3)
    assert metrics['interval_availability']==.75 and metrics['unbounded_interval_share']==.25
    assert metrics['finite_mean_interval_width']==1.75 and np.isnan(metrics['mean_interval_width'])
    assert np.isnan(metrics['interval_width_in_sd'])
    scope = R.RecordScope(O.TOXOPLASMA,'conformal_values','fixture',SC.T_VALUES,'fixed',17,SPLIT.identity,
        'outer_test','fixture','synthetic_control','gene','observed fixture only')
    card = R.aggregate(rows,scope,parameters={'quantity_unit':'fixture unit','nominal_interval_coverage':.9,'baseline_name':'training mean'})
    assert card['extra']['interval_counts']=={'finite':2,'unbounded':1,'unavailable':1}
    assert card['metrics']['promised_coverage']==.9 and card['metrics']['mean_interval_width'] is None
    assert card['extra']['baseline_comparison']['matched_answered_rows']==3
    assert card['extra']['baseline_comparison']['mae_skill']==.75
    assert card['extra']['quantity_unit']=='fixture unit'
    json.dumps(card,allow_nan=False)


def test_finite_and_empty_capacity_and_constant_truth_width_normalization():
    metric = SC.value_intervals([1.,3.],[1.,3.],[0.,2.],[2.,4.],['finite','finite'])
    assert metric['mean_interval_width']==2. and metric['interval_width_in_sd']==2.
    constant = SC.value_intervals([1.,1.],[1.,1.],[0.,0.],[2.,2.],['finite','finite'])
    assert constant['mean_interval_width']==2. and np.isnan(constant['interval_width_in_sd'])
    empty = SC.value_intervals([np.nan],[1.],[np.nan],[np.nan],['unavailable'])
    assert empty['interval_coverage']==0. and empty['interval_availability']==0.
    assert np.isnan(empty['interval_coverage_answered']) and np.isnan(empty['mean_interval_width'])


@pytest.mark.parametrize('lower,upper,status',[
    ([2.],[1.],['finite']),([np.nan],[1.],['finite']),([0.],[1.],['unbounded']),
    ([np.inf],[np.inf],['unbounded']),([np.nan],[np.nan],['wrong']),([-1e308],[1e308],['finite'])])
def test_contradictory_or_nonfinite_bounds_are_refused(lower,upper,status):
    with pytest.raises(ValueError):
        SC.value_intervals([0.],[1.],lower,upper,status)
