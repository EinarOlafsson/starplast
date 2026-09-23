"""Prediction evaluation never fits on held-out families, labels or distributions."""
from dataclasses import replace
import json
import numpy as np
import pandas as pd
import pytest

from starplast import prediction as P


def fixture(kind='classification', n=150):
    rng=np.random.default_rng(11)
    X=rng.normal(size=(n,4))
    data=pd.DataFrame(X,columns=['x1','x2','x3','x4'])
    data['gene_id']=[f'g{i}' for i in range(n)]
    data['orthogroup']=[f'f{i//3}' for i in range(n)]
    outcome=X[:,0]+rng.normal(size=n)
    data['target']=np.where(outcome>0,'A','B') if kind=='classification' else outcome
    data.loc[n-12:,'target']=None
    return data


def test_preprocessing_does_not_fit_the_test_distribution():
    training=pd.DataFrame({'x':[1.,2.,3.], 'unmeasured':[np.nan]*3})
    fitted=P.FoldMatrix().fit(training)
    before=fitted.median_.copy()
    transformed=fitted.transform(pd.DataFrame({'x':[1.e8], 'unmeasured':[100.]}))
    np.testing.assert_array_equal(fitted.median_,before)
    assert fitted.columns_==['x']
    assert np.isfinite(transformed).all()


def test_family_training_calibration_and_test_sets_are_disjoint():
    nodes=fixture(n=240)
    result=P.run(nodes,P.TaskSpec('target'),log=lambda _:None)
    groups=nodes.set_index('gene_id').orthogroup
    for fold in result.folds:
        sets=[set(groups.loc[fold[name]]) for name in ['training_ids','calibration_ids','test_ids']]
        assert not sets[0]&sets[1] and not sets[0]&sets[2] and not sets[1]&sets[2]
        assert 'target' not in fold['features']
        assert 'orthogroup' not in fold['features']
    assert (result.predictions.loc[nodes.target.notna(),'fold']>=0).all()
    assert (result.predictions.loc[nodes.target.isna(),'role']=='unlabelled_candidate').all()
    assert all(status.startswith('held_out') for status in result.predictions.calibration_status)


def test_all_missing_unknown_gene_abstains():
    nodes=fixture();nodes.loc[len(nodes)-1,['x1','x2','x3','x4']]=np.nan
    result=P.run(nodes,P.TaskSpec('target'),log=lambda _:None)
    row=result.predictions.iloc[-1]
    assert not row.supported and pd.isna(row.prediction)
    assert row.abstention_reason=='insufficient_measured_features'


def test_fixed_classes_do_not_remap_wrong_predictions():
    y=np.array([0,0,1,1]);wrong=1-y
    metrics,_=P._evaluate(y,wrong,np.eye(2)[wrong],np.ones(4,dtype=bool),np.ones(4,dtype=bool),['A','B'],'classification')
    assert metrics['accuracy']==0 and metrics['macro_f1']==0


def test_missing_groups_are_independent_singletons():
    n=fixture();n['orthogroup']=None
    groups=P.group_ids(n,'orthogroup')
    assert len(set(groups))==len(n)
    with pytest.raises(ValueError,match='absent'): P.group_ids(n,'missing_column')


def test_masked_factors_accept_missing_views_without_learning_from_test():
    X=np.array([[1,2,np.nan],[2,4,1],[3,np.nan,2],[4,8,3]],dtype=float)
    model=P.MaskedFactors(components=2).fit(X)
    basis=model.components_.copy()
    result=model.transform(np.array([[np.nan,np.nan,np.nan],[100,np.nan,3]]))
    np.testing.assert_array_equal(model.components_,basis)
    assert np.isfinite(result).all() and np.all(result[0]==0)


@pytest.mark.parametrize('method',['prior','linear','boosted','neighbors','pca','multiview'])
def test_regression_preserves_continuous_outcomes_and_reports_intervals(method):
    nodes=fixture('regression')
    result=P.run(nodes,P.TaskSpec('target',kind='regression',method=method,dimensions=3),log=lambda _:None)
    assert 'rmse' in result.metrics
    assert result.predictions.interval_lower.notna().all()
    assert (result.predictions.interval_lower<=result.predictions.interval_upper).all()
    assert result.predictions.prediction.dtype.kind=='f'


def test_multilabel_unknowns_do_not_become_negatives():
    n=fixture();n['trait1']=np.where(n.x2>0,1.,0.);n['trait2']=np.where(n.x3>0,1.,0.)
    n.loc[130:,'trait1']=np.nan;n.loc[140:,'trait2']=np.nan
    results=P.run_multilabel(n,['trait1','trait2'],method='prior')
    assert (results['trait1'].predictions.role=='unlabelled_candidate').sum()==20
    assert (results['trait2'].predictions.role=='unlabelled_candidate').sum()==10
    for result in results.values():
        assert all(not {'trait1','trait2'}&set(f['features']) for f in result.folds)


def test_result_has_roundtrippable_split_and_snapshot_provenance(tmp_path):
    result=P.run(fixture(),P.TaskSpec('target'),log=lambda _:None)
    result.save(tmp_path)
    record=json.loads((tmp_path/'run.json').read_text())
    assert len(record['provenance']['data_sha256'])==64
    assert record['provenance']['executed_method']=='linear'
    assert len(record['folds'])==3
