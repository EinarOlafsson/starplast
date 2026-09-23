"""Generic explanation and budget utilities preserve explicit assumptions."""
import numpy as np
import pandas as pd
import pytest

from starplast.prioritization import prioritize,compare_screen,explain_candidates
from starplast.prediction import TaskSpec,run


def test_budget_and_group_constraints_are_respected():
    candidates=pd.DataFrame({'gene_id':['a','b','c','d'],'score':[.9,.8,.7,.6],
                             'cost':[2.,2.,3.,1.],'benefit':[1.]*4,'family':['f1','f1','f2','f3'],
                             'calibration_status':['held_out_temperature_scaling']*4})
    selected=prioritize(candidates,4.,group_column='family')
    assert selected.cost.sum()<=4 and not selected.family.duplicated().any()
    assert selected.gene_id.tolist()==['d','a']
    assert selected.cost_assumption.eq('user supplied').all()


def test_prioritization_refuses_uncalibrated_scores_and_invalid_costs():
    frame=pd.DataFrame({'gene_id':['a'],'score':[.9],'calibration_status':['not_calibrated']})
    with pytest.raises(ValueError,match='calibrated'): prioritize(frame,10)
    frame.calibration_status='held_out_temperature_scaling';frame['cost']=-1.
    with pytest.raises(ValueError,match='positive'): prioritize(frame,10)


def test_comparison_refuses_implicit_duplicate_aggregation():
    nodes=pd.DataFrame({'gene_id':['g1','g2']})
    screen=pd.DataFrame({'gene_id':['g1','g1'],'effect':[1,2]})
    with pytest.raises(ValueError,match='duplicate'): compare_screen(nodes,screen,'effect')


def test_explanations_only_use_unknown_genes_and_fitting_neighbors():
    rng=np.random.default_rng(42);x=rng.normal(size=180)
    nodes=pd.DataFrame({'gene_id':[f'g{i}' for i in range(180)],'x':x,'y':rng.normal(size=180),
                        'orthogroup':[f'f{i//2}' for i in range(180)],'target':np.where(x>0,'A','B')})
    nodes.loc[170:,'target']=None
    result=run(nodes,TaskSpec('target'),log=lambda _:None)
    explanations=explain_candidates(result,nodes)
    assert explanations and {r['gene_id'] for r in explanations} <= set(nodes.gene_id.iloc[170:])
    for record in explanations:
        assert {n['gene_id'] for n in record['training_neighbors']} <= set(result.provenance['final_training_ids'])
        assert all(r['feature']!='target' for r in record['linear_contributions'])
