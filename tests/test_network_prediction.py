"""Weighted graphs preserve family holdouts and abstain for disconnected genes."""
import numpy as np
import pandas as pd
import pytest
from starplast.network_prediction import run_network
from starplast.prediction import TaskSpec


def fixture():
    n=100
    nodes=pd.DataFrame({'gene_id':[f'g{i}' for i in range(n)],'orthogroup':[f'o{i}' for i in range(n)],
                        'target':['A' if i%2 else 'B' for i in range(n)]})
    nodes.loc[90:,'target']=None
    # Two target-independent synthetic connected components and one isolated gene.
    a=np.arange(97); b=a+2
    graph={'gene_ids':nodes.gene_id.to_numpy(),'signal__a':a,'signal__b':b,'signal__w':np.ones(len(a))}
    return nodes,graph


def test_disconnected_gene_abstains_and_inner_groups_do_not_leak():
    nodes,graph=fixture();graph['signal__w'][-1]=0
    result=run_network(nodes,TaskSpec('target'),graph,['signal'],{'signal':['independent_assay']},log=lambda _:None)
    assert not result.predictions.iloc[-1].supported
    assert result.provenance['evaluation_mode']=='transductive_group_holdout'
    for fold in result.folds:
        keys=['training_ids','weight_validation_ids','calibration_ids','test_ids']
        for i,key in enumerate(keys):
            for other in keys[i+1:]: assert not set(fold[key])&set(fold[other])


def test_reordered_graph_is_refused():
    nodes,graph=fixture();graph['gene_ids']=graph['gene_ids'][::-1]
    with pytest.raises(ValueError,match='ordered gene IDs'):
        run_network(nodes,TaskSpec('target'),graph,['signal'],{'signal':[]})


def test_target_derived_or_unsourced_graph_is_refused():
    nodes,graph=fixture()
    with pytest.raises(ValueError,match='provenance'):
        run_network(nodes,TaskSpec('target'),graph,['signal'],{})
    with pytest.raises(ValueError,match='all network layers'):
        run_network(nodes,TaskSpec('target'),graph,['signal'],{'signal':['target']})


def test_bundled_graphs_have_explicit_sources_and_correct_identity():
    from starplast.network_prediction import bundled_network
    from starplast.paths import cache_file
    for prefix in ('','pf_'):
        nodes=pd.read_parquet(cache_file(prefix+'nodes.parquet'))
        graph,sources=bundled_network(nodes)
        np.testing.assert_array_equal(graph['gene_ids'],nodes.gene_id)
        assert set(sources)=={'domain','coexpression','struct'}
        assert all(sources.values())
