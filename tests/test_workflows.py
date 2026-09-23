"""The guided UI keeps evidence, evaluation and measured imports distinct."""
import numpy as np
import pandas as pd

from starplast.workflows import WorkflowDialog
from starplast.prediction import TaskSpec, run


def nodes():
    rng=np.random.default_rng(12)
    frame=pd.DataFrame({'gene_id':[f'g{i}' for i in range(150)],
                        'orthogroup':[f'family{i//2}' for i in range(150)],
                        'product':['example protein']*150,
                        'x1':rng.normal(size=150),'x2':rng.normal(size=150),
                        'target':['A','B']*75})
    frame.loc[140:,'target']=None
    return frame


def test_explore_and_screen_comparison_keep_missing_data(qtbot,tmp_path):
    frame=nodes(); dialog=WorkflowDialog(frame); qtbot.addWidget(dialog)
    selected=[]; dialog.gene_selected.connect(selected.append)
    dialog.gene.setText('g12'); dialog._explore()
    assert selected[0]=='g12' and dialog.evidence.rowCount()>0
    path=tmp_path/'screen.csv'
    pd.DataFrame({'gene_id':['g1','absent','g2'],'effect':[0.,4.,np.nan]}).to_csv(path,index=False)
    dialog.load_screen(path); dialog._compare()
    assert len(dialog.comparison)==3
    assert dialog.comparison.mapping_status.tolist()==['mapped','unresolved','mapped']
    assert dialog.comparison.screen_effect.iloc[0]==0
    assert pd.isna(dialog.comparison.screen_effect.iloc[-1])
    assert len(dialog.screen_source['sha256'])==64


def test_background_prediction_finishes_with_exportable_result(qtbot):
    dialog=WorkflowDialog(nodes()); qtbot.addWidget(dialog)
    dialog.target.setCurrentText('target')
    dialog._run()
    assert not dialog.run_button.isEnabled()
    qtbot.waitUntil(lambda:dialog.run_button.isEnabled(),timeout=30000)
    assert dialog.result is not None
    assert dialog.export_button.isEnabled() and dialog.predictions.rowCount()==10
    assert 'Held-out evaluation' in dialog.summary.text()
    assert dialog.result.provenance['evaluation_mode']=='inductive_group_holdout'


def test_failure_clears_previous_result(qtbot):
    dialog=WorkflowDialog(nodes()); qtbot.addWidget(dialog)
    dialog.target.setCurrentText('target'); dialog._run()
    qtbot.waitUntil(lambda:dialog.run_button.isEnabled(),timeout=30000)
    dialog.target.setCurrentText('product'); dialog._run()
    qtbot.waitUntil(lambda:dialog.run_button.isEnabled(),timeout=30000)
    assert dialog.result is None and not dialog.export_button.isEnabled()
    assert dialog.predictions.rowCount()==0
    assert 'at least two observed classes' in dialog.summary.text()
