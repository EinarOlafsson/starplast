"""Replay frozen conformal outcomes and record a training-class-vocabulary control."""
from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def verify(pilot):
    """Run exact native parity and freeze the broad-set control in an executed notebook."""
    from notebook_runner import ExecutedNotebook
    pilot = Path(pilot)
    if (pilot/'verification.ipynb').exists() or (pilot/'all_class_control').exists():
        raise ValueError('Verification/control outputs already exist; preserve the snapshot')
    nb = ExecutedNotebook('Native conformal parity, model boundaries and broad-set control')
    nb.md('Check every input/output identity, then compare against the unchanged native conformal function on the identical frozen partitions. Its partition chooser is temporarily replaced for this diagnostic and restored; its inference/calibration code is unchanged. Rare classes remain in the full eligible benchmark population. This is prediction-grade recovery, not biological accuracy or an established exchangeability guarantee.')
    nb.code('import json, hashlib, ast', 'from pathlib import Path', 'from dataclasses import replace',
        'import numpy as np', 'import pandas as pd',
        'from starplast import artifacts as A, conformal_records as F, organisms as O, strategies as S, strategy_learning as N, record_scorecards as R, baselines as B',
        'from starplast.splits import read_split', 'from scripts.freeze_conformal_pilot import frozen_rank_matrix',
        f'out=Path({str(pilot)!r})', "m=json.loads((out/'manifest.json').read_text())",
        'sha=lambda path:hashlib.sha256(Path(path).read_bytes()).hexdigest()',
        "for path,digest in m['input_sha256'].items():\n    assert sha(path)==digest,path",
        "for path,digest in m['outputs'].items():\n    assert sha(out/path)==digest,path",
        "split,exclusions=read_split(out/'split.json')",
        "raw=json.loads((out/'held_out/manifest.json').read_text())", "spec=A.spec_from_dict(raw['spec'])",
        "artifact=A.read_artifact(out/'held_out',expected=spec,split=split)",
        "base_spec=A.spec_from_dict(json.loads((out/'base_model/manifest.json').read_text())['spec'])",
        "base=A.read_artifact(out/'base_model',expected=base_spec,split=split)",
        "state=base.payloads['model_state.json']", "ids=tuple(a.entity for a in split.assignments)",
        "frame=pd.read_parquet(O.nodes_path(O.TOXOPLASMA)).set_index('gene_id')",
        "matrix=frozen_rank_matrix(frame.loc[list(ids),state['input_columns']],state)",
        "where={entity:i for i,entity in enumerate(ids)}",
        "positions={role:np.array([where[e] for e in split.entities(role)]) for role in ('train','calibration','test')}",
        "np.testing.assert_array_equal(matrix[positions['train']],state['training_vectors'])",
        "visible=pd.Series([None]*len(ids),dtype=object)",
        "for role in ('train','calibration'):\n    visible.iloc[positions[role]]=frame.loc[list(split.entities(role)),'compartment'].astype(str).to_numpy()",
        "chooser=N._group_split", "N._group_split=lambda *args:(positions['train'],positions['calibration'])",
        "p={'model':'kNN','k':state['k'],'C':.5,'alpha':.1,'thresholds':'per class'}",
        "try:\n    native,size,sets,thresholds=N._conformal(S.Context(pd.DataFrame({'gene_id':ids}),graph={},organism=O.TOXOPLASMA),matrix,visible,positions['test'],p)\nfinally:\n    N._group_split=chooser",
        "rows=pd.DataFrame(artifact.payloads['rows.json'],columns=artifact.payloads['row_columns.json'])",
        "expected=native.iloc[positions['test']].astype(object).where(native.iloc[positions['test']].notna(),None).tolist()",
        "assert rows.prediction.astype(object).where(rows.prediction.notna(),None).tolist()==expected",
        "assert rows.native_set_text.tolist()==sets.iloc[positions['test']].tolist()",
        "assert rows.set_size.tolist()==size.iloc[positions['test']].tolist()",
        "score_payload=artifact.payloads['class_scores.json']",
        "scores=pd.DataFrame(score_payload['scores'],columns=score_payload['classes'])",
        "np.testing.assert_array_equal(scores,S.class_scores_of(native).iloc[positions['test']])",
        "calibration=artifact.payloads['calibration.json']",
        "for label,value in thresholds.items():\n    assert calibration['thresholds'][label]=={'status':'finite' if np.isfinite(value) else 'unbounded','value':value if np.isfinite(value) else None}",
        "assert hashlib.sha256(A.canonical_object(calibration).encode()).hexdigest()==m['summary']['calibration_identity']",
        "cal_scores=pd.DataFrame(calibration['scores']['scores'],index=calibration['scores']['entities'],columns=calibration['scores']['classes'])",
        "test_scores=scores.copy(); test_scores.index=rows.entity.tolist()",
        "cal_labels=pd.Series(calibration['calibration_labels'],index=calibration['calibration_entities'])",
        "batch=F.conformal_calls(cal_scores,test_scores,cal_labels,split=split,base_model=base,alpha=.1,per_class=True)",
        "assert batch.rows.astype(object).where(batch.rows.notna(),None).to_dict('records')==[{key:row[key] for key in batch.rows.columns} for row in artifact.payloads['rows.json']]",
        "scope_data=dict(artifact.payloads['card.json']['scope']); scope_data['gaps']=tuple(scope_data['gaps'])",
        "scope=R.RecordScope(**scope_data)", "normalize=lambda value:json.loads(json.dumps(value,allow_nan=False))",
        "assert normalize(R.aggregate(rows,scope,parameters={'class_scores':scores}))==artifact.payloads['card.json']",
        "assert normalize(R.class_cards(rows,scope,parameters={'class_scores':scores}))==artifact.payloads['class_cards.json']",
        "labels=pd.Series(state['training_labels'],index=state['fit_entities'])", "baselines=B.label_baselines(labels,split,seed=split.seed)",
        "for column in ('majority','prevalence_call'):\n    simple=rows.drop(columns=['prediction','correct','abstained','prediction_set','native_set_text','set_size']).assign(prediction=baselines[column].to_numpy())\n    assert normalize(R.aggregate(simple,scope))==artifact.payloads['baseline_cards.json'][column]",
        "for path in ('starplast/conformal_records.py','scripts/freeze_conformal_pilot.py','starplast/record_scorecards.py'):\n    ast.parse(Path(path).read_text(),feature_version=(3,10))",
        "{'native_calls_sets_scores_quantiles':'exact','test_rows':len(rows),'independent_biological_admission':False,'metrics':artifact.payloads['card.json']['metrics']}")
    nb.md('The all-training-classes control uses only the fitted training class vocabulary. Its broad sets need no calibration truth. Every test gene is retained, and unsupported truth classes would count as set misses. Compare set size/efficiency alongside coverage; a nearly full set is weak inference capacity. The control is not a calibrated deployment promise.')
    nb.code("classes=sorted(set(labels),key=str)", "assert classes==state['native_classes']",
        "control=rows.copy()", "members=classes",
        "control['prediction_set']=[list(members) for _ in range(len(control))]", "control['set_size']=len(members)",
        "control['native_set_text']=' | '.join(members)",
        "control['prediction']=pd.Series([members[0] if len(members)==1 else None]*len(control),dtype=object)",
        "control['abstained']=len(members)!=1", "control['correct']=pd.array([v==t if isinstance(v,str) else None for v,t in zip(control.prediction,control.truth)],dtype='boolean')",
        "settings=A.canonical_object({'control':'all_training_classes','selection':'complete frozen train vocabulary'})",
        "control_scope=replace(scope,settings=settings)",
        "control_card=R.aggregate(control,control_scope,parameters={'prediction_classes':classes})",
        f"verification_code=Path({str(Path(__file__).resolve())!r})",
        "dependencies=base_spec.dependencies+(A.Dependency('model','fixed_training_classes',base.identity),A.Dependency('code','control_verifier',sha(verification_code)))",
        "control_spec=replace(spec,role='scorecard',partition_id='control:all_training_classes:'+split.identity,dependencies=dependencies,settings_json=settings,confidence_kind='none',calibration_scope='',gaps=tuple(sorted(set((*spec.gaps,'All-training-classes control; no calibrated promise')))))",
        "payload={'rows.json':control.astype(object).where(control.notna(),None).to_dict('records'),'row_columns.json':control.columns.tolist(),'card.json':control_card,'class_cards.json':R.class_cards(control,control_scope,parameters={'prediction_classes':classes})}",
        "payload=normalize(payload)", "control_identity=A.write_artifact(out/'all_class_control',control_spec,payload,split=split)",
        "assert A.read_artifact(out/'all_class_control',expected=control_spec,split=split).payloads==payload",
        "(out/'code'/verification_code.name).write_bytes(verification_code.read_bytes())",
        "{'control_identity':control_identity,'control_metrics':control_card['metrics'],'class_vocabulary_source':'training only'}")
    nb.write(str(pilot/'verification.ipynb'))
    print('Native sets/calls/scores/quantiles and all-class control verified for 560 test genes')


def main():
    """Verify one immutable conformal pilot without overwriting existing evidence."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pilot', type=Path, required=True)
    args = parser.parse_args()
    verify(args.pilot)


if __name__ == '__main__':
    main()
