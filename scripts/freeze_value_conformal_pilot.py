"""Execute and independently replay numeric intervals from the frozen ridge reference."""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from starplast import artifacts as A, baselines as B, capabilities as C, organisms as O, record_scorecards as R, scorecard as SC, strategies as S, strategy_catalog as N, strategy_learning as L, value_conformal_records as F  # noqa: E402
from starplast.splits import read_split, write_split  # noqa: E402
from scripts.freeze_conformal_pilot import frozen_rank_matrix  # noqa: E402

PRIOR = ROOT/'results/ridge_value_pilot_2026_10_08_v2'


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _normal(value):
    return json.loads(json.dumps(value,allow_nan=False))


def _reference():
    split,exclusions = read_split(PRIOR/'split.json')
    spec = A.spec_from_dict(json.loads((PRIOR/'held_out/manifest.json').read_text())['spec'])
    old = A.read_artifact(PRIOR/'held_out',expected=spec,split=split)
    table_path = Path(O.nodes_path(O.TOXOPLASMA))
    if next(d.sha256 for d in spec.dependencies if d.kind=='table')!=_sha(table_path):
        raise ValueError('Frozen model source is stale')
    state = old.payloads['model_state.json']
    ids = [a.entity for a in split.assignments]
    table = pd.read_parquet(table_path).set_index('gene_id',drop=False)
    if state['fit_entities']!=list(split.entities('train')):
        raise ValueError('Frozen ridge fit population differs')
    np.testing.assert_array_equal(table.loc[state['fit_entities'],spec.target],state['training_values'])
    matrix = frozen_rank_matrix(table.loc[ids,state['input_columns']],state)
    positions = {entity:i for i,entity in enumerate(ids)}
    at = {role:np.array([positions[e] for e in split.entities(role)]) for role in ('train','calibration','test')}
    np.testing.assert_array_equal(matrix[at['train']],state['training_vectors'])
    visible = pd.Series(np.nan,index=range(len(ids)))
    visible.iloc[at['train']] = state['training_values']
    predictions = N._fit_predict(matrix,visible,'ridge')
    prior_rows = pd.DataFrame(old.payloads['rows.json'],columns=old.payloads['row_columns.json'])
    np.testing.assert_array_equal(predictions[at['test']],prior_rows.prediction.to_numpy(dtype=float))
    return split,exclusions,old,table,state,matrix,at,visible,predictions


def freeze(output):
    """Pin the verified train-only model reference and calibrate on separate genes."""
    output = Path(output)
    if output.exists():
        raise ValueError('Use a new immutable numeric interval directory')
    split,exclusions,old,table,state,matrix,at,visible,predictions = _reference()
    old_scope = json.loads(old.spec.evaluation_scope_json)
    settings = {'model':'ridge','ridge_alpha':1.,'alpha':.1,'selection':'fixed before final evaluation',
        'preprocessing':'verified prior training ranks and zero imputation','threshold':'native absolute residual quantile'}
    code = [Path(__file__),ROOT/'scripts/freeze_conformal_pilot.py',*(ROOT/'starplast'/name for name in (
        'value_conformal_records.py','strategy_learning.py','strategy_catalog.py','scorecard.py','record_scorecards.py',
        'baselines.py','splits.py','artifacts.py','capabilities.py','organisms.py'))]
    dependencies = [A.Dependency('code',p.stem,_sha(p)) for p in code]
    dependencies += [d for d in old.spec.dependencies if d.kind in {'table','truth','split','exclusions'}]
    dependencies += [A.Dependency('artifact','verified_ridge_state',old.identity)]
    version = subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()+'+numeric-interval-worktree'
    model_scope = dict(old_scope,eligible_population=len(split.entities('train')))
    model_spec = replace(old.spec,entity_order=split.entities('train'),role='fitted_model',
        partition_id='train:'+split.identity,dependencies=tuple(dependencies),evaluation_partition='none',
        evaluation_scope_json=A.canonical_object(model_scope),code_version=version)
    # The full model state stays in its immutable upstream artifact; this typed
    # reference is verified/replayed here and is not a standalone deployment bundle.
    model_payload = {'model_state.json':{'fit_entities':state['fit_entities'],'model':state['model'],'alpha':state['alpha'],
        'upstream_artifact_identity':old.identity,'upstream_path':str(PRIOR.relative_to(ROOT))},
        'reference_scope.json':{'full_model_payload':'model_state.json in verified upstream artifact',
            'standalone_deployment_bundle':False}}
    output.mkdir(parents=True)
    write_split(output/'split.json',split,exclusions=exclusions)
    A.write_artifact(output/'base_model',model_spec,model_payload,split=split)
    base = A.read_artifact(output/'base_model',expected=model_spec,split=split)
    cal = pd.Series(predictions[at['calibration']],index=split.entities('calibration'))
    test = pd.Series(predictions[at['test']],index=split.entities('test'))
    cal_truth = table.loc[list(cal.index),old.spec.target].astype(float)
    batch = F.conformal_values(cal,test,cal_truth,split=split,base_model=base,alpha=.1)
    rows = batch.rows.copy()
    # Outer truth is joined only after the model and calibration state freeze.
    rows['truth'] = table.loc[list(rows.entity),old.spec.target].to_numpy(dtype=float)
    rows['group'] = rows.entity.map({a.entity:a.group for a in split.assignments})
    gaps = tuple(sorted(set((*old.spec.gaps,*batch.gaps))))
    scope = R.RecordScope(O.TOXOPLASMA,'conformal_values',old.spec.target,SC.T_VALUES,A.canonical_object(settings),
        split.seed,split.identity,'outer_test',split.benchmark_id,old_scope['truth_grade'],'gene',old_scope['negative_semantics'],gaps)
    labels = pd.Series(state['training_values'],index=state['fit_entities'])
    baselines = B.value_baselines(labels,split)
    parameters = {'quantity_unit':'unresolved','nominal_interval_coverage':.9}
    cards = {name:R.aggregate(rows.assign(baseline_prediction=baselines[name].to_numpy()),scope,
        parameters=dict(parameters,baseline_name='training '+name)) for name in ('mean','median')}
    rows['baseline_prediction'] = baselines['mean'].to_numpy()
    calibration = batch.model_state
    calibration_hash = hashlib.sha256(A.canonical_object(calibration).encode()).hexdigest()
    dependencies += [A.Dependency('model','base',base.identity),A.Dependency('calibration','absolute_residuals',calibration_hash)]
    spec = A.ArtifactSpec(old.spec.query_json,'conformal_values',SC.T_VALUES,
        next(o for o in C.get('conformal_values').outputs if o.kind=='numeric_intervals'),old.spec.target,
        tuple(rows.entity),'held_out','outer:test:'+split.identity,tuple(dependencies),A.canonical_object(settings),
        A.canonical_object(dict(old_scope,limitations=list(gaps))),split.seed,version,
        fit_entities=split.entities('train'),fit_role='train',benchmark_id=split.benchmark_id,
        confidence_kind='prediction_interval',calibration_scope=A.canonical_object({'split':split.identity,
            'calibration_entities':list(cal.index),'nominal_coverage':.9,'exchangeability':'unresolved','truth_grade':old_scope['truth_grade']}),
        evaluation_partition='test',gaps=gaps)
    payload = _normal({'rows.json':rows.astype(object).where(rows.notna(),None).to_dict('records'),
        'row_columns.json':rows.columns.tolist(),'calibration.json':calibration,'card.json':cards['mean'],
        'baseline_comparisons.json':cards,'base_reference.json':model_payload})
    identity = A.write_artifact(output/'held_out',spec,payload,split=split)
    assert A.read_artifact(output/'held_out',expected=spec,split=split).payloads==payload
    summary = {'strategy':'conformal_values','model':'ridge','target':spec.target,'truth_grade':old_scope['truth_grade'],
        'benchmark_admitted':False,'quantity_unit':'unresolved','train':len(labels),'calibration':len(cal),
        'retained_test_rows':len(rows),'half_width':calibration['half_width'],'metrics':cards['mean']['metrics'],
        'counts':cards['mean']['counts'],'interval_counts':cards['mean']['extra']['interval_counts'],
        'baseline_comparisons':{k:v['extra']['baseline_comparison'] for k,v in cards.items()},
        'artifact_identity':identity,'base_model_identity':base.identity,'upstream_identity':old.identity,
        'interpretation':'fixed stored-score interval recovery; source/units/independence/exchangeability unresolved'}
    (output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    (output/'code').mkdir()
    for path in code:
        (output/'code'/path.name).write_bytes(path.read_bytes())
    inputs = [Path(O.nodes_path(O.TOXOPLASMA)),PRIOR/'split.json',PRIOR/'held_out/manifest.json',
        PRIOR/'held_out/model_state.json',*code]
    (output/'manifest.json').write_text(json.dumps({'summary':summary,'input_sha256':{str(p):_sha(p) for p in inputs},
        'outputs':{str(p.relative_to(output)):_sha(p) for p in output.rglob('*.json')}},indent=2)+'\n')
    return summary


def verify(output):
    """Replay every interval against native fitting on identical frozen populations."""
    output = Path(output)
    manifest = json.loads((output/'manifest.json').read_text())
    for path,digest in manifest['input_sha256'].items():
        assert _sha(path)==digest,path
    for path,digest in manifest['outputs'].items():
        assert _sha(output/path)==digest,path
    split,exclusions,old,table,state,matrix,at,visible,predictions = _reference()
    saved_split,saved_exclusions = read_split(output/'split.json')
    assert (saved_split,saved_exclusions)==(split,exclusions)
    spec = A.spec_from_dict(json.loads((output/'held_out/manifest.json').read_text())['spec'])
    artifact = A.read_artifact(output/'held_out',expected=spec,split=split)
    model_spec = A.spec_from_dict(json.loads((output/'base_model/manifest.json').read_text())['spec'])
    base = A.read_artifact(output/'base_model',expected=model_spec,split=split)
    assert artifact.payloads['base_reference.json']==base.payloads
    assert base.payloads['model_state.json']['upstream_artifact_identity']==old.identity
    cal = pd.Series(predictions[at['calibration']],index=split.entities('calibration'))
    test = pd.Series(predictions[at['test']],index=split.entities('test'))
    cal_truth = table.loc[list(cal.index),spec.target].astype(float)
    batch = F.conformal_values(cal,test,cal_truth,split=split,base_model=base,alpha=.1)
    assert batch.model_state==artifact.payloads['calibration.json']
    cal_hash = hashlib.sha256(A.canonical_object(batch.model_state).encode()).hexdigest()
    assert cal_hash==next(d.sha256 for d in spec.dependencies if d.kind=='calibration')
    rows = pd.DataFrame(artifact.payloads['rows.json'],columns=artifact.payloads['row_columns.json'])
    for column in ('prediction','lower','upper'):
        np.testing.assert_array_equal(batch.rows[column].to_numpy(dtype=float),rows[column].to_numpy(dtype=float))
    assert batch.rows.interval_status.tolist()==rows.interval_status.tolist()
    assert batch.rows.abstained.tolist()==rows.abstained.tolist()
    np.testing.assert_array_equal(rows.truth,table.loc[list(rows.entity),spec.target])
    # Preserve the released native algorithm; replace only its partition chooser
    # while comparing the same already-frozen train/calibration/test roles.
    chooser = L._group_split
    visible.iloc[at['calibration']] = cal_truth.to_numpy()
    try:
        L._group_split = lambda *args:(at['train'],at['calibration'])
        context = S.Context(pd.DataFrame({'gene_id':[a.entity for a in split.assignments]}),graph={},organism=O.TOXOPLASMA)
        pred,lower,upper,half = L._intervals(context,matrix,visible,'ridge',.1)
    finally:
        L._group_split = chooser
    for column,native in (('prediction',pred),('lower',lower),('upper',upper)):
        np.testing.assert_array_equal(rows[column].to_numpy(dtype=float),native[at['test']])
    assert batch.model_state['half_width']=={'status':'finite','value':half}
    scope_data = dict(artifact.payloads['card.json']['scope']);scope_data['gaps']=tuple(scope_data['gaps'])
    scope = R.RecordScope(**scope_data)
    baselines = B.value_baselines(pd.Series(state['training_values'],index=state['fit_entities']),split)
    for name in ('mean','median'):
        card = R.aggregate(rows.assign(baseline_prediction=baselines[name].to_numpy()),scope,
            parameters={'quantity_unit':'unresolved','nominal_interval_coverage':.9,'baseline_name':'training '+name})
        assert _normal(card)==artifact.payloads['baseline_comparisons.json'][name]
    assert artifact.payloads['card.json']==artifact.payloads['baseline_comparisons.json']['mean']
    return {'native_points_bounds_quantile':'exact','calibration_and_cards':'exact replay',
        'retained_test_rows':len(rows),'source_and_code_hashes':'verified','biological_admission':False,'artifact_identity':artifact.identity}


def main():
    """Write an executed notebook of frozen calibration and independent native replay."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    nb = ExecutedNotebook('Frozen numeric ridge intervals: empirical coverage and capacity')
    nb.md('64.11 controlled interval partition: fixed native ridge alpha=1 and conformal alpha=0.1. Reuse verified train-only state and the same 4,015 training / 1,103 calibration / 1,101 test genes. Test truth does not fit the model or threshold. Original source, units, homology, independence and exchangeability gaps remain; nominal coverage is not per-gene certainty. Full model state remains referenced in the immutable upstream artifact; the reference is not a standalone deployment bundle.')
    nb.code('from scripts.freeze_value_conformal_pilot import freeze, verify',f'summary=freeze({str(args.out)!r})','summary')
    nb.md('Verify input/output hashes and replay every native prediction, interval bound and absolute-residual quantile on identical fixed populations. Recompute all interval and matched training mean/median cards. Unbounded or missing intervals must remain explicit; no independent biological benchmark is admitted.')
    nb.code(f'verification=verify({str(args.out)!r})','verification')
    nb.write(str(args.out/'pilot.ipynb'))
    print(nb.ns['summary'])
    print(nb.ns['verification'])


if __name__=='__main__':
    main()
