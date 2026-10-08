"""Freeze and independently replay native ridge on the existing HFF fitness candidate."""
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
sys.path.insert(0, str(ROOT))
from starplast import artifacts as A, baselines as B, capabilities as C, organisms as O, record_scorecards as R, scorecard as SC, strategies as S, strategy_catalog as N, value_records as V  # noqa: E402
from starplast.ground_truth import read_registry  # noqa: E402
from starplast.query import Query  # noqa: E402
from starplast.splits import make_split, make_exclusions, read_split, write_split  # noqa: E402


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _normal(value):
    return json.loads(json.dumps(value, allow_nan=False))


def freeze(output):
    """Fit without held-out truth and pin every source, split and outcome dependency."""
    output = Path(output)
    if output.exists():
        raise ValueError('Use a new immutable ridge-pilot directory')
    registry = ROOT/'results/ground_truth_registry_2026_10_07_v2/registry.json'
    entries, _ = read_registry(registry, relocate_snapshots_to=registry.parent)
    entry = next(e for e in entries if e.organism == O.TOXOPLASMA and e.target == 'fit_invitro_hff' and e.task == SC.T_VALUES)
    universe = json.loads(Path(entry.universe_file.path).read_text())
    source = Path(O.nodes_path(entry.organism))
    table = pd.read_parquet(source).set_index('gene_id', drop=False)
    if table.index.tolist() != universe:
        raise ValueError('Installed and frozen entity order differ')
    mask = entry.mask(universe)
    truth = pd.read_parquet(entry.truth_file.path).set_index('gene_id')[entry.target]
    eligible = mask.index[mask].tolist()
    np.testing.assert_array_equal(table.loc[eligible, entry.target], truth.loc[eligible])
    context = S.Context(table.reset_index(drop=True), graph={}, organism=entry.organism)
    groups = context.groups()[mask.to_numpy()].tolist()
    split = make_split(eligible, groups, organism=entry.organism, benchmark_id=entry.benchmark_id, seed=17, feature_access='inductive')
    train = table.loc[list(split.entities('train'))].reset_index(drop=True)
    train_context = S.Context(train, graph={}, organism=entry.organism)
    exclusions = make_exclusions(train, entry.target, benchmark_id=entry.benchmark_id, split=split)
    own_kind = train_context.same_kind(entry.target)
    exclusions = replace(exclusions, columns=tuple(sorted(set(exclusions.columns)|own_kind)))
    columns = train_context.numeric_columns(exclude=exclusions.columns)
    features = table.loc[eligible, columns]
    labels = truth.loc[list(split.entities('train'))]
    batch = V.trait_ridge(features, labels, split=split, exclusions=exclusions)
    rows = batch.rows.copy()
    # Evaluation truth first enters after feature/model fitting and prediction.
    rows['truth'] = truth.loc[rows.entity].to_numpy(dtype=float)
    rows['group'] = rows.entity.map(dict(zip(eligible, groups)))
    settings = {'model': 'ridge', 'alpha': 1., 'own_kind': 'leave out',
        'selection': 'fixed before evaluation', 'feature_selection': 'training observation and variance',
        'scaling': 'training average rank centered', 'imputation': 'native centered-rank zero'}
    unknown_groups = sum(str(g).startswith('__') for g in groups)
    gaps = tuple(sorted(set((*entry.gaps, *batch.gaps,
        'Stored fitness source requires original assay/mapping/unit review before independent biological admission',
        f'{unknown_groups} eligible genes have identity-only unresolved homology groups'))))
    scope = R.RecordScope(entry.organism, 'trait_regression', entry.target, SC.T_VALUES,
        A.canonical_object(settings), 17, split.identity, 'outer_test', entry.benchmark_id,
        entry.evidence_grade, 'gene', entry.negative_semantics, gaps)
    baseline = B.value_baselines(labels, split)
    cards = {}
    for name in ('mean', 'median'):
        cards[name] = R.aggregate(rows.assign(baseline_prediction=baseline[name].to_numpy()), scope,
            parameters={'quantity_unit': 'unresolved', 'baseline_name': 'training '+name})
    rows['baseline_prediction'] = baseline['mean'].to_numpy()
    code = [Path(__file__), ROOT/'scripts/freeze_conformal_pilot.py', *(ROOT/'starplast'/name for name in (
        'value_records.py','strategy_catalog.py','strategies.py','scorecard.py','record_scorecards.py',
        'baselines.py','splits.py','artifacts.py','capabilities.py','ground_truth.py','search.py'))]
    dependencies = [A.Dependency('code', p.stem, _sha(p)) for p in code]
    dependencies += [A.Dependency('table','installed',_sha(source)), A.Dependency('truth',entry.benchmark_id,entry.truth_file.sha256),
        A.Dependency('split','nested',split.identity),
        A.Dependency('exclusions','training',hashlib.sha256(A.canonical_object(asdict(exclusions)).encode()).hexdigest())]
    eval_scope = {'unit':'gene','eligible_population':len(rows),'truth_grade':entry.evidence_grade,
        'negative_semantics':entry.negative_semantics,'context':{'registered':entry.contexts},
        'limitations':list(gaps),'benchmark_status':entry.status,'source_ids':entry.source_ids,'quantity_unit':'unresolved'}
    version = subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()+'+ridge-adapter-worktree'
    spec = A.ArtifactSpec(Query(entry.organism,'trait',target=entry.target).to_json(), 'trait_regression',SC.T_VALUES,
        next(o for o in C.get('trait_regression').outputs if o.kind=='numeric_estimates'), entry.target,
        tuple(rows.entity),'held_out','outer:test:'+split.identity,tuple(dependencies),A.canonical_object(settings),
        A.canonical_object(eval_scope),17,version,fit_entities=split.entities('train'),fit_role='train',
        benchmark_id=entry.benchmark_id,evaluation_partition='test',gaps=gaps)
    payload = {'rows.json':rows.astype(object).where(rows.notna(),None).to_dict('records'),
        'row_columns.json':rows.columns.tolist(),'model_state.json':batch.model_state,
        'card.json':cards['mean'],'baseline_comparisons.json':cards,
        'baseline_state.json':{'fit_entities':labels.index.tolist(),'fit_role':'train','values':labels.tolist(),
            'mean':float(labels.mean()),'median':float(labels.median())},'own_kind_exclusions.json':sorted(own_kind)}
    payload = _normal(payload)
    output.mkdir(parents=True)
    write_split(output/'split.json',split,exclusions=exclusions)
    identity = A.write_artifact(output/'held_out',spec,payload,split=split)
    assert A.read_artifact(output/'held_out',expected=spec,split=split).payloads==payload
    summary = {'organism':entry.organism,'target':entry.target,'truth_grade':entry.evidence_grade,'benchmark_admitted':False,
        'quantity_unit':'unresolved','eligible':len(eligible),'train':len(labels),'tune':len(split.entities('tune')),
        'calibration':len(split.entities('calibration')),'test':len(rows),'training_selected_columns':len(columns),
        'kept_model_columns':len(batch.model_state['kept_columns']),'counts':cards['mean']['counts'],
        'metrics':cards['mean']['metrics'],'baseline_comparisons':{name:card['extra']['baseline_comparison'] for name,card in cards.items()},
        'artifact_identity':identity,'runtime_changes':'none','interpretation':'stored experimental-score recovery; biological admission unresolved'}
    (output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    inputs = [registry,source,Path(entry.truth_file.path),Path(entry.universe_file.path),*code]
    (output/'manifest.json').write_text(json.dumps({'summary':summary,'input_sha256':{str(p):_sha(p) for p in inputs},
        'outputs':{str(p.relative_to(output)):_sha(p) for p in output.rglob('*.json')}},indent=2)+'\n')
    (output/'code').mkdir()
    for path in code:
        (output/'code'/path.name).write_bytes(path.read_bytes())
    return summary


def verify(output):
    """Replay native training ranks, every prediction, coefficients and both baselines."""
    output = Path(output)
    manifest = json.loads((output/'manifest.json').read_text())
    for path, digest in manifest['input_sha256'].items():
        assert _sha(path)==digest,path
    for path, digest in manifest['outputs'].items():
        assert _sha(output/path)==digest,path
    split, exclusions = read_split(output/'split.json')
    spec = A.spec_from_dict(json.loads((output/'held_out/manifest.json').read_text())['spec'])
    artifact = A.read_artifact(output/'held_out',expected=spec,split=split)
    payload = artifact.payloads
    state = payload['model_state.json']
    ids = tuple(a.entity for a in split.assignments)
    table = pd.read_parquet(O.nodes_path(O.TOXOPLASMA)).set_index('gene_id',drop=False)
    training = table.loc[list(split.entities('train'))].reset_index(drop=True)
    ctx = S.Context(training,graph={},organism=O.TOXOPLASMA)
    assert ctx.numeric_columns(exclude=exclusions.columns)==state['input_columns']
    assert sorted(ctx.same_kind(spec.target))==payload['own_kind_exclusions.json']
    np.testing.assert_array_equal(ctx.matrix(state['kept_columns']),state['training_vectors'])
    labels = pd.Series(state['training_values'],index=state['fit_entities'])
    batch = V.trait_ridge(table.loc[list(ids),state['input_columns']],labels,split=split,exclusions=exclusions)
    assert batch.model_state==state
    rows = pd.DataFrame(payload['rows.json'],columns=payload['row_columns.json'])
    np.testing.assert_array_equal(batch.rows.prediction.to_numpy(dtype=float),rows.prediction.to_numpy(dtype=float))
    assert batch.rows.abstained.tolist()==rows.abstained.tolist()
    # Independent native prediction replay, using only the saved training distributions.
    from scripts.freeze_conformal_pilot import frozen_rank_matrix
    matrix = frozen_rank_matrix(table.loc[list(ids),state['input_columns']],state)
    where = {entity:i for i,entity in enumerate(ids)}
    visible = pd.Series(np.nan,index=range(len(ids)))
    visible.iloc[[where[e] for e in labels.index]] = labels.to_numpy()
    native = N._fit_predict(matrix,visible,'ridge')
    np.testing.assert_array_equal(rows.prediction.to_numpy(dtype=float),native[[where[e] for e in split.entities('test')]])
    np.testing.assert_array_equal(rows.truth,table.loc[list(rows.entity),spec.target])
    scope_data = dict(payload['card.json']['scope']);scope_data['gaps']=tuple(scope_data['gaps'])
    scope = R.RecordScope(**scope_data)
    baselines = B.value_baselines(labels,split)
    for name in ('mean','median'):
        card = R.aggregate(rows.assign(baseline_prediction=baselines[name].to_numpy()),scope,
            parameters={'quantity_unit':'unresolved','baseline_name':'training '+name})
        assert _normal(card)==payload['baseline_comparisons.json'][name]
    assert _normal(R.aggregate(rows,scope,parameters={'quantity_unit':'unresolved','baseline_name':'training mean'}))==payload['card.json']
    assert payload['baseline_state.json']['values']==labels.tolist()
    return {'native_training_matrix':'exact','native_predictions':'exact','retained_test_rows':len(rows),
        'training_only_model_and_baselines':'replayed','biological_admission':False,'artifact_identity':artifact.identity}


def main():
    """Execute one fixed fit and its separate verification as an annotated notebook."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    nb = ExecutedNotebook('Frozen native ridge numeric outcomes and matched-baseline cards')
    nb.md('64.11 first numeric adapter: existing HFF competitive-fitness candidate, seed 17, nested homology groups, training-only exclusions/detection/ranks/imputation/model, fixed native ridge alpha=1. Original measurement units and source/mapping admission gaps remain unresolved; no independent biology is admitted. Tune/calibration/test truth never fits the model or baselines.')
    nb.code('from scripts.freeze_ridge_value_pilot import freeze, verify',f'summary=freeze({str(args.out)!r})','summary')
    nb.md('Independently replay all hashes, exact native training ranks, native ridge predictions and matched training mean/median cards. Every held-out row remains present; no categorical correctness probability or calibrated interval is invented.')
    nb.code(f'verification=verify({str(args.out)!r})','verification')
    nb.write(str(args.out/'pilot.ipynb'))
    print(nb.ns['summary'])
    print(nb.ns['verification'])


if __name__=='__main__':
    main()
