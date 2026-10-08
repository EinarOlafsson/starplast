"""Execute native kNN/logistic conformal sets on frozen calibration/test cohorts."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from starplast import artifacts as A, baselines as B, capabilities as C, conformal_records as F, organisms as O, record_scorecards as R, scorecard as SC, strategies as S, strategy_catalog as N  # noqa: E402
from starplast.ground_truth import read_registry  # noqa: E402
from starplast.query import Query  # noqa: E402
from starplast.splits import read_split, write_split  # noqa: E402


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _records(frame):
    return frame.astype(object).where(frame.notna(), None).to_dict('records')


def _scores(frame):
    return {'entities': frame.index.tolist(), 'classes': frame.columns.tolist(),
        'scores': frame.astype(object).where(frame.notna(), None).to_numpy().tolist()}


def frozen_rank_matrix(features, state):
    """Apply stored training distributions without fitting on held-out features."""
    matrix = np.zeros((len(features), len(state['kept_columns'])))
    for j, column in enumerate(state['kept_columns']):
        ordered = np.asarray(state['training_distributions'][column], dtype=float)
        values = features[column].to_numpy(dtype=float)
        if not len(ordered) or not np.isfinite(ordered).all() or np.isinf(values).any() or not np.array_equal(ordered, np.sort(ordered)):
            raise ValueError('Invalid frozen training rank distribution/features')
        known = np.isfinite(values)
        left = np.searchsorted(ordered, values[known], side='left')
        right = np.searchsorted(ordered, values[known], side='right')
        matrix[known, j] = np.where(right > left, (left+right+1)/2, right)/len(ordered)-.5
    return matrix


def freeze(output, model='kNN'):
    """Reuse fixed training state, calibrate separately, and retain every test set."""
    output = Path(output)
    if output.exists():
        raise ValueError('Use a new conformal-pilot directory')
    if model not in {'kNN', 'logistic'}:
        raise ValueError('Use an explicitly supported native base model')
    prior = ROOT/'results/label_knn_pilot_2026_10_08_v5'
    split, exclusions = read_split(prior/'split.json')
    old_spec = A.spec_from_dict(json.loads((prior/'held_out/manifest.json').read_text())['spec'])
    old = A.read_artifact(prior/'held_out', expected=old_spec, split=split)
    state = old.payloads['model_state.json']
    registry = ROOT/'results/ground_truth_registry_2026_10_07_v2/registry.json'
    entries, _ = read_registry(registry, relocate_snapshots_to=registry.parent)
    entry = next(e for e in entries if e.benchmark_id == split.benchmark_id)
    source = Path(O.nodes_path(entry.organism))
    if next(d.sha256 for d in old_spec.dependencies if d.kind == 'table') != _sha(source):
        raise ValueError('Prior training state has a stale installed source')
    universe = json.loads(Path(entry.universe_file.path).read_text())
    table = pd.read_parquet(source).set_index('gene_id', drop=False)
    if table.index.tolist() != universe:
        raise ValueError('Frozen and installed entity order differ')
    truth = pd.read_parquet(entry.truth_file.path).set_index('gene_id')[entry.target]
    labels = truth.loc[list(split.entities('train'))].astype(str)
    if labels.tolist() != state['training_labels'] or state['fit_entities'] != list(split.entities('train')):
        raise ValueError('Frozen training labels differ')
    ids = tuple(a.entity for a in split.assignments)
    features = table.loc[list(ids), state['input_columns']]
    exclusions.guard_inputs(columns=features.columns, benchmark_id=entry.benchmark_id)
    matrix = frozen_rank_matrix(features, state)
    where = {entity: i for i, entity in enumerate(ids)}
    positions = {role: [where[e] for e in split.entities(role)] for role in ('train', 'calibration', 'test')}
    np.testing.assert_array_equal(matrix[positions['train']], state['training_vectors'])
    visible = pd.Series([None]*len(ids), dtype=object)
    visible.iloc[positions['train']] = labels.to_numpy()
    query_positions = np.union1d(positions['calibration'], positions['test'])
    base_strategy = 'feature_knn' if model == 'kNN' else 'supervised_classifier'
    if model == 'kNN':
        native, _ = S.knn_vote(matrix, visible, state['k'], query=query_positions)
        base_settings = old_spec.settings_json
    else:
        split.guard_fit('model_fit', labels.index)
        native, _, fitted = N._logistic(matrix, visible, query_positions, C=.5)
        if fitted is None:
            raise ValueError('Frozen logistic pilot has no supported native model')
        state = dict(state, strategy=base_strategy, C=.5, max_iter=500, class_weight='balanced',
            coefficients=fitted.coef_.tolist(), intercept=fitted.intercept_.tolist(),
            iterations=fitted.n_iter_.tolist(), native_classes=fitted.classes_.tolist())
        base_settings = A.canonical_object({'C': .5, 'model': 'native balanced logistic', 'selection': 'fixed'})
    scores = S.class_scores_of(native)
    scores.index = list(ids)
    cal_scores, test_scores = (scores.loc[list(split.entities(role))] for role in ('calibration', 'test'))
    if model == 'kNN':
        np.testing.assert_array_equal(test_scores.to_numpy(), np.asarray(old.payloads['class_scores.json']['scores']))
    code = [Path(__file__), *(ROOT/'starplast'/name for name in ('conformal_records.py', 'strategy_learning.py',
        'strategies.py', 'strategy_catalog.py', 'record_scorecards.py', 'scorecard.py', 'baselines.py', 'artifacts.py', 'splits.py', 'capabilities.py', 'ground_truth.py'))]
    dependencies = [A.Dependency('code', p.stem, _sha(p)) for p in code]
    dependencies += [A.Dependency('table', 'installed', _sha(source)), A.Dependency('truth', entry.benchmark_id, entry.truth_file.sha256),
        A.Dependency('split', 'nested', split.identity), A.Dependency('exclusions', 'training', hashlib.sha256(A.canonical_object(asdict(exclusions)).encode()).hexdigest()),
        A.Dependency('artifact', 'fixed_prior_training_state', old.identity)]
    version = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()+'+conformal-adapter-worktree'
    query = Query(entry.organism, 'label', target=entry.target)
    model_scope = dict(json.loads(old_spec.evaluation_scope_json), eligible_population=len(labels))
    model_spec = A.ArtifactSpec(query.to_json(), base_strategy, SC.T_LABEL, C.get(base_strategy).outputs[0], entry.target,
        split.entities('train'), 'fitted_model', 'train:'+split.identity, tuple(dependencies), base_settings,
        A.canonical_object(model_scope), split.seed, version, fit_entities=split.entities('train'), fit_role='train',
        benchmark_id=entry.benchmark_id, confidence_kind='method_support', gaps=('Fixed training state; no selection from outer outcomes',))
    output.mkdir(parents=True)
    write_split(output/'split.json', split, exclusions=exclusions)
    A.write_artifact(output/'base_model', model_spec, {'model_state.json': state}, split=split)
    base = A.read_artifact(output/'base_model', expected=model_spec, split=split)
    cal_labels = truth.loc[list(split.entities('calibration'))].astype(str)
    batch = F.conformal_calls(cal_scores, test_scores, cal_labels, split=split, base_model=base, alpha=.1, per_class=True)
    rows = batch.rows.copy()
    # Hidden truth enters after model and threshold fitting have finished.
    rows['truth'] = truth.loc[rows.entity].astype(str).to_numpy()
    rows['correct'] = pd.array([p == t if isinstance(p, str) else None for p, t in zip(rows.prediction, rows.truth)], dtype='boolean')
    rows['group'] = rows.entity.map({a.entity: a.group for a in split.assignments})
    gaps = tuple(sorted(set((*entry.gaps, *batch.gaps, 'Stored predictions are surrogate truth; independent biological accuracy is unavailable',
        'Single frozen base/target pilot; full outer coverage remains open'))))
    settings = {'model': model, 'k': state['k'], 'C': .5, 'alpha': .1, 'thresholds': 'per class', 'selection': 'fixed without outer-test tuning'}
    scope = R.RecordScope(entry.organism, 'conformal_calls', entry.target, SC.T_LABEL, A.canonical_object(settings), split.seed,
        split.identity, 'outer_test', entry.benchmark_id, entry.evidence_grade, 'gene', entry.negative_semantics, gaps)
    params = {'class_scores': test_scores.reset_index(drop=True)}
    card, classes = R.aggregate(rows, scope, parameters=params), R.class_cards(rows, scope, parameters=params)
    baselines = B.label_baselines(labels, split, seed=split.seed)
    baseline_cards = {}
    for column in ('majority', 'prevalence_call'):
        simple = rows.drop(columns=['prediction', 'correct', 'abstained', 'prediction_set', 'native_set_text', 'set_size']).assign(prediction=baselines[column].to_numpy())
        baseline_cards[column] = R.aggregate(simple, scope)
    calibration = dict(batch.model_state, scores=_scores(cal_scores))
    calibration_hash = hashlib.sha256(A.canonical_object(calibration).encode()).hexdigest()
    dependencies += [A.Dependency('model', 'base', base.identity), A.Dependency('calibration', 'nested_quantiles', calibration_hash)]
    eval_scope = dict(model_scope, eligible_population=len(rows), limitations=list(gaps))
    spec = A.ArtifactSpec(query.to_json(), 'conformal_calls', SC.T_LABEL,
        next(o for o in C.get('conformal_calls').outputs if o.kind == 'label_sets'), entry.target, tuple(rows.entity), 'held_out',
        'outer:test:'+split.identity, tuple(dependencies), A.canonical_object(settings), A.canonical_object(eval_scope), split.seed, version,
        fit_entities=split.entities('train'), fit_role='train', benchmark_id=entry.benchmark_id, confidence_kind='prediction_set',
        calibration_scope=A.canonical_object({'split': split.identity, 'calibration_entities': cal_labels.index.tolist(),
            'nominal_coverage': .9, 'exchangeability': 'unresolved', 'truth_grade': entry.evidence_grade}), evaluation_partition='test', gaps=gaps)
    payloads = {'rows.json': _records(rows), 'row_columns.json': rows.columns.tolist(), 'class_scores.json': _scores(test_scores),
        'calibration.json': calibration, 'card.json': card, 'class_cards.json': classes, 'baseline_cards.json': baseline_cards}
    payloads = json.loads(json.dumps(payloads, allow_nan=False))
    identity = A.write_artifact(output/'held_out', spec, payloads, split=split)
    assert A.read_artifact(output/'held_out', expected=spec, split=split).payloads == payloads
    summary = {'strategy': 'conformal_calls', 'model': model, 'truth_grade': entry.evidence_grade, 'benchmark_admitted': False,
        'train': len(labels), 'calibration': len(cal_labels), 'retained_test_rows': len(rows), 'native_class_columns': len(test_scores.columns),
        'rare_class_overall_fallback': batch.model_state['rare_class_overall_fallback'],
        'unsupported_calibration_classes': batch.model_state['unsupported_calibration_classes'],
        'unsupported_test_classes': sorted(set(rows.truth)-set(test_scores.columns)), 'counts': card['counts'], 'metrics': card['metrics'],
        'nominal_set_coverage': .9, 'baseline_metrics': baseline_cards['majority']['metrics'], 'base_model_identity': base.identity,
        'calibration_identity': calibration_hash, 'artifact_identity': identity, 'interpretation': 'prediction-grade empirical set recovery', 'runtime_changes': 'none'}
    (output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    inputs = [registry, prior/'split.json', prior/'held_out/manifest.json', source, Path(entry.truth_file.path), Path(entry.universe_file.path), *code]
    (output/'manifest.json').write_text(json.dumps({'summary': summary, 'input_sha256': {str(p): _sha(p) for p in inputs},
        'outputs': {str(p.relative_to(output)): _sha(p) for p in output.rglob('*.json')}}, indent=2)+'\n')
    (output/'code').mkdir()
    for path in code:
        (output/'code'/path.name).write_bytes(path.read_bytes())
    return summary


def main():
    """Execute the fixed pilot and save its annotated notebook with actual outputs."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--model', choices=('kNN', 'logistic'), default='kNN')
    args = parser.parse_args()
    nb = ExecutedNotebook('Frozen native '+args.model+' conformal sets and outer-test scorecards')
    nb.md('Fixed seed-17 Toxoplasma compartment candidate: prior training rank state, separate 569-gene calibration and 560-gene test cohorts. Alpha=0.1/per-class thresholds are fixed without outer-outcome selection. Native quantiles and rare-class overall fallback are preserved. Stored labels have prediction grade; exchangeability and independent biological accuracy remain unresolved.')
    nb.code('from scripts.freeze_conformal_pilot import freeze', f'summary=freeze({str(args.out)!r}, model={args.model!r})', 'summary')
    nb.md('Every test set and native score row is retained. Empty and multi-label sets abstain from singleton calling. Scorecards separate both call denominators from empirical set coverage, size, singleton/empty shares and native efficiency. Unbounded thresholds have explicit status/null value. Full-outer/independent-truth work stays open. Installed data, runtime strategies and calibration are unchanged.')
    nb.write(str(args.out/'pilot.ipynb'))
    print(nb.ns['summary'])


if __name__ == '__main__':
    main()
