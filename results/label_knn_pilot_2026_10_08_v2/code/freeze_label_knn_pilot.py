"""Freeze a training-only native-kNN label-record pilot on stored localization calls."""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from starplast import artifacts as A, baselines as B, capabilities as C, label_records as L, organisms as O, record_scorecards as R, scorecard as SC, strategies as S  # noqa: E402
from starplast.ground_truth import read_registry  # noqa: E402
from starplast.query import Query  # noqa: E402
from starplast.splits import make_exclusions, make_split, write_split  # noqa: E402


def _records(frame):
    return json.loads(frame.to_json(orient='records'))


def freeze(output):
    """Fit the frozen pilot without admitting surrogate labels as independent biology."""
    output = Path(output)
    if output.exists():
        raise ValueError('Use a new label-pilot directory')
    registry = ROOT / 'results/ground_truth_registry_2026_10_07_v2/registry.json'
    entries, _ = read_registry(registry, relocate_snapshots_to=registry.parent)
    entry = next(e for e in entries if e.organism == O.TOXOPLASMA and e.target == 'compartment' and e.task == SC.T_LABEL)
    universe = json.loads(Path(entry.universe_file.path).read_text())
    mask = entry.mask(universe)
    source = Path(O.nodes_path(O.TOXOPLASMA))
    frame = pd.read_parquet(source).set_index('gene_id', drop=False)
    if frame.index.tolist() != universe:
        raise ValueError('Installed source and frozen target universe differ')
    truth = pd.read_parquet(entry.truth_file.path).set_index('gene_id')[entry.target]
    eligible = list(mask.index[mask])
    context = S.Context(frame.reset_index(drop=True), graph={}, organism=O.TOXOPLASMA)
    groups = context.groups()[mask.to_numpy()].tolist()
    split = make_split(eligible, groups, organism=O.TOXOPLASMA, benchmark_id=entry.benchmark_id, seed=17, feature_access='inductive')
    training = frame.loc[list(split.entities('train'))].reset_index(drop=True)
    exclusions = make_exclusions(training, entry.target, benchmark_id=entry.benchmark_id, split=split)
    # The column schema is known, but all detection/variance selection uses train rows only.
    train_context = S.Context(training, graph={}, organism=O.TOXOPLASMA)
    columns = train_context.numeric_columns(exclude=exclusions.columns)
    features = frame.loc[eligible, columns]
    training_labels = truth.loc[list(split.entities('train'))].astype(str)
    batch = L.feature_knn(features, training_labels, split=split, exclusions=exclusions, k=15, min_share=.3)
    # Hidden truth is joined only after fitting and prediction finish.
    rows = batch.rows.copy()
    rows['truth'] = truth.loc[rows.entity].astype(str).to_numpy()
    rows['correct'] = pd.array([p == t if isinstance(p, str) else None for p, t in zip(rows.prediction, rows.truth)], dtype='boolean')
    group_index = dict(zip(eligible, groups))
    rows['group'] = rows.entity.map(group_index)
    scores = batch.class_scores.reset_index(drop=True)
    settings = {'k': 15, 'min_share': .3, 'imputation': 'training_median', 'scaling': 'training_mean_std', 'feature_selection': 'training_observation_and_variance'}
    scope = R.RecordScope(O.TOXOPLASMA, 'feature_knn', entry.target, SC.T_LABEL, A.canonical_object(settings), 17,
        split.identity, 'outer_test', entry.benchmark_id, entry.evidence_grade, 'gene', entry.negative_semantics, entry.gaps)
    card = R.aggregate(rows, scope, parameters={'class_scores': scores})
    classes = R.class_cards(rows, scope)
    baseline = B.label_baselines(training_labels, split, seed=17)
    baseline_cards = {}
    for column in ('majority', 'prevalence_call'):
        baseline_rows = rows.drop(columns=['prediction', 'correct', 'abstained']).assign(prediction=baseline[column].to_numpy())
        baseline_cards[column] = R.aggregate(baseline_rows, scope)
    unsupported = rows[~rows.truth.isin(batch.class_scores.columns)].truth.value_counts().to_dict()
    unknown_groups = sum(str(g).startswith('__') for g in groups)
    gaps = tuple(sorted(set((*entry.gaps, *batch.gaps,
        'Candidate target is stored localization prediction; this is surrogate recovery, not independent biological accuracy',
        f'{unknown_groups} eligible genes lack a resolved homology group; identity-only groups retained'))))
    output.mkdir(parents=True)
    write_split(output / 'split.json', split, exclusions=exclusions)
    inputs = [Path(__file__), registry, Path(entry.truth_file.path), Path(entry.universe_file.path), source,
        *(ROOT / 'starplast' / name for name in ('label_records.py', 'artifacts.py', 'record_scorecards.py', 'splits.py', 'strategies.py', 'baselines.py', 'capabilities.py'))]
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    exclusion_hash = hashlib.sha256(A.canonical_object(asdict(exclusions)).encode()).hexdigest()
    dependencies = tuple(A.Dependency(kind, name, digest) for kind, name, digest in (
        ('table', 'installed_gene_universe', sha(source)), ('truth', entry.benchmark_id, entry.truth_file.sha256),
        ('code', 'label_adapter', sha(ROOT / 'starplast/label_records.py')), ('code', 'native_vote', sha(ROOT / 'starplast/strategies.py')),
        ('split', 'nested', split.identity), ('exclusions', 'training_source_closure', exclusion_hash)))
    query = Query(O.TOXOPLASMA, 'label', target=entry.target)
    eval_scope = {'unit': 'gene', 'eligible_population': len(rows), 'truth_grade': entry.evidence_grade,
        'negative_semantics': entry.negative_semantics, 'context': {'registered': entry.contexts},
        'limitations': list(gaps), 'benchmark_status': entry.status, 'source_ids': entry.source_ids}
    code_version = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip() + '+knn-adapter-worktree'
    spec = A.ArtifactSpec(query.to_json(), 'feature_knn', SC.T_LABEL, C.get('feature_knn').outputs[0], entry.target,
        tuple(rows.entity), 'held_out', 'outer:test:' + split.identity, dependencies, A.canonical_object(settings),
        A.canonical_object(eval_scope), 17, code_version, fit_entities=split.entities('train'), fit_role='train',
        benchmark_id=entry.benchmark_id, confidence_kind='method_support', evaluation_partition='test', gaps=gaps)
    payloads = {'rows.json': _records(rows), 'class_scores.json': {'classes': scores.columns.tolist(), 'entities': rows.entity.tolist(), 'scores': json.loads(scores.to_json(orient='values'))},
        'model_state.json': batch.model_state, 'card.json': card, 'class_cards.json': classes, 'baseline_cards.json': baseline_cards}
    # JSON transport preserves list semantics, not Python tuple container types.
    payloads = json.loads(json.dumps(payloads, allow_nan=False))
    identity = A.write_artifact(output / 'held_out', spec, payloads, split=split)
    restored = A.read_artifact(output / 'held_out', expected=spec, split=split)
    assert restored.payloads == payloads
    summary = {'organism': O.TOXOPLASMA, 'target': entry.target, 'truth_grade': entry.evidence_grade,
        'benchmark_admitted': False, 'eligible_candidate_genes': len(eligible),
        'train': len(split.entities('train')), 'tune': len(split.entities('tune')), 'calibration': len(split.entities('calibration')),
        'test': len(rows), 'retained_test_rows': len(restored.payloads['rows.json']), 'selected_training_columns': len(columns),
        'model_columns': len(batch.model_state['kept_columns']), 'class_score_columns': len(scores.columns),
        'unsupported_test_classes': unsupported, 'counts': card['counts'], 'metrics': card['metrics'],
        'majority_baseline_metrics': baseline_cards['majority']['metrics'], 'artifact_identity': identity,
        'interpretation': 'surrogate stored-prediction recovery; not independent biological accuracy',
        'native_runtime_changes': 'none'}
    (output / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    (output / 'manifest.json').write_text(json.dumps({'created_utc': datetime.now(timezone.utc).isoformat(), 'summary': summary,
        'input_sha256': {str(p): sha(p) for p in inputs},
        'outputs': {str(p.relative_to(output)): sha(p) for p in output.rglob('*.json')}}, indent=2) + '\n')
    return summary


def main():
    """Execute one controlled real-data surrogate pilot with an annotated notebook."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    nb = ExecutedNotebook('Native kNN categorical records with training-only preprocessing')
    nb.md('Frozen pilot: the registered Toxoplasma compartment candidate, seed 17, nested homology groups, k=15/min_share=.3, training-only source exclusions, feature detection/variance, median imputation and mean/std scaling. Stored compartment labels are classifier predictions; recovery of them is surrogate validation, not independent biological accuracy. No runtime strategy or installed data is changed.')
    nb.code('from scripts.freeze_label_knn_pilot import freeze', f'summary = freeze({str(args.out)!r})', 'summary')
    nb.md('Every outer-test entity is retained, including abstentions and unsupported classes. Native class scores remain method support. Training-majority/prevalence baselines use the identical outer cohort. Independent biological truth, calibration, full outer-fold coverage and the remaining categorical adapters stay open under 64.10/64.16.')
    nb.write(str(args.out / 'pilot.ipynb'))
    print(nb.ns['summary'])


if __name__ == '__main__':
    main()
