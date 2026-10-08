"""Freeze train-only controls for the verified complete-Pfam profile target.

This partition consumes an immutable preparation packet and retains its exact
ordered test cohort, unsupported profiles and unknown-state semantics. Majority
and seeded prevalence controls estimate only training annotation frequencies.
Their reference-recovery scores are control arithmetic, not independent function
accuracy or calibrated gene confidence. No classifier, feature fitting, source
acquisition, source promotion, biological admission or setting selection occurs.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
SNAPSHOT = ROOT / 'results/functional_profile_targets_2026_10_08_v3'
SNAPSHOT_SHA256 = '3e1d93cf0c587a7a8566e9f271f7cd2bdaab5421c75415e62350b74afa47c7ac'
CONTROL_ID = 'FN-DOM-CONTROL-01-Tg-pfam-complete-profile'
SEED = 20261008


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _write(path, payload):
    Path(path).write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')


def _identity(payload):
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def load_snapshot(path, *, expected_manifest_sha256=SNAPSHOT_SHA256):
    """Verify every packet byte against a pinned manifest and restore typed target data."""
    from starplast.functional_profile_targets import ProfileRow, ProfileSource, ProfileTargetContract
    from starplast.splits import Assignment, SplitManifest

    path = Path(path)
    manifest = path / 'SHA256SUMS.txt'
    if _sha(manifest) != expected_manifest_sha256:
        raise ValueError('Prepared-target manifest differs from the explicitly pinned snapshot')
    receipts = {str(manifest): expected_manifest_sha256}
    for line in manifest.read_text().splitlines():
        digest, relative = line.split('  ', 1)
        member = path / relative
        if not member.resolve().is_relative_to(path.resolve()) or member.is_symlink() or _sha(member) != digest:
            raise ValueError('Prepared-target packet member mismatch: ' + relative)
        receipts[str(member)] = digest
    for name in ('normalized_target.json', 'preparation_scope.json', 'split.json'):
        if str(path / name) not in receipts:
            raise ValueError('Prepared-target manifest omits required member: ' + name)
    saved = json.loads((path / 'normalized_target.json').read_text())
    split = saved.pop('split')
    split['assignments'] = tuple(Assignment(**row) for row in split['assignments'])
    source = saved.pop('source')
    for key in ('source_ids', 'nomenclature_releases'):
        source[key] = tuple(source[key])
    rows = saved.pop('rows')
    for row in rows:
        for key in ('recorded_identifiers', 'malformed_tokens'):
            row[key] = tuple(row[key])
    contract = ProfileTargetContract(**saved, rows=tuple(ProfileRow(**row) for row in rows),
                                     source=ProfileSource(**source), split=SplitManifest(**split))
    scope = json.loads((path / 'preparation_scope.json').read_text())
    for key, expected in (('target_identity', contract.identity), ('cohort_identity', contract.cohort_identity),
                          ('split_identity', contract.split.identity)):
        if scope[key] != expected:
            raise ValueError('Prepared-target semantic identity mismatch: ' + key)
    if _identity(json.loads((path / 'split.json').read_text())) != _identity(asdict(contract.split)):
        raise ValueError('Prepared split file differs from target split')
    return contract, scope, receipts


def build_controls(contract, *, seed=SEED):
    """Return exact train-only predictions, full source-class scores and control cards.

    The complete recorded class universe is reporting scope, not training-class
    selection. Classes absent from training have zero prior and remain in test
    metrics. Unknown gene profiles never become test negatives or target classes.
    """
    import numpy as np
    import pandas as pd
    from starplast import baselines, scorecard as SC
    from starplast.functional_profile_targets import ProfileTargetContract
    from starplast.scorecard_view import build_scorecard_view

    if not isinstance(contract, ProfileTargetContract) or contract.organism != 'Tg' or contract.source_target != 'pfam_id':
        raise ValueError('This control partition requires a typed Tg Pfam profile target')
    if type(seed) is not int:
        raise ValueError('The prevalence control seed must be a frozen integer')
    train = contract.training_labels()
    native = baselines.label_baselines(train, contract.split, role='test', seed=seed)
    counts = Counter(train)
    train_classes = sorted(counts)
    probabilities = np.array([counts[label] / len(train) for label in train_classes])
    assert sum(counts.values()) == len(train)
    np.testing.assert_array_equal(probabilities, np.array([native.attrs['prevalence'][label] for label in train_classes]))
    majority = min(label for label in train_classes if counts[label] == max(counts.values()))
    assert native.majority.tolist() == [majority] * len(native)
    assert native.prevalence_call.tolist() == np.random.default_rng(seed).choice(
        train_classes, size=len(native), p=probabilities).tolist()
    assert tuple(native.index) == contract.split.entities('test')
    # Predictions are frozen above before any test truth is accessed.
    by_gene = {row.gene_id: row for row in contract.rows}
    classes = sorted({row.profile for row in contract.rows if row.eligible})
    rows = pd.DataFrame({'entity': native.index, 'protected_group': [by_gene[gene].protected_group for gene in native.index],
                         'truth': [by_gene[gene].profile for gene in native.index],
                         'training_supported': [by_gene[gene].profile in counts for gene in native.index],
                         'majority': native.majority.to_numpy(), 'prevalence_call': native.prevalence_call.to_numpy()})
    assert all(isinstance(truth, str) for truth in rows.truth)
    control_cards, class_cards, class_scores = {}, {}, {}
    for name, column in (('training_majority', 'majority'), ('training_prevalence', 'prevalence_call')):
        prior = [float(label == majority) if name == 'training_majority' else counts[label] / len(train) for label in classes]
        scores = pd.DataFrame(np.tile(prior, (len(rows), 1)), index=native.index, columns=classes)
        scores.index.name = 'entity'
        # Independently reconstruct from the shared baseline's frozen output;
        # floating row sums need not equal the mathematical value one exactly.
        expected_scores = np.zeros((len(rows), len(classes)), dtype=float)
        for position, label in enumerate(classes):
            expected_scores[:, position] = (float(label == native.majority.iloc[0]) if name == 'training_majority'
                                            else native.attrs['prevalence'].get(label, 0.0))
        np.testing.assert_array_equal(scores.to_numpy(), expected_scores)
        unseen = set(classes) - set(train_classes)
        assert not unseen or scores[list(unseen)].eq(0).all().all()
        metrics = SC.label_calls(rows[column], rows.truth, np.arange(len(rows)), class_scores=scores)
        metrics = {key: float(value) if value is not None and np.isfinite(value) else None for key, value in metrics.items()}
        confusion = [{'truth': truth, 'prediction': prediction, 'count': int(count)}
                     for (truth, prediction), count in rows.groupby(['truth', column], sort=True).size().items()]
        per_class = []
        source_counts = Counter(row.profile for row in contract.rows if row.eligible)
        for label in classes:
            actual, predicted = rows.truth.eq(label), rows[column].eq(label)
            tp, fp, fn = int((actual & predicted).sum()), int((~actual & predicted).sum()), int((actual & ~predicted).sum())
            precision = tp / (tp + fp) if tp + fp else 0.0 if tp + fn else None
            recall = tp / (tp + fn) if tp + fn else None
            f1 = None if precision is None or recall is None else 2 * precision * recall / (precision + recall) if precision + recall else 0.0
            per_class.append({'profile': label, 'source_genes': source_counts[label], 'test_truth_genes': tp + fn,
                              'predicted_genes': tp + fp, 'training_genes': counts[label],
                              'true_positive': tp, 'false_positive': fp, 'false_negative': fn,
                              'precision': precision, 'recall': recall, 'f1': f1,
                              'evaluation_status': 'Reference-profile control' if tp + fn else 'No test truth for this class',
                              'biological_accuracy': None})
        by_class = {item['profile']: item for item in per_class}
        # Match label_calls' declared native truth-class order and reduction;
        # source classes without test truth retain unavailable class recall.
        evaluated = [by_class[label] for label in sorted(set(rows.truth), key=str)]
        assert metrics['accuracy'] == float(rows[column].eq(rows.truth).mean())
        assert metrics['coverage'] == 1.0
        for metric, field in (('macro_precision', 'precision'), ('macro_recall', 'recall'), ('macro_f1', 'f1')):
            assert metrics[metric] == float(np.mean([item[field] for item in evaluated]))
        assert sum(item['count'] for item in confusion) == len(rows)
        correct = int(rows[column].eq(rows.truth).sum())
        card = {'scope': {'organism': contract.organism, 'target': contract.target, 'strategy': name,
                          'task': SC.T_LABEL, 'unit': 'gene', 'seed': seed, 'partition': 'test',
                          'benchmark_id': CONTROL_ID, 'truth_grade': contract.source.evidence_grade,
                          'negative_semantics': contract.source.negative_semantics},
                'counts': {'eligible': len(rows), 'answered': len(rows), 'correct': correct,
                           'abstained': 0, 'wrong': len(rows) - correct}, 'metrics': metrics,
                'extra': {'role': 'training_only_baseline_control', 'confusion': confusion,
                          'class_universe': classes, 'class_scores_semantics': 'Training prior / majority one-hot control; not calibrated gene confidence',
                          'unsupported_test_genes': int((~rows.training_supported).sum()),
                          'target_identity': contract.identity, 'cohort_identity': contract.cohort_identity,
                          'split_identity': contract.split.identity, 'whole_universe': len(contract.rows),
                          'unknown_states': contract.capacity()['unknown_states'], 'biological_accuracy': None,
                          'classifier_or_feature_fitting_performed': False, 'baseline_parameters_estimated_from_training_only': True,
                          'test_support_used_for_selection': False}}
        view = build_scorecard_view(card, details={'source': {'name': list(contract.source.source_ids),
                    'grade': contract.source.evidence_grade, 'version': contract.source.source_release,
                    'negative_semantics': contract.source.negative_semantics},
                    'gaps': ['Source annotation recovery control only; independent function accuracy not admitted',
                             'Source assignment release and homology/source independence remain unresolved']})
        control_cards[name] = {'card': card, 'view_snapshot': json.loads(view.snapshot_json)}
        class_cards[name], class_scores[name] = per_class, scores
    return rows, class_scores, control_cards, class_cards


def prepare(output):
    """Execute and replay controls from verified V3 into a new immutable notebook packet."""
    import numpy as np
    import pandas as pd

    output = Path(output)
    contract, source_scope, receipts = load_snapshot(SNAPSHOT)
    assert len(contract.rows) == 8140 and sum(contract.eligibility) == 4310
    assert len(contract.split.entities('test')) == 635
    source_path = Path(source_scope['source_census_path'])
    assert _sha(source_path) == source_scope['source_census_sha256']
    receipts[str(source_path)] = source_scope['source_census_sha256']
    original_inputs = json.loads((SNAPSHOT / 'input_manifest.json').read_text())
    node_paths = [path for path, digest in original_inputs.items() if digest == source_scope['installed_source_sha256']]
    assert len(node_paths) == 1 and _sha(node_paths[0]) == source_scope['installed_source_sha256']
    receipts[node_paths[0]] = source_scope['installed_source_sha256']
    scope = {'control_id': CONTROL_ID, 'prepared_target_manifest_sha256': SNAPSHOT_SHA256,
             'source_scope': source_scope, 'prevalence_seed': SEED, 'role': 'test',
             'controls': ['training_majority', 'training_prevalence'],
             'selection_rule': 'Predeclared shared controls; no test outcomes/support select settings',
             'classifier_or_feature_fitting_performed': False, 'biological_admission': False}
    _write(output / 'control_scope.json', scope)
    rows, scores, cards, classes = build_controls(contract, seed=SEED)
    assert len(rows) == 635 and int((~rows.training_supported).sum()) == 333
    rows.to_parquet(output / 'test_rows.parquet', index=False)
    _write(output / 'control_cards.json', cards)
    _write(output / 'class_cards.json', classes)
    pd.testing.assert_frame_equal(pd.read_parquet(output / 'test_rows.parquet'), rows, check_exact=True)
    score_hashes = {}
    for name, frame in scores.items():
        path = output / (name + '_class_scores.parquet')
        frame.to_parquet(path)
        pd.testing.assert_frame_equal(pd.read_parquet(path), frame, check_exact=True)
        score_hashes[name] = _sha(path)
    assert json.loads((output / 'control_cards.json').read_text()) == cards
    assert json.loads((output / 'class_cards.json').read_text()) == classes
    code_files = [Path(__file__), ROOT / 'tests/test_functional_profile_controls.py', ROOT / 'scripts/notebook_runner.py', *[
        ROOT / 'starplast' / name for name in ('baselines.py', 'scorecard.py', 'scorecard_view.py',
        'functional_profile_targets.py', 'functional_domain_profiles.py', 'ground_truth.py', 'splits.py', 'provenance.py', 'organisms.py')]]
    (output / 'code').mkdir()
    for path in code_files:
        receipts[str(path)] = _sha(path)
        shutil.copyfile(path, output / 'code' / path.name)
    for path, digest in receipts.items():
        assert _sha(path) == digest, 'Changed input: ' + path
    _write(output / 'input_manifest.json', receipts)
    summary = {'control_id': CONTROL_ID, 'whole_universe': 8140, 'eligible_profiles': 4310,
               'unknown_genes': 3830, 'test_genes': len(rows), 'unsupported_test_genes': 333,
               'unsupported_profiles_retained': True, 'target_identity': contract.identity,
               'split_identity': contract.split.identity, 'cohort_identity': contract.cohort_identity,
               'raw_source_sha256': source_scope['installed_source_sha256'], 'raw_cell_census_sha256': source_scope['source_census_sha256'],
               'control_scope_identity': _identity(scope), 'class_score_file_sha256': score_hashes,
               'controls': {name: {'counts': bundle['card']['counts'], 'metrics': bundle['card']['metrics']} for name, bundle in cards.items()},
               'all_record_and_score_replays_exact': True, 'biological_accuracy': None, 'biological_admission': False,
               'classifier_or_feature_fitting_performed': False, 'baseline_parameters_estimated_from_training_only': True,
               'software_versions': {'python': sys.version.split()[0], 'numpy': np.__version__, 'pandas': pd.__version__}}
    _write(output / 'summary.json', summary)
    (output / 'README.md').write_text(
        '# Complete-Pfam profile baseline controls\n\n'
        'This packet consumes verified prepared-target V3 without changing its source, seed, protected split '
        'or complete recorded profiles. All 635 test genes and 333 unsupported-profile genes remain. '
        'All 3,830 unannotated genes remain unknown outside the eligible split. Training-majority and '
        'seeded training-prevalence controls estimate only the original training annotation frequencies.\n\n'
        'Scores cover the frozen recorded class universe; profiles absent in training have zero prior '
        'and are included in macro metrics and confusion. Class scores are control priors/one-hot vectors, '
        'not calibrated gene probabilities. Classes without test truth retain explicit unavailable recall. '
        'These are annotation-recovery control measurements, not an inference strategy fit, new functional '
        'claim or independent biological accuracy. Assignment release, evidence grade and source/homology '
        'independence gaps remain unresolved. Raw source hashes remain separate from semantic identities.\n\n'
        'Every saved prediction, score, card and count was exactly replayed. Code/input bytes are pinned; '
        'the notebook records actual execution. Earlier immutable source/preparation packets remain unchanged.\n')
    return summary


def main():
    """Execute controls once into a new immutable result directory, preserving failures."""
    from scripts.notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    nb = ExecutedNotebook('FN-DOM-CONTROL-01 complete-Pfam profile baseline controls')
    nb.ns.update(output=args.out, prepare=prepare)
    nb.md('Predeclared training-majority and seeded prevalence controls on the unchanged verified V3 target.',
          'Reference-annotation control arithmetic only; no classifier, source acquisition or independent biological admission.')
    try:
        nb.code('summary = prepare(output)', 'summary')
        print(json.dumps(nb.ns['summary'], indent=2, allow_nan=False))
    except Exception as exc:
        _write(args.out / 'failure.json', {'error': type(exc).__name__, 'detail': str(exc)})
        raise
    finally:
        nb.write(str(args.out / 'executed.ipynb'))
        (args.out / 'SHA256SUMS.txt').write_text('\n'.join(
            f'{_sha(path)}  {path.relative_to(args.out)}' for path in sorted(args.out.rglob('*'))
            if path.is_file() and path.name != 'SHA256SUMS.txt') + '\n')


if __name__ == '__main__':
    main()
