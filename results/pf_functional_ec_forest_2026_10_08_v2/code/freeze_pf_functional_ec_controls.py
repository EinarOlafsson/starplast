"""Freeze direct Pf EC-major profile preparation and train-only baseline controls.

This bounded packet replays existing pinned ENZYME nomenclature and installed
annotations without acquiring sources or fitting a classifier. Complete direct
profiles preserve every resolved major class; unresolved and unannotated genes
stay unknown. Orthology-transferred annotations form a separate census only.
Native grouping, frozen roles, annotation-recovery arithmetic and source gaps
remain explicit. Neither control performance nor nomenclature resolution admits
independent biological accuracy, calibrated gene confidence or deployment scope.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import freeze_functional_profile_controls as P

REVIEW = ROOT / 'results/functional_source_review_2026_10_08_v2'
REVIEW_SHA = '5fb85c26e4939c4737a6969c474cd63395ea2faaa1fb111a71ebd241e0aaa81c'
BENCHMARK = 'FN-EC-PF-CONTROL-01'
TARGET = 'ec_direct_complete_major_profile'
SEED = 20261008
FRACTIONS = (0.55, 0.15, 0.15, 0.15)
NEGATIVE_SEMANTICS = 'unannotated_is_unknown; recovery of complete recorded profiles only'
GAPS = ('Direct EC individual assignment evidence grades/source release remain unresolved',
        'ENZYME nomenclature is not independent measured gene activity',
        'Native orthogroups do not establish complete homology/source independence',
        'Orthology transfer remains a separate census; never combined into direct truth',
        'Calibration and deployment applicability unknown')
SCOPE = {'benchmark_id': BENCHMARK, 'organism': 'Pf', 'source_target': 'ec_number', 'target': TARGET,
         'seed': SEED, 'fractions': list(FRACTIONS), 'controls': ['training_majority', 'training_prevalence'],
         'selection': 'Fixed before source preparation; no test/tune support or metrics select settings',
         'group_policy': 'Exact native Context.groups on whole table positions; collision refused',
         'truth_grade': 'unresolved', 'negative_semantics': NEGATIVE_SEMANTICS,
         'orthology_role': 'separate_census_only', 'classifier_fitting_performed': False,
         'biological_admission': False, 'review_manifest_sha256': REVIEW_SHA}


def _plain(value):
    import numpy as np
    import pandas as pd
    if isinstance(value, dict):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_plain(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    return None if value is None or value is pd.NA or isinstance(value, float) and pd.isna(value) else value


def _records(frame):
    return _plain(frame.to_dict('records'))


def verified_inputs():
    """Verify the archived review, receipt and current Pf nodes before projected reads."""
    import pandas as pd
    import pyarrow.parquet as pq
    from starplast import functional_ontology as F, organisms as O
    manifest_path = REVIEW / 'manifest.json'
    if P._sha(manifest_path) != REVIEW_SHA:
        raise ValueError('Pinned functional source review manifest changed')
    manifest = json.loads(manifest_path.read_text())
    receipts = {str(manifest_path): REVIEW_SHA}
    for name, digest in manifest['output_sha256'].items():
        path = REVIEW / name
        if path.is_symlink() or not path.resolve().is_relative_to(REVIEW.resolve()) or P._sha(path) != digest:
            raise ValueError('Functional review output mismatch: ' + name)
        receipts[str(path)] = digest
    enzyme_receipts = [row for row in json.loads((REVIEW / 'receipts.json').read_text())
                       if row['name'] == 'enzyme.dat' and row['status'] == 'available']
    if len(enzyme_receipts) != 1:
        raise ValueError('One authoritative ENZYME receipt required')
    enzyme_receipt = enzyme_receipts[0]
    enzyme = Path(enzyme_receipt['path'])
    nodes = Path(O.nodes_path(O.FALCIPARUM))
    for path, digest in ((enzyme, enzyme_receipt['sha256']), (nodes, manifest['input_sha256'][str(nodes)])):
        if P._sha(path) != digest:
            raise ValueError('Pinned original source changed: ' + str(path))
        receipts[str(path)] = digest
    if manifest['input_sha256'][str(enzyme)] != enzyme_receipt['sha256']:
        raise ValueError('ENZYME receipt and original manifest differ')
    schema = pq.read_schema(nodes)
    columns = [name for name in ('gene_id', 'ec_number', 'ec_number_orthology', 'orthogroup') if name in schema.names]
    frame = pd.read_parquet(nodes, columns=columns)
    entries = F.parse_enzyme(enzyme.read_text())
    return frame, entries, enzyme_receipt, receipts


def prepare_profiles(nodes, entries):
    """Prepare direct profiles, separate transfer census and the exact native frozen split."""
    import numpy as np
    import pandas as pd
    from starplast import functional_exclusions as E, functional_ontology as F, strategies as S
    from starplast.embedding import as_text
    from starplast.ground_truth import cohort_digest
    from starplast.splits import make_split
    if ('gene_id' not in nodes or 'ec_number' not in nodes or nodes.columns.has_duplicates or
            nodes.gene_id.duplicated().any() or not all(isinstance(gene, str) and gene for gene in nodes.gene_id)):
        raise ValueError('Exact unique Pf original gene/source fields required')
    ctx = S.Context(nodes.copy(deep=True), graph={}, organism='Pf')
    direct, resolutions = F.ec_profiles(ctx, 'ec_number', entries)
    transfer, transfer_resolutions = F.ec_profiles(ctx, 'ec_number_orthology', entries) if 'ec_number_orthology' in nodes else (pd.DataFrame(), pd.DataFrame())
    groups = ctx.groups()
    original = as_text(nodes.orthogroup).to_numpy() if 'orthogroup' in nodes else np.array([''] * len(nodes))
    expected = np.array([group if group else f'__{position}' for position, group in enumerate(original)])
    np.testing.assert_array_equal(groups, expected)
    fallback = original == ''
    collisions = sorted(set(groups[fallback]) & set(groups[~fallback]))
    if collisions:
        raise ValueError('Native missing-group fallback collides with recorded orthogroup: ' + ', '.join(collisions))
    mask = direct.eligible.to_numpy(dtype=bool)
    split = make_split(tuple(direct.gene_id[mask]), tuple(groups[mask]), organism='Pf', benchmark_id=BENCHMARK,
                       group_kind='homology', seed=SEED, fractions=FRACTIONS)
    indexed = nodes.set_index('gene_id', drop=False)
    train = indexed.loc[list(split.entities('train'))].copy()
    labels = direct.set_index('gene_id').profile.loc[list(split.entities('train'))]
    train[TARGET] = labels
    excluded = E.make_exclusions(train.reset_index(drop=True), TARGET, source_targets=('ec_number',),
                                 benchmark_id=BENCHMARK, split=split)
    roles = {row.entity: row.role for row in split.assignments}
    group_frame = pd.DataFrame({'gene_id': nodes.gene_id, 'group': groups, 'fallback': fallback,
        'eligible': mask, 'role': [roles.get(gene) for gene in nodes.gene_id]})
    group_frame['source_group'] = nodes.orthogroup if 'orthogroup' in nodes else None
    return {'direct': direct, 'direct_resolutions': resolutions, 'orthology': transfer,
            'orthology_resolutions': transfer_resolutions, 'groups': group_frame, 'split': split,
            'exclusions': excluded, 'labels': labels, 'cohort_identity': cohort_digest(tuple(nodes.gene_id), mask),
            'gaps': list(GAPS) + [str(int(fallback.sum())) + ' whole genes use unresolved native singleton groups']}


def control_records(prepared):
    """Freeze train-only majority/prevalence and retain complete matched test records/cards."""
    import numpy as np
    import pandas as pd
    from starplast import baselines as B, functional_ontology as F, scorecard as SC
    from starplast.scorecard_view import build_scorecard_view
    split, labels = prepared['split'], prepared['labels']
    native = B.label_baselines(labels, split, role='test', seed=SEED)
    assert tuple(native.index) == split.entities('test')
    counts = Counter(labels)
    classes = sorted(counts)
    majority = min(label for label in classes if counts[label] == max(counts.values()))
    probabilities = np.array([counts[label] / len(labels) for label in classes])
    assert native.majority.tolist() == [majority] * len(native)
    assert native.prevalence_call.tolist() == np.random.default_rng(SEED).choice(classes, len(native), p=probabilities).tolist()
    source_classes = sorted(set(prepared['direct'].profile.dropna()))
    by_gene = prepared['direct'].set_index('gene_id')
    groups = prepared['groups'].set_index('gene_id')
    rows = pd.DataFrame({'entity': native.index, 'truth': by_gene.loc[list(native.index), 'profile'].to_numpy(),
        'group': groups.loc[list(native.index), 'group'].to_numpy(),
        'training_supported': by_gene.loc[list(native.index), 'profile'].isin(labels).to_numpy(),
        'majority': native.majority.to_numpy(), 'prevalence_call': native.prevalence_call.to_numpy()})
    cards, scores = {}, {}
    for name, column in (('training_majority', 'majority'), ('training_prevalence', 'prevalence_call')):
        prior = [float(label == majority) if column == 'majority' else counts[label] / len(labels) for label in source_classes]
        score = pd.DataFrame(np.tile(prior, (len(rows), 1)), columns=source_classes)
        expected = np.zeros(score.shape)
        for position, label in enumerate(source_classes):
            expected[:, position] = float(label == native.majority.iloc[0]) if column == 'majority' else native.attrs['prevalence'].get(label, 0.0)
        np.testing.assert_array_equal(score.to_numpy(), expected)
        outcomes = rows.drop(columns=['majority', 'prevalence_call']).assign(prediction=rows[column], abstained=False)
        metrics = _plain(SC.label_calls(outcomes.prediction, outcomes.truth, np.arange(len(rows)), class_scores=score))
        correct = int(outcomes.prediction.eq(outcomes.truth).sum())
        confusion = [{'truth': truth, 'prediction': prediction, 'count': int(count)}
            for (truth, prediction), count in outcomes.groupby(['truth', 'prediction'], sort=True).size().items()]
        card = {'scope': {'organism': 'Pf', 'strategy': name, 'target': TARGET, 'task': SC.T_LABEL,
            'unit': 'gene', 'settings': {'role': 'training_only_baseline_control', 'control_name': name,
            'selection': SCOPE['selection'], 'classifier_fitting_performed': False}, 'seed': SEED,
            'protocol': split.identity, 'partition': 'outer_test', 'benchmark_id': BENCHMARK,
            'truth_grade': 'unresolved', 'negative_semantics': NEGATIVE_SEMANTICS},
            'counts': {'eligible': len(rows), 'answered': len(rows), 'abstained': 0, 'correct': correct, 'wrong': len(rows) - correct},
            'metrics': metrics, 'extra': {'control_name': name, 'role': 'training_only_baseline_control',
            'classifier_fitting_performed': False, 'baseline_parameters_estimated_from_training_only': True,
            'test_support_used_for_selection': False, 'confusion': confusion, 'class_universe': source_classes,
            'unsupported_test_genes': int((~rows.training_supported).sum()), 'biological_accuracy': None,
            'biological_admission': False, 'unknown_states': prepared['direct'].status[~prepared['direct'].eligible].value_counts().to_dict(),
            'class_score_semantics': 'Training prior/majority one-hot control; not calibrated gene confidence'}}
        profile_cards = []
        for label in source_classes:
            actual, predicted = outcomes.truth.eq(label), outcomes.prediction.eq(label)
            tp, fp, fn = int((actual & predicted).sum()), int((~actual & predicted).sum()), int((actual & ~predicted).sum())
            precision = tp / (tp + fp) if tp + fp else 0.0 if tp + fn else None
            recall = tp / (tp + fn) if tp + fn else None
            f1 = None if precision is None or recall is None else 2 * precision * recall / (precision + recall) if precision + recall else 0.0
            profile_cards.append({'class': label, 'source_genes': int(prepared['direct'].profile.eq(label).sum()),
                'training_genes': counts[label], 'test_truth_genes': tp + fn, 'predicted_genes': tp + fp,
                'true_positive': tp, 'false_positive': fp, 'false_negative': fn,
                'precision': precision, 'recall': recall, 'f1': f1, 'biological_accuracy': None,
                'status': 'reference_profile_control' if tp + fn else 'unavailable_no_test_truth'})
        by_class = {item['class']: item for item in profile_cards}
        evaluated = [by_class[label] for label in sorted(set(rows.truth), key=str)]
        for metric, field in (('macro_precision', 'precision'), ('macro_recall', 'recall'), ('macro_f1', 'f1')):
            assert metrics[metric] == float(np.mean([item[field] for item in evaluated]))
        member_cards = F.member_class_cards(outcomes)
        view = build_scorecard_view(card, title=name + ': reference-profile control', details={'source': {
            'grade': 'unresolved', 'name': 'pf_enzyme_classification; direct ec_number', 'version': 'unresolved',
            'negative_semantics': NEGATIVE_SEMANTICS}, 'gaps': prepared['gaps']})
        cards[name] = {'card': _plain(card), 'view_snapshot': json.loads(view.snapshot_json),
                       'profile_class_cards': _plain(profile_cards), 'major_class_cards': member_cards}
        scores[name] = score
        assert sum(row['count'] for row in card['extra']['confusion']) == len(rows)
    return rows, scores, cards


def freeze(output):
    """Execute source replay and exact train-only controls once in an immutable notebook packet."""
    import pandas as pd
    from starplast.splits import write_split, read_split
    output = Path(output)
    P._write(output / 'predeclared_scope.json', SCOPE)
    nodes, entries, enzyme_receipt, receipts = verified_inputs()
    prepared = prepare_profiles(nodes, entries)
    for target, key in (('ec_number', 'direct'), ('ec_number_orthology', 'orthology')):
        if not prepared[key].empty:
            archived = pd.read_parquet(REVIEW / ('Pf_' + target + '_profiles.parquet'))
            assert _records(prepared[key]) == _records(archived), 'Source profile replay: ' + target
    rows, scores, cards = control_records(prepared)
    replay = prepare_profiles(nodes, entries)
    repeated_rows, repeated_scores, repeated_cards = control_records(replay)
    assert cards == repeated_cards and prepared['split'] == replay['split'] and prepared['exclusions'] == replay['exclusions']
    pd.testing.assert_frame_equal(prepared['groups'], replay['groups'], check_exact=True)
    pd.testing.assert_frame_equal(rows, repeated_rows, check_exact=True)
    frames = {'source_cells': nodes, 'direct_profiles': prepared['direct'], 'direct_resolutions': prepared['direct_resolutions'],
        'orthology_profiles': prepared['orthology'], 'orthology_resolutions': prepared['orthology_resolutions'],
        'group_roles': prepared['groups'], 'control_rows': rows, **scores}
    for name, frame in frames.items():
        frame.to_parquet(output / (name + '.parquet'), index=False)
        restored = pd.read_parquet(output / (name + '.parquet'))
        assert _records(frame) == _records(restored), 'Exact stored records: ' + name
        assert frame.columns.tolist() == restored.columns.tolist()
    for name in scores:
        pd.testing.assert_frame_equal(scores[name], repeated_scores[name], check_exact=True)
    write_split(output / 'split.json', prepared['split'], prepared['exclusions'])
    assert read_split(output / 'split.json') == (prepared['split'], prepared['exclusions'])
    P._write(output / 'control_cards.json', cards)
    assert json.loads((output / 'control_cards.json').read_text()) == cards
    normalized = {'profiles': _records(prepared['direct']), 'source_target': 'ec_number',
        'source_release': 'unresolved', 'evidence_grade': 'unresolved', 'enzyme_sha256': enzyme_receipt['sha256'],
        'groups': _records(prepared['groups']), 'split': asdict(prepared['split']), 'negative_semantics': NEGATIVE_SEMANTICS}
    P._write(output / 'normalized_target.json', normalized)
    assert json.loads((output / 'normalized_target.json').read_text()) == _plain(normalized)
    P._write(output / 'enzyme_receipt.json', enzyme_receipt)
    code = [Path(__file__), ROOT / 'tests/test_pf_functional_ec_controls.py', Path(P.__file__), ROOT / 'scripts/notebook_runner.py',
        *[ROOT / 'starplast' / name for name in ('functional_ontology.py', 'functional_exclusions.py', 'discovery_labels.py',
        'strategies.py', 'datasets.py', 'slots.py', 'search.py', 'embedding.py', 'baselines.py', 'record_scorecards.py',
        'scorecard.py', 'scorecard_view.py', 'artifacts.py', 'capabilities.py', 'splits.py', 'ground_truth.py', 'organisms.py', 'provenance.py')]]
    (output / 'code').mkdir()
    for path in code:
        receipts[str(path)] = P._sha(path)
        shutil.copyfile(path, output / 'code' / path.name)
    for path, digest in receipts.items():
        assert P._sha(path) == digest, 'Changed input: ' + path
    P._write(output / 'input_manifest.json', receipts)
    summary = {'benchmark_id': BENCHMARK, 'whole_genes': len(nodes), 'direct_status': prepared['direct'].status.value_counts().to_dict(),
        'direct_eligible': int(prepared['direct'].eligible.sum()), 'direct_profile_classes': int(prepared['direct'].profile.nunique()),
        'orthology_status': prepared['orthology'].status.value_counts().to_dict() if not prepared['orthology'].empty else None,
        'orthology_role': 'separate_census_only', 'roles': {role: len(prepared['split'].entities(role)) for role in ('train', 'tune', 'calibration', 'test')},
        'missing_groups': int(prepared['groups'].fallback.sum()), 'fallback_collision_count': 0,
        'unsupported_test_genes': int((~rows.training_supported).sum()),
        'unsupported_test_profiles': rows.truth[~rows.training_supported].value_counts().to_dict(),
        'cohort_identity': prepared['cohort_identity'], 'normalized_target_identity': P._identity(normalized),
        'split_identity': prepared['split'].identity, 'predeclared_scope_identity': P._identity(SCOPE),
        'original_node_sha256': receipts[str(__import__('starplast.organisms', fromlist=['nodes_path']).nodes_path('Pf'))],
        'review_manifest_sha256': REVIEW_SHA, 'enzyme_sha256': enzyme_receipt['sha256'],
        'input_receipt_count': len(receipts), 'exact_source_groups_controls_cards_and_serialization_replay': True,
        'controls': {name: {'counts': bundle['card']['counts'], 'metrics': bundle['card']['metrics']} for name, bundle in cards.items()},
        'classifier_fitting_performed': False, 'biological_accuracy': None, 'biological_admission': False,
        'calibration': 'unavailable', 'deployment': 'unknown', 'gaps': prepared['gaps']}
    P._write(output / 'summary.json', summary)
    (output / 'README.md').write_text('Predeclared direct Pf EC-major complete-profile annotation-recovery controls. '
        'ENZYME 02-Sep-2026 (SIB, CC BY 4.0) nomenclature resolution does not admit independent gene activity. '
        'Original source cells/order, every unresolved/unannotated profile and native missing-group fallback remain explicit. '
        'Transferred ec_number_orthology is a separate census, excluded from target inputs. Majority/prevalence parameters '
        'use training profiles only. Each card declares its control name; no classifier '
        'inference was run. Cards include full cohort/confusion, complete-profile and overlapping major-class '
        'metrics. Scores are control priors/one-hot vectors, not calibrated gene probabilities. All exact replays and source/code '
        'receipts are recorded; source/homology independence, biological accuracy and deployment applicability remain unknown.\n')
    return summary


def main():
    """Execute controls once under a notebook namespace and preserve any failed packet."""
    from scripts.notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    nb = ExecutedNotebook(BENCHMARK)
    nb.ns.update(freeze=freeze, output=args.out)
    nb.md('Predeclared seed/fractions/direct-target/controls. Orthology transfer is separate census only. No classifier, source acquisition or independent biological admission.')
    try:
        nb.code('summary = freeze(output)', 'summary')
        print(json.dumps(nb.ns['summary'], indent=2, allow_nan=False))
    except Exception as exc:
        P._write(args.out / 'failure.json', {'error': type(exc).__name__, 'detail': str(exc)})
        nb.md('Failed: ' + type(exc).__name__ + ': ' + str(exc))
        raise
    finally:
        nb.write(str(args.out / 'executed.ipynb'))
        (args.out / 'SHA256SUMS.txt').write_text('\n'.join(f'{P._sha(path)}  {path.relative_to(args.out)}'
            for path in sorted(args.out.rglob('*')) if path.is_file() and path.name != 'SHA256SUMS.txt') + '\n')


if __name__ == '__main__':
    main()
