"""Freeze a serial native forest on the unchanged prepared Pf direct EC profiles.

This bounded candidate reuses the accepted complete annotation target and protected
split without source acquisition, truth admission, tuning or deployment. Numeric
selection and rank transforms use training only. The native forest statistics are
unchanged; its constructor uses one job solely to make tree probability reduction
exactly reproducible. Every unsupported profile and unknown source state remains
visible, and method support is never represented as calibrated biological confidence.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts import freeze_pf_functional_ec_knn as K
from scripts import freeze_pf_functional_ec_controls as C
from scripts import freeze_functional_profile_controls as P, freeze_functional_profile_knn as N

PILOT = 'FN-EC-PF-RF-01'
PREPARATION, PREPARATION_SHA = K.PREPARATION, K.PREPARATION_SHA
SETTINGS = {'pilot_id': PILOT, 'trees': 300, 'min_leaf': 2, 'execution_n_jobs': 1, 'native_n_jobs': 4,
    'max_features': 'sqrt', 'class_weight': 'balanced_subsample',
    'execution_difference': 'Single job for exact serial tree probability reduction; statistical defaults unchanged',
    'selection': 'Fixed before execution; no test/tune/calibration selection',
    'feature_selection': 'Registered numeric columns; native training coverage/variation floor',
    'feature_transform': 'Native training average ranks; frozen query ECDF', 'graphs_or_maps_used': False}


def fit_candidate(features, labels, packet, exclusions):
    """Fit fixed native statistics and score the entire original outer-test cohort."""
    import numpy as np
    import pandas as pd
    from starplast import artifacts as A, forest_records as RF, functional_ontology as F, record_scorecards as R, scorecard as SC
    from starplast.scorecard_view import build_scorecard_view
    split = packet['split']
    batch = RF.random_forest(features, labels, split=split, exclusions=exclusions,
        trees=SETTINGS['trees'], min_leaf=SETTINGS['min_leaf'], seed=split.seed)
    pd.testing.assert_frame_equal(batch.class_scores, RF.replay_scores(features, batch.model_state), check_exact=True)
    rows = batch.rows.copy()
    reference, groups = packet['profiles'].set_index('gene_id'), packet['groups'].set_index('gene_id')
    rows['truth'] = reference.loc[list(rows.entity), 'profile'].to_numpy()
    rows['group'] = groups.loc[list(rows.entity), 'group'].to_numpy()
    rows['training_supported'] = rows.truth.isin(labels)
    for column in ('entity', 'truth', 'group', 'training_supported'):
        assert rows[column].tolist() == packet['control_rows'][column].tolist(), 'Frozen control scope differs: ' + column
    classes = sorted(set(packet['profiles'].profile.dropna()))
    full_scores = batch.class_scores.reindex(columns=classes, fill_value=0.0) if len(batch.class_scores.columns) else pd.DataFrame(
        np.nan, index=batch.class_scores.index, columns=classes)
    gaps = tuple((*packet['gaps'], *batch.gaps, 'Biological accuracy/calibration/deployment unknown; reference-profile recovery only'))
    scope = R.RecordScope('Pf', 'random_forest', C.TARGET, SC.T_LABEL, A.canonical_object(SETTINGS), split.seed,
        split.identity, 'outer_test', split.benchmark_id, 'unresolved', 'gene', C.NEGATIVE_SEMANTICS, gaps)
    card = R.aggregate(rows, scope, parameters={'class_scores': full_scores.reset_index(drop=True)})
    available = bool(batch.model_state['trees'])
    card['extra'].update(pilot_id=PILOT, target_identity=packet['target_identity'], prepared_cohort_identity=packet['cohort_identity'],
        unsupported_test_genes=int((~rows.training_supported).sum()),
        unknown_states=packet['profiles'].status[~packet['profiles'].eligible].value_counts().to_dict(),
        biological_accuracy=None, biological_admission=False, calibration='unavailable', deployment='unknown',
        native_n_jobs=4, execution_n_jobs=1,
        model_status='available_reference_recovery_candidate' if available else 'unavailable_native_forest')
    view = build_scorecard_view(card, details={'source': {'grade': 'unresolved', 'name': 'pf_enzyme_classification; direct ec_number',
        'version': 'unresolved', 'negative_semantics': C.NEGATIVE_SEMANTICS}, 'gaps': list(gaps)})
    comparisons = {}
    for name, bundle in packet['controls'].items():
        control = bundle['card']
        if control['scope']['strategy'] != name or control['scope']['target'] != C.TARGET or control['scope']['organism'] != 'Pf' or control['scope']['protocol'] != split.identity:
            raise ValueError('Frozen control card scope differs')
        comparisons[name] = {'source_scope': control['scope'], 'counts': control['counts'], 'metrics': control['metrics'],
            'same_entity_order_truth_groups_and_support': True,
            'forest_minus_control': {metric: card['metrics'][metric] - value if card['metrics'][metric] is not None and value is not None else None
                                    for metric, value in control['metrics'].items()}}
    return {'batch': batch, 'rows': rows, 'class_scores': full_scores, 'card': C._plain(card),
        'view_snapshot': json.loads(view.snapshot_json),
        'profile_class_cards': C._plain(R.class_cards(rows, scope, parameters={'class_scores': full_scores.reset_index(drop=True)})),
        'major_class_cards': F.member_class_cards(rows), 'baseline_comparisons': comparisons, 'model_available': available}


def artifact_spec(packet, dependencies, *, gaps=(), model_available=True):
    """Declare exact population, native execution difference and unresolved evidence limits."""
    from starplast import artifacts as A, capabilities as Cap, scorecard as SC
    from starplast.query import Query
    split = packet['split']
    limitations = tuple(dict.fromkeys((*packet['gaps'], *gaps, 'Reference-profile recovery; independent biological accuracy unknown',
        'Native n_jobs4 replaced with1 for exact serial probability reduction; statistical defaults unchanged')))
    evaluation = {'pilot_id': PILOT, 'unit': 'gene', 'eligible_population': len(split.entities('test')),
        'truth_grade': 'unresolved', 'negative_semantics': C.NEGATIVE_SEMANTICS,
        'target_identity': packet['target_identity'], 'cohort_identity': packet['cohort_identity'],
        'prepared_packet_manifest_sha256': packet.get('manifest_hash', PREPARATION_SHA),
        'context': {'organism': 'Pf', 'assay': 'direct complete EC-major annotation', 'stage': 'unresolved', 'strain': 'unresolved'},
        'limitations': list(limitations), 'biological_admission': False, 'calibration': 'unavailable', 'deployment': 'unknown',
        'native_n_jobs': 4, 'execution_n_jobs': 1}
    return A.ArtifactSpec(Query('Pf', 'label', target=C.TARGET).to_json(), 'random_forest', SC.T_LABEL,
        Cap.get('random_forest').outputs[0], C.TARGET, split.entities('test'), 'held_out', 'outer:test:' + split.identity,
        tuple(dependencies), A.canonical_object(SETTINGS), A.canonical_object(evaluation), split.seed,
        PILOT + '-pinned-code', fit_entities=split.entities('train'), fit_role='train', benchmark_id=split.benchmark_id,
        confidence_kind='method_support', evaluation_partition='test', gaps=limitations, status='ok' if model_available else 'unavailable')


def freeze(output):
    """Freeze source, training state, exact refit and typed artifact in one new packet."""
    import gc
    import pandas as pd
    from starplast import artifacts as A, strategies as S
    from starplast.splits import read_split, write_split
    output = Path(output)
    P._write(output / 'predeclared_scope.json', {'pilot_id': PILOT, 'target': C.TARGET, 'settings': SETTINGS,
        'preparation_manifest_sha256': PREPARATION_SHA, 'biological_admission': False})
    code = [Path(__file__), ROOT / 'tests/test_pf_functional_ec_forest.py', ROOT / 'tests/test_forest_records.py',
        Path(K.__file__), Path(N.__file__), Path(C.__file__), Path(P.__file__), ROOT / 'scripts/notebook_runner.py',
        *[ROOT / 'starplast' / name for name in ('forest_records.py', 'strategy_learning.py', 'functional_exclusions.py', 'functional_ontology.py',
        'label_records.py', 'datasets.py', 'strategies.py', 'search.py', 'embedding.py', 'slots.py', 'record_scorecards.py', 'scorecard.py',
        'scorecard_view.py', 'artifacts.py', 'capabilities.py', 'splits.py', 'query.py', 'ground_truth.py', 'provenance.py', 'organisms.py')]]
    receipts = {str(path): P._sha(path) for path in code}
    (output / 'code').mkdir()
    for path in code:
        shutil.copyfile(path, output / 'code' / path.name)
    P._write(output / 'input_manifest.json', receipts)
    packet = K.load_packet()
    receipts.update(packet['receipts'])
    P._write(output / 'input_manifest.json', receipts)
    nodes = pd.read_parquet(packet['source_path'])
    whole_genes = len(nodes)
    assert P._sha(packet['source_path']) == packet['source_hash']
    cells = pd.read_parquet(PREPARATION / 'source_cells.parquet')
    assert C._records(nodes[list(cells.columns)]) == C._records(cells)
    assert S.Context(nodes, graph={}, organism='Pf').groups().tolist() == packet['groups'].group.tolist()
    features, labels, excluded, provenance = K.source_closed_inputs(nodes, packet)
    del cells, nodes
    gc.collect()
    assert len(features) == 1050 and len(labels) == 600
    P._write(output / 'feature_scope.json', {'provenance': provenance, 'exclusions': asdict(excluded),
        'target_identity': packet['target_identity'], 'cohort_identity': packet['cohort_identity'], 'split_identity': packet['split'].identity})
    result = fit_candidate(features, labels, packet, excluded)
    assert len(result['rows']) == 152 and int((~result['rows'].training_supported).sum()) == 4
    repeated = fit_candidate(features, labels, packet, excluded)
    assert result['batch'].model_state == repeated['batch'].model_state
    for key in ('card', 'view_snapshot', 'profile_class_cards', 'major_class_cards', 'baseline_comparisons'):
        assert result[key] == repeated[key], key
    for key in ('rows', 'class_scores'):
        pd.testing.assert_frame_equal(result[key], repeated[key], check_exact=True)
    pd.testing.assert_frame_equal(result['batch'].class_scores, repeated['batch'].class_scores, check_exact=True)
    del repeated
    gc.collect()
    for name, frame in {'features': features, 'rows': result['rows'], 'class_scores': result['class_scores'],
                        'native_class_scores': result['batch'].class_scores}.items():
        frame.to_parquet(output / (name + '.parquet'))
        pd.testing.assert_frame_equal(frame, pd.read_parquet(output / (name + '.parquet')), check_exact=True)
    write_split(output / 'split.json', packet['split'], excluded)
    assert read_split(output / 'split.json') == (packet['split'], excluded)
    dependencies = [A.Dependency('code', str(path.relative_to(ROOT)), receipts[str(path)]) for path in code]
    dependencies.extend((A.Dependency('table', 'installed_nodes', packet['source_hash']),
        A.Dependency('source', 'prepared_EC_packet', PREPARATION_SHA),
        A.Dependency('truth', 'complete_direct_profiles', P._sha(PREPARATION / 'direct_profiles.parquet')),
        A.Dependency('split', 'prepared_protected', packet['split'].identity),
        A.Dependency('exclusions', 'training_full_source_closure', P._identity(asdict(excluded))),
        A.Dependency('model', 'native_training_state', P._identity(N._json(result['batch'].model_state)))))
    spec = artifact_spec(packet, dependencies, gaps=provenance['gaps'], model_available=result['model_available'])
    payload = {key + '.json': result[key] for key in ('card', 'view_snapshot', 'profile_class_cards', 'major_class_cards', 'baseline_comparisons')}
    payload.update({'model_state.json': result['batch'].model_state, 'exclusions.json': asdict(excluded), 'feature_provenance.json': provenance,
        'rows.json': C._records(result['rows']), 'native_class_scores.json': N._score_payload(result['batch'].class_scores),
        'class_scores.json': N._score_payload(result['class_scores']), 'baseline_cards.json': packet['controls']})
    payload = C._plain(payload)
    identity = A.write_artifact(output / 'held_out', spec, payload, split=packet['split'])
    assert A.read_artifact(output / 'held_out', expected=spec, split=packet['split']).payloads == payload
    for path, digest in receipts.items():
        assert P._sha(path) == digest, 'Changed input: ' + path
    summary = {'pilot_id': PILOT, 'organism': 'Pf', 'target': C.TARGET, 'whole_genes': whole_genes,
        'eligible_profiles': len(features), 'train': len(labels), 'tune': len(packet['split'].entities('tune')),
        'calibration': len(packet['split'].entities('calibration')), 'test': len(result['rows']),
        'source_annotated_genes': int(packet['profiles'].status.ne('unannotated').sum()),
        'unknown_genes': int(packet['profiles'].status.eq('unannotated').sum()),
        'unresolved_annotated_genes': int(packet['profiles'].status.eq('unresolved').sum()),
        'unsupported_test_genes': int((~result['rows'].training_supported).sum()),
        'unknown_statuses': result['card']['extra']['unknown_states'],
        'selected_features': features.columns.tolist(), 'kept_features': result['batch'].model_state['kept_columns'],
        'native_classes': len(result['batch'].class_scores.columns), 'recorded_source_classes': len(result['class_scores'].columns),
        'counts': result['card']['counts'], 'metrics': result['card']['metrics'], 'baseline_comparisons': result['baseline_comparisons'],
        'artifact_identity': identity, 'artifact_status': spec.status, 'target_identity': packet['target_identity'],
        'cohort_identity': packet['cohort_identity'], 'split_identity': packet['split'].identity,
        'source_sha256': packet['source_hash'], 'preparation_manifest_sha256': PREPARATION_SHA,
        'input_receipt_count': len(receipts), 'exact_source_native_refit_and_serialization_replay': True,
        'native_n_jobs': 4, 'execution_n_jobs': 1, 'tree_count': len(result['batch'].model_state['trees']),
        'biological_accuracy': None, 'biological_admission': False, 'calibration_status': 'unavailable', 'deployment': 'unknown',
        'gaps': list(dict.fromkeys((*provenance['gaps'], *result['batch'].gaps)))}
    P._write(output / 'summary.json', summary)
    (output / 'README.md').write_text('Fixed native random forest,300trees/min_leaf2, serial n_jobs1 replacing native4 solely '
        'for exact tree probability reduction. All statistical defaults, prepared direct EC-major target, protected roles and '
        'original training-only controls retained. Complete native tree states, training-only ranks, scores/calls/support, '
        'source closure and exact source/refit/tree-JSON/parquet/artifact replay recorded. Unsupported profiles and source '
        'unknowns remain visible. Source assignment lineage, independent biological accuracy, calibration and deployment '
        'remain unknown. No graphs/maps/tuning/permutation importance or biological admission.\n')
    return summary


def main():
    """Execute once under the externally enforced cap and retain failures without retries."""
    from scripts.notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    notebook = ExecutedNotebook(PILOT)
    notebook.ns.update(freeze=freeze, output=args.out)
    notebook.md('Fixed native300trees/min_leaf2, serial n_jobs1; unchanged accepted direct Pf EC target/protected split. '
        'Method support only; no biological admission, tuning or deployment.')
    try:
        notebook.code('summary = freeze(output)', 'summary')
        print(json.dumps(notebook.ns['summary'], indent=2, allow_nan=False))
    except Exception as exc:
        P._write(args.out / 'failure.json', {'error': type(exc).__name__, 'detail': str(exc)})
        notebook.md('Failure: ' + type(exc).__name__ + ': ' + str(exc))
        raise
    finally:
        notebook.write(str(args.out / 'executed.ipynb'))
        (args.out / 'SHA256SUMS.txt').write_text('\n'.join(f'{P._sha(path)}  {path.relative_to(args.out)}'
            for path in sorted(args.out.rglob('*')) if path.is_file() and path.name != 'SHA256SUMS.txt') + '\n')


if __name__ == '__main__':
    main()
