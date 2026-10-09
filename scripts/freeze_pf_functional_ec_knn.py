"""Freeze native Pf EC-major annotation-profile recovery on the accepted control split.

This candidate consumes an immutable direct-annotation preparation packet without
changing its complete profiles, eligibility, gene order or protected partitions.
Registered numeric input selection and source closure use training rows only.
Native distance votes remain separate from full-source reporting scores and from
calibrated gene confidence. Direct and orthology-transferred annotations cannot
seed features. No graph, map, tuning, deployment or independent biological truth
admission is performed; unresolved source lineage remains an explicit limitation.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
import json
from pathlib import Path
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from starplast import organisms as _ORGANISMS
from scripts import freeze_functional_profile_controls as P, freeze_functional_profile_knn as N
from scripts import freeze_pf_functional_ec_controls as C

PREPARATION = ROOT / 'results/pf_functional_ec_controls_2026_10_08'
PREPARATION_SHA = 'a29e416488ba7b0722afca63764a27afcaa6f1fa521a5539a247837a9a7d38c7'
PILOT = 'FN-EC-PF-KNN-01'
SETTINGS = {'pilot_id': PILOT, 'k': 15, 'min_share': 0.3,
    'selection': 'Fixed before execution; no test/tune/calibration selection',
    'feature_selection': 'Registered numeric columns; native training coverage/variation floor',
    'feature_transform': 'Native training average ranks; frozen query ECDF', 'graphs_or_maps_used': False}


def load_packet(path=PREPARATION, *, manifest_sha256=PREPARATION_SHA):
    """Verify the accepted preparation bytes and restore its unchanged direct EC split."""
    import pandas as pd
    from starplast.ground_truth import cohort_digest
    from starplast.splits import read_split
    path = Path(path)
    checksum = path / 'SHA256SUMS.txt'
    if P._sha(checksum) != manifest_sha256:
        raise ValueError('Preparation checksum manifest differs from pinned source')
    receipts = {str(checksum): manifest_sha256}
    for line in checksum.read_text().splitlines():
        digest, relative = line.split('  ', 1)
        member = path / relative
        if member.is_symlink() or not member.resolve().is_relative_to(path.resolve()) or P._sha(member) != digest:
            raise ValueError('Prepared EC packet member mismatch: ' + relative)
        receipts[str(member)] = digest
    names = ('direct_profiles.parquet', 'group_roles.parquet', 'control_rows.parquet', 'control_cards.json',
             'split.json', 'normalized_target.json', 'summary.json', 'source_cells.parquet', 'predeclared_scope.json')
    if any(str(path / name) not in receipts for name in names):
        raise ValueError('Prepared EC manifest omits required data')
    profiles = pd.read_parquet(path / 'direct_profiles.parquet')
    groups = pd.read_parquet(path / 'group_roles.parquet')
    split, _ = read_split(path / 'split.json')
    summary = json.loads((path / 'summary.json').read_text())
    normalized = json.loads((path / 'normalized_target.json').read_text())
    scope = json.loads((path / 'predeclared_scope.json').read_text())
    if scope['organism'] != _ORGANISMS.FALCIPARUM or scope['source_target'] != 'ec_number' or scope['target'] != C.TARGET or scope['biological_admission'] is not False:
        raise ValueError('Only the accepted direct Pf EC preparation scope is supported')
    assert profiles.gene_id.tolist() == groups.gene_id.tolist()
    assert C._records(profiles) == normalized['profiles'] and C._records(groups) == normalized['groups']
    assert P._identity(normalized) == summary['normalized_target_identity']
    assert cohort_digest(tuple(profiles.gene_id), profiles.eligible.to_numpy(dtype=bool)) == summary['cohort_identity']
    assert split.identity == summary['split_identity'] and tuple(profiles.gene_id[profiles.eligible]) == tuple(a.entity for a in split.assignments)
    original = json.loads((path / 'input_manifest.json').read_text())
    source_paths = [name for name, digest in original.items() if digest == summary['original_node_sha256']]
    assert len(source_paths) == 1 and P._sha(source_paths[0]) == summary['original_node_sha256']
    receipts[source_paths[0]] = summary['original_node_sha256']
    return {'profiles': profiles, 'groups': groups, 'split': split,
        'control_rows': pd.read_parquet(path / 'control_rows.parquet'),
        'controls': json.loads((path / 'control_cards.json').read_text()),
        'target_identity': summary['normalized_target_identity'], 'cohort_identity': summary['cohort_identity'],
        'source_path': source_paths[0], 'source_hash': summary['original_node_sha256'],
        'manifest_hash': manifest_sha256, 'receipts': receipts, 'gaps': summary['gaps']}


def source_closed_inputs(nodes, packet, *, derived_inputs=()):
    """Rebuild full-schema functional exclusions and close every late source veto."""
    import pandas as pd
    from starplast import datasets as D, functional_exclusions as E, strategies as S
    from starplast.splits import SplitManifest
    split, profiles = packet['split'], packet['profiles']
    if (not isinstance(split, SplitManifest) or split.organism != _ORGANISMS.FALCIPARUM or 'gene_id' not in nodes or 'ec_number' not in nodes or profiles.organism.ne(_ORGANISMS.FALCIPARUM).any() or
            profiles.source_target.ne('ec_number').any() or nodes.columns.has_duplicates or C.TARGET in nodes or
            tuple(nodes.gene_id) != tuple(profiles.gene_id) or tuple(profiles.gene_id[profiles.eligible]) != tuple(a.entity for a in split.assignments)):
        raise ValueError('Exact unchanged Pf direct EC source universe and typed split required')
    indexed = nodes.set_index('gene_id', drop=False)
    labels = profiles.set_index('gene_id').profile.loc[list(split.entities('train'))].copy()
    train = indexed.loc[list(labels.index)].copy()
    train[C.TARGET] = labels
    excluded = E.make_exclusions(train.reset_index(drop=True), C.TARGET, source_targets=('ec_number',),
        benchmark_id=split.benchmark_id, split=split, derived_inputs=derived_inputs)
    metadata, unregistered, unresolved = {}, [], []
    for column in nodes:
        source = D.provenance(column, _ORGANISMS.FALCIPARUM)
        if source is None:
            unregistered.append(column)
        elif any(parent not in nodes for parent in source.derived_from):
            unresolved.append(column)
        else:
            metadata[column] = {'source_id': source.key, 'parents': list(source.derived_from),
                'evidence_grade': 'unresolved', 'source_release': 'unresolved', 'status': 'Registry assertion; assignment independence unverified'}
    banned = set(excluded.columns) | set(unregistered) | set(unresolved)
    parents = {column: set(source['parents']) for column, source in metadata.items()}
    for item in derived_inputs:
        if item.kind == 'column':
            parents.setdefault(item.name, set()).update(name for kind, name in item.parents if kind == 'column')
    reasons = {}
    while True:
        added = set(D.derived_dependents(banned, _ORGANISMS.FALCIPARUM)) - banned
        for column in added:
            reasons.setdefault(column, {'reason': 'Pf registered dependency closure',
                'vetoed_parents': sorted(set(D.derived_sources(column, _ORGANISMS.FALCIPARUM)) & banned)})
        for column, source_parents in parents.items():
            vetoed = sorted(source_parents & banned)
            if column not in banned and vetoed:
                added.add(column)
                reasons[column] = {'reason': 'Declared descendant of vetoed input', 'vetoed_parents': vetoed}
        if not added:
            break
        banned.update(added)
    excluded = replace(excluded, columns=tuple(sorted(banned)))
    columns = S.Context(train.reset_index(drop=True), graph={}, organism=_ORGANISMS.FALCIPARUM).numeric_columns(exclude=excluded.columns)
    features = indexed.loc[[a.entity for a in split.assignments], columns].copy()
    provenance = {'policy': E.POLICY_VERSION, 'selected_columns': columns, 'fit_entities': list(labels.index),
        'sources': {column: metadata[column] for column in columns}, 'unregistered_withheld': unregistered,
        'unknown_parent_withheld': unresolved, 'dependency_closure_withheld': reasons,
        'gaps': list(packet['gaps']) + ['Registry/empirical closure is not complete source-independence proof',
                                     'No graph/map inputs; feature assignment evidence grades remain unresolved']}
    return features, labels, excluded, provenance


def fit_candidate(features, labels, packet, exclusions):
    """Fit one fixed native adapter, preserving complete test truth and matched controls."""
    import numpy as np
    import pandas as pd
    from starplast import artifacts as A, functional_ontology as F, label_records as L, record_scorecards as R, scorecard as SC
    from starplast.scorecard_view import build_scorecard_view
    split = packet['split']
    batch = L.feature_knn(features, labels, split=split, exclusions=exclusions,
                          k=SETTINGS['k'], min_share=SETTINGS['min_share'])
    calls, support, native_scores, votes = N.native_replay(features, labels, split, batch.model_state)
    assert batch.rows.prediction.astype(object).where(batch.rows.prediction.notna(), None).tolist() == calls.tolist()
    np.testing.assert_array_equal(batch.rows.support.to_numpy(dtype=float), support.to_numpy(dtype=float))
    pd.testing.assert_frame_equal(batch.class_scores, native_scores, check_exact=True)
    rows = batch.rows.copy()
    reference = packet['profiles'].set_index('gene_id')
    groups = packet['groups'].set_index('gene_id')
    rows['truth'] = reference.loc[list(rows.entity), 'profile'].to_numpy()
    rows['group'] = groups.loc[list(rows.entity), 'group'].to_numpy()
    rows['training_supported'] = rows.truth.isin(labels)
    for column in ('entity', 'truth', 'group', 'training_supported'):
        assert rows[column].tolist() == packet['control_rows'][column].tolist(), 'Frozen control scope differs: ' + column
    classes = sorted(set(packet['profiles'].profile.dropna()))
    full_scores = native_scores.reindex(columns=classes, fill_value=0.0) if len(native_scores.columns) else pd.DataFrame(
        np.nan, index=native_scores.index, columns=classes)
    gaps = tuple((*packet['gaps'], *batch.gaps, 'Biological accuracy/calibration/deployment unknown; reference-profile recovery only'))
    scope = R.RecordScope(_ORGANISMS.FALCIPARUM, 'feature_knn', C.TARGET, SC.T_LABEL, A.canonical_object(SETTINGS), split.seed,
        split.identity, 'outer_test', split.benchmark_id, 'unresolved', 'gene', C.NEGATIVE_SEMANTICS, gaps)
    card = R.aggregate(rows, scope, parameters={'class_scores': full_scores.reset_index(drop=True)})
    available = bool(batch.model_state['kept_columns']) and len(labels) >= 2
    card['extra'].update(pilot_id=PILOT, target_identity=packet['target_identity'], prepared_cohort_identity=packet['cohort_identity'],
        unsupported_test_genes=int((~rows.training_supported).sum()),
        unknown_states=packet['profiles'].status[~packet['profiles'].eligible].value_counts().to_dict(),
        biological_accuracy=None, biological_admission=False, calibration='unavailable', deployment='unknown',
        model_status='available_reference_recovery_candidate' if available else 'unavailable_no_supported_native_vote')
    view = build_scorecard_view(card, details={'source': {'grade': 'unresolved', 'name': 'pf_enzyme_classification; direct ec_number',
        'version': 'unresolved', 'negative_semantics': C.NEGATIVE_SEMANTICS}, 'gaps': list(gaps)})
    comparisons = {}
    for name, bundle in packet['controls'].items():
        control = bundle['card']
        if control['scope']['strategy'] != name or control['scope']['target'] != C.TARGET or control['scope']['organism'] != _ORGANISMS.FALCIPARUM or control['scope']['protocol'] != split.identity:
            raise ValueError('Frozen control card scope differs')
        comparisons[name] = {'source_scope': control['scope'], 'counts': control['counts'], 'metrics': control['metrics'],
            'same_entity_order_truth_groups_and_support': True,
            'knn_minus_control': {metric: card['metrics'][metric] - value if card['metrics'][metric] is not None and value is not None else None
                                 for metric, value in control['metrics'].items()}}
    return {'batch': batch, 'rows': rows, 'class_scores': full_scores, 'raw_votes': votes, 'card': C._plain(card),
        'view_snapshot': json.loads(view.snapshot_json), 'profile_class_cards': C._plain(R.class_cards(rows, scope, parameters={'class_scores': full_scores.reset_index(drop=True)})),
        'major_class_cards': F.member_class_cards(rows), 'baseline_comparisons': comparisons, 'model_available': available}


def artifact_spec(packet, dependencies, *, gaps=(), model_available=True):
    """Declare exact held-out population/unit, source limits and unavailable model state."""
    from starplast import artifacts as A, capabilities as Cap, scorecard as SC
    from starplast.query import Query
    split = packet['split']
    limitations = tuple(dict.fromkeys((*packet['gaps'], *gaps, 'Reference-profile recovery; independent biological accuracy unknown')))
    evaluation = {'pilot_id': PILOT, 'unit': 'gene', 'eligible_population': len(split.entities('test')),
        'truth_grade': 'unresolved', 'negative_semantics': C.NEGATIVE_SEMANTICS,
        'target_identity': packet['target_identity'], 'cohort_identity': packet['cohort_identity'],
        'prepared_packet_manifest_sha256': packet.get('manifest_hash', PREPARATION_SHA),
        'context': {'organism': _ORGANISMS.FALCIPARUM, 'assay': 'direct complete EC-major annotation', 'stage': 'unresolved', 'strain': 'unresolved'},
        'limitations': list(limitations), 'biological_admission': False, 'calibration': 'unavailable', 'deployment': 'unknown'}
    return A.ArtifactSpec(Query(_ORGANISMS.FALCIPARUM, 'label', target=C.TARGET).to_json(), 'feature_knn', SC.T_LABEL,
        Cap.get('feature_knn').outputs[0], C.TARGET, split.entities('test'), 'held_out', 'outer:test:' + split.identity,
        tuple(dependencies), A.canonical_object(SETTINGS), A.canonical_object(evaluation), split.seed,
        'FN-EC-PF-KNN-01-pinned-code', fit_entities=split.entities('train'), fit_role='train', benchmark_id=split.benchmark_id,
        confidence_kind='method_support', evaluation_partition='test', gaps=limitations, status='ok' if model_available else 'unavailable')


def freeze(output):
    """Execute one candidate with exact source/input/native/artifact serialization replay."""
    import pandas as pd
    from starplast import artifacts as A, strategies as S
    from starplast.splits import read_split, write_split
    output = Path(output)
    P._write(output / 'predeclared_scope.json', {'pilot_id': PILOT, 'target': C.TARGET, 'settings': SETTINGS,
        'preparation_manifest_sha256': PREPARATION_SHA, 'biological_admission': False})
    packet = load_packet()
    nodes = pd.read_parquet(packet['source_path'])
    assert P._sha(packet['source_path']) == packet['source_hash']
    assert C._records(nodes[list(pd.read_parquet(PREPARATION / 'source_cells.parquet').columns)]) == C._records(pd.read_parquet(PREPARATION / 'source_cells.parquet'))
    native_groups = S.Context(nodes, graph={}, organism=_ORGANISMS.FALCIPARUM).groups()
    assert native_groups.tolist() == packet['groups'].group.tolist()
    features, labels, excluded, provenance = source_closed_inputs(nodes, packet)
    assert len(features) == 1050 and len(labels) == 600
    P._write(output / 'feature_scope.json', {'provenance': provenance, 'exclusions': asdict(excluded),
        'target_identity': packet['target_identity'], 'cohort_identity': packet['cohort_identity'], 'split_identity': packet['split'].identity})
    result = fit_candidate(features, labels, packet, excluded)
    assert len(result['rows']) == 152 and int((~result['rows'].training_supported).sum()) == 4
    repeated = fit_candidate(features, labels, packet, excluded)
    assert N._json(result['batch'].model_state) == N._json(repeated['batch'].model_state)
    for key in ('raw_votes', 'card', 'view_snapshot', 'profile_class_cards', 'major_class_cards', 'baseline_comparisons'):
        assert result[key] == repeated[key], key
    for name, frame in {'features': features, 'rows': result['rows'], 'class_scores': result['class_scores'],
                        'native_class_scores': result['batch'].class_scores}.items():
        frame.to_parquet(output / (name + '.parquet'))
        pd.testing.assert_frame_equal(frame, pd.read_parquet(output / (name + '.parquet')), check_exact=True)
    for key in ('rows', 'class_scores'):
        pd.testing.assert_frame_equal(result[key], repeated[key], check_exact=True)
    pd.testing.assert_frame_equal(result['batch'].class_scores, repeated['batch'].class_scores, check_exact=True)
    write_split(output / 'split.json', packet['split'], excluded)
    assert read_split(output / 'split.json') == (packet['split'], excluded)
    code = [Path(__file__), ROOT / 'tests/test_pf_functional_ec_knn.py', Path(N.__file__), Path(C.__file__), Path(P.__file__),
        ROOT / 'scripts/notebook_runner.py', *[ROOT / 'starplast' / name for name in ('functional_exclusions.py', 'functional_ontology.py',
        'label_records.py', 'datasets.py', 'strategies.py', 'search.py', 'embedding.py', 'slots.py', 'record_scorecards.py', 'scorecard.py',
        'scorecard_view.py', 'artifacts.py', 'capabilities.py', 'splits.py', 'query.py', 'ground_truth.py', 'provenance.py', 'organisms.py')]]
    dependencies = [A.Dependency('code', str(path.relative_to(ROOT)), P._sha(path)) for path in code]
    dependencies.extend((A.Dependency('table', 'installed_nodes', packet['source_hash']),
        A.Dependency('source', 'prepared_EC_packet', PREPARATION_SHA),
        A.Dependency('truth', 'complete_direct_profiles', P._sha(PREPARATION / 'direct_profiles.parquet')),
        A.Dependency('split', 'prepared_protected', packet['split'].identity),
        A.Dependency('exclusions', 'training_full_source_closure', P._identity(asdict(excluded))),
        A.Dependency('model', 'native_training_state', P._identity(N._json(result['batch'].model_state)))))
    spec = artifact_spec(packet, dependencies, gaps=provenance['gaps'], model_available=result['model_available'])
    payload = {key + '.json': result[key] for key in ('raw_votes', 'card', 'view_snapshot', 'profile_class_cards', 'major_class_cards', 'baseline_comparisons')}
    payload.update({'model_state.json': N._json(result['batch'].model_state), 'exclusions.json': asdict(excluded), 'feature_provenance.json': provenance,
        'rows.json': result['rows'].astype(object).where(result['rows'].notna(), None).to_dict('records'),
        'native_class_scores.json': N._score_payload(result['batch'].class_scores), 'class_scores.json': N._score_payload(result['class_scores']),
        'baseline_cards.json': packet['controls']})
    payload = C._plain(payload)
    identity = A.write_artifact(output / 'held_out', spec, payload, split=packet['split'])
    assert A.read_artifact(output / 'held_out', expected=spec, split=packet['split']).payloads == payload
    receipts = dict(packet['receipts'])
    (output / 'code').mkdir()
    for path in code:
        receipts[str(path)] = P._sha(path)
        shutil.copyfile(path, output / 'code' / path.name)
    for path, digest in receipts.items():
        assert P._sha(path) == digest, 'Changed input: ' + path
    P._write(output / 'input_manifest.json', receipts)
    summary = {'pilot_id': PILOT, 'organism': _ORGANISMS.FALCIPARUM, 'target': C.TARGET, 'whole_genes': len(nodes),
        'eligible_profiles': len(features), 'train': len(labels), 'tune': len(packet['split'].entities('tune')),
        'calibration': len(packet['split'].entities('calibration')), 'test': len(result['rows']),
        'source_annotated_genes': int(packet['profiles'].status.ne('unannotated').sum()),
        'unknown_genes': int(packet['profiles'].status.eq('unannotated').sum()),
        'unresolved_annotated_genes': int(packet['profiles'].status.eq('unresolved').sum()),
        'unsupported_test_genes': 4, 'unknown_statuses': result['card']['extra']['unknown_states'],
        'selected_features': features.columns.tolist(), 'kept_features': result['batch'].model_state['kept_columns'],
        'native_classes': len(result['batch'].class_scores.columns), 'recorded_source_classes': len(result['class_scores'].columns),
        'counts': result['card']['counts'], 'metrics': result['card']['metrics'], 'baseline_comparisons': result['baseline_comparisons'],
        'artifact_identity': identity, 'artifact_status': spec.status, 'target_identity': packet['target_identity'],
        'cohort_identity': packet['cohort_identity'], 'split_identity': packet['split'].identity,
        'source_sha256': packet['source_hash'], 'preparation_manifest_sha256': PREPARATION_SHA,
        'input_receipt_count': len(receipts), 'exact_source_native_refit_and_serialization_replay': True,
        'biological_accuracy': None, 'biological_admission': False, 'calibration_status': 'unavailable', 'deployment': 'unknown', 'gaps': provenance['gaps']}
    P._write(output / 'summary.json', summary)
    (output / 'README.md').write_text('Fixed native Pf direct EC-major complete-profile recovery candidate. Accepted source/target/split '
        'and original controls remain unchanged; every unsupported profile and abstention stays in the test cohort. Numeric selection '
        'and native ranks use training only, with full-schema functional/orthology/attention closure and transitive late vetoes. '
        'Native scores/votes are uncalibrated method support, separate from full-source reporting. Exact source, native, refit, '
        'serialization and typed-artifact replay is recorded. No graph/map/tuning/calibration/deployment or independent biological '
        'truth admission occurs; source/feature assignment lineage remains unresolved.\n')
    return summary


def main():
    """Execute once in a new immutable notebook packet and retain every failure."""
    from scripts.notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    notebook = ExecutedNotebook(PILOT)
    notebook.ns.update(freeze=freeze, output=args.out)
    notebook.md('Fixed k15/min_share0.3 on unchanged accepted direct Pf EC target/protected split. Source recovery only; no biological admission.')
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
