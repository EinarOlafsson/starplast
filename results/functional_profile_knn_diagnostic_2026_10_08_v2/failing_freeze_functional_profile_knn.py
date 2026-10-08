"""Freeze one source-closed native kNN annotation-profile recovery candidate.

The verified complete-Pfam preparation fixes every gene, protected group and
split before fitting. Numeric feature eligibility, distributions and ranks use
training rows only. Native complete-profile votes remain atomic, unsupported
test profiles remain in denominators, and raw support remains uncalibrated.
This candidate neither admits biological function accuracy nor supplies a
deployment model. Source closure is conservative declared/empirical protection,
not proof that incomplete annotation and feature lineage are independent.
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
from scripts import freeze_functional_profile_controls as P

PILOT_ID = 'FN-DOM-KNN-01-Tg-pfam-complete-profile'
SETTINGS = {'k': 15, 'min_share': 0.3, 'selection': 'Predeclared; no tune/test selection',
            'feature_transform': 'native training average ranks; frozen query ECDF',
            'feature_selection': 'registered numeric columns; native training coverage/variation floor',
            'graphs_or_maps_used': False, 'pilot_id': PILOT_ID}


def _json(value):
    import numpy as np
    return json.loads(json.dumps(value, allow_nan=False, default=lambda item: item.item()
                               if isinstance(item, np.generic) else _unsupported(item)))


def _unsupported(value):
    raise ValueError('Unsupported finite result scalar: ' + type(value).__name__)


def _score_payload(frame):
    return {'classes': frame.columns.tolist(),
            'scores': frame.astype(object).where(frame.notna(), None).to_numpy().tolist()}


def prepare_inputs(nodes, contract, *, derived_inputs=()):
    """Freeze registered numeric inputs and source closure from exact training rows."""
    import pandas as pd
    from starplast import datasets as D, functional_exclusions as E, strategies as S
    from starplast.functional_profile_targets import ProfileTargetContract

    if not isinstance(contract, ProfileTargetContract) or contract.organism != 'Tg' or contract.source_target != 'pfam_id':
        raise ValueError('A typed Tg complete-Pfam preparation is required')
    if (not isinstance(nodes, pd.DataFrame) or nodes.columns.has_duplicates or 'gene_id' not in nodes or
            tuple(nodes.gene_id) != contract.universe or contract.target in nodes):
        raise ValueError('Exact original whole-gene table/order and an unused encoded target are required')
    if contract.source_target not in nodes:
        raise ValueError('Original functional source column is missing')
    indexed = nodes.set_index('gene_id', drop=False)
    train = indexed.loc[list(contract.split.entities('train'))].copy()
    train[contract.target] = contract.training_labels()
    excluded = E.make_exclusions(train.reset_index(drop=True), contract.target,
        source_targets=(contract.source_target,), benchmark_id=contract.split.benchmark_id,
        split=contract.split, derived_inputs=derived_inputs)
    # Unknown registered derivation parents cannot be interpreted as independent.
    unregistered, unresolved_derived, metadata = [], [], {}
    for column in nodes:
        source = D.provenance(column, contract.organism)
        if source is None:
            unregistered.append(column)
        elif any(parent not in nodes for parent in source.derived_from):
            unresolved_derived.append(column)
        else:
            metadata[column] = {'source_id': source.key, 'declared_derivation': list(source.derived_from),
                                'source_release': 'unresolved', 'evidence_grade': 'unresolved',
                                'lineage_status': 'Registry assertion; independent assignment lineage unverified'}
    banned = set(excluded.columns) | set(unregistered) | set(unresolved_derived)
    dependency_reasons = {}
    parents = {column: tuple(source['declared_derivation']) for column, source in metadata.items()}
    for item in derived_inputs:
        if item.kind == 'column':
            parents[item.name] = tuple(set(parents.get(item.name, ())) |
                {name for kind, name in item.parents if kind == 'column'})
    while True:
        added = set(D.derived_dependents(banned, contract.organism)) - banned
        for column in added:
            dependency_reasons.setdefault(column, {'reason': 'Organism-scoped registered dependency closure',
                'vetoed_parents': sorted(set(D.derived_sources(column, contract.organism)) & banned)})
        for column, dependencies in parents.items():
            vetoed = sorted(set(dependencies) & banned)
            if column not in banned and vetoed:
                added.add(column)
                dependency_reasons[column] = {'reason': 'Declared descendant of vetoed input', 'vetoed_parents': vetoed}
        if not added:
            break
        banned.update(added)
    excluded = replace(excluded, columns=tuple(sorted(banned)))
    columns = S.Context(train.reset_index(drop=True), graph={}, organism=contract.organism).numeric_columns(exclude=excluded.columns)
    features = indexed.loc[[row.entity for row in contract.split.assignments], columns].copy()
    provenance = {'policy': E.POLICY_VERSION, 'selection_fit_entities': list(contract.split.entities('train')),
                  'selected_columns': columns, 'sources': {column: metadata[column] for column in columns},
                  'unregistered_columns_withheld': unregistered, 'unresolved_derivations_withheld': unresolved_derived,
                  'dependency_closure_withheld': dependency_reasons,
                  'gaps': ['Registry and empirical closure do not prove complete source/homology independence',
                           'Gene annotation assignment release and feature evidence grades remain unresolved']}
    return features, excluded, provenance


def native_replay(features, labels, split, model):
    """Independently reconstruct training ranks, native votes and raw neighbor weights."""
    import numpy as np
    import pandas as pd
    from sklearn.neighbors import NearestNeighbors
    from starplast import strategies as S

    ids = tuple(a.entity for a in split.assignments)
    if tuple(features.index) != ids or tuple(labels.index) != split.entities('train'):
        raise ValueError('Replay requires exact split order and training-only labels')
    matrix = np.zeros((len(ids), len(model['kept_columns'])), dtype=float)
    for position, column in enumerate(model['kept_columns']):
        ordered = np.asarray(model['training_distributions'][column], dtype=float)
        for row, value in enumerate(features[column].to_numpy(dtype=float)):
            if np.isfinite(value):
                equal = np.flatnonzero(ordered == value)
                rank = (equal[0] + equal[-1] + 2) / 2 if len(equal) else int((ordered <= value).sum())
                matrix[row, position] = rank / len(ordered) - 0.5
    at = {entity: position for position, entity in enumerate(ids)}
    train = np.asarray([at[entity] for entity in labels.index], dtype=int)
    query = np.asarray([at[entity] for entity in split.entities('test')], dtype=int)
    expected = (features.loc[list(labels.index), model['kept_columns']].rank(pct=True) - 0.5).fillna(0).to_numpy()
    np.testing.assert_array_equal(matrix[train], expected)
    np.testing.assert_array_equal(matrix[train], model['training_vectors'])
    visible = pd.Series([None] * len(ids), dtype=object)
    visible.iloc[train] = labels.to_numpy()
    raw_calls, support = S.knn_vote(matrix, visible, model['k'], query=query)
    native_scores = S.class_scores_of(raw_calls)
    native_scores = native_scores.iloc[query].copy() if native_scores is not None else pd.DataFrame(index=query)
    native_scores.index = list(split.entities('test'))
    raw_votes = []
    if len(train) >= 2 and matrix.shape[1] and len(query):
        neighbors = NearestNeighbors(n_neighbors=min(model['k'] + 1, len(train))).fit(matrix[train])
        distances, positions = neighbors.kneighbors(matrix[query])
        for row, (entity, dist, pos) in enumerate(zip(split.entities('test'), distances, positions)):
            dist, pos = dist[:model['k']], pos[:model['k']]
            weights = 1.0 / (dist + 1e-6)
            votes = {}
            for label, weight in zip(labels.to_numpy()[pos], weights):
                votes[label] = votes.get(label, 0.0) + weight
            winner, total = max(votes, key=votes.get), sum(votes.values())
            assert raw_calls.iloc[query[row]] == winner and support.iloc[query[row]] == votes[winner] / total
            np.testing.assert_array_equal(native_scores.iloc[row].to_numpy(),
                np.array([votes.get(label, 0.0) / total for label in native_scores.columns]))
            raw_votes.append({'entity': entity, 'neighbors': [str(labels.index[position]) for position in pos],
                              'distances': dist.tolist(), 'weights': weights.tolist(), 'votes': votes,
                              'total_weight': total, 'raw_call': winner, 'raw_support': votes[winner] / total})
    else:
        raw_votes = [{'entity': entity, 'neighbors': [], 'distances': [], 'weights': [], 'votes': {},
                      'total_weight': None, 'raw_call': None, 'raw_support': None} for entity in split.entities('test')]
    calls = raw_calls.where(support >= model['min_share']).iloc[query].astype(object)
    calls = calls.where(calls.notna(), None)
    return calls, support.iloc[query], native_scores, _json(raw_votes)


def fit_recovery(features, contract, exclusions):
    """Fit fixed native kNN and retain every test profile, abstention and class score."""
    import numpy as np
    import pandas as pd
    from starplast import artifacts as A, label_records as L, record_scorecards as R, scorecard as SC
    from starplast.scorecard_view import build_scorecard_view

    labels = contract.training_labels()
    batch = L.feature_knn(features, labels, split=contract.split, exclusions=exclusions,
                          k=SETTINGS['k'], min_share=SETTINGS['min_share'])
    calls, support, scores, votes = native_replay(features, labels, contract.split, batch.model_state)
    assert batch.rows.prediction.tolist() == calls.tolist()
    np.testing.assert_array_equal(batch.rows.support.to_numpy(dtype=float), support.to_numpy(dtype=float))
    pd.testing.assert_frame_equal(batch.class_scores, scores, check_exact=True)
    rows = batch.rows.copy()
    by_gene = {row.gene_id: row for row in contract.rows}
    rows['truth'] = [by_gene[entity].profile for entity in rows.entity]
    rows['group'] = [by_gene[entity].protected_group for entity in rows.entity]
    rows['training_supported'] = rows.truth.isin(labels)
    # Native classes stay untouched in batch; reporting explicitly supplies zero
    # vote support for source profiles never seen in training, including all test truth.
    classes = sorted({row.profile for row in contract.rows if row.eligible})
    full_scores = scores.reindex(columns=classes, fill_value=0.0) if len(scores.columns) else pd.DataFrame(
        np.nan, index=scores.index, columns=classes)
    gaps = tuple((*batch.gaps, 'Recorded profile recovery only; independent function accuracy unknown',
                  'Assignment release, empirical/source closure completeness and homology independence unresolved',
                  'Calibration and deployment applicability are not established'))
    scope = R.RecordScope(contract.organism, 'feature_knn', contract.target, SC.T_LABEL,
        A.canonical_object(SETTINGS), contract.split.seed, contract.split.identity, 'outer_test',
        contract.split.benchmark_id, contract.source.evidence_grade, 'gene', contract.source.negative_semantics, gaps)
    card = R.aggregate(rows, scope, parameters={'class_scores': full_scores.reset_index(drop=True)})
    model_available = bool(batch.model_state['kept_columns']) and len(labels) >= 2
    card['extra'].update(target_identity=contract.identity, prepared_cohort_identity=contract.cohort_identity,
        unsupported_test_genes=int((~rows.training_supported).sum()), unknown_states=contract.capacity()['unknown_states'],
        biological_accuracy=None, biological_admission=False, calibrated=False, deployment_applicability='unknown',
        model_status='available_reference_recovery_candidate' if model_available else 'unavailable_no_supported_native_vote')
    view = build_scorecard_view(card, details={'source': {'grade': contract.source.evidence_grade,
        'name': list(contract.source.source_ids), 'version': contract.source.source_release,
        'negative_semantics': contract.source.negative_semantics}, 'gaps': list(gaps)})
    class_cards = R.class_cards(rows, scope, parameters={'class_scores': full_scores.reset_index(drop=True)})
    baseline_rows, baseline_scores, baseline_cards, _ = P.build_controls(contract, seed=P.SEED)
    assert rows.entity.tolist() == baseline_rows.entity.tolist() and rows.truth.tolist() == baseline_rows.truth.tolist()
    assert rows.training_supported.tolist() == baseline_rows.training_supported.tolist()
    comparisons = {name: {'same_entity_order': True, 'same_truth': True, 'same_class_order': classes == baseline_scores[name].columns.tolist(),
        'metrics': bundle['card']['metrics'], 'difference_from_knn': {key: card['metrics'][key] - value
        if card['metrics'][key] is not None and value is not None else None for key, value in bundle['card']['metrics'].items()}}
        for name, bundle in baseline_cards.items()}
    return {'batch': batch, 'rows': rows, 'class_scores': full_scores, 'raw_votes': votes, 'model_available': model_available,
            'card': _json(card), 'view_snapshot': json.loads(view.snapshot_json), 'class_cards': _json(class_cards),
            'baseline_cards': baseline_cards, 'baseline_comparisons': comparisons}


def freeze(output):
    """Execute one fixed candidate into a new result packet with byte-pinned inputs."""
    import numpy as np
    import pandas as pd
    from starplast import artifacts as A, capabilities as C, scorecard as SC
    from starplast.query import Query
    from starplast.splits import write_split, read_split

    output = Path(output)
    contract, source_scope, receipts = P.load_snapshot(P.SNAPSHOT)
    inputs = json.loads((P.SNAPSHOT / 'input_manifest.json').read_text())
    node_paths = [path for path, digest in inputs.items() if digest == source_scope['installed_source_sha256']]
    assert len(node_paths) == 1 and P._sha(node_paths[0]) == source_scope['installed_source_sha256']
    receipts[node_paths[0]] = source_scope['installed_source_sha256']
    census = Path(source_scope['source_census_path'])
    assert P._sha(census) == source_scope['source_census_sha256']
    receipts[str(census)] = source_scope['source_census_sha256']
    nodes = pd.read_parquet(node_paths[0])
    features, exclusions, provenance = prepare_inputs(nodes, contract)
    assert len(features) == 4310 and len(contract.training_labels()) == 2381
    P._write(output / 'candidate_scope.json', {'pilot_id': PILOT_ID, 'settings': SETTINGS,
        'target_identity': contract.identity, 'cohort_identity': contract.cohort_identity,
        'split_identity': contract.split.identity, 'prepared_manifest_sha256': P.SNAPSHOT_SHA256,
        'source_scope': source_scope, 'feature_provenance': provenance, 'exclusions': asdict(exclusions)})
    result = fit_recovery(features, contract, exclusions)
    assert len(result['rows']) == 635 and int((~result['rows'].training_supported).sum()) == 333
    replay = fit_recovery(features, contract, exclusions)
    for key in ('card', 'class_cards', 'raw_votes', 'baseline_cards', 'baseline_comparisons', 'view_snapshot'):
        assert result[key] == replay[key], key
    assert _json(result['batch'].model_state) == _json(replay['batch'].model_state)
    for key in ('rows', 'class_scores'):
        pd.testing.assert_frame_equal(result[key], replay[key], check_exact=True)
        result[key].to_parquet(output / (key + '.parquet'))
        pd.testing.assert_frame_equal(result[key], pd.read_parquet(output / (key + '.parquet')), check_exact=True)
    features.to_parquet(output / 'features.parquet')
    pd.testing.assert_frame_equal(features, pd.read_parquet(output / 'features.parquet'), check_exact=True)
    result['batch'].class_scores.to_parquet(output / 'native_class_scores.parquet')
    pd.testing.assert_frame_equal(result['batch'].class_scores, pd.read_parquet(output / 'native_class_scores.parquet'), check_exact=True)
    write_split(output / 'split.json', contract.split, exclusions)
    assert read_split(output / 'split.json') == (contract.split, exclusions)
    code = [Path(__file__), ROOT / 'tests/test_functional_profile_knn.py', Path(P.__file__), ROOT / 'scripts/notebook_runner.py',
            *[ROOT / 'starplast' / name for name in ('label_records.py', 'functional_exclusions.py', 'functional_profile_targets.py',
            'functional_domain_profiles.py', 'strategies.py', 'datasets.py', 'search.py', 'slots.py', 'embedding.py',
            'scorecard.py', 'record_scorecards.py', 'scorecard_view.py', 'baselines.py', 'splits.py', 'artifacts.py',
            'capabilities.py', 'query.py', 'ground_truth.py', 'provenance.py', 'organisms.py')]]
    dependencies = [A.Dependency('code', str(path.relative_to(ROOT)), P._sha(path)) for path in code]
    dependencies.extend((A.Dependency('table', 'installed_nodes', receipts[node_paths[0]]),
        A.Dependency('truth', 'normalized_complete_profiles', P._sha(P.SNAPSHOT / 'normalized_target.json')),
        A.Dependency('source', 'original_cell_census', receipts[str(census)]),
        A.Dependency('source', 'prepared_target_manifest', P.SNAPSHOT_SHA256),
        A.Dependency('split', 'prepared_protected', contract.split.identity),
        A.Dependency('exclusions', 'training_source_closure', P._identity(asdict(exclusions))),
        A.Dependency('model', 'native_training_state', P._identity(_json(result['batch'].model_state)))))
    evaluation = {'pilot_id': PILOT_ID, 'truth_grade': contract.source.evidence_grade,
        'target_identity': contract.identity, 'cohort_identity': contract.cohort_identity,
        'negative_semantics': contract.source.negative_semantics, 'biological_admission': False,
        'context': 'Tg ME49 recorded complete-Pfam annotation; assignment assay/stage unresolved',
        'limitations': list(provenance['gaps']), 'calibration': 'unavailable', 'deployment': 'unknown'}
    spec = A.ArtifactSpec(Query(contract.organism, 'label', target=contract.target).to_json(), 'feature_knn',
        SC.T_LABEL, C.get('feature_knn').outputs[0], contract.target, tuple(result['rows'].entity), 'held_out',
        'outer:test:' + contract.split.identity, tuple(dependencies), A.canonical_object(SETTINGS),
        A.canonical_object(evaluation), contract.split.seed, 'FN-DOM-KNN-01-pinned-code',
        fit_entities=contract.split.entities('train'), fit_role='train', benchmark_id=contract.split.benchmark_id,
        confidence_kind='method_support', evaluation_partition='test', gaps=tuple(provenance['gaps']),
        status='ok' if result['model_available'] else 'unavailable')
    payloads = {key + '.json': result[key] for key in ('raw_votes', 'card', 'class_cards', 'view_snapshot', 'baseline_cards', 'baseline_comparisons')}
    payloads.update({'model_state.json': _json(result['batch'].model_state), 'exclusions.json': asdict(exclusions),
        'feature_provenance.json': provenance, 'rows.json': result['rows'].astype(object).where(result['rows'].notna(), None).to_dict('records'),
        'native_class_scores.json': _score_payload(result['batch'].class_scores),
        'class_scores.json': _score_payload(result['class_scores'])})
    payloads = _json(payloads)
    artifact_identity = A.write_artifact(output / 'held_out', spec, payloads, split=contract.split)
    assert A.read_artifact(output / 'held_out', expected=spec, split=contract.split).payloads == payloads
    (output / 'code').mkdir()
    for path in code:
        receipts[str(path)] = P._sha(path)
        shutil.copyfile(path, output / 'code' / (('scripts_' if path.parent.name == 'scripts' else '') + path.name))
    for path, digest in receipts.items():
        assert P._sha(path) == digest, 'Changed input: ' + path
    P._write(output / 'input_manifest.json', receipts)
    summary = {'pilot_id': PILOT_ID, 'whole_genes': len(contract.rows), 'eligible_genes': len(features),
        'train_genes': 2381, 'test_genes': 635, 'unsupported_test_genes': 333, 'unknown_genes': 3830,
        'selected_feature_columns': features.columns.tolist(), 'kept_feature_columns': result['batch'].model_state['kept_columns'],
        'native_training_classes': len(result['batch'].class_scores.columns), 'recorded_source_classes': len(result['class_scores'].columns),
        'metrics': result['card']['metrics'], 'counts': result['card']['counts'], 'baseline_comparisons': result['baseline_comparisons'],
        'artifact_identity': artifact_identity, 'artifact_status': spec.status,
        'target_identity': contract.identity, 'cohort_identity': contract.cohort_identity,
        'split_identity': contract.split.identity, 'source_sha256': source_scope['installed_source_sha256'],
        'source_cell_census_sha256': source_scope['source_census_sha256'], 'exact_native_and_serialization_replay': True,
        'biological_accuracy': None, 'biological_admission': False, 'calibration': 'unavailable', 'deployment': 'unknown',
        'software_versions': {'python': sys.version.split()[0], 'numpy': np.__version__, 'pandas': pd.__version__}}
    P._write(output / 'summary.json', summary)
    (output / 'README.md').write_text('Complete-Pfam annotation-profile recovery candidate, fixed k=15/min_share=0.3. '
        'All original eligible/test profiles, unsupported profiles and unknown states remain. Native ranks and feature selection '
        'are training-only; graphs/maps and unregistered or unresolved derived inputs are withheld. Native class scores and '
        'raw neighbor weights remain method support, with separate full-source-class reporting scores. Exact native and '
        'matched baseline replay is recorded. Independent biological accuracy, calibration, source/homology independence '
        'and deployment applicability are unknown; registry source assertions are not verified measurement lineage.\n')
    return summary


def main():
    """Execute the fixed candidate as an immutable notebook packet; retain failures."""
    from scripts.notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    notebook = ExecutedNotebook(PILOT_ID)
    notebook.ns.update(output=args.out, freeze=freeze)
    notebook.md('Predeclared complete-Pfam reference-profile recovery candidate; no tune/test selection or biological admission.')
    try:
        notebook.code('summary = freeze(output)', 'summary')
        print(json.dumps(notebook.ns['summary'], indent=2, allow_nan=False))
    except Exception as exc:
        P._write(args.out / 'failure.json', {'error': type(exc).__name__, 'detail': str(exc)})
        notebook.md('Execution failed: ' + type(exc).__name__ + ': ' + str(exc))
        raise
    finally:
        notebook.write(str(args.out / 'executed.ipynb'))
        (args.out / 'SHA256SUMS.txt').write_text('\n'.join(f'{P._sha(path)}  {path.relative_to(args.out)}'
            for path in sorted(args.out.rglob('*')) if path.is_file() and path.name != 'SHA256SUMS.txt') + '\n')


if __name__ == '__main__':
    main()
