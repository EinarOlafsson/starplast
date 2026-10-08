"""Freeze data-only artifacts for six task meanings and distinct evaluation roles."""
from __future__ import annotations

import argparse
from dataclasses import asdict, replace
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from starplast import artifacts as A, capabilities as C, organisms as O  # noqa: E402
from starplast.query import EntityRef, Query  # noqa: E402
from starplast.splits import Assignment, SplitManifest  # noqa: E402


def freeze(output):
    """Create an immutable synthetic round-trip/invalidation pilot, without fitting."""
    output = Path(output)
    if output.exists():
        raise ValueError('Use a new artifact-pilot directory')
    output.mkdir(parents=True)
    ids = tuple(f'TGME49_{i:06d}' for i in range(100001, 100009))
    split = SplitManifest(O.TOXOPLASMA, 'synthetic_artifact_fixture', 'entity', tuple(
        Assignment(i, i, role) for role, cohort in zip(('train', 'tune', 'calibration', 'test'), (ids[:2], ids[2:4], ids[4:6], ids[6:])) for i in cohort), 17)
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    fixture_hash = hashlib.sha256(b'synthetic_artifact_inputs_v1').hexdigest()
    code_version = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip() + '+artifact-contract-worktree'
    code_hash = sha(ROOT / 'starplast/artifacts.py')
    split_payload = asdict(split)
    (output / 'split.json').write_text(json.dumps(split_payload, indent=2) + '\n')
    loaded, reports = [], []
    cases = (('feature_knn', 'label', 'label_calls'), ('trait_regression', 'trait', 'numeric_estimates'),
             ('positive_unlabeled', 'gene_set', 'gene_ranking'), ('geneset_hunt', 'gene_set', 'gene_set'),
             ('holdout_search', 'label', 'map_modules'), ('blind_battery', 'gene', 'replication_findings'))
    for key, query_kind, output_kind in cases:
        cap = C.get(key)
        target = 'synthetic_value' if query_kind == 'trait' else 'synthetic_target'
        entities = tuple(EntityRef(O.TOXOPLASMA, 'gene', i) for i in ids[:5]) if query_kind == 'gene_set' else (EntityRef(O.TOXOPLASMA, 'gene', ids[0]),) if query_kind == 'gene' else ()
        query = Query(O.TOXOPLASMA, query_kind, entities, target=target)
        scope = {'unit': cap.benchmark_unit, 'eligible_population': 2, 'truth_grade': 'synthetic_control',
                 'negative_semantics': 'declared fixture outcomes only; not biological absence', 'context': {},
                 'limitations': ['schema control only; not biological accuracy'], 'quantity_unit': 'synthetic unit' if query_kind == 'trait' else 'not applicable'}
        deps = tuple(A.Dependency(kind, name, digest) for kind, name, digest in (
            ('table', 'synthetic_universe', fixture_hash), ('code', 'artifact_adapter', code_hash),
            ('truth', target, fixture_hash), ('split', 'nested', split.identity), ('exclusions', 'explicit_fixture_closure', fixture_hash)))
        spec = A.ArtifactSpec(query.to_json(), key, cap.benchmark_tasks[0], next(o for o in cap.outputs if o.kind == output_kind),
            target, split.entities('test'), 'held_out', 'outer:test:' + split.identity,
            deps, A.canonical_object({'fixture': 'no runner invoked'}), A.canonical_object(scope), 17, code_version,
            fit_entities=split.entities('train'), fit_role='train', benchmark_id=split.benchmark_id,
            confidence_kind='method_support', evaluation_partition='test', gaps=('No independent biological benchmark admitted',))
        payloads = {'rows.json': [{'entity': ids[6], 'status': 'answered', 'support': 0.6, 'calibrated_confidence': None},
                                 {'entity': ids[7], 'status': 'abstained', 'support': None, 'calibrated_confidence': None}],
                    'interpretation.json': {'fixture_only': True, 'output': asdict(spec.output), 'fit_performed': False}}
        identity = A.write_artifact(output / key, spec, payloads, split=split)
        restored = A.read_artifact(output / key, expected=spec, split=split)
        assert restored.payloads == payloads and restored.identity == identity and restored.spec == spec
        loaded.append(restored)
        reports.append({'strategy': key, 'task': spec.task, 'output_kind': output_kind, 'role': spec.role, 'identity': identity, 'roundtrip': 'exact'})
    baseline = loaded[0].spec
    base = replace(baseline, role='reusable_base', evaluation_partition='none', fit_entities=(), fit_role='none',
        dependencies=tuple(d for d in baseline.dependencies if d.kind not in {'split', 'truth'}), benchmark_id='', partition_id='reusable:fixture')
    A.write_artifact(output / 'reusable_base', base, {'data.json': {'fixture_only': True}})
    base_artifact = A.read_artifact(output / 'reusable_base', expected=base)
    model = replace(baseline, role='fitted_model', evaluation_partition='none', entity_order=split.entities('train'),
        partition_id='fit:train:' + split.identity, dependencies=baseline.dependencies + (A.Dependency('artifact', 'base', base_artifact.identity),))
    A.write_artifact(output / 'model_metadata', model, {'model.json': {'fixture_only': True, 'fit_performed': False}}, split=split)
    model_artifact = A.read_artifact(output / 'model_metadata', expected=model, split=split)
    deployment = replace(baseline, role='deployment', evaluation_partition='deployment', entity_order=('synthetic_unknown_gene',),
        partition_id='deployment:fixture', dependencies=baseline.dependencies + (A.Dependency('model', 'fixture_model', model_artifact.identity),))
    A.write_artifact(output / 'deployment', deployment, {'rows.json': [{'entity': 'synthetic_unknown_gene', 'support': None, 'status': 'unavailable'}]}, split=split)
    deployment_artifact = A.read_artifact(output / 'deployment', expected=deployment, split=split)
    card = replace(baseline, role='scorecard', fit_entities=(), fit_role='none',
        dependencies=baseline.dependencies + (A.Dependency('artifact', 'held_out_rows', loaded[0].identity),))
    A.write_artifact(output / 'scorecard', card, {'card.json': {'fixture_only': True, 'accuracy': None, 'limitations': ['metrics aggregation follows under 64.09']}}, split=split)
    card_artifact = A.read_artifact(output / 'scorecard', expected=card, split=split)
    graph = [base_artifact, model_artifact, deployment_artifact, card_artifact, *loaded]
    changed = hashlib.sha256(b'changed fixture').hexdigest()
    stale = A.invalidated(graph, {('code', 'artifact_adapter'): changed})
    assert len(stale) == len(graph)
    unchanged = A.invalidated(graph, {('code', 'artifact_adapter'): code_hash})
    assert not unchanged
    refusals = []
    for name, kwargs in (('entity_order', {'entity_order': tuple(reversed(baseline.entity_order))}),
                         ('settings', {'settings_json': A.canonical_object({'fixture': 'changed'})}),
                         ('partition', {'partition_id': 'another:test'}), ('model_population', {'fit_entities': (ids[0],)})):
        try:
            A.read_artifact(output / 'feature_knn', expected=replace(baseline, **kwargs), split=split)
        except ValueError as error:
            refusals.append({'change': name, 'refused': str(error)})
        else:
            raise AssertionError('Stale artifact accepted: ' + name)
    summary = {'scorecard_tasks': len(reports), 'roundtrips': len(graph), 'distinct_roles': sorted({a.spec.role for a in graph}),
               'changed_code_invalidates': len(stale), 'unchanged_inputs_invalidate': len(unchanged),
               'deliberate_scope_refusals': refusals, 'biological_accuracy': 'not measured', 'fitted_models': 0}
    (output / 'report.json').write_text(json.dumps({'summary': summary, 'tasks': reports}, indent=2) + '\n')
    inputs = [Path(__file__), ROOT / 'starplast/artifacts.py', ROOT / 'starplast/capabilities.py', ROOT / 'starplast/splits.py', ROOT / 'starplast/query.py']
    files = [p for p in output.rglob('*.json')]
    (output / 'manifest.json').write_text(json.dumps({'created_utc': datetime.now(timezone.utc).isoformat(), 'summary': summary,
        'input_sha256': {str(p): sha(p) for p in inputs}, 'outputs': {str(p.relative_to(output)): sha(p) for p in files}}, indent=2) + '\n')
    return summary


def main():
    """Execute an annotated data-only artifact pilot without altering runtime outputs."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    nb = ExecutedNotebook('Artifact roles, exact round trips and cache invalidation')
    nb.md('Frozen scope: six scorecard tasks, held-out/deployment/model-metadata/base/scorecard roles, finite JSON payloads, ordered populations, nested split guards and declared dependency invalidation. All inputs and outputs are synthetic schema fixtures; no model is fitted and no biological accuracy is measured.')
    nb.code('from scripts.freeze_artifact_pilot import freeze', f'summary = freeze({str(args.out)!r})', 'summary')
    nb.md('Method support and calibrated confidence remain separate. A model metadata placeholder does not represent a fitted estimator. Unavailable biological validation and later task-specific row adapters remain explicit.')
    nb.write(str(args.out / 'pilot.ipynb'))
    print(nb.ns['summary'])


if __name__ == '__main__':
    main()
