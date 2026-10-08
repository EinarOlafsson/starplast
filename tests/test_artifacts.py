"""Artifact role separation, stale-input refusals and dependency propagation."""
from dataclasses import replace
import hashlib
import json

import pytest

from starplast import artifacts as A, capabilities as C, organisms as O
from starplast.query import Query
from starplast.splits import Assignment, SplitManifest

HASH = hashlib.sha256(b'fixture').hexdigest()
OTHER = hashlib.sha256(b'changed').hexdigest()
IDS = tuple(f'TGME49_{i:06d}' for i in range(100001, 100009))
SPLIT = SplitManifest(O.TOXOPLASMA, 'synthetic_fixture', 'entity', tuple(
    Assignment(i, i, role) for role, cohort in zip(('train', 'tune', 'calibration', 'test'), (IDS[:2], IDS[2:4], IDS[4:6], IDS[6:])) for i in cohort), 3)


def spec(key='feature_knn', *, role='held_out', target='location'):
    cap = C.get(key)
    query_kind = 'gene_set' if 'gene_set' in cap.query_kinds and 'label' not in cap.query_kinds else 'trait' if cap.benchmark_tasks[0] == 'values' else 'label' if 'label' in cap.query_kinds else 'gene'
    from starplast.query import EntityRef
    entities = tuple(EntityRef(O.TOXOPLASMA, 'gene', i) for i in IDS[:5]) if query_kind == 'gene_set' else (EntityRef(O.TOXOPLASMA, 'gene', IDS[0]),) if query_kind == 'gene' else ()
    query = Query(O.TOXOPLASMA, query_kind, entities, target=target)
    deps = [A.Dependency('code', 'adapter', HASH), A.Dependency('table', 'entities', HASH), A.Dependency('exclusions', 'explicit_empty_fixture', HASH)]
    if role != 'reusable_base':
        deps.append(A.Dependency('split', 'nested', SPLIT.identity))
    if role in {'held_out', 'scorecard'}:
        deps.append(A.Dependency('truth', target, HASH))
    if role == 'deployment':
        deps.append(A.Dependency('model', key, HASH))
    scope = {'unit': cap.benchmark_unit, 'eligible_population': 2, 'truth_grade': 'synthetic_control',
             'negative_semantics': 'fixture labels only', 'context': {}, 'limitations': ['not biological validation']}
    return A.ArtifactSpec(query.to_json(), key, cap.benchmark_tasks[0], cap.outputs[0], target,
        SPLIT.entities('test') if role in {'held_out', 'scorecard'} else ('TGME49_999999',) if role == 'deployment' else IDS[:2], role, 'fixture:' + role,
        tuple(deps), A.canonical_object({'k': 15}), A.canonical_object(scope), 4, 'synthetic-fixture-v1',
        fit_entities=IDS[:2] if role in {'fitted_model', 'deployment', 'held_out'} else (),
        fit_role='train' if role in {'fitted_model', 'deployment', 'held_out'} else 'none',
        benchmark_id='synthetic_fixture' if role != 'reusable_base' else '', confidence_kind='method_support',
        evaluation_partition='test' if role in {'held_out', 'scorecard'} else 'deployment' if role == 'deployment' else 'none')


@pytest.mark.parametrize('key', ['feature_knn', 'trait_regression', 'positive_unlabeled', 'geneset_hunt', 'holdout_search', 'blind_battery'])
def test_every_scorecard_task_roundtrips_without_losing_scope(tmp_path, key):
    current = spec(key)
    payload = {'rows.json': [{'entity': IDS[6], 'support': 0.6, 'calibrated_confidence': None},
                             {'entity': IDS[7], 'support': None, 'status': 'abstained'}]}
    identity = A.write_artifact(tmp_path / key, current, payload, split=SPLIT)
    restored = A.read_artifact(tmp_path / key, expected=current, split=SPLIT)
    assert restored.spec == current and restored.payloads == payload and restored.identity == identity
    assert restored.spec.role == 'held_out' and restored.spec.evaluation_partition == 'test'
    assert restored.spec.confidence_kind == 'method_support'
    with pytest.raises(FileExistsError):
        A.write_artifact(tmp_path / key, current, payload, split=SPLIT)


@pytest.mark.parametrize('field,value', [('entity_order', tuple(reversed(SPLIT.entities('test')))),
    ('settings_json', A.canonical_object({'k': 20})), ('seed', 9), ('partition_id', 'another:test'),
    ('target', 'different'), ('fit_entities', IDS[:1]), ('gaps', ('new uncertainty',))])
def test_changed_scope_never_reuses_cached_result(tmp_path, field, value):
    current = spec()
    A.write_artifact(tmp_path / 'result', current, {'rows.json': []}, split=SPLIT)
    try:
        changed = replace(current, **{field: value})
        with pytest.raises(ValueError):
            A.read_artifact(tmp_path / 'result', expected=changed, split=SPLIT)
    except ValueError:
        assert field == 'target'


@pytest.mark.parametrize('kind', ['table', 'code', 'truth', 'split', 'exclusions', 'model', 'graph', 'calibration', 'mapping', 'source'])
def test_every_input_hash_participates_in_identity(tmp_path, kind):
    current = spec()
    deps = tuple(d for d in current.dependencies if d.kind != kind) + (A.Dependency(kind, 'identity', HASH if kind != 'split' else SPLIT.identity),)
    current = replace(current, dependencies=deps)
    A.write_artifact(tmp_path / 'result', current, {'rows.json': []}, split=SPLIT)
    changed = replace(current, dependencies=tuple(replace(d, sha256=OTHER) if d.kind == kind else d for d in current.dependencies))
    with pytest.raises(ValueError):
        A.read_artifact(tmp_path / 'result', expected=changed, split=SPLIT)


def test_model_deployment_scorecard_and_held_out_roles_stay_distinct(tmp_path):
    for role in ('fitted_model', 'deployment', 'held_out', 'scorecard', 'reusable_base'):
        current = spec(role=role)
        A.write_artifact(tmp_path / role, current, {'data.json': {'role': role}}, split=None if role == 'reusable_base' else SPLIT)
    assert len({spec(role=r).identity for r in A.ROLES}) == len(A.ROLES)
    with pytest.raises(ValueError, match='must be distinct'):
        replace(spec(role='deployment'), evaluation_partition='test')
    with pytest.raises(ValueError, match='cannot be model fitting'):
        replace(spec(), fit_entities=SPLIT.entities('test'))
    with pytest.raises(ValueError, match='matching frozen split'):
        A.write_artifact(tmp_path / 'no_split', spec(), {'data.json': []})
    with pytest.raises(ValueError, match='Forbidden model_fit'):
        A.write_artifact(tmp_path / 'leak', replace(spec(role='fitted_model'), fit_entities=SPLIT.entities('test')), {'model.json': {}}, split=SPLIT)


def test_confidence_scope_and_unavailable_outputs_are_explicit():
    with pytest.raises(ValueError, match='calibration identity'):
        replace(spec(), confidence_kind='calibrated_probability')
    deps = spec().dependencies + (A.Dependency('calibration', 'probability_calibration', HASH),)
    calibrated = replace(spec(), dependencies=deps, confidence_kind='calibrated_probability', calibration_scope='fixture only')
    assert calibrated.identity != spec().identity
    with pytest.raises(ValueError, match='reinterpret geometry'):
        replace(spec('holdout_search'), dependencies=deps, confidence_kind='calibrated_probability', calibration_scope='fixture only')
    with pytest.raises(ValueError, match='explicit gaps'):
        replace(spec(), status='unavailable')
    assert replace(spec(), status='unavailable', gaps=('no truth',)).status == 'unavailable'


def test_payload_corruption_manifest_corruption_and_symlinks_refused(tmp_path):
    current = spec()
    for case in ('payload', 'manifest', 'symlink'):
        folder = tmp_path / case
        A.write_artifact(folder, current, {'data.json': {'value': 1}}, split=SPLIT)
        if case == 'payload':
            (folder / 'data.json').write_text('{"value":2}')
        elif case == 'manifest':
            raw = json.loads((folder / 'manifest.json').read_text())
            raw['spec']['seed'] += 1
            (folder / 'manifest.json').write_text(json.dumps(raw))
        else:
            data = (folder / 'data.json').read_bytes()
            (tmp_path / 'external.json').write_bytes(data)
            (folder / 'data.json').unlink()
            (folder / 'data.json').symlink_to(tmp_path / 'external.json')
        with pytest.raises(ValueError):
            A.read_artifact(folder, expected=current, split=SPLIT)
    with pytest.raises(ValueError, match='safe named'):
        A.write_artifact(tmp_path / 'escape', current, {'../unsafe.json': []}, split=SPLIT)
    with pytest.raises(ValueError):
        A.write_artifact(tmp_path / 'nonfinite', current, {'data.json': [float('nan')]}, split=SPLIT)
    assert not (tmp_path / 'nonfinite').exists()


def test_invalidation_propagates_only_through_declared_dependencies(tmp_path):
    base = spec(role='reusable_base')
    A.write_artifact(tmp_path / 'base', base, {'data.json': []})
    first = A.read_artifact(tmp_path / 'base', expected=base)
    child = replace(spec(), dependencies=spec().dependencies + (A.Dependency('artifact', 'base', first.identity),))
    A.write_artifact(tmp_path / 'child', child, {'data.json': []}, split=SPLIT)
    second = A.read_artifact(tmp_path / 'child', expected=child, split=SPLIT)
    detached = replace(base, dependencies=tuple(replace(d, name='unrelated_code') if d.kind == 'code' else d for d in base.dependencies))
    A.write_artifact(tmp_path / 'detached', detached, {'data.json': []})
    third = A.read_artifact(tmp_path / 'detached', expected=detached)
    stale = A.invalidated((first, second, third), {('code', 'adapter'): OTHER})
    assert stale == {first.identity, second.identity}
    assert not A.invalidated((first, second, third), {('code', 'adapter'): HASH})


def test_known_benchmark_entities_cannot_be_published_as_unknown_deployment(tmp_path):
    with pytest.raises(ValueError, match='Known benchmark entities'):
        A.write_artifact(tmp_path / 'false_deployment', replace(spec(role='deployment'), entity_order=IDS[:2]), {'rows.json': []}, split=SPLIT)


@pytest.mark.parametrize('dependency_kind', ['model', 'calibration'])
def test_model_and_calibration_artifact_descendants_invalidate(tmp_path, dependency_kind):
    base = spec(role='reusable_base')
    A.write_artifact(tmp_path / 'base', base, {'data.json': []})
    first = A.read_artifact(tmp_path / 'base', expected=base)
    child = replace(spec(), dependencies=spec().dependencies + (A.Dependency(dependency_kind, 'upstream', first.identity),))
    A.write_artifact(tmp_path / 'child', child, {'data.json': []}, split=SPLIT)
    second = A.read_artifact(tmp_path / 'child', expected=child, split=SPLIT)
    assert A.invalidated((first, second), {('code', 'adapter'): OTHER}) == {first.identity, second.identity}
