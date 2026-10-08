"""Tiny software fixtures; no fitting, installed data, Qt or biological claims."""
from dataclasses import replace
import fcntl
import hashlib
import json
from types import SimpleNamespace

import pytest

from starplast import artifacts as A, capabilities as C, precompute_jobs as P
from starplast.query import Query
from starplast.splits import Assignment, SplitManifest

HASH = hashlib.sha256(b'fixture').hexdigest()
OTHER = hashlib.sha256(b'changed').hexdigest()
LIMIT = P.Budget(20, 30, 1000000, 1000000, 400 * 1024 * 1024)
SPLIT = SplitManifest('Tg', 'software_fixture', 'entity', tuple(
    Assignment(str(i), str(i), role) for i, role in enumerate(('train', 'tune', 'calibration', 'test'))), 0)


@pytest.fixture(autouse=True)
def controlled_rss(monkeypatch):
    """Exercise budget logic without depending on earlier Qt/test process work.

    The actual focused test process still runs under an external MemoryMax400M
    scope. The budget rejection test explicitly replaces this simulated value.
    Return the real reader for its separate simulated /proc contract check.
    """
    reader = P._rss
    monkeypatch.setattr(P, '_rss', lambda: 50 * 1024 * 1024)
    return reader


def spec(*, role='reusable_base', table='entities'):
    cap = C.get('feature_knn')
    dependencies = (A.Dependency('code', 'fixture_builder', HASH), A.Dependency('table', table, HASH),
                    A.Dependency('exclusions', 'explicit_empty_fixture', HASH))
    if role != 'reusable_base':
        dependencies += (A.Dependency('split', 'nested', SPLIT.identity),)
    if role in {'held_out', 'scorecard'}:
        dependencies += (A.Dependency('truth', 'fixture_target', HASH),)
    if role == 'deployment':
        dependencies += (A.Dependency('model', 'fitted', HASH),)
    scope = {'unit': cap.benchmark_unit, 'eligible_population': 1, 'truth_grade': 'synthetic_control',
             'negative_semantics': 'fixture only', 'context': {}, 'limitations': ['software test only']}
    return A.ArtifactSpec(Query('Tg', 'label', target='fixture_target').to_json(), 'feature_knn',
        cap.benchmark_tasks[0], cap.outputs[0], 'fixture_target', ('3',) if role in {'held_out', 'scorecard'} else ('unknown',),
        role, 'fixture:' + role, dependencies, A.canonical_object({'setting': 1}), A.canonical_object(scope),
        0, 'software-fixture-v1', fit_entities=('0',) if role in {'fitted_model', 'deployment'} else (),
        fit_role='train' if role in {'fitted_model', 'deployment'} else 'none',
        benchmark_id='software_fixture' if role != 'reusable_base' else '',
        evaluation_partition='test' if role in {'held_out', 'scorecard'} else 'deployment' if role == 'deployment' else 'none')


def child(name='child', parent='base', *, base=None, kind='artifact', address='upstream'):
    current = base or spec()
    current = replace(current, dependencies=current.dependencies + (A.Dependency(kind, address, HASH),))
    return P.Job(name, current, (P.Upstream(parent, kind, address),), None if current.role == 'reusable_base' else SPLIT)


def graph():
    return P.Plan((child(), P.Job('base', spec()), child('grandchild', 'child'), P.Job('independent', spec(table='other'))))


def builders(plan, calls):
    def build(context):
        calls.append(context.job.name)
        context.checkpoint()
        assert all(a.identity in {d.sha256 for d in context.expected.dependencies} for a in context.parents.values())
        return {'rows.json': {'job': context.job.name, 'parent_ids': sorted(a.identity for a in context.parents.values())}}
    return {j.name: build for j in plan.jobs}


def journal(root):
    return json.loads((root / 'journal.json').read_text())['journal']


def test_topology_and_verified_content_ids_resume(tmp_path):
    plan, calls = graph(), []
    result = P.run(plan, builders(plan, calls), tmp_path, LIMIT)
    assert calls == ['base', 'independent', 'child', 'grandchild']
    assert set(result.states.values()) == {'succeeded'}
    for job in plan.jobs:
        output = result.artifacts[job.name]
        output.verify_contents()
        for upstream in job.upstream:
            assert A.Dependency(upstream.kind, upstream.name, result.artifacts[upstream.job].identity) in output.spec.dependencies
    resumed = P.run(plan, {}, tmp_path, LIMIT)
    assert set(resumed.states.values()) == {'succeeded'}
    assert {k: a.identity for k, a in resumed.artifacts.items()} == {k: a.identity for k, a in result.artifacts.items()}
    log = journal(tmp_path)
    assert len(log['runs']) == 2
    assert all(r['resumed'] for r in log['runs'][-1]['jobs'].values())
    assert log['runs'][0]['package_bytes'] == log['runs'][0]['artifact_storage_bytes']
    assert log['runs'][0]['process_peak_rss_bytes'] > 0
    assert all(h['attempts'][0]['runtime_seconds'] >= 0 for h in log['records'].values())


@pytest.mark.parametrize('change', ['table', 'settings', 'target'])
def test_only_declared_descendants_rebuild(tmp_path, change):
    original, calls = graph(), []
    first = P.run(original, builders(original, calls), tmp_path, LIMIT)
    base = next(j for j in original.jobs if j.name == 'base')
    if change == 'table':
        updated = replace(base.spec, dependencies=tuple(replace(d, sha256=OTHER) if d.kind == 'table' else d for d in base.spec.dependencies))
    elif change == 'settings':
        updated = replace(base.spec, settings_json=A.canonical_object({'setting': 2}))
    else:
        updated = replace(base.spec, target='another', query_json=Query('Tg', 'label', target='another').to_json())
    plan = replace(original, jobs=tuple(replace(j, spec=updated) if j.name == 'base' else j for j in original.jobs))
    assert plan.fingerprints()['independent'] == original.fingerprints()['independent']
    assert all(plan.fingerprints()[name] != original.fingerprints()[name] for name in ('base', 'child', 'grandchild'))
    calls.clear()
    second = P.run(plan, builders(plan, calls), tmp_path, LIMIT)
    assert calls == ['base', 'child', 'grandchild']
    assert second.artifacts['independent'].identity == first.artifacts['independent'].identity


@pytest.mark.parametrize('case', ['payload', 'manifest', 'missing', 'symlink'])
def test_invalid_completed_partition_rebuilds_and_audits(tmp_path, case):
    plan, calls = graph(), []
    first = P.run(plan, builders(plan, calls), tmp_path, LIMIT)
    path = tmp_path / 'artifacts' / plan.fingerprints()['base'] / '1'
    if case == 'payload':
        (path / 'rows.json').write_text('{}')
    elif case == 'manifest':
        (path / 'manifest.json').write_text('{}')
    elif case == 'missing':
        (path / 'manifest.json').unlink()
    else:
        (path / 'rows.json').unlink()
        (path / 'rows.json').symlink_to(tmp_path / 'external.json')
    calls.clear()
    second = P.run(plan, builders(plan, calls), tmp_path, LIMIT)
    assert calls == ['base']  # Deterministic repaired content keeps descendants valid.
    assert second.artifacts['base'].identity == first.artifacts['base'].identity
    history = journal(tmp_path)['records'][plan.fingerprints()['base']]['attempts']
    assert len(history) == 2 and history[0]['resume_refusals']
    assert path.exists()


def test_failure_blocks_descendants_but_independent_job_completes(tmp_path):
    plan, calls = graph(), []
    mapping = builders(plan, calls)
    def fail(context):
        raise RuntimeError('fixture failure')
    mapping['base'] = fail
    result = P.run(plan, mapping, tmp_path, LIMIT)
    assert dict(result.states) == {'base': 'failed', 'independent': 'succeeded', 'child': 'blocked', 'grandchild': 'blocked'}
    assert calls == ['independent']
    assert journal(tmp_path)['runs'][0]['jobs']['base']['error'] == {'type': 'RuntimeError', 'message': 'fixture failure'}
    recovered = P.run(plan, builders(plan, calls), tmp_path, LIMIT)
    assert set(recovered.states.values()) == {'succeeded'}


def test_cancellation_preserves_completed_partitions(tmp_path):
    plan = P.Plan((P.Job('base', spec()), child()))
    calls, cancel = [], [False]
    mapping = builders(plan, calls)
    build = mapping['child']
    def stop(context):
        cancel[0] = True
        context.checkpoint()
        return build(context)
    mapping['child'] = stop
    result = P.run(plan, mapping, tmp_path, LIMIT, cancelled=lambda: cancel[0])
    assert dict(result.states) == {'base': 'succeeded', 'child': 'pending'}
    assert result.stop_reason == 'Cancellation requested'
    calls.clear()
    resumed = P.run(plan, builders(plan, calls), tmp_path, LIMIT)
    assert calls == ['child'] and set(resumed.states.values()) == {'succeeded'}
    attempts = journal(tmp_path)['records'][plan.fingerprints()['child']]['attempts']
    assert [r['status'] for r in attempts] == ['cancelled', 'succeeded']


def test_attempt_budget_resumes_completed_partitions(tmp_path):
    plan, calls = graph(), []
    result = P.run(plan, builders(plan, calls), tmp_path, replace(LIMIT, max_jobs=1))
    assert calls == ['base'] and result.states['base'] == 'succeeded'
    assert result.stop_reason == 'Run job attempt budget exhausted'
    calls.clear()
    P.run(plan, builders(plan, calls), tmp_path, LIMIT)
    assert calls == ['independent', 'child', 'grandchild']


@pytest.mark.parametrize('limit', ['max_storage_bytes', 'max_package_bytes'])
def test_output_budget_refuses_before_writing(tmp_path, limit):
    plan = P.Plan((P.Job('base', spec()),))
    result = P.run(plan, {'base': lambda context: {'data.json': 'x' * 5000}}, tmp_path, replace(LIMIT, **{limit: 1000}))
    assert result.states['base'] == 'pending' and 'budget exhausted' in result.stop_reason
    assert not list((tmp_path / 'artifacts').rglob('manifest.json'))
    assert not list((tmp_path / 'artifacts').rglob('data.json'))


def test_memory_budget_is_reported_without_starting_builder(tmp_path, monkeypatch):
    monkeypatch.setattr(P, '_rss', lambda: LIMIT.memory_limit_bytes + 1)
    plan = P.Plan((P.Job('base', spec()),))
    result = P.run(plan, {}, tmp_path, LIMIT)
    assert result.states['base'] == 'pending' and 'Current process RSS' in result.stop_reason


def test_linux_memory_reader_uses_current_rss_and_separate_process_peak(monkeypatch, controlled_rss):
    monkeypatch.setattr(P.sys, 'platform', 'linux')
    monkeypatch.setattr(P, 'Path', lambda path: SimpleNamespace(
        read_text=lambda: 'Name:\tfixture\nVmHWM:\t123456 kB\nVmRSS:\t12000 kB\n'))
    def inherited_peak_is_not_current_usage(*args):
        raise AssertionError('Inherited launcher peak must not replace /proc values')
    monkeypatch.setattr(P.resource, 'getrusage', inherited_peak_is_not_current_usage)
    assert controlled_rss() == 12000 * 1024
    assert P._peak_rss() == 123456 * 1024


def test_time_budget_cooperative_stop(tmp_path, monkeypatch):
    clock = [0.0]
    monkeypatch.setattr(P.time, 'monotonic', lambda: clock[0])
    plan = P.Plan((P.Job('base', spec()),))
    def build(context):
        clock[0] = LIMIT.max_seconds
        context.checkpoint()
    result = P.run(plan, {'base': build}, tmp_path, LIMIT)
    assert result.states['base'] == 'pending' and 'time budget' in result.stop_reason


def test_snapshot_changes_during_builder_block_output_and_descendants(tmp_path):
    source = tmp_path / 'input.txt'
    source.write_bytes(b'fixture')
    plan = replace(graph(), snapshots=(P.Snapshot(A.Dependency('table', 'entities', HASH), str(source)),))
    calls = []
    mapping = builders(plan, calls)
    def mutate(context):
        source.write_bytes(b'changed')
        return {'data.json': []}
    mapping['base'] = mutate
    result = P.run(plan, mapping, tmp_path / 'run', LIMIT)
    assert result.states['base'] == 'failed'
    assert result.states['child'] == result.states['grandchild'] == 'blocked'
    assert result.states['independent'] == 'succeeded'
    assert not list((tmp_path / 'run' / 'artifacts' / plan.fingerprints()['base']).rglob('manifest.json'))


def test_missing_snapshot_fails_and_never_resumes_stale_partition(tmp_path):
    source = tmp_path / 'input.txt'
    source.write_bytes(b'fixture')
    plan = P.Plan((P.Job('base', spec()),), (P.Snapshot(A.Dependency('table', 'entities', HASH), str(source)),))
    P.run(plan, builders(plan, []), tmp_path / 'run', LIMIT)
    source.unlink()
    result = P.run(plan, {}, tmp_path / 'run', LIMIT)
    assert result.states['base'] == 'failed'


def test_unavailable_partition_records_gaps_and_blocks_child(tmp_path):
    base = replace(spec(), status='unavailable', gaps=('fixture missing input',))
    plan = P.Plan((P.Job('base', base), child()))
    result = P.run(plan, {'base': lambda context: {'gap.json': base.gaps}}, tmp_path, LIMIT)
    assert dict(result.states) == {'base': 'unavailable', 'child': 'blocked'}
    resumed = P.run(plan, {}, tmp_path, LIMIT)
    assert dict(resumed.states) == dict(result.states)


@pytest.mark.parametrize('role', ['reusable_base', 'fitted_model', 'held_out', 'deployment', 'scorecard'])
def test_artifact_and_split_roles_remain_distinct(tmp_path, role):
    job = P.Job('partition', spec(role=role), split=None if role == 'reusable_base' else SPLIT)
    result = P.run(P.Plan((job,)), {'partition': lambda context: {'data.json': {'role': role}}}, tmp_path, LIMIT)
    assert result.states['partition'] == 'succeeded'
    assert result.artifacts['partition'].spec.role == role
    if role != 'reusable_base':
        with pytest.raises(ValueError):
            replace(job, split=replace(SPLIT, seed=3))


@pytest.mark.parametrize('kind', ['model', 'calibration'])
def test_model_and_calibration_bind_complete_parent_id(tmp_path, kind):
    plan = P.Plan((P.Job('base', spec()), child(kind=kind)))
    result = P.run(plan, builders(plan, []), tmp_path, LIMIT)
    assert A.Dependency(kind, 'upstream', result.artifacts['base'].identity) in result.artifacts['child'].spec.dependencies


def test_immutable_graph_validation():
    with pytest.raises(ValueError, match='Unknown upstream'):
        P.Plan((child(),))
    with pytest.raises(ValueError, match='cycle'):
        P.Plan((child('base', 'child'), child()))
    with pytest.raises(ValueError, match='Duplicate job'):
        P.Plan((P.Job('base', spec()), P.Job('base', spec())))
    with pytest.raises(ValueError, match='absent'):
        P.Job('base', spec(), (P.Upstream('other', 'artifact', 'undeclared'),))
    with pytest.raises(ValueError, match='immutable'):
        P.Plan([P.Job('base', spec())])
    with pytest.raises(ValueError, match='positive'):
        replace(LIMIT, max_jobs=0)
    with pytest.raises(ValueError, match='finite'):
        replace(LIMIT, max_seconds=float('nan'))


def test_atomic_journal_corruption_is_refused(tmp_path, monkeypatch):
    plan = P.Plan((P.Job('base', spec()),))
    P.run(plan, builders(plan, []), tmp_path, LIMIT)
    original = (tmp_path / 'journal.json').read_bytes()
    def fail_replace(source, target):
        raise OSError('fixture replace failed')
    monkeypatch.setattr(P.os, 'replace', fail_replace)
    with pytest.raises(OSError, match='replace failed'):
        P.run(plan, {}, tmp_path, LIMIT)
    assert (tmp_path / 'journal.json').read_bytes() == original
    assert not list(tmp_path.glob('.journal-*'))
    monkeypatch.undo()
    payload = json.loads(original)
    payload['journal']['runs'][0]['plan'] = OTHER
    (tmp_path / 'journal.json').write_text(json.dumps(payload))
    with pytest.raises(ValueError, match='corrupt'):
        P.run(plan, {}, tmp_path, LIMIT)


def test_crash_orphan_is_preserved_and_rebuilt(tmp_path):
    plan = P.Plan((P.Job('base', spec()),))
    key = plan.fingerprints()['base']
    orphan = tmp_path / 'artifacts' / key / '1'
    orphan.mkdir(parents=True)
    (orphan / 'rows.json').write_text('partial')
    result = P.run(plan, builders(plan, []), tmp_path, LIMIT)
    assert result.states['base'] == 'succeeded'
    assert (orphan / 'rows.json').read_text() == 'partial'
    assert (tmp_path / 'artifacts' / key / '2' / 'manifest.json').is_file()


def test_single_writer_lock(tmp_path):
    plan = P.Plan((P.Job('base', spec()),))
    with (tmp_path / '.writer.lock').open('a') as stream:
        fcntl.flock(stream, fcntl.LOCK_EX | fcntl.LOCK_NB)
        with pytest.raises(RuntimeError, match='Another runner'):
            P.run(plan, {}, tmp_path, LIMIT)


def test_completed_artifact_recovers_after_interrupted_journal_write(tmp_path, monkeypatch):
    plan = P.Plan((P.Job('base', spec()),))
    original_read = A.read_artifact
    def interrupt_after_write(*args, **kwargs):
        raise KeyboardInterrupt('fixture interrupt after completion manifest')
    monkeypatch.setattr(A, 'read_artifact', interrupt_after_write)
    with pytest.raises(KeyboardInterrupt):
        P.run(plan, builders(plan, []), tmp_path, LIMIT)
    monkeypatch.setattr(A, 'read_artifact', original_read)
    recovered = P.run(plan, {}, tmp_path, LIMIT)
    assert recovered.states['base'] == 'succeeded'
    attempts = journal(tmp_path)['records'][plan.fingerprints()['base']]['attempts']
    assert len(attempts) == 1 and attempts[0]['recovered_after_interruption']
    assert journal(tmp_path)['runs'][0]['status'] == 'interrupted'


def test_incomplete_interrupted_partition_is_preserved_and_retried(tmp_path):
    plan = P.Plan((P.Job('base', spec()),))
    def interrupt(context):
        partial = tmp_path / 'artifacts' / plan.fingerprints()['base'] / '1'
        partial.mkdir(parents=True)
        (partial / 'rows.json').write_text('partial')
        raise KeyboardInterrupt('fixture interrupted builder')
    with pytest.raises(KeyboardInterrupt):
        P.run(plan, {'base': interrupt}, tmp_path, LIMIT)
    result = P.run(plan, builders(plan, []), tmp_path, LIMIT)
    assert result.states['base'] == 'succeeded'
    attempts = journal(tmp_path)['records'][plan.fingerprints()['base']]['attempts']
    assert [r['status'] for r in attempts] == ['interrupted', 'succeeded']
    assert attempts[0]['resume_refusals']


@pytest.mark.parametrize('mutation', ['memory', 'disk'])
def test_builder_parent_mutation_fails_before_child_publication(tmp_path, mutation):
    plan = P.Plan((P.Job('base', spec()), child()))
    mapping = builders(plan, [])
    def corrupt(context):
        if mutation == 'memory':
            context.parents['base'].payloads['rows.json']['job'] = 'changed'
        else:
            path = tmp_path / 'artifacts' / plan.fingerprints()['base'] / '1' / 'rows.json'
            path.write_text('{}')
        context.checkpoint()
        return {'data.json': []}
    mapping['child'] = corrupt
    result = P.run(plan, mapping, tmp_path, LIMIT)
    assert result.states['child'] == 'failed'
    result.artifacts['base'].verify_contents()
    assert not list((tmp_path / 'artifacts' / plan.fingerprints()['child']).rglob('manifest.json'))
