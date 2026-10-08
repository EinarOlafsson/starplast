"""Serial, resumable artifact partitions with explicit dependency and budget records.

Each job is one independently reusable partition. Its ArtifactSpec is an immutable
template: Upstream bindings replace only declared artifact/model/calibration hashes
with verified complete parent content identities before the builder is called.
Fresh input identities remain the caller's responsibility, as in artifacts.py;
optional Snapshot records additionally verify local input bytes at every checkpoint.

Builders receive a cooperative checkpoint, not a preemptive sandbox. Run them under
an OS memory/time cap when required. This module supplies no scientific builders,
source acquisition, fitting, or evidence that a fixture is biologically accurate.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import sys
import tempfile
import time
from types import MappingProxyType

from . import artifacts as A

SCHEMA_VERSION = 1


def _hash(value):
    return hashlib.sha256(A.canonical_object(value).encode('utf-8')).hexdigest()


@dataclass(frozen=True)
class Upstream:
    """Bind a declared dependency address to a named job's complete output ID."""

    job: str
    kind: str
    name: str

    def __post_init__(self):
        A.Dependency(self.kind, self.name, '0' * 64)
        A._name(self.job)
        if self.kind not in {'artifact', 'model', 'calibration'}:
            raise ValueError('Upstream jobs require artifact, model or calibration dependencies')


@dataclass(frozen=True)
class Snapshot:
    """A local immutable file whose byte digest is the declared dependency ID."""

    dependency: A.Dependency
    path: str

    def __post_init__(self):
        if not isinstance(self.dependency, A.Dependency) or not isinstance(self.path, str) or not self.path:
            raise ValueError('Snapshot requires a typed dependency and explicit file path')

    def verify(self):
        """Rehash the local file, rejecting missing, symlinked or changed inputs."""
        path = Path(self.path)
        if path.is_symlink() or not path.is_file():
            raise ValueError('Snapshot must remain a regular local file')
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(1024 * 1024), b''):
                digest.update(block)
        if digest.hexdigest() != self.dependency.sha256:
            raise ValueError('Input snapshot changed: ' + self.dependency.name)


@dataclass(frozen=True)
class Job:
    """One bounded output partition; fitting/evaluation roles come from its spec."""

    name: str
    spec: A.ArtifactSpec
    upstream: tuple[Upstream, ...] = ()
    split: object = None

    def __post_init__(self):
        A._name(self.name)
        if not isinstance(self.spec, A.ArtifactSpec):
            raise ValueError('Job requires a validated ArtifactSpec')
        if not isinstance(self.upstream, tuple) or not all(isinstance(u, Upstream) for u in self.upstream):
            raise ValueError('Upstream bindings must be immutable typed tuples')
        addresses = [(u.kind, u.name) for u in self.upstream]
        if len(set(addresses)) != len(addresses):
            raise ValueError('Duplicate upstream dependency address')
        declared = {(d.kind, d.name) for d in self.spec.dependencies}
        if not set(addresses) <= declared:
            raise ValueError('Upstream address is absent from the artifact specification')
        if any(d.kind == 'split' for d in self.spec.dependencies):
            self.spec.validate_split(self.split)
        elif self.split is not None:
            raise ValueError('Split must be declared in the artifact dependencies')


@dataclass(frozen=True)
class Plan:
    """An acyclic named job graph; unrelated job changes retain their cache keys."""

    jobs: tuple[Job, ...]
    snapshots: tuple[Snapshot, ...] = ()

    def __post_init__(self):
        if not isinstance(self.jobs, tuple) or not self.jobs or not all(isinstance(j, Job) for j in self.jobs):
            raise ValueError('Plan requires a nonempty immutable tuple of jobs')
        if len({j.name for j in self.jobs}) != len(self.jobs):
            raise ValueError('Duplicate job name')
        if not isinstance(self.snapshots, tuple) or not all(isinstance(s, Snapshot) for s in self.snapshots):
            raise ValueError('Snapshots must be immutable typed tuples')
        addresses = [(s.dependency.kind, s.dependency.name) for s in self.snapshots]
        if len(set(addresses)) != len(addresses):
            raise ValueError('Duplicate snapshot address')
        for snapshot in self.snapshots:
            users = [d for j in self.jobs for d in j.spec.dependencies
                     if (d.kind, d.name) == (snapshot.dependency.kind, snapshot.dependency.name)]
            if not users or any(d != snapshot.dependency for d in users):
                raise ValueError('Snapshot identity must match every declared use')
            if any((u.kind, u.name) == (snapshot.dependency.kind, snapshot.dependency.name)
                   for j in self.jobs for u in j.upstream):
                raise ValueError('A dependency cannot bind both a snapshot and an upstream job')
        self.ordered()

    def ordered(self):
        """Return stable topological order, rejecting missing parents and cycles."""
        pending = list(self.jobs)
        known = {j.name for j in pending}
        if any(u.job not in known for j in pending for u in j.upstream):
            raise ValueError('Unknown upstream job')
        result, done = [], set()
        while pending:
            ready = [j for j in pending if all(u.job in done for u in j.upstream)]
            if not ready:
                raise ValueError('Job dependency cycle')
            for job in ready:
                result.append(job)
                done.add(job.name)
                pending.remove(job)
        return tuple(result)

    def fingerprints(self):
        """Transitive immutable recipe keys, distinct from artifact content IDs."""
        keys = {}
        for job in self.ordered():
            keys[job.name] = _hash({'name': job.name, 'spec': asdict(job.spec),
                'upstream': [{**asdict(u), 'job_key': keys[u.job]} for u in job.upstream]})
        return MappingProxyType(keys)

    @property
    def identity(self):
        return _hash(dict(self.fingerprints()))


@dataclass(frozen=True)
class Budget:
    """Per-run attempt/time limits and total disk/current-package byte limits."""

    max_jobs: int
    max_seconds: float
    max_storage_bytes: int
    max_package_bytes: int
    memory_limit_bytes: int

    def __post_init__(self):
        for name in ('max_jobs', 'max_storage_bytes', 'max_package_bytes', 'memory_limit_bytes'):
            if type(getattr(self, name)) is not int or getattr(self, name) <= 0:
                raise ValueError('Budget limits must be explicit positive integers')
        if isinstance(self.max_seconds, bool) or not isinstance(self.max_seconds, (float, int)) or not math.isfinite(self.max_seconds) or self.max_seconds <= 0:
            raise ValueError('Time budget must be positive and finite')


class Cancelled(Exception):
    """A cooperative cancellation; completed earlier partitions remain valid."""


class BudgetExceeded(Exception):
    """A cooperative budget stop; builders may be retried under a new budget."""


def _process_memory(field):
    # Linux getrusage may retain a launcher's pre-exec high-water mark. VmHWM
    # records this executed process instead; do not reject a tiny new worker
    # because its launcher once used more RAM.
    if sys.platform.startswith('linux'):
        for row in Path('/proc/self/status').read_text().splitlines():
            if row.startswith(field + ':'):
                return int(row.split()[1]) * 1024
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * (1 if sys.platform == 'darwin' else 1024))


def _rss():
    return _process_memory('VmRSS')


def _peak_rss():
    # Process lifetime high-water mark, including preceding partitions. This
    # is deliberately not presented as an isolated per-job memory measurement.
    return _process_memory('VmHWM')


def _size(path):
    return sum(p.stat().st_size for p in path.rglob('*') if p.is_file() and not p.is_symlink())


def _artifact_size(spec, payloads):
    """Preflight exact writer bytes so a budget refusal leaves no new payloads."""
    if not isinstance(payloads, dict) or not payloads:
        raise ValueError('Builder must return nonempty named JSON payloads')
    files = {}
    for name, value in payloads.items():
        A._payload_name(name)
        data = (A._json(value) + '\n').encode('utf-8')
        files[name] = {'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data)}
    body = {'schema_version': A.SCHEMA_VERSION, 'key': spec.identity, 'spec': asdict(spec), 'files': files}
    return sum(r['bytes'] for r in files.values()) + len((A._json({**body, 'identity': A._digest(body)}) + '\n').encode('utf-8'))


def _atomic_json(path, value):
    """Fsync data, atomically replace the journal, then fsync its directory."""
    data = (A.canonical_object(value) + '\n').encode('utf-8')
    descriptor, temporary = tempfile.mkstemp(prefix='.journal-', dir=path.parent)
    try:
        with os.fdopen(descriptor, 'wb') as stream:
            stream.write(data)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
        directory = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(directory)
        finally:
            os.close(directory)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


@dataclass(frozen=True)
class Context:
    """Builder view with verified parent artifacts and a cooperative checkpoint."""

    job: Job
    expected: A.ArtifactSpec
    parents: object
    checkpoint: object


@dataclass(frozen=True)
class RunResult:
    """Current job outcomes; journal retains all attempts and past plan records."""

    states: object
    artifacts: object
    stop_reason: str
    journal: Path


def run(plan, builders, directory, budget, *, cancelled=lambda: False):
    """Build/resume serial partitions, preserving invalid/partial outputs for audit.

    Builders map names to ``builder(Context) -> {safe_json_name: payload}``.
    Successful resume checks current spec, split, payload bytes and journal content
    identity; a failed check rebuilds in a new attempt directory. Failures block
    declared descendants and independent jobs continue. Cancellation/budget stops
    leave remaining jobs pending. One writer owns a directory at a time; an OS
    advisory lock automatically releases after process death.
    """
    import fcntl

    if not isinstance(plan, Plan) or not isinstance(budget, Budget):
        raise TypeError('A validated Plan and Budget are required')
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    for path in (root, root / 'journal.json', root / '.writer.lock', root / 'artifacts'):
        if path.is_symlink():
            raise ValueError('Run paths cannot be symlinks')
    with (root / '.writer.lock').open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError('Another runner owns this journal') from error
        return _run_locked(plan, builders, root, budget, cancelled)


def _run_locked(plan, builders, root, budget, cancelled):
    journal_path = root / 'journal.json'
    if journal_path.exists():
        envelope = json.loads(journal_path.read_text())
        journal = envelope.get('journal')
        if (set(envelope) != {'journal', 'sha256'} or not isinstance(journal, dict)
                or _hash(journal) != envelope['sha256'] or journal.get('schema_version') != SCHEMA_VERSION
                or not isinstance(journal.get('records'), dict) or not isinstance(journal.get('runs'), list)):
            raise ValueError('Journal is corrupt or unsupported')
    else:
        journal = {'schema_version': SCHEMA_VERSION, 'records': {}, 'runs': []}
    for history in journal['records'].values():
        if not isinstance(history, dict) or not isinstance(history.get('attempts'), list):
            raise ValueError('Journal job record is invalid')
        for record in history['attempts']:
            if record.get('status') == 'running':
                record.update(status='interrupted', error={'type': 'Interrupted',
                    'message': 'Previous runner ended without recording completion'})
    for previous in journal['runs']:
        if previous.get('status') == 'running':
            previous.update(status='interrupted', stop_reason='Previous runner ended without recording completion')
    artifacts_root = root / 'artifacts'
    artifacts_root.mkdir(exist_ok=True)
    started = time.monotonic()
    run_record = {'plan': plan.identity, 'budget': asdict(budget), 'started_unix': time.time(),
        'status': 'running', 'stop_reason': '', 'jobs': {}}
    journal['runs'].append(run_record)
    states, results, used_package, attempts = {}, {}, 0, 0
    artifact_paths, builder_parents = {}, {}
    fingerprints = plan.fingerprints()
    ordered = plan.ordered()

    def save():
        _atomic_json(journal_path, {'journal': journal, 'sha256': _hash(journal)})

    def checkpoint(job):
        if cancelled():
            raise Cancelled('Cancellation requested')
        if time.monotonic() - started >= budget.max_seconds:
            raise BudgetExceeded('Run time budget exhausted')
        if _rss() > budget.memory_limit_bytes:
            raise BudgetExceeded('Current process RSS exceeds the declared memory budget')
        if _size(artifacts_root) > budget.max_storage_bytes:
            raise BudgetExceeded('Artifact storage budget exhausted')
        addresses = {(d.kind, d.name) for d in job.spec.dependencies}
        for snapshot in plan.snapshots:
            if (snapshot.dependency.kind, snapshot.dependency.name) in addresses:
                snapshot.verify()
        for upstream in job.upstream:
            parent = results[upstream.job]
            parent.verify_contents()
            verified = A.read_artifact(artifact_paths[upstream.job], expected=parent.spec,
                split=next(j.split for j in ordered if j.name == upstream.job))
            if verified.identity != parent.identity:
                raise ValueError('Upstream artifact content identity changed')
        for parent in builder_parents.get(job.name, {}).values():
            parent.verify_contents()

    save()
    stop_reason = ''
    for job in ordered:
        key = fingerprints[job.name]
        history = journal['records'].setdefault(key, {'name': job.name, 'attempts': []})
        if history.get('name') != job.name or not isinstance(history.get('attempts'), list):
            raise ValueError('Journal job record is invalid')
        blocked = [u.job for u in job.upstream if states[u.job] in {'failed', 'blocked', 'unavailable'}]
        if blocked:
            states[job.name] = 'blocked'
            run_record['jobs'][job.name] = {'status': 'blocked', 'parents': sorted(set(blocked)), 'key': key}
            save()
            continue
        binding = {(u.kind, u.name): results[u.job].identity for u in job.upstream}
        expected = replace(job.spec, dependencies=tuple(
            replace(d, sha256=binding.get((d.kind, d.name), d.sha256)) for d in job.spec.dependencies))
        current = None
        job_started = time.monotonic()
        try:
            checkpoint(job)
            for record in reversed(history['attempts']):
                if record.get('status') not in {'succeeded', 'interrupted'} or record.get('spec_key') != expected.identity:
                    continue
                try:
                    path = artifacts_root / key / str(record['attempt'])
                    if path.is_symlink() or path.parent.is_symlink():
                        raise ValueError('Artifact directory cannot be a symlink')
                    artifact = A.read_artifact(path, expected=expected, split=job.split)
                    if artifact.identity != record.get('content_identity', artifact.identity):
                        raise ValueError('Journal and artifact content identities differ')
                    size = _size(path)
                except (OSError, ValueError, KeyError, TypeError) as error:
                    record.setdefault('resume_refusals', []).append({'type': type(error).__name__, 'message': str(error)})
                    save()
                    continue
                if used_package + size > budget.max_package_bytes:
                    raise BudgetExceeded('Current package byte budget exhausted')
                checkpoint(job)
                if record['status'] == 'interrupted':
                    record.update(status='succeeded', content_identity=artifact.identity,
                        storage_bytes=size, recovered_after_interruption=True)
                results[job.name] = artifact
                artifact_paths[job.name] = path
                used_package += size
                states[job.name] = 'succeeded' if artifact.spec.status == 'ok' else 'unavailable'
                run_record['jobs'][job.name] = {'status': states[job.name], 'resumed': True, 'key': key,
                    'content_identity': artifact.identity, 'storage_bytes': size}
                break
            if job.name in results:
                save()
                continue
            if attempts >= budget.max_jobs:
                raise BudgetExceeded('Run job attempt budget exhausted')
            attempts += 1
            number = max([r['attempt'] for r in history['attempts']] + [0]) + 1
            path = artifacts_root / key / str(number)
            # Crash-orphaned and partially written directories are never overwritten.
            while path.exists():
                number += 1
                path = artifacts_root / key / str(number)
            current = {'attempt': number, 'status': 'running', 'spec_key': expected.identity,
                'started_unix': time.time()}
            history['attempts'].append(current)
            save()
            builder_parents[job.name] = {u.job: A.read_artifact(artifact_paths[u.job],
                expected=results[u.job].spec, split=next(j.split for j in ordered if j.name == u.job)) for u in job.upstream}
            context = Context(job, expected, MappingProxyType(builder_parents[job.name]), lambda: checkpoint(job))
            payloads = builders[job.name](context)
            checkpoint(job)
            planned_size = _artifact_size(expected, payloads)
            if _size(artifacts_root) + planned_size > budget.max_storage_bytes:
                raise BudgetExceeded('Artifact storage budget exhausted')
            if used_package + planned_size > budget.max_package_bytes:
                raise BudgetExceeded('Current package byte budget exhausted')
            identity = A.write_artifact(path, expected, payloads, split=job.split)
            checkpoint(job)
            size = _size(path)
            if used_package + size > budget.max_package_bytes:
                raise BudgetExceeded('Current package byte budget exhausted')
            artifact = A.read_artifact(path, expected=expected, split=job.split)
            current.update(status='succeeded', content_identity=identity, storage_bytes=size)
            results[job.name] = artifact
            artifact_paths[job.name] = path
            states[job.name] = 'succeeded' if artifact.spec.status == 'ok' else 'unavailable'
            used_package += size
        except (Cancelled, BudgetExceeded) as error:
            stop_reason = str(error)
            if current is not None:
                current.update(status='cancelled' if isinstance(error, Cancelled) else 'budget_exhausted',
                    error={'type': type(error).__name__, 'message': str(error)})
            states[job.name] = 'pending'
            run_record['jobs'][job.name] = {'status': 'pending', 'key': key, 'reason': stop_reason}
            break
        except Exception as error:
            states[job.name] = 'failed'
            failure = {'type': type(error).__name__, 'message': str(error)}
            if current is not None:
                current.update(status='failed', error=failure)
            run_record['jobs'][job.name] = {'status': 'failed', 'key': key, 'error': failure}
        finally:
            if current is not None:
                current.update(runtime_seconds=time.monotonic() - job_started,
                    process_peak_rss_bytes=_peak_rss(), process_rss_bytes=_rss())
            if job.name in results:
                run_record['jobs'][job.name] = {**run_record['jobs'].get(job.name, {}),
                    'status': states[job.name], 'key': key, 'content_identity': results[job.name].identity}
            save()
    for job in ordered:
        if job.name not in states:
            states[job.name] = 'pending'
            run_record['jobs'][job.name] = {'status': 'pending', 'key': fingerprints[job.name], 'reason': stop_reason}
    run_record.update(stop_reason=stop_reason, runtime_seconds=time.monotonic() - started,
        status='stopped' if stop_reason else 'finished', process_peak_rss_bytes=_peak_rss(),
        process_rss_bytes=_rss(), package_bytes=used_package, artifact_storage_bytes=_size(artifacts_root))
    save()
    return RunResult(MappingProxyType(states), MappingProxyType(results), stop_reason, journal_path)
