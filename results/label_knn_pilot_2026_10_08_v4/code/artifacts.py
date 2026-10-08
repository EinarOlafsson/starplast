"""Immutable inference artifacts with explicit evaluation roles and cache lineage.

This data-only contract complements legacy CSV exports. Callers must supply
fresh, verified input identities and expected scope when loading; a cache never
infers them from a filename. It does not certify biological truth or fitted code.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path
import re

from . import capabilities as C, scorecard as SC
from .query import Query

SCHEMA_VERSION = 1
DEPENDENCY_KINDS = frozenset({'table', 'graph', 'truth', 'code', 'split', 'exclusions',
                              'model', 'calibration', 'artifact', 'mapping', 'source'})
ROLES = frozenset({'held_out', 'deployment', 'fitted_model', 'scorecard', 'reusable_base'})
STATUS = frozenset({'ok', 'unavailable', 'outside_scope', 'failed'})


def _json(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False, allow_nan=False, separators=(',', ':'))


def _digest(value):
    return hashlib.sha256(_json(value).encode('utf-8')).hexdigest()


def _hash(value):
    if not isinstance(value, str) or not re.fullmatch(r'[a-f0-9]{64}', value):
        raise ValueError('A SHA-256 content identity is required')


def _name(value):
    if not isinstance(value, str) or not value or value != value.strip():
        raise ValueError('An explicit trimmed identity is required')


@dataclass(frozen=True, order=True)
class Dependency:
    """One named input identity; upstream artifacts pin their complete content ID."""

    kind: str
    name: str
    sha256: str

    def __post_init__(self):
        if self.kind not in DEPENDENCY_KINDS:
            raise ValueError('Unknown dependency kind')
        _name(self.name)
        _hash(self.sha256)


@dataclass(frozen=True)
class ArtifactSpec:
    """Complete expected scope for an inference, held-out row or scorecard snapshot.

    Entity order is part of the key. Partition, fit population, model/calibration
    identity and transductive access remain distinct. Settings and evaluation
    scope are canonical JSON objects, containing no functions or executable data.
    """

    query_json: str
    strategy: str
    task: str
    output: C.Output
    target: str
    entity_order: tuple[str, ...]
    role: str
    partition_id: str
    dependencies: tuple[Dependency, ...]
    settings_json: str
    evaluation_scope_json: str
    seed: int
    code_version: str
    fit_entities: tuple[str, ...] = ()
    fit_role: str = 'none'
    feature_access: str = 'inductive'
    benchmark_id: str = ''
    confidence_kind: str = 'none'
    calibration_scope: str = ''
    evaluation_partition: str = 'none'
    status: str = 'ok'
    gaps: tuple[str, ...] = ()
    schema_version: int = SCHEMA_VERSION

    def __post_init__(self):
        query = Query.from_json(self.query_json)
        cap = C.get(self.strategy)
        if self.task not in cap.benchmark_tasks or self.task not in SC.TASKS:
            raise ValueError('Task is not declared for this strategy')
        cap.validate_output(self.output)
        if query.organism not in cap.organisms:
            raise ValueError('Organism adapter is not declared')
        if query.kind not in cap.query_kinds:
            raise ValueError('Query kind is not declared for this strategy')
        if query.target and query.target != self.target:
            raise ValueError('Query and artifact target differ')
        if self.role not in ROLES or self.status not in STATUS or self.schema_version != SCHEMA_VERSION:
            raise ValueError('Unsupported artifact role, status or schema')
        _name(self.target)
        _name(self.partition_id)
        _name(self.code_version)
        for values in (self.entity_order, self.fit_entities):
            if not isinstance(values, tuple) or len(set(values)) != len(values):
                raise ValueError('Entity populations must be ordered unique tuples')
            for value in values:
                _name(value)
        if not self.entity_order and self.status == 'ok':
            raise ValueError('Successful artifacts require an explicit entity population')
        if type(self.seed) is not int or self.feature_access not in {'inductive', 'transductive'}:
            raise ValueError('Declare seed and feature-access regime')
        if not isinstance(self.dependencies, tuple) or not all(isinstance(d, Dependency) for d in self.dependencies):
            raise ValueError('Dependencies must be immutable typed identities')
        if len({(d.kind, d.name) for d in self.dependencies}) != len(self.dependencies):
            raise ValueError('Duplicate dependency address')
        kinds = {d.kind for d in self.dependencies}
        if not {'code', 'table', 'exclusions'} <= kinds:
            raise ValueError('Artifacts require code, entity-table and explicit exclusion identities')
        if self.role in {'held_out', 'scorecard'} and (not {'truth', 'split'} <= kinds or not self.benchmark_id):
            raise ValueError('Held-out outputs and scorecards require truth, split and benchmark identity')
        if self.role == 'deployment' and 'model' not in kinds:
            raise ValueError('Deployment outputs must identify their fitted model')
        expected_partition = 'test' if self.role in {'held_out', 'scorecard'} else 'deployment' if self.role == 'deployment' else 'none'
        if self.evaluation_partition != expected_partition:
            raise ValueError('Held-out test and unknown-entity deployment partitions must be distinct')
        if self.role in {'fitted_model', 'deployment'} and (not self.fit_entities or 'split' not in kinds):
            raise ValueError('Model/deployment artifacts require fit population and split identity')
        if self.fit_role not in {'none', 'train', 'train+tune', 'calibration'}:
            raise ValueError('Unsupported fit role')
        if bool(self.fit_entities) != (self.fit_role != 'none'):
            raise ValueError('Fit entities and fit role must agree')
        if self.role == 'held_out' and set(self.entity_order) & set(self.fit_entities):
            raise ValueError('Held-out entities cannot be model fitting entities')
        if self.confidence_kind not in {'none', 'method_support', 'calibrated_probability', 'prediction_set', 'prediction_interval'}:
            raise ValueError('Unknown confidence meaning')
        if self.confidence_kind in {'calibrated_probability', 'prediction_set', 'prediction_interval'} and (
                'calibration' not in kinds or not self.calibration_scope):
            raise ValueError('Calibrated confidence requires calibration identity and applicability scope')
        if self.confidence_kind == 'calibrated_probability' and self.output.kind != 'label_calls':
            raise ValueError('This artifact cannot reinterpret geometry or ranking as label probability')
        if self.confidence_kind == 'prediction_set' and self.output.kind != 'label_sets':
            raise ValueError('Set confidence needs a label-set output')
        if self.confidence_kind == 'prediction_interval' and self.output.kind != 'numeric_intervals':
            raise ValueError('Interval confidence needs a numerical interval output')
        if self.status != 'ok' and not self.gaps:
            raise ValueError('Unavailable or failed outputs require explicit gaps')
        for value in (self.settings_json, self.evaluation_scope_json):
            parsed = json.loads(value, parse_constant=lambda s: _nonfinite(s))
            if not isinstance(parsed, dict) or _json(parsed) != value:
                raise ValueError('Settings and scope must be canonical finite JSON objects')
        scope = json.loads(self.evaluation_scope_json)
        if not {'unit', 'eligible_population', 'truth_grade', 'negative_semantics', 'context', 'limitations'} <= set(scope):
            raise ValueError('Evaluation scope must retain population, truth and interpretation limits')
        if scope['unit'] != cap.benchmark_unit:
            raise ValueError('Evaluation unit differs from the strategy contract')
        if type(scope['eligible_population']) is not int or scope['eligible_population'] < 0:
            raise ValueError('Eligible population must be an explicit nonnegative count')
        if self.role == 'held_out' and scope['eligible_population'] != len(self.entity_order):
            raise ValueError('Held-out cohort must match the frozen eligible population')

    @property
    def organism(self):
        """The organism explicitly recorded by the typed query."""
        return Query.from_json(self.query_json).organism

    @property
    def identity(self):
        """Canonical cache key including ordered entities and every declared input."""
        return _digest(asdict(self))

    def validate_split(self, split):
        """Check the recorded frozen split, fitting access and exact test cohort."""
        from .splits import SplitManifest
        if not isinstance(split, SplitManifest) or split.organism != self.organism:
            raise ValueError('A matching frozen split is required')
        identities = [d.sha256 for d in self.dependencies if d.kind == 'split']
        if identities != [split.identity] or self.feature_access != split.feature_access:
            raise ValueError('Split identity or feature-access regime differs')
        if self.benchmark_id and self.benchmark_id != split.benchmark_id:
            raise ValueError('Benchmark and split identity differ')
        if self.role in {'held_out', 'scorecard'} and self.entity_order != split.entities('test'):
            raise ValueError('Artifact does not retain the complete ordered outer-test cohort')
        if self.role == 'deployment' and set(self.entity_order) & {a.entity for a in split.assignments}:
            raise ValueError('Known benchmark entities require held-out outputs, not unknown-entity deployment')
        if self.fit_entities:
            operation = 'refit' if self.fit_role == 'train+tune' else 'calibration' if self.fit_role == 'calibration' else 'model_fit'
            split.guard_fit(operation, self.fit_entities)


def _nonfinite(value):
    raise ValueError('Nonfinite JSON number is not a recorded measurement: ' + value)


def canonical_object(value):
    """Encode finite JSON metadata deterministically; unavailable values use null."""
    if not isinstance(value, dict):
        raise ValueError('Metadata requires a JSON object')
    return _json(value)


def spec_from_dict(value):
    """Reconstruct an exact schema record without accepting extra/missing fields."""
    if not isinstance(value, dict) or set(value) != set(ArtifactSpec.__dataclass_fields__):
        raise ValueError('Artifact specification fields differ from this schema')
    value = dict(value)
    value['output'] = C.Output(**value['output'])
    value['dependencies'] = tuple(Dependency(**d) for d in value['dependencies'])
    for name in ('entity_order', 'fit_entities', 'gaps'):
        value[name] = tuple(value[name])
    return ArtifactSpec(**value)


@dataclass(frozen=True)
class Artifact:
    """A verified specification, data-only payloads and complete content identity."""

    spec: ArtifactSpec
    payloads: dict
    identity: str


def _payload_name(name):
    if not isinstance(name, str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*\.json', name) or name == 'manifest.json':
        raise ValueError('Artifact payloads need safe named JSON files')


def write_artifact(path, spec, payloads, *, split=None):
    """Write an immutable data-only artifact; an existing directory is never reused.

    The completion manifest is written last. A interrupted partial directory is
    refused by readers and preserved for inspection, never silently resumed.
    """
    if not isinstance(spec, ArtifactSpec) or not isinstance(payloads, dict) or not payloads:
        raise ValueError('A validated spec and nonempty named payloads are required')
    if any(d.kind == 'split' for d in spec.dependencies):
        spec.validate_split(split)
    encoded = {}
    for name, value in payloads.items():
        _payload_name(name)
        encoded[name] = (_json(value) + '\n').encode('utf-8')
    files = {name: {'sha256': hashlib.sha256(data).hexdigest(), 'bytes': len(data)} for name, data in encoded.items()}
    body = {'schema_version': SCHEMA_VERSION, 'key': spec.identity, 'spec': asdict(spec), 'files': files}
    identity = _digest(body)
    path = Path(path)
    path.mkdir(parents=True, exist_ok=False)
    for name, data in encoded.items():
        with (path / name).open('xb') as stream:
            stream.write(data)
    with (path / 'manifest.json').open('x') as stream:
        stream.write(_json({**body, 'identity': identity}) + '\n')
    return identity


def read_artifact(path, *, expected, split=None):
    """Read only the exact current scope, checking manifest and all payload bytes.

    Supply a freshly constructed expected spec from verified inputs. Reordered
    genes, changed sources/settings/exclusions or another fit/partition/model
    identity refuse reuse. Paths are relative; payload symlinks are refused.
    """
    if not isinstance(expected, ArtifactSpec):
        raise TypeError('A validated current expected specification is required')
    if any(d.kind == 'split' for d in expected.dependencies):
        expected.validate_split(split)
    path = Path(path)
    manifest = path / 'manifest.json'
    if manifest.is_symlink():
        raise ValueError('Manifest cannot be an external symlink')
    raw = json.loads(manifest.read_text(), parse_constant=lambda s: _nonfinite(s))
    if set(raw) != {'schema_version', 'key', 'spec', 'files', 'identity'} or raw['schema_version'] != SCHEMA_VERSION:
        raise ValueError('Unsupported artifact envelope')
    spec = spec_from_dict(raw['spec'])
    body = {k: raw[k] for k in ('schema_version', 'key', 'spec', 'files')}
    if _digest(body) != raw['identity'] or raw['key'] != spec.identity:
        raise ValueError('Artifact manifest identity changed')
    if spec.identity != expected.identity:
        raise ValueError('Artifact input, entity order or evaluation scope is stale or mismatched')
    if not isinstance(raw['files'], dict) or not raw['files']:
        raise ValueError('Artifact has no data-only payloads')
    payloads = {}
    for name, record in raw['files'].items():
        _payload_name(name)
        member = path / name
        if member.is_symlink() or not member.is_file():
            raise ValueError('Artifact payload is not a regular local file')
        data = member.read_bytes()
        if set(record) != {'sha256', 'bytes'} or len(data) != record['bytes'] or hashlib.sha256(data).hexdigest() != record['sha256']:
            raise ValueError('Artifact payload changed or was truncated')
        payloads[name] = json.loads(data, parse_constant=lambda s: _nonfinite(s))
    return Artifact(spec, payloads, raw['identity'])


def invalidated(artifacts, changed_inputs):
    """Find stale artifacts and declared descendants after named input hash changes.

    ``changed_inputs`` maps (dependency kind, name) to current SHA-256. Unlisted
    inputs are assumed unchanged; callers must inspect all relevant live inputs.
    Descendants pin complete upstream content IDs, not only a settings key.
    """
    artifacts = tuple(artifacts)
    if len({a.identity for a in artifacts}) != len(artifacts):
        raise ValueError('Duplicate artifact content identity in dependency graph')
    for (kind, name), digest in changed_inputs.items():
        Dependency(kind, name, digest)
    stale = {a.identity for a in artifacts if any(
        (d.kind, d.name) in changed_inputs and changed_inputs[(d.kind, d.name)] != d.sha256
        for d in a.spec.dependencies)}
    while True:
        descendants = {a.identity for a in artifacts if any(
            d.kind in {'artifact', 'model', 'calibration'} and d.sha256 in stale
            for d in a.spec.dependencies)}
        if descendants <= stale:
            return frozenset(stale)
        stale |= descendants
