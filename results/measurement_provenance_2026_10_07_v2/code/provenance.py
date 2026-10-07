"""Trace organism-qualified outputs to source files and explicit processing steps.

Legacy matrices do not contain every assay unit, mapping audit or license. The
contract keeps those fields unresolved instead of inferring them from a column
name or turning an installed cache into a raw experiment. Mapping audits count
source entities, ambiguity and reverse cardinality before any values are joined.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass, field
import hashlib
import json
from pathlib import Path
import re

from . import organisms

SCHEMA_VERSION = 2
EVIDENCE_GRADES = frozenset({"direct_experiment", "curation", "orthology_transfer",
                             "prediction", "derived_quantity", "synthetic_control", "unresolved"})
FILE_ROLES = frozenset({"raw_input", "processed_input", "installed_cache", "mapping_reference"})


def _organism(value):
    if value not in organisms.SPACES and value not in organisms.HOST_TABLES:
        raise ValueError("Unknown explicitly declared organism")


def _hash(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


@dataclass(frozen=True)
class SourceFile:
    """Content identity of a located input, distinct from its paper association."""

    path: str
    sha256: str
    bytes: int
    role: str
    url: str = ""
    association: str = "registry_asserted_not_independently_verified"

    def __post_init__(self):
        if not self.path or not re.fullmatch(r"[a-f0-9]{64}", self.sha256):
            raise ValueError("A source file needs a path and SHA-256 content identity")
        if type(self.bytes) is not int or self.bytes < 0 or self.role not in FILE_ROLES:
            raise ValueError("Invalid byte count or file role")
        if not self.association:
            raise ValueError("State the source-file association status")

    @classmethod
    def inspect(cls, path, role, url="", association="registry_asserted_not_independently_verified"):
        """Hash a regular file; directories are not source-file identities."""
        path = Path(path)
        if not path.is_file():
            raise ValueError("A regular input file is required")
        before = path.stat()
        digest = _hash(path)
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError("Source file changed while hashing")
        return cls(str(path.resolve()), digest, after.st_size, role, url, association)

    def verify(self):
        """Refuse a moved, truncated or altered input rather than silently reusing it."""
        path = Path(self.path)
        if not path.is_file() or path.stat().st_size != self.bytes or _hash(path) != self.sha256:
            raise ValueError("Source content no longer matches its frozen identity")


@dataclass(frozen=True)
class MappingAudit:
    """A measured mapping census for explicit source/target namespaces and versions.

    Counts use unique nonempty source entities. Ambiguous sources are withheld;
    several unambiguous sources mapping to one target remain visible, because a
    later aggregation needs its own policy. A mapping does not establish assay QC.
    """

    source_organism: str
    target_organism: str
    source_namespace: str
    target_namespace: str
    source_version: str
    target_version: str
    source_entities: int
    mapped_entities: int
    ambiguous_entities: int
    unmapped_entities: int
    target_entities: int
    many_to_one_targets: int
    source_identifier_missing_rows: int
    reference: str

    def __post_init__(self):
        _organism(self.source_organism)
        _organism(self.target_organism)
        for name in ("source_namespace", "target_namespace", "source_version", "target_version", "reference"):
            if not getattr(self, name):
                raise ValueError("Mapping namespaces, versions and reference must be explicit")
        for name in ("source_entities", "mapped_entities", "ambiguous_entities", "unmapped_entities",
                     "target_entities", "many_to_one_targets", "source_identifier_missing_rows"):
            if type(getattr(self, name)) is not int or getattr(self, name) < 0:
                raise ValueError("Mapping counts must be nonnegative integers")
        if self.mapped_entities + self.ambiguous_entities + self.unmapped_entities != self.source_entities:
            raise ValueError("Mapping losses do not reconcile to source entities")
        if self.target_entities > self.mapped_entities or self.many_to_one_targets > self.target_entities:
            raise ValueError("Mapping target cardinality is inconsistent")


def audit_mapping(source_ids, pairs, *, source_organism, target_organism, source_namespace,
                  target_namespace, source_version, target_version, reference) -> MappingAudit:
    """Count candidate mappings, omitting duplicate edges and refusing ambiguous joins.

    ``pairs`` contains explicit (source_id, target_id) pairs from a verified
    reference, never accession strings manufactured from a suffix or symbol.
    Missing source IDs are counted as rows and do not become a shared entity.
    """
    ids, missing = set(), 0
    for value in source_ids:
        if value is None or not isinstance(value, str) or not value.strip():
            missing += 1
        else:
            ids.add(value)
    candidates = {value: set() for value in ids}
    for source, target in pairs:
        if source in candidates and isinstance(target, str) and target.strip():
            candidates[source].add(target)
    chosen = [next(iter(values)) for values in candidates.values() if len(values) == 1]
    counts = Counter(chosen)
    return MappingAudit(source_organism, target_organism, source_namespace, target_namespace,
                        source_version, target_version, len(ids), len(chosen),
                        sum(len(values) > 1 for values in candidates.values()),
                        sum(not values for values in candidates.values()), len(counts),
                        sum(count > 1 for count in counts.values()), missing, reference)


@dataclass(frozen=True)
class Transform:
    """One recorded processing step and its exact implementation identity."""

    name: str
    implementation: str
    code_sha256: str
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]
    parameters: dict = field(default_factory=dict)

    def __post_init__(self):
        if not self.name or not self.implementation or not re.fullmatch(r"[a-f0-9]{64}", self.code_sha256):
            raise ValueError("A transform requires a named implementation and code hash")
        if not self.inputs or not self.outputs or any(not value for value in (*self.inputs, *self.outputs)):
            raise ValueError("A transform needs input and output addresses")
        json.dumps(self.parameters, allow_nan=False)


@dataclass(frozen=True)
class MeasurementTrace:
    """One source/organism/output address with explicit unknowns and derived ancestry.

    ``storage_unit`` is gene/protein/pair/metabolite, whereas ``quantity_unit``
    describes the measurement (TPM, log2FC, etc.). Unknown is not dimensionless.
    Missing license/redistribution or mapping evidence stays an unresolved gap.
    """

    source_id: str
    registry_organism: str
    output_organism: str
    storage_unit: str
    column: str
    evidence_grade: str = "unresolved"
    quantity_unit: str = "unresolved"
    context: tuple[str, ...] = ()
    source_files: tuple[SourceFile, ...] = ()
    transforms: tuple[Transform, ...] = ()
    mappings: tuple[MappingAudit, ...] = ()
    upstream_sources: tuple[str, ...] = ()
    publication_identity: str = ""
    publication_status: str = "unresolved"
    license: str = "unresolved"
    redistribution: str = "unresolved"
    gaps: tuple[str, ...] = ()
    source_species: tuple[str, ...] = ()
    source_species_status: str = "unresolved"
    source_species_reference: str = "unresolved"
    schema_version: int = SCHEMA_VERSION

    def __post_init__(self):
        _organism(self.registry_organism)
        _organism(self.output_organism)
        if not isinstance(self.source_species, tuple) or any(not isinstance(s, str) or not s.strip() for s in self.source_species):
            raise ValueError("Source species must be an explicit tuple of scientific names")
        if self.source_species_status not in {"unresolved", "source_file_verified", "primary_metadata_verified"}:
            raise ValueError("Invalid source species association status")
        if self.source_species_status != "unresolved" and (not self.source_species or self.source_species_reference == "unresolved"):
            raise ValueError("Verified source species requires a source reference")
        if not self.source_species_reference:
            raise ValueError("Use an explicit unresolved species reference")
        if not self.source_id or not self.column or self.evidence_grade not in EVIDENCE_GRADES:
            raise ValueError("Invalid source/output address or evidence grade")
        if self.storage_unit not in {"gene", "protein", "pair", "metabolite", "source_records", "unresolved_host"}:
            raise ValueError("Invalid output storage unit")
        if type(self.schema_version) is not int or self.schema_version != SCHEMA_VERSION:
            raise ValueError("Unsupported provenance schema version")
        if not self.quantity_unit or not self.license or not self.redistribution:
            raise ValueError("Use explicit unresolved fields rather than omitted semantics")
        for values, kind in ((self.source_files, SourceFile), (self.transforms, Transform), (self.mappings, MappingAudit)):
            if not isinstance(values, tuple) or not all(isinstance(value, kind) for value in values):
                raise ValueError("Provenance components must be typed tuples")
        if self.evidence_grade == "orthology_transfer" and not self.mappings:
            if "transfer_mapping_unresolved" not in self.gaps:
                raise ValueError("Transferred outputs require mapping evidence or an explicit transfer gap")

    @property
    def trace_id(self):
        """Content identity of a trace, including gaps and all processing versions."""
        return hashlib.sha256(json.dumps(asdict(self), sort_keys=True, allow_nan=False).encode()).hexdigest()


def trace_sources(source_id, sources) -> dict:
    """Walk declared upstream columns without dropping unknown or ambiguous parents.

    Ancestors are resolved within the registry's explicit organism address; a
    same-named column in another species cannot silently become an input. Cycles
    are rejected, and unknown/multiply declared inputs remain graph gaps.
    """
    sources = list(sources)
    by_id = {source.key: source for source in sources}
    if len(by_id) != len(sources):
        raise ValueError("Source keys must be unique")
    if source_id not in by_id:
        raise ValueError("Source not registered")
    visited, active, gaps, edges = set(), set(), [], []
    def walk(key):
        if key in active:
            raise ValueError("Cycle in declared source derivations")
        if key in visited:
            return
        active.add(key)
        source = by_id[key]
        for column in source.derived_from:
            parents = [item.key for item in sources if item.organism == source.organism and column in item.columns and item.key != key]
            if len(parents) != 1:
                gaps.append({"source_id": key, "column": column, "status": "ambiguous" if parents else "unresolved",
                             "candidates": sorted(parents)})
            else:
                edges.append({"source_id": key, "input_column": column, "upstream_source": parents[0]})
                walk(parents[0])
        active.remove(key)
        visited.add(key)
    walk(source_id)
    return {"source_id": source_id, "sources": sorted(visited), "edges": edges, "gaps": gaps}


def write_traces(traces, path):
    """Write unique output addresses and content identities without lossy flattening."""
    traces = list(traces)
    addresses = [(t.source_id, t.output_organism, t.storage_unit, t.column) for t in traces]
    if len(set(addresses)) != len(addresses):
        raise ValueError("Duplicate measurement trace address")
    records = [{**asdict(trace), "trace_id": trace.trace_id} for trace in traces]
    path = Path(path)
    if path.exists():
        raise ValueError("Use a new provenance snapshot")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(records, indent=2, allow_nan=False) + "\n")
    return _hash(path)


def read_traces(path) -> list[MeasurementTrace]:
    """Validate stored trace identities and reconstruct typed source/mapping steps."""
    result = []
    for record in json.loads(Path(path).read_text()):
        identity = record.pop("trace_id")
        for name in ("context", "upstream_sources", "gaps", "source_species"):
            record[name] = tuple(record[name])
        record["source_files"] = tuple(SourceFile(**item) for item in record["source_files"])
        record["mappings"] = tuple(MappingAudit(**item) for item in record["mappings"])
        record["transforms"] = tuple(Transform(**{**item, "inputs": tuple(item["inputs"]), "outputs": tuple(item["outputs"])})
                                     for item in record["transforms"])
        trace = MeasurementTrace(**record)
        if trace.trace_id != identity:
            raise ValueError("Measurement trace content identity changed")
        result.append(trace)
    addresses = [(t.source_id, t.output_organism, t.storage_unit, t.column) for t in result]
    if len(set(addresses)) != len(addresses):
        raise ValueError("Duplicate measurement trace address")
    return result
