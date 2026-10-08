"""Versioned biological benchmark candidates, eligibility and explicit validation gaps.

An installed value can be a classifier output, transfer or computed annotation.
Candidate eligibility records observed values; admission for independent biology
also requires verified lineage, an assayed population and a validation protocol.
This module does not fit models or alter the shipped evaluation behavior.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass, replace
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from . import organisms as O, scorecard
from .provenance import EVIDENCE_GRADES, SourceFile

SCHEMA_VERSION = 1
TASKS = (scorecard.T_LABEL, scorecard.T_RANK, scorecard.T_SET,
         scorecard.T_CLUSTER, scorecard.T_VALUES, scorecard.T_REPL)
ABSENT = frozenset({"unassigned", "unknown", "", "nan", "none", "unlabelled", "unlabeled", "<na>"})


def eligibility_mask(values, kind, *, excluded_ids=()) -> pd.Series:
    """Select observed finite values without manufacturing negatives or filling gaps.

    The index is the explicitly named entity universe. Numeric zero and boolean
    False remain observations, but do not establish biological absence. Excluded
    entities (for example ambiguous mappings) cannot enter the candidate cohort.
    """
    if kind not in {"categorical", "numeric"}:
        raise ValueError("Declare categorical or numeric eligibility")
    values = pd.Series(values)
    if values.index.has_duplicates or values.index.isna().any():
        raise ValueError("Eligibility requires unique nonmissing entity identifiers")
    if kind == "numeric":
        valid = pd.Series(np.isfinite(pd.to_numeric(values, errors="coerce").to_numpy(dtype=float, na_value=np.nan)),
                          index=values.index)
    else:
        valid = values.notna() & ~values.astype(str).str.strip().str.lower().isin(ABSENT)
    return (valid & ~values.index.isin(excluded_ids)).astype(bool)


def cohort_digest(ids, mask) -> str:
    """Pin the ordered entity universe and eligibility together, including exclusions."""
    ids, mask = list(ids), list(mask)
    if len(ids) != len(mask) or len(set(ids)) != len(ids):
        raise ValueError("Cohort must have a unique ordered universe and matching mask")
    if any(not isinstance(value, str) or not value for value in ids):
        raise ValueError("Cohort identifiers must be nonempty strings")
    if any(type(value) not in (bool, np.bool_) for value in mask):
        raise ValueError("Eligibility mask must be boolean")
    return hashlib.sha256(json.dumps(list(zip(ids, map(bool, mask))), separators=(",", ":")).encode()).hexdigest()


@dataclass(frozen=True)
class BenchmarkEntry:
    """One organism/target/task truth candidate with reproducible cohort semantics."""

    benchmark_id: str
    organism: str
    target: str
    task: str
    evaluation_unit: str
    evidence_grade: str
    grade_basis: str
    source_ids: tuple[str, ...]
    truth_file: SourceFile
    universe_file: SourceFile
    eligibility_hex: str
    cohort_sha256: str
    stored_population: int
    eligible_population: int
    measured_population: int | None
    contexts: tuple[str, ...]
    quantity_unit: str
    exclusions: tuple[str, ...]
    negative_semantics: str
    negative_values: tuple[str, ...] = ()
    status: str = "candidate"
    gaps: tuple[str, ...] = ()

    def __post_init__(self):
        if self.organism not in O.SPACES and self.organism not in O.HOST_TABLES:
            raise ValueError("Unknown explicitly declared organism")
        if self.task not in TASKS or self.evidence_grade not in EVIDENCE_GRADES:
            raise ValueError("Unknown task or evidence grade")
        if self.evaluation_unit not in {"gene", "protein", "pair", "module", "finding"}:
            raise ValueError("Declare the biological evaluation unit")
        for name in ("benchmark_id", "target", "grade_basis", "quantity_unit", "negative_semantics"):
            if not getattr(self, name):
                raise ValueError("Truth semantics must be explicit: " + name)
        if self.status not in {"candidate", "admitted"}:
            raise ValueError("Unknown benchmark admission status")
        for count in (self.stored_population, self.eligible_population):
            if type(count) is not int or count < 0:
                raise ValueError("Population counts must be nonnegative integers")
        if self.eligible_population > self.stored_population:
            raise ValueError("Eligible population exceeds stored universe")
        if self.measured_population is not None and (type(self.measured_population) is not int or self.measured_population < 0):
            raise ValueError("Assayed population must be measured or unknown")
        if not self.contexts or not self.exclusions:
            raise ValueError("Declare context and exclusions, including unresolved scope")
        try:
            packed = bytes.fromhex(self.eligibility_hex)
        except ValueError as error:
            raise ValueError("Invalid packed eligibility") from error
        if len(packed) != (self.stored_population + 7) // 8:
            raise ValueError("Packed eligibility has the wrong universe length")
        bits = np.unpackbits(np.frombuffer(packed, dtype=np.uint8))
        if int(bits[:self.stored_population].sum()) != self.eligible_population or bits[self.stored_population:].any():
            raise ValueError("Packed eligibility does not reconcile to counts")
        if len(self.cohort_sha256) != 64 or any(c not in "0123456789abcdef" for c in self.cohort_sha256):
            raise ValueError("Pin the entity-qualified cohort identity")
        if self.status == "admitted" and (self.gaps or not self.source_ids or self.measured_population is None
                                         or self.evidence_grade not in {"direct_experiment", "curation"}):
            raise ValueError("Independent biological admission requires verified truth and no unresolved gaps")

    def mask(self, ids) -> pd.Series:
        """Restore eligibility only for the original entity universe and order."""
        ids = list(ids)
        bits = np.unpackbits(np.frombuffer(bytes.fromhex(self.eligibility_hex), dtype=np.uint8))
        mask = bits[:self.stored_population].astype(bool)
        if len(ids) != self.stored_population or cohort_digest(ids, mask) != self.cohort_sha256:
            raise ValueError("Entity universe or order differs from the frozen cohort")
        return pd.Series(mask, index=ids)


@dataclass(frozen=True)
class StrategyBenchmark:
    """Task-specific references and unresolved needs for one strategy and species."""

    organism: str
    strategy: str
    task: str
    benchmark_ids: tuple[str, ...]
    gaps: tuple[str, ...]

    def __post_init__(self):
        if self.organism not in O.SPACES and self.organism not in O.HOST_TABLES:
            raise ValueError("Unknown strategy benchmark organism")
        if not self.strategy or self.task not in TASKS or not (self.benchmark_ids or self.gaps):
            raise ValueError("Every strategy task needs references or explicit gaps")


def write_registry(path, entries, strategies):
    """Write a new registry, rejecting duplicate addresses and dangling references."""
    path, entries, strategies = Path(path), tuple(entries), tuple(strategies)
    _validate(entries, strategies)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump({"schema_version": SCHEMA_VERSION, "entries": [asdict(e) for e in entries],
                   "strategies": [asdict(s) for s in strategies]}, stream, indent=2)
        stream.write("\n")


def _validate(entries, strategies):
    ids = {e.benchmark_id: e for e in entries}
    if len(ids) != len(entries) or len({(e.organism, e.target, e.task) for e in entries}) != len(entries):
        raise ValueError("Duplicate benchmark identity or address")
    if len({(s.organism, s.strategy, s.task) for s in strategies}) != len(strategies):
        raise ValueError("Duplicate strategy/task address")
    for strategy in strategies:
        for ref in strategy.benchmark_ids:
            if ref not in ids or ids[ref].organism != strategy.organism or ids[ref].task != strategy.task:
                raise ValueError("Benchmark reference must match strategy organism and task")


def read_registry(path, *, verify_files=True, relocate_snapshots_to=None):
    """Read a registry, optionally relocating explicitly copied flat snapshot files.

    Relocation requires checksum verification and retains the content/cohort
    identities. It never searches for files or changes the recorded registry.
    Conflicting basenames are refused instead of selecting an arbitrary source.
    """
    if relocate_snapshots_to is not None and not verify_files:
        raise ValueError("Snapshot relocation requires content verification")
    payload = json.loads(Path(path).read_text())
    if payload.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("Unsupported ground-truth registry version")
    entries = []
    for row in payload["entries"]:
        row = dict(row)
        for name in ("source_ids", "contexts", "exclusions", "negative_values", "gaps"):
            row[name] = tuple(row[name])
        for name in ("truth_file", "universe_file"):
            row[name] = SourceFile(**row[name])
            if relocate_snapshots_to is not None:
                root = Path(relocate_snapshots_to).resolve()
                if not root.is_dir():
                    raise ValueError("Explicit snapshot relocation directory is unavailable")
                row[name] = replace(row[name], path=str(root / Path(row[name].path).name))
        entries.append(BenchmarkEntry(**row))
    strategies = [StrategyBenchmark(row["organism"], row["strategy"], row["task"],
                                    tuple(row["benchmark_ids"]), tuple(row["gaps"])) for row in payload["strategies"]]
    _validate(entries, strategies)
    if verify_files:
        files = {}
        for entry in entries:
            for item in (entry.truth_file, entry.universe_file):
                if item.path in files and files[item.path] != item:
                    raise ValueError("Conflicting content identity for a shared snapshot")
                files[item.path] = item
        for item in files.values():
            item.verify()
        universes = {}
        for entry in entries:
            if entry.universe_file.path not in universes:
                universes[entry.universe_file.path] = json.loads(Path(entry.universe_file.path).read_text())
            entry.mask(universes[entry.universe_file.path])
    return tuple(entries), tuple(strategies)
