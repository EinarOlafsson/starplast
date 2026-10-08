"""Immutable complete-profile targets for future frozen categorical adapters.

This prepares supplied recorded annotations; it does not fit, select settings,
admit biological truth, or convert missing annotation into a negative class.
Capacity describes whether a complete reference profile occurs in training,
not whether a model will recover it or whether the underlying function is true.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
import hashlib
import json

import numpy as np
import pandas as pd

from . import functional_domain_profiles as D
from .ground_truth import cohort_digest
from .provenance import EVIDENCE_GRADES
from .splits import Assignment, ROLES, SplitManifest


def _strings(values, name):
    if not isinstance(values, (tuple, list, np.ndarray)) or any(not isinstance(value, str) or not value for value in values):
        raise ValueError(name + ' requires an explicit string collection')
    return tuple(values)


@dataclass(frozen=True)
class ProfileRow:
    """One original gene and protected group, including unavailable truth states."""

    gene_id: str
    protected_group: str
    status: str
    profile: str | None
    recorded_identifiers: tuple[str, ...]
    malformed_tokens: tuple[str, ...]

    @property
    def eligible(self) -> bool:
        return self.profile is not None


@dataclass(frozen=True)
class ProfileSource:
    """Original assignment lineage, distinct from current nomenclature versions."""

    source_ids: tuple[str, ...]
    source_release: str
    evidence_grade: str
    nomenclature_releases: tuple[str, ...]
    negative_semantics: str = D.NEGATIVE_SEMANTICS
    biological_admission: bool = False


@dataclass(frozen=True)
class ProfileTargetContract:
    """Whole-universe truth and an exact protected split, copied into immutable data."""

    organism: str
    source_target: str
    target: str
    rows: tuple[ProfileRow, ...]
    source: ProfileSource
    split: SplitManifest

    def __post_init__(self):
        if (self.source_target not in D.TARGET_PATTERNS or not isinstance(self.target, str) or not self.target or
                not isinstance(self.split, SplitManifest) or self.split.organism != self.organism):
            raise ValueError('Explicit matching source/target and typed organism split are required')
        if not isinstance(self.split.assignments, tuple) or not all(isinstance(a, Assignment) for a in self.split.assignments):
            raise ValueError('Split assignments must be immutable typed records')
        if not isinstance(self.rows, tuple) or not self.rows or not all(isinstance(row, ProfileRow) for row in self.rows):
            raise ValueError('A nonempty immutable whole-gene universe is required')
        if not isinstance(self.source, ProfileSource):
            raise ValueError('Immutable source lineage is required')
        source = self.source
        if (not isinstance(source.source_ids, tuple) or len(set(source.source_ids)) != len(source.source_ids) or
                not isinstance(source.nomenclature_releases, tuple) or not isinstance(source.source_release, str) or
                not source.source_release or source.evidence_grade not in EVIDENCE_GRADES or
                source.negative_semantics != D.NEGATIVE_SEMANTICS or source.biological_admission is not False):
            raise ValueError('Source lineage must retain reference-profile and unknown-negative semantics')
        _strings(source.source_ids, 'Source IDs')
        _strings(source.nomenclature_releases, 'Nomenclature releases')
        genes = set()
        for row in self.rows:
            if (not isinstance(row.gene_id, str) or not row.gene_id.strip() or row.gene_id in genes or
                    not isinstance(row.protected_group, str) or not row.protected_group.strip()):
                raise ValueError('Unique explicit genes and protected groups are required')
            genes.add(row.gene_id)
            if not isinstance(row.recorded_identifiers, tuple) or not isinstance(row.malformed_tokens, tuple):
                raise ValueError('Recorded and malformed identifiers require immutable collections')
            _strings(row.recorded_identifiers, 'Recorded identifiers')
            _strings(row.malformed_tokens, 'Malformed source tokens')
            if list(row.recorded_identifiers) != sorted(set(row.recorded_identifiers)):
                raise ValueError('Recorded source identifiers must be sorted and unique')
            if any(not D.TARGET_PATTERNS[self.source_target].fullmatch(term) for term in row.recorded_identifiers):
                raise ValueError('Recorded identifier does not match the source namespace')
            if row.status == 'recorded_complete':
                members = D.profile_members(row.profile, self.source_target)
                canonical = json.dumps(sorted(members), separators=(',', ':'))
                if row.profile != canonical or members != set(row.recorded_identifiers) or row.malformed_tokens:
                    raise ValueError('A complete target must retain all source identifiers canonically')
            elif row.status == 'unannotated':
                if row.profile is not None or row.recorded_identifiers or row.malformed_tokens:
                    raise ValueError('Unannotated genes cannot receive target classes')
            elif row.status == 'malformed':
                if (row.profile is not None or not row.malformed_tokens or
                        any(D.TARGET_PATTERNS[self.source_target].fullmatch(token) for token in row.malformed_tokens)):
                    raise ValueError('Malformed source profiles cannot become partial truth')
            else:
                raise ValueError('Unknown complete-profile source state')
        eligible = tuple(row for row in self.rows if row.eligible)
        if tuple(row.gene_id for row in eligible) != tuple(a.entity for a in self.split.assignments):
            raise ValueError('Split must cover exactly the eligible genes in original universe order')
        if any(row.protected_group != assignment.group for row, assignment in zip(eligible, self.split.assignments)):
            raise ValueError('Split protected groups differ from the supplied whole-universe grouping')

    @property
    def universe(self) -> tuple[str, ...]:
        return tuple(row.gene_id for row in self.rows)

    @property
    def eligibility(self) -> tuple[bool, ...]:
        return tuple(row.eligible for row in self.rows)

    @property
    def cohort_identity(self) -> str:
        return cohort_digest(self.universe, self.eligibility)

    @property
    def identity(self) -> str:
        payload = json.dumps(asdict(self), sort_keys=True, separators=(',', ':'), allow_nan=False)
        return hashlib.sha256(payload.encode()).hexdigest()

    def training_labels(self) -> pd.Series:
        """Return a fresh exact training-only label series for typed fitting guards."""
        ids = self.split.entities('train')
        self.split.guard_fit('model_fit', ids)
        by_gene = {row.gene_id: row.profile for row in self.rows}
        return pd.Series([by_gene[gene] for gene in ids], index=ids, dtype=object, name=self.target)

    def capacity(self, role: str = 'test') -> dict:
        """Describe fixed training-profile support without fitting or dropping test genes.

        These truth-dependent counts are reporting data. They must not select
        settings, features, thresholds or eligible populations after a split.
        Unsupported profiles remain in all future evaluation denominators.
        """
        if role not in ROLES:
            raise ValueError('Declare a frozen split role')
        ids = self.split.entities(role)
        by_gene = {row.gene_id: row for row in self.rows}
        train_counts = Counter(by_gene[gene].profile for gene in self.split.entities('train'))
        counts = Counter(by_gene[gene].profile for gene in ids)
        unsupported = {profile: count for profile, count in sorted(counts.items()) if profile not in train_counts}
        n = len(ids)
        unsupported_genes = sum(unsupported.values())
        return {'organism': self.organism, 'source_target': self.source_target, 'target': self.target,
                'benchmark_id': self.split.benchmark_id, 'split_identity': self.split.identity,
                'target_identity': self.identity, 'cohort_identity': self.cohort_identity,
                'role': role, 'feature_access': self.split.feature_access, 'group_kind': self.split.group_kind,
                'whole_universe': len(self.rows), 'eligible_whole_universe': sum(self.eligibility),
                'unknown_states': dict(sorted(Counter(row.status for row in self.rows if not row.eligible).items())),
                'role_genes': n, 'role_protected_groups': len({by_gene[gene].protected_group for gene in ids}),
                'training_profile_counts': dict(sorted(train_counts.items())), 'role_profile_counts': dict(sorted(counts.items())),
                'supported_genes': n - unsupported_genes, 'unsupported_genes': unsupported_genes,
                'unsupported_profile_counts': unsupported,
                'training_profile_support_fraction': (n - unsupported_genes) / n if n else None,
                'model_accuracy': None, 'biological_accuracy': None, 'biological_admission': False,
                'negative_semantics': self.source.negative_semantics,
                'source': asdict(self.source),
                'interpretation': 'Training support for whole recorded profiles only; no model fit or function accuracy measured'}


def freeze_target(profiles, *, organism, source_target, universe, protected_groups,
                  split, target=None) -> ProfileTargetContract:
    """Copy supplied recorded profiles with exact universe, source and split validation."""
    required = {'organism', 'gene_id', 'source_target', 'recorded_identifiers', 'profile', 'eligible',
                'status', 'malformed_tokens', 'source_ids', 'source_release', 'evidence_grade',
                'negative_semantics', 'biological_admission'}
    if not isinstance(profiles, pd.DataFrame) or profiles.columns.has_duplicates or not required <= set(profiles.columns):
        raise ValueError('A complete recorded-profile dataframe is required')
    universe = _strings(universe, 'Whole-gene universe')
    protected_groups = _strings(protected_groups, 'Protected groups')
    if len(universe) != len(protected_groups) or tuple(profiles.gene_id) != universe:
        raise ValueError('Supplied whole-gene universe/order and grouping must match every source row')
    if not isinstance(split, SplitManifest):
        raise TypeError('A typed frozen split is required')
    source_keys, rows, nomenclature = set(), [], set()
    for record, group in zip(profiles.to_dict('records'), protected_groups):
        if record['organism'] != organism or record['source_target'] != source_target:
            raise ValueError('Source organism/target address differs')
        eligible = record['eligible']
        if not isinstance(eligible, (bool, np.bool_)):
            raise ValueError('Eligibility must be explicitly boolean')
        profile = record['profile']
        if not eligible and (profile is None or pd.isna(profile)):
            profile = None
        if bool(eligible) != (record['status'] == 'recorded_complete') or bool(eligible) != isinstance(profile, str):
            raise ValueError('Eligibility, complete state and profile disagree')
        source_ids = _strings(record['source_ids'], 'Source IDs')
        if len(set(source_ids)) != len(source_ids):
            raise ValueError('Source identifiers cannot repeat')
        source_keys.add((source_ids, record['source_release'], record['evidence_grade'],
                         record['negative_semantics'], record['biological_admission']))
        nomenclature.update(_strings(record.get('nomenclature_releases', ()), 'Nomenclature releases'))
        rows.append(ProfileRow(record['gene_id'], group, record['status'], profile,
                               _strings(record['recorded_identifiers'], 'Recorded identifiers'),
                               _strings(record['malformed_tokens'], 'Malformed source tokens')))
    if len(source_keys) != 1:
        raise ValueError('One explicit source lineage is required for the frozen target')
    source_ids, release, grade, negative, admission = next(iter(source_keys))
    source = ProfileSource(source_ids, release, grade, tuple(sorted(nomenclature)), negative, admission)
    name = source_target + '_complete_profile' if target is None else target
    return ProfileTargetContract(organism, source_target, name, tuple(rows), source, split)
