"""Complete recorded domain profiles for reference-recovery benchmark candidates.

The target is an unordered set of every assigned source identifier, preserving
Pfam accession versions and overlapping memberships. It is not an ordered domain
architecture, a current-active rewrite, or experimental functional truth. Missing
current names and officially retired identifiers stay in the recorded profile.
An invalid nonempty token makes the whole profile unavailable rather than silently
leaving a partial truth label. Missing annotation is never a biological negative.
"""
from __future__ import annotations

import json
import re

import pandas as pd

from .provenance import EVIDENCE_GRADES

TARGET_PATTERNS = {
    'interpro_id': re.compile(r'IPR[0-9]{6}'),
    'interpro_ids': re.compile(r'IPR[0-9]{6}'),
    'pfam_id': re.compile(r'PF[0-9]{5}(?:\.[0-9]+)?'),
    'pfam_ids': re.compile(r'PF[0-9]{5}(?:\.[0-9]+)?'),
}
NEGATIVE_SEMANTICS = ('Nonmembership in a complete recorded profile means reference-profile omission only; '
                      'unannotated/malformed profiles are unknown; biological domain absence is not adjudicated')
PROFILE_COLUMNS = ('organism', 'gene_id', 'source_target', 'source_value', 'source_tokens',
                   'recorded_identifiers', 'profile', 'eligible', 'status', 'malformed_tokens',
                   'ignored_empty_tokens', 'duplicate_identifier_tokens', 'metadata_statuses',
                   'nomenclature_releases', 'source_ids', 'source_release', 'evidence_grade',
                   'negative_semantics', 'biological_admission')
TERM_COLUMNS = ('organism', 'gene_id', 'source_target', 'identifier', 'profile_eligible',
                'metadata_status', 'current_name', 'current_entry_type', 'nomenclature_source',
                'nomenclature_release', 'metadata_license', 'accession_version_status',
                'requested_accession_version', 'current_accession_version')


def parse_source_profile(value: str | None, target: str) -> dict:
    """Validate every semicolon token before assigning a complete recorded profile.

    Empty delimiter tokens are ignored but counted, including leading/trailing
    delimiters. A delimiter-only value is explicitly unannotated. Duplicate IDs
    collapse into set membership; versions remain distinct recorded identifiers.
    No names, embedded descriptions, commas or prefix matches are accepted.
    """
    if target not in TARGET_PATTERNS:
        raise ValueError('Declare one installed InterPro or Pfam source field')
    if value is not None and not isinstance(value, str):
        raise TypeError('Source domain annotation must be text or explicitly missing')
    tokens = [] if value is None else [token.strip() for token in value.split(';')]
    nonempty = [token for token in tokens if token]
    malformed = [token for token in nonempty if not TARGET_PATTERNS[target].fullmatch(token)]
    valid = sorted({token for token in nonempty if TARGET_PATTERNS[target].fullmatch(token)})
    eligible = bool(valid) and not malformed
    return {'source_tokens': tokens, 'recorded_identifiers': valid,
            'profile': json.dumps(valid, separators=(',', ':')) if eligible else None,
            'eligible': eligible, 'status': 'recorded_complete' if eligible else 'malformed' if malformed else 'unannotated',
            'malformed_tokens': malformed, 'ignored_empty_tokens': len(tokens) - len(nonempty),
            'duplicate_identifier_tokens': len(nonempty) - len(malformed) - len(valid)}


def profile_members(value: str, target: str) -> frozenset[str]:
    """Decode only a complete canonical recorded profile, retaining accession versions."""
    if target not in TARGET_PATTERNS:
        raise ValueError('Declare the recorded domain profile source field')
    try:
        identifiers = json.loads(value)
    except (ValueError, TypeError) as exc:
        raise ValueError('Explicit complete recorded domain profile required') from exc
    if (not isinstance(identifiers, list) or not identifiers or
            any(not isinstance(term, str) or not TARGET_PATTERNS[target].fullmatch(term) for term in identifiers) or
            identifiers != sorted(set(identifiers))):
        raise ValueError('Recorded profile must contain sorted unique complete identifiers')
    return frozenset(identifiers)


def recorded_profiles(ctx, target: str, *, domain_lookup=None, source_ids=(),
                      source_release='unresolved', evidence_grade='unresolved') -> tuple[pd.DataFrame, pd.DataFrame]:
    """Retain the entire gene universe and separately disclose nomenclature status.

    Metadata availability/status cannot alter the eligibility or members of a
    valid recorded source profile. Valid tokens inside a malformed source appear
    only as diagnostic term rows with profile_eligible=False. Neither dataframe
    modifies the source table, its descriptions, original IDs or gene membership.
    Source release means the assignment/report release, separately from current
    nomenclature release; unresolved source lineage must remain explicit.
    """
    if target not in TARGET_PATTERNS or target not in ctx.nodes:
        raise ValueError('Declare an installed InterPro or Pfam source field')
    if ('gene_id' not in ctx.nodes or ctx.nodes.gene_id.isna().any() or len(set(ctx.gene_ids)) != ctx.n or
            any(not str(gene).strip() for gene in ctx.gene_ids)):
        raise ValueError('An explicit unique source gene universe is required')
    if isinstance(source_ids, str):
        raise ValueError('Source identifiers require an explicit collection, not one string')
    source_ids = tuple(source_ids)
    if any(not isinstance(source, str) or not source for source in source_ids) or len(set(source_ids)) != len(source_ids):
        raise ValueError('Source identifiers must be explicit unique strings')
    if not isinstance(source_release, str) or not source_release or evidence_grade not in EVIDENCE_GRADES:
        raise ValueError('Declare source assignment release and evidence grade, including unresolved')
    truth = ctx.truth(target)
    profiles, term_rows = [], []
    for position, gene in enumerate(ctx.gene_ids):
        value = str(truth.iloc[position]) if pd.notna(truth.iloc[position]) else None
        parsed = parse_source_profile(value, target)
        metadata_statuses, releases = {}, set()
        for identifier in parsed['recorded_identifiers']:
            metadata = domain_lookup.lookup(identifier) if domain_lookup is not None else {}
            status = metadata.get('status', 'nomenclature_unavailable')
            if status not in {'current_metadata', 'withdrawn', 'missing_current_metadata', 'outside_snapshot', 'nomenclature_unavailable'}:
                raise ValueError('Unexpected nomenclature status for a validated source identifier')
            metadata_statuses[status] = metadata_statuses.get(status, 0) + 1
            if metadata.get('ontology_release'):
                releases.add(metadata['source'] + ':' + metadata['ontology_release'])
            term_rows.append({'organism': ctx.organism, 'gene_id': str(gene), 'source_target': target,
                              'identifier': identifier, 'profile_eligible': parsed['eligible'],
                              'metadata_status': status, 'current_name': metadata.get('name'),
                              'current_entry_type': metadata.get('entry_type'), 'nomenclature_source': metadata.get('source'),
                              'nomenclature_release': metadata.get('ontology_release'),
                              'metadata_license': metadata.get('metadata_license'),
                              'accession_version_status': metadata.get('accession_version_status', 'unknown'),
                              'requested_accession_version': metadata.get('requested_accession_version'),
                              'current_accession_version': metadata.get('current_accession_version')})
        raw_value = ctx.nodes[target].iloc[position]
        profiles.append({'organism': ctx.organism, 'gene_id': str(gene), 'source_target': target,
                         'source_value': str(raw_value) if pd.notna(raw_value) else None, **parsed,
                         'metadata_statuses': json.dumps(metadata_statuses, sort_keys=True, separators=(',', ':')),
                         'nomenclature_releases': sorted(releases), 'source_ids': list(source_ids),
                         'source_release': source_release, 'evidence_grade': evidence_grade,
                         'negative_semantics': NEGATIVE_SEMANTICS, 'biological_admission': False})
    profile_frame = pd.DataFrame(profiles, columns=PROFILE_COLUMNS)
    # pandas 3 infers nullable string columns and otherwise converts None to NaN.
    # Preserve explicit missing truth/raw-cell sentinels for artifact serialization.
    for column in ('profile', 'source_value'):
        profile_frame[column] = pd.Series([row[column] for row in profiles], dtype=object)
    return profile_frame, pd.DataFrame(term_rows, columns=TERM_COLUMNS)
