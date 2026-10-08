"""Pinned InterPro/Pfam names for existing annotations, without inventing assignments.

Names and types come from current public nomenclature. They describe a source
identifier; they do not verify the original gene assignment, its model release,
or biological activity. Missing identifiers stay unknown rather than silently
acquiring a replacement. Pfam accession versions are preserved in lookup results.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import json
from pathlib import Path
import re

IPR_ID = re.compile(r'IPR\d{6}')
PFAM_ID = re.compile(r'(PF\d{5})(?:\.(\d+))?')


@dataclass(frozen=True)
class DomainTerm:
    """One authoritative current entry, independent of any gene membership."""
    identifier: str
    name: str
    entry_type: str
    source: str
    accession_version: str | None = None


def parse_interpro(text: str) -> dict[str, DomainTerm]:
    """Parse InterPro entry.list; refuse malformed rows and duplicate IDs.

    Complete release counts are verified by the acquisition audit separately.
    """
    lines = text.splitlines()
    if not lines or lines[0] != 'ENTRY_AC\tENTRY_TYPE\tENTRY_NAME':
        raise ValueError('Invalid InterPro entry.list header')
    terms = {}
    for line in lines[1:]:
        if not line.strip():
            continue
        fields = line.split('\t')
        if len(fields) != 3:
            raise ValueError('Malformed InterPro metadata row')
        identifier, kind, name = fields
        if not IPR_ID.fullmatch(identifier) or not kind.strip() or not name.strip():
            raise ValueError('Incomplete InterPro metadata row')
        if identifier in terms:
            raise ValueError('Duplicate InterPro identifier')
        terms[identifier] = DomainTerm(identifier, name, kind, 'InterPro')
    return terms


def parse_interpro_names(text: str) -> dict[str, str]:
    """Parse the companion names.dat, including its literal escaped-tab delimiter."""
    names = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = re.split(r'\t|\\t', line, maxsplit=1)
        if len(fields) != 2 or not IPR_ID.fullmatch(fields[0]) or not fields[1].strip():
            raise ValueError('Malformed InterPro names row')
        if fields[0] in names:
            raise ValueError('Duplicate InterPro name identifier')
        names[fields[0]] = fields[1]
    return names


def parse_pfam(text: str) -> dict[str, DomainTerm]:
    """Read Pfam descriptor records only, retaining accessions, names and types.

    The metadata file uses Stockholm headers but contains no alignment sequences
    or fitted HMM parameters. A partial final record is refused.
    """
    terms, record = {}, None
    for line in text.splitlines():
        if line == '# STOCKHOLM 1.0':
            if record is not None:
                raise ValueError('Unterminated Pfam descriptor')
            record = {}
        elif line == '//':
            if record is None or not {'AC', 'DE', 'TP'} <= set(record):
                raise ValueError('Incomplete Pfam descriptor')
            match = PFAM_ID.fullmatch(record['AC'])
            if not match or not record['DE'].strip() or not record['TP'].strip():
                raise ValueError('Malformed Pfam descriptor')
            identifier, version = match.groups()
            if identifier in terms:
                raise ValueError('Duplicate Pfam accession')
            terms[identifier] = DomainTerm(identifier, record['DE'], record['TP'], 'Pfam', version)
            record = None
        elif line.startswith('#=GF '):
            if record is None:
                raise ValueError('Pfam metadata outside descriptor')
            fields = line.split(maxsplit=2)
            if len(fields) != 3:
                raise ValueError('Incomplete Pfam field')
            _, key, value = fields
            if key in {'AC', 'DE', 'TP'} and key in record:
                raise ValueError('Duplicate Pfam descriptor field')
            record[key] = value
        elif line.strip() and not line.startswith('#'):
            raise ValueError('Unexpected sequence/model content in Pfam metadata')
    if record is not None:
        raise ValueError('Truncated Pfam descriptor')
    return terms


def normalize_identifier(identifier: str) -> tuple[str | None, str | None]:
    """Normalize a valid Pfam versioned ID; never guess malformed identifiers."""
    if IPR_ID.fullmatch(identifier):
        return identifier, None
    match = PFAM_ID.fullmatch(identifier)
    return match.groups() if match else (None, None)


def parse_retired_pfam(text: str) -> dict[str, dict]:
    """Read retired Pfam metadata; retain history without applying forwarding IDs."""
    retired, record = {}, None
    for line in text.splitlines():
        if line == '# STOCKHOLM 1.0':
            if record is not None:
                raise ValueError('Unterminated retired Pfam descriptor')
            record = {}
        elif line == '//':
            if record is None or not {'AC', 'KL'} <= set(record):
                raise ValueError('Incomplete retired Pfam descriptor')
            identifier, version = normalize_identifier(record['AC'])
            if identifier is None or not identifier.startswith('PF') or identifier in retired:
                raise ValueError('Invalid or duplicate retired Pfam accession')
            forwards = sorted({match[0] for match in PFAM_ID.findall(record.get('FW', ''))})
            retired[identifier] = {'retired_short_name': record.get('ID'),
                                   'retirement_note': record['KL'], 'comment': record.get('CC', ''),
                                   'forwarding_identifiers': forwards}
            record = None
        elif line.startswith('#=GF '):
            if record is None:
                raise ValueError('Retired Pfam field outside descriptor')
            fields = line.split(maxsplit=2)
            if len(fields) < 2:
                raise ValueError('Malformed retired Pfam field')
            key, value = fields[1], fields[2] if len(fields) == 3 else ''
            if key in {'AC', 'ID', 'KL'} and key in record:
                raise ValueError('Duplicate retired Pfam field')
            record[key] = (record.get(key, '') + ' ' + value).strip()
        elif line.strip() and not line.startswith('#'):
            raise ValueError('Unexpected content in retired Pfam metadata')
    if record is not None:
        raise ValueError('Truncated retired Pfam descriptor')
    return retired


def make_snapshot(identifiers, interpro, pfam, *, sources: dict, retired_pfam=None) -> dict:
    """Freeze only requested identifiers with explicit current-entry/unknown status."""
    records = {}
    retired_pfam = {} if retired_pfam is None else retired_pfam
    if set(pfam) & set(retired_pfam):
        raise ValueError('Current and withdrawn Pfam identifiers overlap')
    for requested in sorted(set(identifiers)):
        identifier, version = normalize_identifier(requested)
        if identifier is None:
            raise ValueError('Snapshot identifier is malformed: ' + requested)
        term = (interpro if identifier.startswith('IPR') else pfam).get(identifier)
        records[requested] = {
            'requested_identifier': requested, 'identifier': identifier,
            'requested_accession_version': version,
            'name': term.name if term else None, 'entry_type': term.entry_type if term else None,
            'source': 'InterPro' if identifier.startswith('IPR') else 'Pfam',
            'current_accession_version': term.accession_version if term else None,
            'status': 'current_metadata' if term else 'missing_current_metadata',
            'accession_version_status': ('different_from_current' if term and version and
                                        version != term.accession_version else 'matches_current' if
                                        term and version else 'original_version_not_recorded'),
        }
        if identifier in retired_pfam:
            records[requested].update(status='withdrawn', retirement_history=retired_pfam[identifier])
    return {'schema': 1, 'sources': sources, 'terms': records,
            'interpretation': 'Current nomenclature only; original gene assignments and source versions unchanged. '
                              'Missing-current identifiers have unknown retirement/replacement status. '
                              'No new gene functions, GO mapping or accuracy inferred.'}


class DomainLookup:
    """Offline installed-ID name lookup with bounded, explicit missing coverage."""

    def __init__(self, snapshot: dict):
        if snapshot.get('schema') != 1 or not isinstance(snapshot.get('terms'), dict):
            raise ValueError('Unsupported domain metadata snapshot')
        self.sources = snapshot['sources']
        self._terms = snapshot['terms']
        self._canonical = {r['identifier']: r for r in self._terms.values()}

    @classmethod
    def read(cls, path: str | Path) -> 'DomainLookup':
        """Read a pinned snapshot without downloading or changing gene annotations."""
        return cls(json.loads(Path(path).read_text()))

    def lookup(self, requested: str) -> dict:
        """Return metadata, preserving a requested accession version and unknowns."""
        identifier, version = normalize_identifier(requested)
        record = self._terms.get(requested) or self._terms.get(identifier)
        # Snapshot may contain only a versioned source ID; allow canonical name
        # lookup without asserting that a different source version was assigned.
        if record is None and identifier:
            record = self._canonical.get(identifier)
        if record is None:
            return {'requested_identifier': requested, 'identifier': identifier,
                    'requested_accession_version': version, 'name': None, 'entry_type': None,
                    'source': 'InterPro' if identifier and identifier.startswith('IPR') else
                              'Pfam' if identifier else None, 'current_accession_version': None,
                    'accession_version_status': 'unknown',
                    'status': 'outside_snapshot' if identifier else 'invalid_identifier'}
        result = dict(record, requested_identifier=requested, requested_accession_version=version)
        if version and result['status'] == 'current_metadata':
            result['accession_version_status'] = ('matches_current' if version == result['current_accession_version']
                                                 else 'different_from_current')
        elif not version:
            result['accession_version_status'] = 'original_version_not_recorded'
        source = self.sources[result['source']]
        result.update(ontology_release=source['release'], metadata_license=source['license'],
                      metadata_license_url=source['license_url'])
        return result


@lru_cache(maxsize=1)
def shipped() -> DomainLookup:
    """Load the packaged nomenclature subset once, without network access."""
    return DomainLookup.read(Path(__file__).parent / 'data' / 'functional_domain_names.json')
