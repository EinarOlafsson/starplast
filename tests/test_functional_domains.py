"""Domain metadata names preserve version, missingness and assignment boundaries."""
import json
import hashlib
from pathlib import Path

import pytest

from starplast import functional_domains as F

IPR = 'ENTRY_AC\tENTRY_TYPE\tENTRY_NAME\nIPR000001\tDomain\tKringle\n'
PFAM = '# STOCKHOLM 1.0\n#=GF ID   Kinase\n#=GF AC   PF00001.5\n#=GF DE   Kinase family\n#=GF TP   Family\n//\n'
DEAD = '# STOCKHOLM 1.0\n#=GF ID   old_name\n#=GF AC   PF00002\n#=GF KL   This family has been killed\n#=GF FW   PF00001\n#=GF CC   Split families require new assignment\n//\n'


def _snapshot():
    return F.make_snapshot(['IPR000001', 'IPR000002', 'PF00001.3', 'PF00002'],
        F.parse_interpro(IPR), F.parse_pfam(PFAM), retired_pfam=F.parse_retired_pfam(DEAD),
        sources={key: {'release': 'test-release', 'license': 'CC0-1.0', 'license_url': 'official-license'}
                 for key in ('InterPro', 'Pfam')})


def test_interpro_names_pair_by_identifier_and_exact_text():
    entries = F.parse_interpro(IPR)
    assert entries['IPR000001'].name == 'Kringle'
    assert entries['IPR000001'].entry_type == 'Domain'
    assert F.parse_interpro_names('IPR000001\\tKringle\n') == {'IPR000001': 'Kringle'}
    assert F.parse_interpro_names('IPR000001\tKringle\n') == {'IPR000001': 'Kringle'}


@pytest.mark.parametrize('text', ['', IPR.replace('ENTRY_AC', 'accession'), IPR + IPR.splitlines()[1] + '\n',
                                  IPR.replace('\tKringle', '\t'), IPR + 'broken\n'])
def test_invalid_interpro_metadata_is_refused(text):
    with pytest.raises(ValueError):
        F.parse_interpro(text)


@pytest.mark.parametrize('text', ['IPR000001\\t\n', 'IPR00001\\tInvalid\n',
                                  'IPR000001\\tKringle\nIPR000001\\tOther\n'])
def test_invalid_names_metadata_is_refused(text):
    with pytest.raises(ValueError):
        F.parse_interpro_names(text)


@pytest.mark.parametrize('text', [PFAM.replace('//\n', ''), PFAM + PFAM,
                                  PFAM.replace('#=GF TP   Family\n', ''),
                                  PFAM.replace('PF00001.5', 'bad'),
                                  PFAM.replace('//', 'protein ACDE\n//')])
def test_incomplete_or_sequence_pfam_content_is_refused(text):
    with pytest.raises(ValueError):
        F.parse_pfam(text)


def test_pfam_names_types_and_accession_versions_are_distinct():
    term = F.parse_pfam(PFAM)['PF00001']
    assert (term.name, term.entry_type, term.accession_version) == ('Kinase family', 'Family', '5')
    lookup = F.DomainLookup(_snapshot())
    old = lookup.lookup('PF00001.3')
    assert old['requested_identifier'] == 'PF00001.3'
    assert old['accession_version_status'] == 'different_from_current'
    assert lookup.lookup('PF00001.5')['accession_version_status'] == 'matches_current'
    assert lookup.lookup('PF00001')['accession_version_status'] == 'original_version_not_recorded'


def test_withdrawn_forwarding_is_history_not_new_membership_or_name():
    snap = _snapshot()
    old = snap['terms']['PF00002']
    assert old['status'] == 'withdrawn'
    assert old['name'] is None and old['entry_type'] is None
    assert old['retirement_history']['forwarding_identifiers'] == ['PF00001']
    assert set(snap['terms']) == {'IPR000001', 'IPR000002', 'PF00001.3', 'PF00002'}
    assert 'PF00001' not in snap['terms']


def test_missing_current_and_outside_snapshot_are_different_unknowns():
    lookup = F.DomainLookup(_snapshot())
    missing = lookup.lookup('IPR000002')
    assert missing['status'] == 'missing_current_metadata' and missing['name'] is None
    assert lookup.lookup('IPR999999')['status'] == 'outside_snapshot'
    assert lookup.lookup('IPR000001.4')['status'] == 'invalid_identifier'
    assert lookup.lookup('pf00001')['status'] == 'invalid_identifier'


@pytest.mark.parametrize('text', [DEAD.replace('//\n', ''), DEAD + DEAD,
                                  DEAD.replace('#=GF KL   This family has been killed\n', '')])
def test_retired_descriptors_require_complete_unique_records(text):
    with pytest.raises(ValueError):
        F.parse_retired_pfam(text)


def test_retired_current_collision_is_refused():
    with pytest.raises(ValueError, match='overlap'):
        F.make_snapshot(['PF00001'], {}, F.parse_pfam(PFAM), sources={}, retired_pfam={'PF00001': {}})


def test_shipped_snapshot_is_offline_and_keeps_unknown_names():
    lookup = F.shipped()
    snap = json.loads((Path(F.__file__).parent / 'data/functional_domain_names.json').read_text())
    assert lookup.lookup('IPR000001')['name'] == 'Kringle'
    assert lookup.lookup('IPR000001')['ontology_release'] == '110.0'
    assert all(lookup.lookup(key)['metadata_license'] == 'CC0-1.0' for key in snap['terms'])
    for key, term in snap['terms'].items():
        metadata = lookup.lookup(key)
        assert metadata['identifier'] == term['identifier']
        if term['status'] != 'current_metadata':
            assert metadata['name'] is None and metadata['entry_type'] is None


@pytest.mark.parametrize('change', ['missing_sources', 'missing_release', 'bad_identifier', 'bad_source',
                                  'bad_status', 'missing_name', 'invented_unknown_name', 'no_retirement_history'])
def test_malformed_snapshot_records_are_refused(change):
    snapshot = _snapshot()
    if change == 'missing_sources':
        del snapshot['sources']
    elif change == 'missing_release':
        del snapshot['sources']['InterPro']['release']
    elif change == 'bad_identifier':
        snapshot['terms']['IPR000001']['identifier'] = 'IPR000002'
    elif change == 'bad_source':
        snapshot['terms']['IPR000001']['source'] = 'Pfam'
    elif change == 'bad_status':
        snapshot['terms']['IPR000001']['status'] = 'verified_function'
    elif change == 'missing_name':
        snapshot['terms']['IPR000001']['name'] = None
    elif change == 'invented_unknown_name':
        snapshot['terms']['IPR000002']['name'] = 'Guessed kinase'
    else:
        del snapshot['terms']['PF00002']['retirement_history']
    with pytest.raises(ValueError):
        F.DomainLookup(snapshot)


@pytest.mark.parametrize('snapshot', [None, [], {'schema': 2, 'terms': {}, 'sources': {}},
                                     {'schema': 1, 'terms': []}])
def test_incompatible_snapshot_schema_is_refused(snapshot):
    with pytest.raises(ValueError):
        F.DomainLookup(snapshot)


def test_packaged_snapshot_integrity_is_verified(monkeypatch):
    F.shipped.cache_clear()
    monkeypatch.setattr(Path, 'read_bytes', lambda path: b'{"changed":true}')
    try:
        with pytest.raises(ValueError, match='SHA256 mismatch'):
            F.shipped()
    finally:
        F.shipped.cache_clear()


def test_missing_packaged_snapshot_stays_unavailable(monkeypatch):
    F.shipped.cache_clear()

    def missing(path):
        raise FileNotFoundError('packaged nomenclature is missing')

    monkeypatch.setattr(Path, 'read_bytes', missing)
    try:
        with pytest.raises(FileNotFoundError, match='nomenclature is missing'):
            F.shipped()
    finally:
        F.shipped.cache_clear()


def test_integrity_does_not_substitute_for_valid_snapshot_schema(monkeypatch):
    F.shipped.cache_clear()
    payload = b'{"schema":1,"terms":{}}'
    monkeypatch.setattr(Path, 'read_bytes', lambda path: payload)
    monkeypatch.setattr(F, 'SNAPSHOT_SHA256', hashlib.sha256(payload).hexdigest())
    try:
        with pytest.raises(ValueError, match='metadata sources'):
            F.shipped()
    finally:
        F.shipped.cache_clear()


def test_invalid_json_remains_unavailable_even_with_matching_hash(monkeypatch):
    F.shipped.cache_clear()
    payload = b'not-json'
    monkeypatch.setattr(Path, 'read_bytes', lambda path: payload)
    monkeypatch.setattr(F, 'SNAPSHOT_SHA256', hashlib.sha256(payload).hexdigest())
    try:
        with pytest.raises(ValueError):
            F.shipped()
    finally:
        F.shipped.cache_clear()
