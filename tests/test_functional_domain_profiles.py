"""Whole domain source profiles preserve overlaps, versions and unknown absence."""
from dataclasses import dataclass

import pandas as pd
import pytest

from starplast import functional_domain_profiles as P, organisms as O


@dataclass
class Context:
    nodes: pd.DataFrame
    organism: str = O.TOXOPLASMA

    @property
    def gene_ids(self):
        return self.nodes.gene_id.astype(str).to_numpy()

    @property
    def n(self):
        return len(self.nodes)

    def truth(self, target):
        result = self.nodes[target].astype('string').str.strip()
        return result.mask(result.str.lower().isin({'', 'nan', 'unknown', 'unassigned', 'none', '<na>'}))


class Metadata:
    def lookup(self, identifier):
        status = 'withdrawn' if identifier == 'PF00002' else 'missing_current_metadata' if identifier == 'IPR999999' else 'current_metadata'
        return {'status': status, 'source': 'Pfam' if identifier.startswith('PF') else 'InterPro',
                'ontology_release': 'current-name-release', 'metadata_license': 'CC0-1.0',
                'name': 'known name' if status == 'current_metadata' else None,
                'entry_type': 'Domain' if status == 'current_metadata' else None,
                'requested_accession_version': identifier.split('.')[1] if '.' in identifier else None,
                'current_accession_version': '8' if identifier.startswith('PF') else None,
                'accession_version_status': 'different_from_current' if '.' in identifier else 'original_version_not_recorded'}


def test_whole_profiles_preserve_overlaps_versions_and_canonical_order():
    parsed = P.parse_source_profile(' PF00069.2 ;PF00001;PF00069.12;PF00069.2 ', 'pfam_ids')
    assert parsed['profile'] == '["PF00001","PF00069.12","PF00069.2"]'
    assert parsed['eligible'] and parsed['duplicate_identifier_tokens'] == 1
    assert P.profile_members(parsed['profile'], 'pfam_ids') == {'PF00001', 'PF00069.12', 'PF00069.2'}


@pytest.mark.parametrize('value, empty_tokens', [(None, 0), ('', 1), (';', 2), (';; ;', 4)])
def test_missing_and_delimiter_only_profiles_are_unknown(value, empty_tokens):
    parsed = P.parse_source_profile(value, 'interpro_id')
    assert parsed['status'] == 'unannotated' and not parsed['eligible'] and parsed['profile'] is None
    assert parsed['ignored_empty_tokens'] == empty_tokens


def test_empty_delimiters_do_not_discard_other_assigned_identifiers():
    parsed = P.parse_source_profile(';IPR000002;; IPR000001;', 'interpro_ids')
    assert parsed['profile'] == '["IPR000001","IPR000002"]'
    assert parsed['ignored_empty_tokens'] == 3 and parsed['malformed_tokens'] == []


@pytest.mark.parametrize('value', ['IPR000001;malformed', 'xIPR000001', 'IPR000001.1',
                                  'IPR000001,IPR000002', 'IPR000001 (name)', 'ipr000001',
                                  'IPR000001;PF00001', 'IPR000001suffix', 'IPR 000001', 'IPR０００００１', 'IPR١٢٣٤٥٦'])
def test_invalid_nonempty_token_refuses_entire_profile(value):
    parsed = P.parse_source_profile(value, 'interpro_id')
    assert parsed['status'] == 'malformed' and not parsed['eligible'] and parsed['profile'] is None
    assert parsed['malformed_tokens']


@pytest.mark.parametrize('value', ['PF00001.2x', 'PF0001', 'PF00001.-1', 'PF00001.2.3', 'PF00001,PF00002'])
def test_pfam_version_syntax_is_validated_as_whole_token(value):
    assert P.parse_source_profile(value, 'pfam_id')['status'] == 'malformed'


@pytest.mark.parametrize('value', ['[]', 'null', '["IPR000002","IPR000001"]', '["IPR000001","IPR000001"]',
                                  '["IPR000001","invalid"]', '[1]', '{"IPR000001":true}', None])
def test_incomplete_noncanonical_profiles_cannot_enter_truth(value):
    with pytest.raises(ValueError):
        P.profile_members(value, 'interpro_id')


def test_withdrawn_and_missing_names_remain_recorded_membership():
    ctx = Context(pd.DataFrame({'gene_id': ['gene-a', 'gene-b', 'gene-c', 'gene-d'],
                                'pfam_id': ['PF00001.2;PF00002', None, 'PF00002;bad', 'PF00001'],
                                'pfam_desc': ['old source description', None, 'preserved malformed evidence', 'source']}))
    original = ctx.nodes.copy(deep=True)
    profiles, terms = P.recorded_profiles(ctx, 'pfam_id', domain_lookup=Metadata(), source_ids=('recorded-source',))
    assert profiles.gene_id.tolist() == ctx.gene_ids.tolist()
    assert profiles.profile.tolist() == ['["PF00001.2","PF00002"]', None, None, '["PF00001"]']
    assert profiles.status.tolist() == ['recorded_complete', 'unannotated', 'malformed', 'recorded_complete']
    assert terms[terms.identifier.eq('PF00002')].metadata_status.eq('withdrawn').all()
    assert terms[terms.gene_id.eq('gene-c')].profile_eligible.eq(False).all()
    assert profiles.source_release.eq('unresolved').all() and profiles.evidence_grade.eq('unresolved').all()
    assert profiles.biological_admission.eq(False).all()
    assert profiles.negative_semantics.str.contains('not adjudicated').all()
    pd.testing.assert_frame_equal(ctx.nodes, original, check_exact=True)
    unknown = Context(pd.DataFrame({'gene_id': ['gene-a'], 'interpro_ids': ['IPR000001;IPR999999']}))
    ipr_profiles, ipr_terms = P.recorded_profiles(unknown, 'interpro_ids', domain_lookup=Metadata())
    assert ipr_profiles.profile.iloc[0] == '["IPR000001","IPR999999"]'
    assert ipr_terms.metadata_status.tolist() == ['current_metadata', 'missing_current_metadata']


def test_metadata_availability_never_changes_truth_profiles():
    ctx = Context(pd.DataFrame({'gene_id': ['gene-a', 'gene-b'], 'pfam_id': ['PF00002', None]}))
    unavailable, terms = P.recorded_profiles(ctx, 'pfam_id')
    named, _ = P.recorded_profiles(ctx, 'pfam_id', domain_lookup=Metadata())
    columns = ['gene_id', 'recorded_identifiers', 'profile', 'eligible', 'status', 'source_value']
    pd.testing.assert_frame_equal(unavailable[columns], named[columns], check_exact=True)
    assert terms.metadata_status.tolist() == ['nomenclature_unavailable']
    assert terms.current_name.isna().all()


def test_every_raw_source_cell_is_retained_and_unknown_spellings_normalize_only_truth():
    ctx = Context(pd.DataFrame({'gene_id': ['gene-a', 'gene-b', 'gene-c'],
                                'interpro_id': [' ; IPR000001 ; ', 'unknown', None]}))
    profiles, _ = P.recorded_profiles(ctx, 'interpro_id')
    assert profiles.source_value.tolist() == [' ; IPR000001 ; ', 'unknown', None]
    assert profiles.profile.tolist() == ['["IPR000001"]', None, None]


@pytest.mark.parametrize('nodes', [pd.DataFrame({'gene_id': ['same', 'same'], 'interpro_id': [None, None]}),
                                   pd.DataFrame({'gene_id': [None], 'interpro_id': [None]}),
                                   pd.DataFrame({'gene_id': [' '], 'interpro_id': [None]})])
def test_unique_explicit_source_gene_universe_required(nodes):
    with pytest.raises(ValueError):
        P.recorded_profiles(Context(nodes), 'interpro_id')


@pytest.mark.parametrize('kwargs', [{'source_ids': 'source'}, {'source_ids': ('source', 'source')},
                                    {'source_release': ''}, {'evidence_grade': 'verified_domain_activity'}])
def test_source_scope_must_remain_explicit(kwargs):
    ctx = Context(pd.DataFrame({'gene_id': ['gene-a'], 'interpro_id': ['IPR000001']}))
    with pytest.raises(ValueError):
        P.recorded_profiles(ctx, 'interpro_id', **kwargs)


def test_no_other_annotation_or_numeric_field_is_implicitly_domain_truth():
    with pytest.raises(ValueError):
        P.parse_source_profile('2.7.11.1', 'ec_number')
    with pytest.raises(TypeError):
        P.parse_source_profile(['IPR000001'], 'interpro_id')
