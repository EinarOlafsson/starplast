"""Synthetic complete-profile target preparation preserves split and unknown truth."""
from dataclasses import FrozenInstanceError, replace

import pandas as pd
import pytest

from starplast import functional_domain_profiles as D, functional_profile_targets as F, organisms as O
from starplast.ground_truth import cohort_digest
from starplast.splits import Assignment, SplitManifest


def fixture():
    universe = tuple(f'fixture_{i}' for i in range(10))
    groups = ('train-a', 'train-b', 'tune-a', 'cal-a', 'test-a', 'test-b', 'test-c', 'train-b', 'unknown-a', 'unknown-b')
    source = ('PF00001', 'PF00001;PF00002', 'PF00001', 'PF00002', 'PF00003',
              'PF00003;PF00004', 'PF00001', 'PF00002;PF00001', None, 'PF00005;bad')
    roles = ('train', 'train', 'tune', 'calibration', 'test', 'test', 'test', 'train')
    records = []
    for gene, value in zip(universe, source):
        parsed = D.parse_source_profile(value, 'pfam_id')
        records.append({'organism': O.TOXOPLASMA, 'gene_id': gene, 'source_target': 'pfam_id', **parsed,
                        'source_ids': ['synthetic-fixture'], 'source_release': 'fixture-release',
                        'evidence_grade': 'synthetic_control', 'nomenclature_releases': ['Pfam:fixture-names'],
                        'negative_semantics': D.NEGATIVE_SEMANTICS, 'biological_admission': False})
    profiles = pd.DataFrame(records)
    split = SplitManifest(O.TOXOPLASMA, 'complete_profile_fixture', 'homology',
                          tuple(Assignment(gene, group, role) for gene, group, role in zip(universe, groups, roles)), 7)
    return profiles, universe, groups, split


def freeze(profiles, universe, groups, split):
    return F.freeze_target(profiles, organism=O.TOXOPLASMA, source_target='pfam_id', universe=universe,
                           protected_groups=groups, split=split, target='fixture_complete_domain_profile')


def test_whole_universe_unknown_states_and_exact_group_split_are_retained():
    data, universe, groups, split = fixture()
    target = freeze(data, universe, groups, split)
    assert target.universe == universe
    assert target.eligibility == (True,) * 8 + (False, False)
    assert target.rows[-1].recorded_identifiers == ('PF00005',)
    assert target.rows[-1].profile is None and target.rows[-1].malformed_tokens == ('bad',)
    assert tuple(row.protected_group for row in target.rows) == groups
    assert target.split == split and target.cohort_identity == cohort_digest(universe, target.eligibility)
    assert target.source.source_release == 'fixture-release'
    assert target.source.nomenclature_releases == ('Pfam:fixture-names',)
    assert target.source.biological_admission is False


def test_capacity_uses_complete_profiles_and_keeps_unsupported_test_genes():
    target = freeze(*fixture())
    card = target.capacity()
    assert card['role_genes'] == 3 and card['supported_genes'] == 1 and card['unsupported_genes'] == 2
    assert card['unsupported_profile_counts'] == {'["PF00003"]': 1, '["PF00003","PF00004"]': 1}
    assert card['training_profile_counts'] == {'["PF00001"]': 1, '["PF00001","PF00002"]': 2}
    assert card['training_profile_support_fraction'] == 1 / 3
    assert card['unknown_states'] == {'malformed': 1, 'unannotated': 1}
    assert card['model_accuracy'] is None and card['biological_accuracy'] is None
    assert card['role_protected_groups'] == 3 and card['group_kind'] == 'homology'
    assert target.capacity('calibration')['unsupported_genes'] == 1
    assert target.capacity('train')['unsupported_genes'] == 0
    assert 'absence' in card['negative_semantics']


def test_training_labels_are_ordered_train_only_and_do_not_borrow_hidden_truth():
    target = freeze(*fixture())
    labels = target.training_labels()
    assert tuple(labels.index) == target.split.entities('train')
    assert labels.tolist() == ['["PF00001"]', '["PF00001","PF00002"]', '["PF00001","PF00002"]']
    labels.iloc[0] = '["PF99999"]'
    assert target.training_labels().iloc[0] == '["PF00001"]'
    with pytest.raises(ValueError, match='Forbidden'):
        target.split.guard_fit('model_fit', target.split.entities('test'))


def test_input_and_exported_count_mutation_cannot_change_immutable_target():
    profiles, universe, groups, split = fixture()
    target = freeze(profiles, universe, groups, split)
    identity = target.identity
    profiles.at[0, 'profile'] = '["PF99999"]'
    profiles.at[0, 'recorded_identifiers'].append('PF99999')
    profiles.at[0, 'source_ids'].append('other-source')
    card = target.capacity();card['training_profile_counts'].clear();card['source']['source_ids'] = ('other-source',)
    assert target.identity == identity and target.rows[0].recorded_identifiers == ('PF00001',)
    assert target.source.source_ids == ('synthetic-fixture',)
    with pytest.raises(FrozenInstanceError):
        target.rows[0].profile = '["PF99999"]'


@pytest.mark.parametrize('change', ['reordered', 'omitted_unknown', 'extra_gene', 'duplicate_gene', 'wrong_organism',
                                  'wrong_source_target', 'incomplete_split', 'split_order', 'group_mismatch', 'mutable_split'])
def test_source_universe_and_split_address_mismatches_are_refused(change):
    profiles, universe, groups, split = fixture()
    if change == 'reordered':
        profiles = profiles.iloc[::-1]
    elif change == 'omitted_unknown':
        profiles = profiles.iloc[:-1]
    elif change == 'extra_gene':
        universe += ('extra',);groups += ('extra',)
    elif change == 'duplicate_gene':
        profiles.at[1, 'gene_id'] = universe[0]
        universe = tuple(profiles.gene_id)
    elif change == 'wrong_organism':
        profiles.at[0, 'organism'] = O.FALCIPARUM
    elif change == 'wrong_source_target':
        profiles.at[0, 'source_target'] = 'interpro_id'
    elif change == 'incomplete_split':
        split = replace(split, assignments=split.assignments[:-1])
    elif change == 'split_order':
        split = replace(split, assignments=split.assignments[::-1])
    elif change == 'group_mismatch':
        groups = ('wrong-group',) + groups[1:]
    else:
        split = replace(split, assignments=list(split.assignments))
    with pytest.raises(ValueError):
        freeze(profiles, universe, groups, split)


@pytest.mark.parametrize('change', ['numeric_eligibility', 'partial_malformed', 'unknown_as_class', 'first_domain_only',
                                  'uncanonical_json', 'unresolved_as_negative', 'source_lineage_mix', 'biological_admission'])
def test_partial_classes_or_invented_absence_cannot_enter_target(change):
    profiles, universe, groups, split = fixture()
    if change == 'numeric_eligibility':
        profiles['eligible'] = profiles.eligible.astype(int)
    elif change == 'partial_malformed':
        profiles.at[9, 'profile'] = '["PF00005"]';profiles.at[9, 'eligible'] = True
        profiles.at[9, 'status'] = 'recorded_complete'
    elif change == 'unknown_as_class':
        profiles.at[8, 'profile'] = '[]';profiles.at[8, 'eligible'] = True
        profiles.at[8, 'status'] = 'recorded_complete'
    elif change == 'first_domain_only':
        profiles.at[1, 'profile'] = '["PF00001"]'
    elif change == 'uncanonical_json':
        profiles.at[1, 'profile'] = '["PF00001", "PF00002"]'
    elif change == 'unresolved_as_negative':
        profiles['negative_semantics'] = 'Missing domains are verified negatives'
    elif change == 'source_lineage_mix':
        profiles.at[0, 'source_release'] = 'other-release'
    else:
        profiles['biological_admission'] = True
    with pytest.raises(ValueError):
        freeze(profiles, universe, groups, split)


def test_typed_split_and_complete_frame_are_required():
    data, universe, groups, split = fixture()
    with pytest.raises(TypeError):
        freeze(data, universe, groups, {'identity': split.identity})
    with pytest.raises(ValueError):
        freeze(data.drop(columns='recorded_identifiers'), universe, groups, split)
    with pytest.raises(ValueError):
        freeze(data, 'gene-list', groups, split)
    with pytest.raises(ValueError):
        freeze(data, universe, groups[:-1], split)
    with pytest.raises(ValueError):
        freeze(*fixture()).capacity('held_out')


def test_encoded_target_name_cannot_be_explicitly_empty():
    profiles, universe, groups, split = fixture()
    with pytest.raises(ValueError):
        F.freeze_target(profiles, organism=O.TOXOPLASMA, source_target='pfam_id', universe=universe,
                        protected_groups=groups, split=split, target='')


def test_valid_identifier_cannot_be_declared_malformed_to_hide_a_source_gene():
    profiles, universe, groups, split = fixture()
    profiles.at[9, 'malformed_tokens'] = ['PF00006']
    with pytest.raises(ValueError):
        freeze(profiles, universe, groups, split)
