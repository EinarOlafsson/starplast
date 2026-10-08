"""Synthetic control arithmetic retains complete profiles and verifies frozen packets."""
from dataclasses import asdict, replace
import hashlib
import json

import numpy as np
import pandas as pd
import pytest

from scripts import freeze_functional_profile_controls as C
from starplast import functional_domain_profiles as D, functional_profile_targets as F, organisms as O
from starplast.splits import Assignment, SplitManifest


def fixture():
    genes = tuple(f'fixture_{i}' for i in range(10))
    groups = ('train-a', 'train-b', 'tune-a', 'cal-a', 'test-a', 'test-b', 'test-c', 'train-b', 'unknown-a', 'unknown-b')
    values = ('PF00001', 'PF00001;PF00002', 'PF00001', 'PF00002', 'PF00003',
              'PF00003;PF00004', 'PF00001', 'PF00002;PF00001', None, 'PF00005;bad')
    roles = ('train', 'train', 'tune', 'calibration', 'test', 'test', 'test', 'train')
    profiles = pd.DataFrame([
        {'organism': O.TOXOPLASMA, 'gene_id': gene, 'source_target': 'pfam_id',
         **D.parse_source_profile(value, 'pfam_id'), 'source_ids': ['synthetic-fixture'],
         'source_release': 'fixture-release', 'evidence_grade': 'synthetic_control',
         'nomenclature_releases': ['Pfam:fixture-names'], 'negative_semantics': D.NEGATIVE_SEMANTICS,
         'biological_admission': False} for gene, value in zip(genes, values)])
    split = SplitManifest(O.TOXOPLASMA, 'control-fixture', 'homology',
                          tuple(Assignment(gene, group, role) for gene, group, role in zip(genes, groups, roles)), 7)
    return F.freeze_target(profiles, organism=O.TOXOPLASMA, source_target='pfam_id', universe=genes,
                           protected_groups=groups, split=split, target='fixture_complete_profile')


def packet(path, target, *, change=None):
    scope = {'target_identity': target.identity, 'cohort_identity': target.cohort_identity,
             'split_identity': target.split.identity}
    split = asdict(target.split)
    if change == 'scope':
        scope['target_identity'] = 'changed'
    elif change == 'split':
        split['seed'] += 1
    payloads = {'normalized_target.json': asdict(target), 'preparation_scope.json': scope, 'split.json': split}
    for name, data in payloads.items():
        (path / name).write_text(json.dumps(data))
    names = [name for name in payloads if not (change == 'omit' and name == 'split.json')]
    lines = [f'{C._sha(path / name)}  {name}' for name in names]
    if change == 'escape':
        outside = path.parent / 'outside-packet.json'
        outside.write_text('{}')
        lines.append(f'{C._sha(outside)}  ../outside-packet.json')
    elif change == 'symlink':
        (path / 'linked.json').symlink_to(path / 'split.json')
        lines.append(f'{C._sha(path / "linked.json")}  linked.json')
    (path / 'SHA256SUMS.txt').write_text('\n'.join(lines) + '\n')
    return C._sha(path / 'SHA256SUMS.txt')


def test_complete_profiles_unknowns_and_unsupported_test_rows_remain():
    target = fixture()
    rows, scores, cards, classes = C.build_controls(target, seed=7)
    assert rows.entity.tolist() == ['fixture_4', 'fixture_5', 'fixture_6']
    assert rows.training_supported.tolist() == [False, False, True]
    assert rows.majority.tolist() == ['["PF00001","PF00002"]'] * 3
    universe = sorted({row.profile for row in target.rows if row.eligible})
    for name, frame in scores.items():
        assert frame.columns.tolist() == universe
        assert frame.index.tolist() == rows.entity.tolist()
        assert frame['["PF00003","PF00004"]'].eq(0).all()
        assert frame['["PF00003"]'].eq(0).all()
        expected = np.array([float(label == '["PF00001","PF00002"]') if name == 'training_majority'
                             else {'["PF00001"]': 1 / 3, '["PF00001","PF00002"]': 2 / 3}.get(label, 0.0)
                             for label in universe])
        np.testing.assert_array_equal(frame.to_numpy(), np.broadcast_to(expected, frame.shape))
        raw = cards[name]['card']
        assert raw['extra']['unknown_states'] == {'malformed': 1, 'unannotated': 1}
        assert raw['extra']['unsupported_test_genes'] == 2
        assert raw['extra']['whole_universe'] == 10
        assert raw['counts']['eligible'] == 3 and raw['counts']['abstained'] == 0
        assert sum(cell['count'] for cell in raw['extra']['confusion']) == 3
        assert len(classes[name]) == 5
        assert raw['scope']['truth_grade'] == 'synthetic_control'
        assert raw['extra']['biological_accuracy'] is None
        assert raw['extra']['classifier_or_feature_fitting_performed'] is False
        assert raw['extra']['test_support_used_for_selection'] is False


def test_literal_macro_metrics_include_unsupported_truth_and_unknown_recall():
    rows, scores, cards, classes = C.build_controls(fixture(), seed=7)
    majority = cards['training_majority']['card']
    assert majority['metrics']['accuracy'] == 0
    assert majority['metrics']['coverage'] == 1
    assert majority['metrics']['macro_recall'] == 0
    assert majority['metrics']['macro_auroc'] == 0.5
    assert majority['metrics']['macro_auprc'] == float(np.mean([1 / 3, 1 / 3, 1 / 3]))
    by_class = {item['profile']: item for item in classes['training_majority']}
    assert by_class['["PF00001","PF00002"]']['false_positive'] == 3
    assert by_class['["PF00001","PF00002"]']['recall'] is None
    assert by_class['["PF00003","PF00004"]']['false_negative'] == 1
    assert by_class['["PF00003","PF00004"]']['recall'] == 0
    assert scores['training_prevalence']['["PF00001","PF00002"]'].eq(2 / 3).all()
    json.dumps(cards, allow_nan=False)
    json.dumps(classes, allow_nan=False)


def test_fixed_seed_and_changed_test_truth_cannot_select_predictions_or_priors():
    target = fixture()
    first = C.build_controls(target, seed=7)
    replay = C.build_controls(target, seed=7)
    pd.testing.assert_frame_equal(first[0], replay[0], check_exact=True)
    assert first[2:] == replay[2:]
    expected = np.random.default_rng(7).choice(
        ['["PF00001","PF00002"]', '["PF00001"]'], size=3, p=[2 / 3, 1 / 3])
    assert first[0].prevalence_call.tolist() == expected.tolist()
    changed_rows = list(target.rows)
    changed_rows[4] = replace(changed_rows[4], profile='["PF99999"]', recorded_identifiers=('PF99999',))
    changed = replace(target, rows=tuple(changed_rows))
    assert changed.identity != target.identity and changed.split.identity == target.split.identity
    second = C.build_controls(changed, seed=7)
    pd.testing.assert_frame_equal(first[0][['entity', 'majority', 'prevalence_call']],
                                  second[0][['entity', 'majority', 'prevalence_call']], check_exact=True)
    for name in first[1]:
        common = sorted(set(first[1][name]) & set(second[1][name]))
        pd.testing.assert_frame_equal(first[1][name][common], second[1][name][common], check_exact=True)
        assert second[1][name]['["PF99999"]'].eq(0).all()


@pytest.mark.parametrize('seed', [True, None, 1.0, '7'])
def test_noninteger_control_seeds_are_refused(seed):
    with pytest.raises(ValueError, match='frozen integer'):
        C.build_controls(fixture(), seed=seed)


def test_untyped_or_foreign_organism_targets_are_refused():
    with pytest.raises(ValueError, match='typed Tg Pfam'):
        C.build_controls({'organism': 'Tg'})
    target = fixture()
    other = replace(target, organism=O.FALCIPARUM, split=replace(target.split, organism=O.FALCIPARUM))
    with pytest.raises(ValueError, match='typed Tg Pfam'):
        C.build_controls(other)


def test_hidden_truth_cannot_be_supplied_to_baseline_training(monkeypatch):
    target = fixture()
    original = F.ProfileTargetContract.training_labels

    def contaminated_training(contract):
        return pd.concat([original(contract), pd.Series({'fixture_4': '["PF00003"]'})])

    monkeypatch.setattr(F.ProfileTargetContract, 'training_labels', contaminated_training)
    with pytest.raises(ValueError, match='Forbidden'):
        C.build_controls(target, seed=7)


def test_synthetic_packet_roundtrip_retains_typed_identities(tmp_path):
    target = fixture()
    manifest_sha = packet(tmp_path, target)
    restored, scope, receipts = C.load_snapshot(tmp_path, expected_manifest_sha256=manifest_sha)
    assert restored == target and restored.identity == target.identity
    assert scope['cohort_identity'] == target.cohort_identity
    assert len(receipts) == 4
    assert hashlib.sha256((tmp_path / 'SHA256SUMS.txt').read_bytes()).hexdigest() == manifest_sha


@pytest.mark.parametrize('change, message', [('scope', 'semantic identity'), ('split', 'split file'),
                                           ('omit', 'omits required'), ('escape', 'member mismatch'),
                                           ('symlink', 'member mismatch')])
def test_packet_scope_split_manifest_and_path_failures_are_explicit(tmp_path, change, message):
    manifest_sha = packet(tmp_path, fixture(), change=change)
    with pytest.raises(ValueError, match=message):
        C.load_snapshot(tmp_path, expected_manifest_sha256=manifest_sha)


def test_unpinned_or_changed_packet_bytes_are_refused(tmp_path):
    manifest_sha = packet(tmp_path, fixture())
    with pytest.raises(ValueError, match='pinned snapshot'):
        C.load_snapshot(tmp_path, expected_manifest_sha256='0' * 64)
    (tmp_path / 'normalized_target.json').write_text('{}')
    with pytest.raises(ValueError, match='member mismatch'):
        C.load_snapshot(tmp_path, expected_manifest_sha256=manifest_sha)
