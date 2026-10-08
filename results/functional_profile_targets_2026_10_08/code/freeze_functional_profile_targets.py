"""Execute one bounded complete-profile preparation census without model fitting.

The verified original-cell census remains a separate immutable dependency. This
script freezes the named source, full ordered universe, protected grouping and
eligible split before reporting training support. Normalized target identity is
distinct from raw source, cohort, split and capacity payload identities. Unknown
annotations and unresolved homology remain explicit; no biological admission,
prediction, function accuracy, source acquisition or negative absence is inferred.
"""
from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
CENSUS = ROOT / 'results/functional_domain_profile_census_2026_10_08'
PREPARATION_ID = 'FN-DOM-PREP-01-Tg-pfam-complete-profile'
SEED = 20261008
FRACTIONS = (0.55, 0.15, 0.15, 0.15)


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _identity(payload):
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def _write(path, payload):
    path.write_text(json.dumps(payload, indent=2, allow_nan=False) + '\n')


def _verified_census():
    checksums = CENSUS / 'SHA256SUMS.txt'
    dependencies = {str(checksums): _sha(checksums)}
    for line in checksums.read_text().splitlines():
        digest, relative = line.split('  ', 1)
        path = CENSUS / relative
        if path.resolve().is_relative_to(CENSUS.resolve()) is False or _sha(path) != digest:
            raise ValueError('Source census checksum mismatch: ' + relative)
        dependencies[str(path)] = digest
    return dependencies


def prepare(output):
    """Freeze verified recorded profiles and replay every saved record and capacity count."""
    import numpy as np
    import pandas as pd
    from starplast import organisms as O
    from starplast.embedding import as_text
    from starplast.functional_profile_targets import ProfileRow, ProfileSource, ProfileTargetContract, freeze_target
    from starplast.splits import Assignment, ROLES, SplitManifest, make_split

    output = Path(output)
    dependencies = _verified_census()
    census_inputs = json.loads((CENSUS / 'input_manifest.json').read_text())
    node_path = Path(O.nodes_path(O.TOXOPLASMA))
    metadata_path = ROOT / 'starplast/data/functional_domain_names.json'
    for path in (node_path, metadata_path):
        expected = census_inputs[str(path)]
        if _sha(path) != expected:
            raise ValueError('Installed source changed since census: ' + str(path))
        dependencies[str(path)] = expected
    code_paths = [Path(__file__), ROOT / 'scripts/notebook_runner.py', *[
        ROOT / 'starplast' / name for name in ('functional_profile_targets.py', 'functional_domain_profiles.py',
                                              'splits.py', 'ground_truth.py', 'organisms.py', 'provenance.py',
                                              'embedding.py', 'scorecard.py')]]
    (output / 'code').mkdir()
    for path in code_paths:
        dependencies[str(path)] = _sha(path)
        shutil.copyfile(path, output / 'code' / path.name)

    source_path = CENSUS / 'Tg_pfam_id_profiles.parquet'
    profiles = pd.read_parquet(source_path)
    # Column projection avoids loading measurements, graph edges or full nodes.
    nodes = pd.read_parquet(node_path, columns=['gene_id', 'orthogroup', 'pfam_id'])
    universe = tuple(nodes.gene_id.astype(str))
    assert tuple(profiles.gene_id) == universe and len(universe) == 8140
    assert int(profiles.eligible.sum()) == 4310
    for original, row in zip(nodes.pfam_id, profiles.itertuples(index=False)):
        raw = str(original) if pd.notna(original) else None
        assert (None if pd.isna(row.source_value) else row.source_value) == raw
        # Literal independent whole-token replay, rather than taking a partial
        # extracted term or current metadata availability as target truth.
        normalized = None if raw is None or raw.strip().lower() in {
            'unassigned', 'unknown', '', 'nan', 'none', 'unlabelled', 'unlabeled', '<na>'} else raw.strip()
        tokens = [] if normalized is None else [token.strip() for token in normalized.split(';')]
        terms = sorted({token for token in tokens if re.fullmatch(r'PF[0-9]{5}(?:\.[0-9]+)?', token)})
        invalid = [token for token in tokens if token and not re.fullmatch(r'PF[0-9]{5}(?:\.[0-9]+)?', token)]
        expected = json.dumps(terms, separators=(',', ':')) if terms and not invalid else None
        assert (None if pd.isna(row.profile) else row.profile) == expected
        assert bool(row.eligible) == (expected is not None)
        assert list(row.recorded_identifiers) == terms and list(row.malformed_tokens) == invalid

    # Same grouping convention as Context.groups: unassigned entries receive
    # unique positional singleton keys. This does not resolve missing homology.
    original_groups = as_text(nodes.orthogroup).tolist()
    groups = tuple(group if group != '' else f'__{position}' for position, group in enumerate(original_groups))
    resolved = tuple(group != '' for group in original_groups)
    fallback = {group for group, known in zip(groups, resolved) if not known}
    if fallback & {group for group, known in zip(groups, resolved) if known}:
        raise ValueError('Unresolved singleton group collides with an installed homology identifier')
    eligible = profiles.eligible.to_numpy(dtype=bool)
    split = make_split(tuple(gene for gene, keep in zip(universe, eligible) if keep),
                       tuple(group for group, keep in zip(groups, eligible) if keep),
                       organism=O.TOXOPLASMA, benchmark_id=PREPARATION_ID, group_kind='homology',
                       seed=SEED, fractions=FRACTIONS, feature_access='inductive')
    contract = freeze_target(profiles, organism=O.TOXOPLASMA, source_target='pfam_id', universe=universe,
                             protected_groups=groups, split=split)
    scope = {'preparation_id': PREPARATION_ID, 'organism': O.TOXOPLASMA, 'source_target': 'pfam_id',
             'target': contract.target, 'seed': SEED, 'group_fractions': list(FRACTIONS),
             'split_identity': split.identity, 'target_identity': contract.identity,
             'cohort_identity': contract.cohort_identity, 'source_census_path': str(source_path),
             'source_census_sha256': dependencies[str(source_path)], 'installed_source_sha256': dependencies[str(node_path)],
             'source_census_checksums_sha256': dependencies[str(CENSUS / 'SHA256SUMS.txt')],
             'group_source_column': 'orthogroup', 'group_source_release': 'unresolved',
             'group_policy': 'Context.groups equivalent: original nonempty orthogroup; missing values use unique __position singleton',
             'source_assignment_release': contract.source.source_release,
             'source_evidence_grade': contract.source.evidence_grade,
             'nomenclature_releases': list(contract.source.nomenclature_releases),
             'feature_access': split.feature_access, 'purpose': 'source preparation and training-profile support only',
             'capacity_used_for_selection': False, 'fitting_performed': False, 'biological_admission': False,
             'raw_source_cells_retained_in': str(source_path),
             'target_identity_semantics': 'Normalized complete profiles, diagnostic IDs, source lineage, protected groups and split; excludes raw source cells',
             'cohort_identity_semantics': 'Exact whole-gene universe/order and eligibility only',
             'gap': 'Missing orthogroups remain unresolved singleton protection; homology/source independence and biological function accuracy are not established'}
    _write(output / 'preparation_scope.json', scope)
    _write(output / 'normalized_target.json', asdict(contract))
    _write(output / 'split.json', asdict(split))
    # Scope, seed and all assignments are saved BEFORE held-out truth-dependent
    # capacity is read. No subsequent filtering or setting selection occurs.
    role_by_gene = {assignment.entity: assignment.role for assignment in split.assignments}
    records = [{**asdict(row), 'eligible': row.eligible, 'split_role': role_by_gene.get(row.gene_id)} for row in contract.rows]
    frame = pd.DataFrame(records)
    frame.to_parquet(output / 'target_records.parquet', index=False)
    pd.DataFrame({'gene_id': universe, 'source_orthogroup': original_groups,
                  'protected_group': groups, 'homology_resolved': resolved}).to_parquet(output / 'protected_groups.parquet', index=False)
    capacities = {role: contract.capacity(role) for role in ROLES}
    capacity_ids = {role: _identity(payload) for role, payload in capacities.items()}
    _write(output / 'capacity.json', capacities)

    # Replay saved immutable typed contract and every normalized full-universe
    # record. Unknown rows cannot enter the split or become target classes.
    saved = json.loads((output / 'normalized_target.json').read_text())
    saved_split = saved.pop('split')
    saved_split['assignments'] = tuple(Assignment(**row) for row in saved_split['assignments'])
    saved_source = saved.pop('source')
    for key in ('source_ids', 'nomenclature_releases'):
        saved_source[key] = tuple(saved_source[key])
    saved_rows = saved.pop('rows')
    for row in saved_rows:
        for key in ('recorded_identifiers', 'malformed_tokens'):
            row[key] = tuple(row[key])
    replay = ProfileTargetContract(**saved, rows=tuple(ProfileRow(**row) for row in saved_rows),
                                   source=ProfileSource(**saved_source), split=SplitManifest(**saved_split))
    assert replay == contract and replay.identity == contract.identity
    pd.testing.assert_frame_equal(pd.read_parquet(output / 'target_records.parquet'), frame, check_exact=True)
    assert len(records) == 8140 and sum(row['eligible'] for row in records) == 4310
    assert all(row['profile'] is None and row['split_role'] is None for row in records if not row['eligible'])
    role_groups = {}
    for row in records:
        if row['eligible']:
            assert role_groups.setdefault(row['protected_group'], row['split_role']) == row['split_role']
    saved_capacities = json.loads((output / 'capacity.json').read_text())
    train_counts = Counter(row['profile'] for row in records if row['split_role'] == 'train')
    for role in ROLES:
        counts = Counter(row['profile'] for row in records if row['split_role'] == role)
        unsupported = {profile: count for profile, count in sorted(counts.items()) if profile not in train_counts}
        payload = saved_capacities[role]
        assert payload == replay.capacity(role)
        assert payload['role_profile_counts'] == dict(counts)
        assert payload['training_profile_counts'] == dict(train_counts)
        assert payload['unsupported_profile_counts'] == unsupported
        assert payload['unsupported_genes'] == sum(unsupported.values())
        assert payload['role_genes'] == sum(counts.values())
        assert payload['supported_genes'] + payload['unsupported_genes'] == payload['role_genes']
        assert _identity(payload) == capacity_ids[role]
    assert list(contract.training_labels().index) == list(split.entities('train'))
    for path, digest in dependencies.items():
        assert _sha(path) == digest, 'Dependency changed during preparation: ' + path
    _write(output / 'input_manifest.json', dependencies)
    summary = {'preparation_id': PREPARATION_ID, 'organism': O.TOXOPLASMA, 'source_target': 'pfam_id',
               'whole_universe': len(records), 'eligible_complete_profiles': sum(eligible),
               'unknown_states': dict(Counter(row['status'] for row in records if not row['eligible'])),
               'unique_complete_profiles': len({row['profile'] for row in records if row['eligible']}),
               'recorded_memberships': sum(len(row['recorded_identifiers']) for row in records if row['eligible']),
               'unresolved_homology_whole_universe': sum(not known for known in resolved),
               'unresolved_homology_eligible': sum(keep and not known for keep, known in zip(eligible, resolved)),
               'roles': {role: {key: capacities[role][key] for key in (
                   'role_genes', 'role_protected_groups', 'supported_genes', 'unsupported_genes',
                   'training_profile_support_fraction')} for role in ROLES},
               'identities': {'raw_installed_source': dependencies[str(node_path)],
                              'raw_cell_census': dependencies[str(source_path)],
                              'target': contract.identity, 'cohort': contract.cohort_identity,
                              'split': split.identity, 'capacity_payloads': capacity_ids,
                              'preparation_scope': _identity(scope)},
               'exact_source_cell_replay': True, 'exact_full_universe_record_replay': True,
               'exact_capacity_replay': True, 'protected_group_boundaries_verified': True,
               'unsupported_profiles_retained': True, 'fitting_performed': False,
               'model_accuracy': None, 'biological_accuracy': None, 'biological_admission': False,
               'source_assignment_release': 'unresolved', 'individual_gene_evidence_grade': 'unresolved',
               'negative_semantics': contract.source.negative_semantics,
               'capacity_interpretation': 'Observed training-profile support, not model inference capacity or accuracy',
               'software_versions': {'python': sys.version.split()[0], 'numpy': np.__version__, 'pandas': pd.__version__}}
    # Convert numpy summation scalars without changing literal numeric counts.
    for key in ('eligible_complete_profiles', 'unresolved_homology_eligible'):
        summary[key] = int(summary[key])
    _write(output / 'summary.json', summary)
    (output / 'README.md').write_text(
        '# Frozen complete-profile source preparation\n\n'
        'FN-DOM-PREP-01 freezes the verified Tg Pfam recorded profiles, all original 8,140 genes, '
        'unknown annotation states, original groups and one protected eligible split. Scope, seed and '
        'assignments were saved before reporting held-out training-profile support. No classifier, '
        'model fitting, source download, biological function inference, calibration or admission occurs.\n\n'
        'The prior immutable census retains original source cells. This target contract retains normalized '
        'complete recorded profiles; raw source/census byte hashes, normalized target, cohort, split and '
        'role-specific capacity payload hashes have distinct meanings in summary.json. All unknown genes '
        'remain outside eligible roles. Every unsupported test profile remains in the test denominator. '
        'Training support is not measured model inference capacity or function accuracy and cannot select settings.\n\n'
        'Original assignment releases and evidence grades remain unresolved, separately from current '
        'nomenclature. Nonmembership means recorded-profile omission only, not adjudicated biological absence. '
        'Missing orthogroups use explicit unique singleton groups, leaving homology independence unresolved. '
        'Column-projected installed sources, every saved normalized record, capacities and hashes were exactly replayed. '
        'The executed notebook records the actual run; input_manifest.json pins dependencies/code, and '
        'SHA256SUMS.txt covers all output bytes. These artifacts are preparation evidence, not benchmark outcomes.\n')
    return summary


def main():
    """Execute the preparation once into a new immutable notebook-backed directory."""
    from scripts.notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    nb = ExecutedNotebook('FN-DOM-PREP-01 Tg Pfam complete-profile source preparation')
    nb.ns.update(output=args.out, prepare=prepare)
    nb.md('Freeze scope and protected eligible split before reporting training-profile support.',
          'Complete recorded profile preparation only; unknown absence remains unknown. No fitting or biological admission.')
    try:
        nb.code('summary = prepare(output)', 'summary')
        print(json.dumps(nb.ns['summary'], indent=2, allow_nan=False))
    except Exception as exc:
        _write(args.out / 'failure.json', {'error': type(exc).__name__, 'detail': str(exc)})
        raise
    finally:
        nb.write(str(args.out / 'executed.ipynb'))
        (args.out / 'SHA256SUMS.txt').write_text('\n'.join(
            f'{_sha(path)}  {path.relative_to(args.out)}' for path in sorted(args.out.rglob('*'))
            if path.is_file() and path.name != 'SHA256SUMS.txt') + '\n')


if __name__ == '__main__':
    main()
