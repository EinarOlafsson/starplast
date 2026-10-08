"""Freeze complete recorded-domain profile candidates without fitting or downloads."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def census(output):
    """Pin all original source cells, validated profiles, metadata and evidence gaps."""
    import pandas as pd
    from starplast import datasets as D, functional_domain_profiles as P, functional_domains as FD, organisms as O, strategies as S

    output = Path(output)
    lookup = FD.shipped()
    metadata_path = ROOT / 'starplast/data/functional_domain_names.json'
    node_paths = {organism: Path(O.nodes_path(organism)) for organism in (O.TOXOPLASMA, O.FALCIPARUM)}
    input_hashes = {str(path): _sha(path) for path in (*node_paths.values(), metadata_path)}
    summary = {'benchmark_partition': 'FN-DOM-01', 'fitting_performed': False, 'biological_admission': False,
               'target_semantics': 'Complete original recorded identifier set; not ordered domain architecture',
               'negative_semantics': P.NEGATIVE_SEMANTICS,
               'source_assignment_release': 'unresolved', 'source_gene_evidence_grade': 'unresolved',
               'current_nomenclature': lookup.sources, 'labels': []}
    for organism, node_path in node_paths.items():
        ctx = S.Context.shipped(organism, graph={})
        original_nodes = ctx.nodes.copy(deep=True)
        for target in P.TARGET_PATTERNS:
            if target not in ctx.nodes:
                continue
            owners = [dataset for dataset in D.REGISTRY if dataset.organism == organism and target in dataset.columns]
            profiles, terms = P.recorded_profiles(ctx, target, domain_lookup=lookup,
                                                  source_ids=tuple(dataset.key for dataset in owners))
            assert profiles.gene_id.tolist() == ctx.gene_ids.tolist()
            # Independent whole-token replay validates each original cell and
            # complete profile. It never uses the parser's extracted IDs as truth.
            expected_members = set()
            pattern = r'IPR[0-9]{6}' if target.startswith('interpro') else r'PF[0-9]{5}(?:\.[0-9]+)?'
            truth = ctx.truth(target)
            for position, row in enumerate(profiles.itertuples(index=False)):
                value = str(truth.iloc[position]) if pd.notna(truth.iloc[position]) else None
                tokens = [] if value is None else [part.strip() for part in value.split(';')]
                nonempty = [token for token in tokens if token]
                valid = sorted({token for token in nonempty if re.fullmatch(pattern, token)})
                malformed = [token for token in nonempty if not re.fullmatch(pattern, token)]
                eligible = bool(valid) and not malformed
                expected = json.dumps(valid, separators=(',', ':')) if eligible else None
                assert row.profile == expected and row.eligible == eligible
                assert row.recorded_identifiers == valid and row.malformed_tokens == malformed
                raw = ctx.nodes[target].iloc[position]
                assert row.source_value == (str(raw) if pd.notna(raw) else None)
                expected_members.update((row.gene_id, term, eligible) for term in valid)
                if eligible:
                    assert P.profile_members(row.profile, target) == set(valid)
            assert set(zip(terms.gene_id, terms.identifier, terms.profile_eligible)) == expected_members
            assert not terms.duplicated(['organism', 'gene_id', 'source_target', 'identifier']).any()
            profiles.to_parquet(output / f'{organism}_{target}_profiles.parquet', index=False)
            terms.to_parquet(output / f'{organism}_{target}_terms.parquet', index=False)
            eligible_profiles = profiles[profiles.eligible]
            profile_sizes = eligible_profiles.profile.value_counts()
            assigned_metadata = terms[terms.profile_eligible]
            eligible_groups = ctx.groups()[profiles.eligible.to_numpy(dtype=bool)]
            summary['labels'].append({
                'organism': organism, 'source_target': target, 'genes': ctx.n,
                'statuses': profiles.status.value_counts().to_dict(), 'eligible_profiles': int(profiles.eligible.sum()),
                'recorded_term_memberships': len(assigned_metadata), 'diagnostic_valid_terms_in_malformed_profiles': int((~terms.profile_eligible).sum()),
                'unique_profiles': len(profile_sizes), 'singleton_profiles': int(profile_sizes.eq(1).sum()),
                'multi_identifier_profiles': int(eligible_profiles.recorded_identifiers.map(len).gt(1).sum()),
                'empty_delimiter_genes': int(profiles.ignored_empty_tokens.gt(0).sum()),
                'ignored_empty_tokens': int(profiles.ignored_empty_tokens.sum()),
                'duplicate_identifier_tokens': int(profiles.duplicate_identifier_tokens.sum()),
                'metadata_status_memberships': assigned_metadata.metadata_status.value_counts().to_dict(),
                'unique_recorded_ids': int(assigned_metadata.identifier.nunique()),
                'eligible_homology_groups': len(set(eligible_groups)),
                'eligible_genes_without_resolved_homology': sum(str(group).startswith('__') for group in eligible_groups),
                'source_ids': [dataset.key for dataset in owners],
                'source_records': [{'key': dataset.key, 'url': dataset.url, 'path': dataset.path, 'note': dataset.note} for dataset in owners],
                'source_assignment_release': 'unresolved', 'individual_gene_evidence_grade': 'unresolved',
                'top_recorded_profiles': profile_sizes.head(20).to_dict(),
                'biological_admission': False, 'exact_source_and_profile_replay': True})
        pd.testing.assert_frame_equal(ctx.nodes, original_nodes, check_exact=True)
        assert _sha(node_path) == input_hashes[str(node_path)]
    assert _sha(metadata_path) == input_hashes[str(metadata_path)]
    code_paths = [Path(__file__), ROOT / 'scripts/notebook_runner.py', *[
        ROOT / 'starplast' / name for name in ('functional_domain_profiles.py', 'functional_domains.py',
                                              'strategies.py', 'embedding.py', 'provenance.py', 'organisms.py', 'datasets.py')]]
    (output / 'code').mkdir()
    for path in code_paths:
        input_hashes[str(path)] = _sha(path)
        shutil.copyfile(path, output / 'code' / path.name)
    (output / 'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False) + '\n')
    (output / 'input_manifest.json').write_text(json.dumps(input_hashes, indent=2) + '\n')
    (output / 'README.md').write_text(
        '# Complete recorded-domain profile census\n\n'
        'FN-DOM-01 source preparation only: both installed organisms and all four InterPro/Pfam source fields. '
        'No fitting, gene function inference, biological admission, source downloads or GO mapping. '
        'Each original gene and source cell is retained; strict whole-token validation rejects an entire malformed '
        'profile. Empty delimiters are ignored with explicit counts. Valid overlapping identifiers and Pfam '
        'accession versions remain in a sorted unique recorded set. This is not ordered domain architecture.\n\n'
        'Current nomenclature availability cannot alter a valid recorded profile: missing-current or retired IDs '
        'remain original source membership with separate status. No forwarding memberships are applied. '
        'Missing or malformed source profiles are unknown; omission from a recorded profile is not verified '
        'biological absence. Gene-assignment source releases and individual evidence grades remain unresolved, '
        'separately from pinned InterPro/Pfam nomenclature.\n\n'
        'All source cells, complete profiles and diagnostic terms were independently replayed. Installed node '
        'tables and packaged metadata retain their input hashes. input_manifest.json pins source/code bytes, '
        'the executed notebook contains the actual census, and summary.json retains eligibility, singleton '
        'profiles, unresolved homology counts and exact scope. This census does not establish benchmark accuracy.\n')
    return summary


def main():
    """Execute the source census once into a new immutable result directory."""
    from scripts.notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    nb = ExecutedNotebook('FN-DOM-01 complete recorded-domain source profiles')
    nb.ns.update(output=args.out, census=census)
    nb.md('Freeze source-recovery candidates while retaining every recorded domain identifier and unknown absence.',
          'This census reads installed sources and pinned name metadata only; it does not download, fit or admit biological activity truth.')
    try:
        nb.code('summary = census(output)', 'summary')
    except Exception as exc:
        (args.out / 'failure.json').write_text(json.dumps({'error': type(exc).__name__, 'detail': str(exc)}, indent=2) + '\n')
        raise
    finally:
        nb.write(str(args.out / 'executed.ipynb'))
        (args.out / 'SHA256SUMS.txt').write_text('\n'.join(
            f'{_sha(path)}  {path.relative_to(args.out)}' for path in sorted(args.out.rglob('*'))
            if path.is_file() and path.name != 'SHA256SUMS.txt') + '\n')


if __name__ == '__main__':
    main()
