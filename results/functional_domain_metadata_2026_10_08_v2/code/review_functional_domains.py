"""Pin public InterPro/Pfam metadata and audit installed annotation names offline."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import shutil
import sys
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

DOWNLOADS = {
    'interpro_entry.list': 'https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list',
    'interpro_names.dat': 'https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/names.dat',
    'interpro_release_notes.txt': 'https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/release_notes.txt',
    'Pfam-A.hmm.dat.gz': 'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz',
    'Pfam.version.gz': 'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz',
    'pfam_relnotes.txt': 'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt',
    'interpro_license.html': 'https://interpro-documentation.readthedocs.io/en/latest/license.html',
    'pfam_about.html': 'https://pfam-docs.readthedocs.io/en/latest/pfam.html',
}


def fetch(archive, downloads=None):
    """Acquire bounded metadata only, retaining HTTP, date and hash receipts."""
    archive = Path(archive)
    archive.mkdir(parents=True, exist_ok=False)
    receipts = []
    for name, url in (DOWNLOADS if downloads is None else downloads).items():
        receipt = {'name': name, 'requested_url': url,
                   'retrieved_utc': datetime.now(timezone.utc).isoformat()}
        try:
            request = urllib.request.Request(url, headers={'User-Agent': 'Starplast-domain-metadata/0.54'})
            with urllib.request.urlopen(request, timeout=30) as response:
                payload = response.read(8_000_001)
                if len(payload) > 8_000_000:
                    raise ValueError('Metadata exceeds 8 MB acquisition limit')
                if not name.endswith('.html') and b'<html' in payload[:2000].lower():
                    raise ValueError('Ontology metadata returned HTML')
                path = archive / name
                path.write_bytes(payload)
                receipt.update(status='available', final_url=response.url, http_status=response.status,
                               headers=dict(response.headers), bytes=len(payload), path=str(path.resolve()),
                               sha256=hashlib.sha256(payload).hexdigest())
        except Exception as exc:
            receipt.update(status='unavailable', error=type(exc).__name__ + ': ' + str(exc))
        receipts.append(receipt)
        (archive / 'receipts.json').write_text(json.dumps(receipts, indent=2) + '\n')
    (archive / 'URLS.txt').write_text('\n'.join(f'{r["name"]}\t{r["requested_url"]}' for r in receipts) + '\n')
    (archive / 'SHA256SUMS.txt').write_text('\n'.join(
        f'{r["sha256"]}  {r["name"]}' for r in receipts if r['status'] == 'available') + '\n')
    return receipts


def main():
    """Execute acquisition and analysis in a new immutable evidence notebook."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--acquire-only', action='store_true')
    parser.add_argument('--reuse-archive', action='store_true')
    parser.add_argument('--status-archive', type=Path)
    args = parser.parse_args()
    if args.out.exists():
        raise ValueError('Use a new immutable result directory')
    args.out.mkdir(parents=True)
    from scripts.notebook_runner import ExecutedNotebook
    nb = ExecutedNotebook('Pinned functional domain nomenclature and exact membership audit')
    nb.ns.update(archive=args.archive, output=args.out, fetch=fetch,
                 status_archive=args.status_archive, verify_receipts=verify_receipts)
    nb.md('Metadata only: InterPro current entries/names, Pfam descriptors/version, and official license evidence.',
          'Source nomenclature does not add gene assignments or establish biological inference accuracy.')
    try:
        if args.reuse_archive:
            nb.code('receipts = verify_receipts(archive)', 'receipts')
        else:
            nb.code('receipts = fetch(archive)', 'receipts')
        if args.status_archive:
            if args.status_archive.exists():
                nb.code('status_receipts = verify_receipts(status_archive)', 'status_receipts')
            else:
                nb.code("status_receipts = fetch(status_archive, {'Pfam-A.dead.gz': "
                        "'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.dead.gz'})",
                        'status_receipts')
        if not args.acquire_only:
            nb.ns.update(review=review)
            nb.code('summary = review(output, archive, receipts, status_archive)', 'summary')
    finally:
        nb.write(str(args.out / 'executed.ipynb'))


def verify_receipts(archive):
    """Verify exact bytes before reusing a completed public acquisition."""
    archive = Path(archive)
    receipts = json.loads((archive / 'receipts.json').read_text())
    for record in receipts:
        if record['status'] == 'available':
            path = archive / record['name']
            if hashlib.sha256(path.read_bytes()).hexdigest() != record['sha256']:
                raise ValueError('Pinned metadata source hash changed: ' + record['name'])
    return receipts


def review(output, archive, receipts, status_archive=None):
    """Build and independently replay an installed-ID metadata subset."""
    from starplast import functional_domains as F, organisms as O, strategies as S
    import pandas as pd
    import re
    archive, output = Path(archive), Path(output)
    entries = F.parse_interpro((archive / 'interpro_entry.list').read_text())
    names = F.parse_interpro_names((archive / 'interpro_names.dat').read_text())
    if names != {key: entry.name for key, entry in entries.items()}:
        raise ValueError('InterPro entry.list and names.dat disagree')
    release_notes = (archive / 'interpro_release_notes.txt').read_text()
    version_text = gzip.decompress((archive / 'Pfam.version.gz').read_bytes()).decode()
    interpro_version = re.search(r'^Release (\S+), (.+)$', release_notes, re.M).groups()
    pfam_version = re.search(r'^Pfam release\s*:\s*(\S+)', version_text, re.M).group(1)
    declared_ipr_count = int(re.search(r'contains (\d+) entries', release_notes).group(1))
    declared_pfam_count = int(re.search(r'^Pfam-A families\s*:\s*(\d+)', version_text, re.M).group(1))
    pfam = F.parse_pfam(gzip.decompress((archive / 'Pfam-A.hmm.dat.gz').read_bytes()).decode())
    assert len(entries) == declared_ipr_count
    assert len(pfam) == declared_pfam_count
    license_text = (archive / 'interpro_license.html').read_text()
    assert 'CC0 1.0' in license_text and 'InterPro, Pfam, PRINTS and SFLD downloadable data' in license_text
    sources = {
        'InterPro': {'release': interpro_version[0], 'release_date': interpro_version[1],
                    'license': 'CC0-1.0', 'license_url': DOWNLOADS['interpro_license.html'],
                    'metadata_url': DOWNLOADS['interpro_entry.list'], 'entries': len(entries)},
        'Pfam': {'release': pfam_version, 'version_metadata': version_text,
                 'license': 'CC0-1.0', 'license_url': DOWNLOADS['interpro_license.html'],
                 'metadata_url': DOWNLOADS['Pfam-A.hmm.dat.gz'], 'entries': len(pfam),
                 'version_gap': 'Pfam.version states UniProtKB 2025_03; online about page states 2026_01. '
                                'Descriptor and version bytes pinned; original gene assignment release unresolved.'}}
    members, hashes = [], {}
    for organism in (O.TOXOPLASMA, O.FALCIPARUM):
        path = Path(O.nodes_path(organism))
        hashes[organism] = hashlib.sha256(path.read_bytes()).hexdigest()
        context = S.Context.shipped(organism, graph={})
        for field in ('interpro_id', 'interpro_ids', 'pfam_id', 'pfam_ids'):
            if field in context.nodes:
                members.append(source_members(context, field))
    membership = pd.concat(members, ignore_index=True)
    # Missing IDs remain unresolved unless the separate authoritative retired
    # descriptor confirms withdrawal. Replacement mentions never change members.
    retired_path = Path(status_archive) / 'Pfam-A.dead.gz' if status_archive else None
    retired_text = gzip.decompress(retired_path.read_bytes()).decode() if retired_path and retired_path.exists() else ''
    retired = F.parse_retired_pfam(retired_text)
    snapshot = F.make_snapshot(membership.value, entries, pfam, sources=sources, retired_pfam=retired)
    sources['Pfam']['retired_metadata_entries'] = len(retired)
    (output / 'lookup.json').write_text(json.dumps(snapshot, indent=2, ensure_ascii=False) + '\n')
    lookup = F.DomainLookup.read(output / 'lookup.json')
    # Independent parser replay uses literal source rows/blocks rather than the
    # production metadata parsers. A current name must agree with both name
    # files; known withdrawal must appear in the dedicated retirement source.
    direct_ipr = {line.split('\t')[0]: tuple(line.split('\t')[1:])
                  for line in (archive / 'interpro_entry.list').read_text().splitlines()[1:] if line}
    direct_pfam = {}
    for block in gzip.decompress((archive / 'Pfam-A.hmm.dat.gz').read_bytes()).decode().split('//'):
        fields = {line[5:7]: line[7:].strip() for line in block.splitlines()
                  if line.startswith('#=GF ') and line[5:7] in {'AC', 'DE', 'TP'}}
        if fields:
            direct_pfam[fields['AC'].split('.')[0]] = (fields['DE'], fields['TP'], fields['AC'].split('.')[-1])
    direct_retired = set(re.findall(r'^#=GF AC\s+(PF\d{5})', retired_text, re.M))
    for requested, row in snapshot['terms'].items():
        identifier = row['identifier']
        if row['status'] == 'current_metadata':
            if identifier.startswith('IPR'):
                kind, name = direct_ipr[identifier]
                assert (row['name'], row['entry_type']) == (name, kind)
            else:
                assert (row['name'], row['entry_type'], row['current_accession_version']) == direct_pfam[identifier]
        elif row['status'] == 'withdrawn':
            assert identifier in direct_retired and identifier not in direct_pfam
            assert row['name'] is None and row['entry_type'] is None
        else:
            assert identifier not in direct_ipr and identifier not in direct_pfam and identifier not in direct_retired
            assert row['name'] is None and row['entry_type'] is None
    annotated = membership.copy()
    rows = [lookup.lookup(value) for value in annotated.value]
    for field in ('name', 'entry_type', 'status', 'ontology_release', 'metadata_license', 'metadata_license_url',
                  'current_accession_version', 'accession_version_status'):
        annotated['ontology_' + field if not field.startswith('ontology_') else field] = [r.get(field) for r in rows]
    # This table adds metadata columns only. The entire original membership table
    # and every source description must reproduce bit-for-bit/element-for-element.
    pd.testing.assert_frame_equal(annotated[membership.columns], membership, check_exact=True)
    for organism, digest in hashes.items():
        assert hashlib.sha256(Path(O.nodes_path(organism)).read_bytes()).hexdigest() == digest
    membership.to_parquet(output / 'source_memberships.parquet', index=False)
    annotated.to_parquet(output / 'named_memberships.parquet', index=False)
    summary = {'source_nodes_sha256': hashes, 'source_releases': sources, 'retired_metadata_acquired': bool(retired_text),
               'independent_metadata_replay': True, 'exact_source_membership_replay': True,
               'terms': len(snapshot['terms']), 'membership_rows': len(membership),
               'status_counts': pd.Series([r['status'] for r in snapshot['terms'].values()]).value_counts().to_dict(),
               'labels': []}
    for (organism, target), group in annotated.groupby(['organism', 'target']):
        summary['labels'].append({'organism': organism, 'target': target, 'membership_rows': len(group),
                                 'genes': group.gene_id.nunique(), 'terms': group.value.nunique(),
                                 'named_rows': int(group.ontology_name.notna().sum()),
                                 'status_counts': group.ontology_status.value_counts().to_dict()})
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    shutil.copyfile(archive / 'receipts.json', output / 'source_receipts.json')
    if retired_path and retired_path.exists():
        shutil.copyfile(Path(status_archive) / 'receipts.json', output / 'status_receipts.json')
    (output / 'code').mkdir()
    for source in (Path(__file__), ROOT / 'scripts/notebook_runner.py', ROOT / 'starplast/functional_domains.py'):
        shutil.copyfile(source, output / 'code' / source.name)
    (output / 'SHA256SUMS.txt').write_text('\n'.join(
        f'{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(output)}'
        for path in sorted(output.rglob('*')) if path.is_file()) + '\n')
    return summary


def source_members(context, target):
    """Replay raw installed domain IDs/descriptions independently of UI enrichment."""
    import pandas as pd
    import re
    rows = []
    pattern = r'\bPF\d{5}(?:\.\d+)?\b' if target.startswith('pfam') else r'\bIPR\d{6}\b'
    description_field = 'pfam_desc' if target.startswith('pfam') else 'interpro_desc'
    descriptions = context.nodes.get(description_field, pd.Series('', index=context.nodes.index))
    expected = set()
    for position, text in context.truth(target).dropna().items():
        terms = list(dict.fromkeys(re.findall(pattern, str(text))))
        raw_terms = str(text).split(';')
        raw_names = str(descriptions.iloc[position]).split(';') if pd.notna(descriptions.iloc[position]) else []
        paired = dict(zip((term.strip() for term in raw_terms), raw_names)) if len(raw_terms) == len(raw_names) else {}
        gene = str(context.gene_ids[position])
        expected.update((gene, term) for term in re.findall(pattern, str(text)))
        for term in terms:
            rows.append((context.organism, target, term, paired.get(term, ''), gene))
    members = pd.DataFrame(rows, columns=['organism', 'target', 'value', 'description', 'gene_id'])
    assert set(zip(members.gene_id, members.value)) == expected
    assert not members.duplicated(['organism', 'target', 'value', 'gene_id']).any()
    return members


if __name__ == '__main__':
    main()
