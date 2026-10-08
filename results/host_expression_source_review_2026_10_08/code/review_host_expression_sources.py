"""Review two existing host expression sources and the public GTEx v11 summary.

GT-HOST-TX-01 checks exact stored-input receipts, primary release/assay metadata
and bounded processed-summary downloads. Individual access failures do not stop
other sources. Original inputs remain external; no runtime replacement is made.
"""
from __future__ import annotations

import argparse
import base64
from dataclasses import asdict
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import sys
import urllib.parse
import urllib.request

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'scripts'))
from starplast import datasets as D  # noqa: E402
from starplast.provenance import SourceFile  # noqa: E402
from recover_source_files import download  # noqa: E402

# Exact filename publicly listed by the official GTEx downloads catalogue.
V11_SUMMARY = 'GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz'


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def _metadata(url, target):
    request = urllib.request.Request(url, headers={'User-Agent': 'starplast-source-review/1'})
    with urllib.request.urlopen(request, timeout=25) as response:
        data = response.read((8 << 20)+1)
    if len(data) > 8 << 20:
        raise ValueError('Primary metadata exceeds bounded review size')
    target.write_bytes(data)
    return data


def _gct_profile(path):
    with gzip.open(path, 'rt') as stream:
        version = stream.readline().strip()
        shape = tuple(map(int, stream.readline().strip().split('\t')))
        frame = pd.read_csv(stream, sep='\t')
    if version != '#1.2' or shape != (len(frame), len(frame.columns)-2) or 'Name' not in frame:
        raise ValueError('GCT header/dimensions do not reconcile')
    return {'format': version, 'genes': len(frame), 'tissue_columns': frame.columns[2:].tolist(),
        'duplicate_versioned_gene_rows': int(frame.Name.duplicated(keep=False).sum()),
        'quantity': 'gene-level RNA-seq TPM; median summary; not protein abundance',
        'selected_fields_present': {name: name in frame for name in ('Cells_Cultured_fibroblasts', 'Liver_Hepatocyte')}}


def review(root, output):
    """Verify relocated summaries and retrieve one exact public GTEx version candidate."""
    root, output = Path(root).resolve(), Path(output)
    if not root.is_dir() or output.exists():
        raise ValueError('Use the existing archive and a new host-review snapshot')
    output.mkdir(parents=True)
    binding_path = ROOT/'results/source_recovery_final_2026_10_07_v2/bindings.json'
    inputs = [Path(__file__), binding_path, ROOT/'starplast/datasets.py', ROOT/'scripts/recover_source_files.py']
    records = []
    for key in ('host_gtex_transcriptome', 'host_mouse_tissue_transcriptome'):
        source = D.get(key)
        row = {'source_id': key, 'pmid': source.pmid, 'benchmark_admitted': False,
            'runtime_replacement': False, 'primary_requests': [], 'gaps': []}
        try:
            if key == 'host_gtex_transcriptome':
                binding = json.loads(binding_path.read_text())[key]
                path = Path(binding['path'])
                if _sha(path) != binding['sha256']:
                    raise ValueError('Frozen source-binding hash differs')
                row.update(file=binding, profile=_gct_profile(path))
            else:
                receipt_path = root/'spaces/_acquisition/parts/Mm/expression__FANTOM5_E-MTAB-3579.json'
                receipt = json.loads(receipt_path.read_text())
                if receipt['space'] != source.organism or receipt['identifiers']['arrayexpress'] != source.accession or str(receipt['identifiers']['paper_pmid']) != source.pmid:
                    raise ValueError('FANTOM accession/organism/publication receipt differs')
                inputs.append(receipt_path)
                files = []
                for item in receipt['files']:
                    path = (root/item['path']).resolve()
                    if not path.is_relative_to(root) or _sha(path) != item['sha256'] or path.stat().st_size != item['bytes']:
                        raise ValueError('FANTOM original acquisition receipt/hash differs')
                    files.append(dict(item, absolute_path=str(path)))
                    inputs.append(path)
                frame = pd.read_csv(next(f['absolute_path'] for f in files if f['name'] == 'E-MTAB-3579-mouse-tissue-tpm.tsv'), sep='\t', comment='#')
                row.update(receipt=receipt, verified_files=files, profile={'rows': len(frame),
                    'columns': frame.columns.tolist(), 'assay': 'FANTOM5 CAGE; Atlas gene-level export, not RNA-seq TPM',
                    'scope': 'four adult brain regions and juvenile biceps femoris kept separate'})
                row['gaps'].append('Historical local table filename differs from primary resource tpmss.tsv; verify transformation/association, never rename to imply exact original filename')
        except Exception as error:
            row['gaps'].append('Existing source: '+type(error).__name__+': '+str(error))
        if key == 'host_gtex_transcriptome':
            requests = [('official_download_catalogue', 'https://www.gtexportal.org/home/downloads/adult-gtex#qtl', 'gtex_downloads.html'),
                ('v11_public_object_listing', 'https://storage.googleapis.com/storage/v1/b/adult-gtex/o?'+urllib.parse.urlencode({'prefix':'bulk-gex/v11/rna-seq/','maxResults':1000}), 'gtex_v11_objects.json')]
        else:
            requests = [('primary_resource_index', 'https://www.ebi.ac.uk/gxa/json/experiments/E-MTAB-3579/resources/DATA', 'fantom_resources.json'),
                ('primary_deposit', 'https://www.ebi.ac.uk/biostudies/api/v1/studies/E-MTAB-3579', 'fantom_deposit.json')]
        requests.append(('registered_publication_metadata', 'https://www.ebi.ac.uk/europepmc/webservices/rest/search?'+urllib.parse.urlencode({'query':'EXT_ID:'+source.pmid+' AND SRC:MED', 'format':'json','resultType':'core'}), key+'_publication.json'))
        for kind, url, filename in requests:
            request = {'kind': kind, 'url': url}
            try:
                data = _metadata(url, output/filename)
                request.update(status='retrieved_metadata_not_biology', sha256=hashlib.sha256(data).hexdigest(), bytes=len(data))
                if kind == 'official_download_catalogue' and V11_SUMMARY.encode() not in data:
                    request['display_gap'] = 'JavaScript catalogue does not expose the listed filename in static HTML; primary public object metadata reviewed separately'
                if kind == 'registered_publication_metadata':
                    matches = json.loads(data)['resultList']['result']
                    if len(matches) != 1 or str(matches[0]['id']) != source.pmid:
                        raise ValueError('Publication response identity is unresolved')
                    row['publication_title'] = matches[0]['title']
                    row['publication_first_date'] = matches[0].get('firstPublicationDate')
                    row['publication_scope_gap'] = 'Resource-family citation alone does not identify this release/deposit/sample lineage'
            except Exception as error:
                request.update(status='unavailable_continue', error=type(error).__name__+': '+str(error))
            row['primary_requests'].append(request)
        if key == 'host_gtex_transcriptome':
            try:
                listing = json.loads((output/'gtex_v11_objects.json').read_text())
                if listing.get('nextPageToken'):
                    raise ValueError('Public version listing truncated; no filename association chosen')
                matches = [item for item in listing.get('items', []) if Path(item['name']).name == V11_SUMMARY]
                if len(matches) != 1 or matches[0].get('bucket') != 'adult-gtex' or not matches[0]['name'].startswith('bulk-gex/v11/rna-seq/'):
                    raise ValueError('Exact official summary filename has no unique primary object association')
                item = matches[0]
                if int(item['size']) > 32 << 20 or not item.get('md5Hash'):
                    raise ValueError('Summary exceeds bounded download or lacks provider checksum')
                url = 'https://storage.googleapis.com/adult-gtex/'+urllib.parse.quote(item['name'], safe='/')
                target = root/'spaces/Hs/expression/GTEx_v11_review'/V11_SUMMARY
                if not target.exists():
                    download(url, target, limit=32 << 20)
                if base64.b64encode(hashlib.md5(target.read_bytes()).digest()).decode() != item['md5Hash'] or target.stat().st_size != int(item['size']):
                    raise ValueError('Candidate summary differs from public object checksum/size')
                row['candidate_v11'] = {'file': asdict(SourceFile.inspect(target, 'processed_input', url,
                    'exact_official_catalogue_filename_primary_public_object_MD5_size_verified; eligibility_pending')),
                    'object_metadata': item, 'profile': _gct_profile(target),
                    'decision': 'candidate only; shared-cohort mapping, quantities/context and affected checks pending'}
                inputs.append(target)
            except Exception as error:
                row['gaps'].append('v11 candidate: '+type(error).__name__+': '+str(error))
        records.append(row)
        print(key, 'candidate', 'candidate_v11' in row, 'gaps', len(row['gaps']), flush=True)
    (output/'sources.json').write_text(json.dumps(records, indent=2, allow_nan=False)+'\n')
    summary = {'existing_source_addresses': len(records), 'verified_v11_candidates': sum('candidate_v11' in r for r in records),
        'runtime_changes': 'none', 'biological_benchmarks_admitted': 0}
    (output/'manifest.json').write_text(json.dumps({'created_utc': datetime.now(timezone.utc).isoformat(),
        'summary': summary, 'input_sha256': {str(p.resolve()): _sha(p) for p in inputs},
        'outputs': {p.name: _sha(p) for p in output.iterdir()}}, indent=2)+'\n')
    return summary


def main():
    """Execute public source/summary review and retain every access gap in a notebook."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    nb = ExecutedNotebook('Host expression source review: GTEx v10/v11 and FANTOM5')
    nb.md('GT-HOST-TX-01: two existing registered host expression sources, plus the exact GTEx v11 median-summary filename listed by the official catalogue. Verify original stored acquisition/binding receipts. Query public primary metadata; discover the exact candidate object through the public storage listing and check provider MD5/size. Retain every failure/JavaScript display gap while the next source continues. Download only this bounded summary, not per-sample or instrument raw data. No source replacement or new biological admission.')
    nb.code('from scripts.review_host_expression_sources import review',
        f'summary=review({str(args.root)!r},{str(args.out)!r})', 'summary')
    nb.md('Median gene-level RNA-seq TPM and Atlas CAGE exports are distinct measurements. Release date is not a new publication/citation lineage. GTEx cultured fibroblasts do not establish the exact infection-study donor/culture match; LCM hepatocyte enrichment requires sample/method review. FANTOM adult brain regions and juvenile muscle remain separate contexts. Do not promote a newer release before mapping/common-universe/assay checks and affected validation. Unavailable content does not stop independent work.')
    nb.write(str(args.out/'review.ipynb'))
    print(nb.ns['summary'])


if __name__ == '__main__':
    main()
