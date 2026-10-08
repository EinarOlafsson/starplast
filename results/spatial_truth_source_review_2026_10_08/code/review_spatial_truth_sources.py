"""Recover exact spatial-proteomics truth companions and inspect their schemas.

Scope is two already registered primary papers. No marker/classifier output is
automatically admitted as biological truth; missing content becomes a gap while
the next independent source proceeds. Original files remain in the data archive.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import sys
import urllib.parse
import xml.etree.ElementTree as ET

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from starplast import datasets as D  # noqa: E402
from starplast.provenance import SourceFile  # noqa: E402
from recover_pmc_sources import BASE, _metadata, named_media  # noqa: E402
from recover_source_files import download  # noqa: E402

# Original names are from the linked primary article's supplementary captions.
PAPERS = (('lopit_tgon', 'PMC7670262', ('mmc4.xls',)),
          ('pf_spatial_proteome', 'PMC13369866', ('41467_2026_73664_MOESM4_ESM.xlsx', '41467_2026_73664_MOESM1_ESM.pdf')))


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _profile(path):
    if path.suffix.lower() == '.pdf':
        return {'format': 'pdf', 'review': 'figure-level visual/assay interpretation pending'}
    if path.suffix.lower() == '.csv':
        frame = pd.read_csv(path, nrows=4)
        return {'format': 'processed_csv', 'columns': frame.columns.tolist(),
            'sample_rows': frame.fillna('').astype(str).to_dict('records'), 'scope': 'schema/sample only; not an assay population estimate'}
    with pd.ExcelFile(path) as workbook:
        return {'format': 'spreadsheet', 'sheets': {name: pd.read_excel(workbook, sheet_name=name, header=None, nrows=8)
            .fillna('').astype(str).values.tolist() for name in workbook.sheet_names},
            'scope': 'eight header/sample rows per sheet; field-level grade/mapping still pending'}


def review(root, output):
    """Verify article identities and supplemental MD5s; continue after source-specific failures."""
    root, output = Path(root).resolve(), Path(output)
    if not root.is_dir() or output.exists():
        raise ValueError('Use the existing archive and a new source-review directory')
    output.mkdir(parents=True)
    bindings_path = ROOT/'results/source_recovery_final_2026_10_07_v2/bindings.json'
    bindings = json.loads(bindings_path.read_text())
    records = []
    for key, pmcid, filenames in PAPERS:
        source = D.get(key)
        row = {'source_id': key, 'organism': source.organism, 'pmid': source.pmid, 'pmcid': pmcid,
            'primary_article': 'https://pmc.ncbi.nlm.nih.gov/articles/'+pmcid+'/', 'benchmark_admitted': False,
            'scope': 'GT-SPATIAL-01: two existing spatial papers; mixed marker origins and targeted microscopy reviewed separately',
            'companions': [], 'gaps': []}
        existing = Path(bindings[key]['path'])
        try:
            if _sha(existing) != bindings[key]['sha256']:
                raise ValueError('Recorded processed binding changed')
            row['existing_file'] = bindings[key]
            row['existing_profile'] = _profile(existing)
        except Exception as error:
            row['gaps'].append('Existing source inspection: '+type(error).__name__+': '+str(error))
        try:
            listing = BASE+'?'+urllib.parse.urlencode({'list-type': '2', 'prefix': pmcid+'.', 'delimiter': '/'})
            tree = ET.fromstring(_metadata(listing, output/(pmcid+'_versions.xml')))
            if tree.findtext('{*}IsTruncated') != 'false':
                raise ValueError('Article version listing is incomplete')
            versions = []
            for element in tree.findall('{*}CommonPrefixes/{*}Prefix'):
                prefix = element.text
                if not re.fullmatch(pmcid+r'\.\d+/', prefix or ''):
                    raise ValueError('Unexpected article prefix')
                url = BASE+prefix+prefix.rstrip('/')+'.json'
                metadata = json.loads(_metadata(url, output/(prefix.rstrip('/')+'.json')))
                if str(metadata.get('pmcid')) != pmcid or str(metadata.get('pmid')) != source.pmid:
                    raise ValueError('Primary PMC metadata differs from registered paper identity')
                versions.append(metadata)
            row['published_media_names'] = sorted({Path(urllib.parse.urlparse(url).path).name
                for metadata in versions for url in metadata.get('media_urls', [])})
            for filename in filenames:
                companion = {'declared_primary_filename': filename, 'status': 'unresolved'}
                try:
                    matches = [(metadata, url) for metadata in versions for url in named_media(metadata, pmcid, source.pmid, filename)]
                    if len(matches) != 1:
                        raise ValueError('Exact filename/published-version association has '+str(len(matches))+' matches')
                    metadata, url = matches[0]
                    target = root/'spatial_truth_review_inputs'/key/filename
                    expected_md5 = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)['md5'][0].lower()
                    if target.exists():
                        if hashlib.md5(target.read_bytes()).hexdigest() != expected_md5:
                            raise ValueError('Existing companion differs from primary MD5')
                    else:
                        download(url, target)
                    if hashlib.md5(target.read_bytes()).hexdigest() != expected_md5:
                        raise ValueError('Downloaded companion differs from published MD5')
                    file = SourceFile.inspect(target, 'processed_input', url,
                        'exact_primary_PMCID_PMID_named_published_media_and_MD5; benchmark_field_grade_pending')
                    companion.update(status='retrieved_verified_primary_companion', file=asdict(file),
                        published_md5=expected_md5, license=metadata.get('license_code', 'unresolved'),
                        article_doi=metadata.get('doi'), article_version=metadata.get('version'))
                    try:
                        companion['profile'] = _profile(target)
                    except Exception as error:
                        companion['profile_gap'] = type(error).__name__+': '+str(error)
                except Exception as error:
                    companion.update(status='unavailable_or_association_refused', error=type(error).__name__+': '+str(error))
                row['companions'].append(companion)
        except Exception as error:
            row['gaps'].append('Primary companion discovery: '+type(error).__name__+': '+str(error))
        # Full text is kept externally; do not redistribute the source paper in Git.
        try:
            folder = root/'spatial_truth_review_inputs'/key
            folder.mkdir(parents=True, exist_ok=True)
            fulltext = folder/(pmcid+'_fulltext.xml')
            url = 'https://www.ebi.ac.uk/europepmc/webservices/rest/'+pmcid+'/fullTextXML'
            payload = fulltext.read_bytes() if fulltext.exists() else _metadata(url, fulltext)
            article = ET.fromstring(payload)
            identifiers = {element.get('pub-id-type'): element.text for element in article.findall('.//article-id')}
            if identifiers.get('pmid') != source.pmid:
                raise ValueError('Full-text XML publication identity differs')
            row['fulltext'] = asdict(SourceFile.inspect(fulltext, 'raw_input', url, 'primary_article_context; not assay truth'))
            row['article_identifiers'] = identifiers
            row['supplement_captions'] = [{'hrefs': [e.get('{http://www.w3.org/1999/xlink}href') for e in part.iter() if e.get('{http://www.w3.org/1999/xlink}href')],
                'label': part.findtext('label', ''), 'caption': ''.join(part.find('caption').itertext())[:300] if part.find('caption') is not None else ''}
                for part in article.findall('.//supplementary-material')]
        except Exception as error:
            row['gaps'].append('Primary full-text context: '+type(error).__name__+': '+str(error))
        row['gaps'].extend(('Marker flag alone does not establish experimental truth or feature independence',
            'Per-gene microscopy success, compartment granularity, mapping and selection bias require review'))
        records.append(row)
        print(key, 'primary companions', sum(c['status'].startswith('retrieved') for c in row['companions']), 'gaps', len(row['gaps']), flush=True)
    (output/'sources.json').write_text(json.dumps(records, indent=2, allow_nan=False)+'\n')
    summary = {'reviewed_source_addresses': len(records), 'verified_companions': sum(c['status'].startswith('retrieved') for r in records for c in r['companions']),
        'biological_benchmarks_admitted': 0, 'runtime_data_changes': 'none', 'scope': 'primary acquisition/schema preflight; field-level review pending'}
    inputs = [Path(__file__), bindings_path, ROOT/'starplast/datasets.py', ROOT/'scripts/recover_pmc_sources.py', ROOT/'scripts/recover_source_files.py']
    (output/'manifest.json').write_text(json.dumps({'created_utc': datetime.now(timezone.utc).isoformat(), 'summary': summary,
        'input_sha256': {str(p): _sha(p) for p in inputs}, 'outputs': {p.name: _sha(p) for p in output.iterdir()}}, indent=2)+'\n')
    return summary


def main():
    """Execute bounded source acquisition/schema inspection with actual notebook outputs."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    nb = ExecutedNotebook('Spatial marker and microscopy truth-source preflight')
    nb.md('GT-SPATIAL-01: two already registered primary spatial papers. Verify stored source bindings; discover public PMC article versions; retrieve only the primary named metadata/microscopy companions and retain their original filenames/checksums. Source-specific failures become gaps and the next source continues. No restricted access is bypassed, no source is renamed or promoted.')
    nb.code('from scripts.review_spatial_truth_sources import review', f'summary=review({str(args.root)!r},{str(args.out)!r})', 'summary')
    nb.md('Published marker sets combine literature, functional/sequence inference and profile-based selection; they are not automatically independent measured truth. Review microscopy separately with gene-level success, assay compartment resolution and selection scope. Schema samples are not complete population or performance estimates. No biological benchmark or runtime dataset is admitted by this preflight.')
    nb.write(str(args.out/'preflight.ipynb'))
    print(nb.ns['summary'])


if __name__ == '__main__':
    main()
