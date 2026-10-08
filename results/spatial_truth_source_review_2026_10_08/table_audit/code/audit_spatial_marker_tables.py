"""Audit verified spatial marker tables and recover primary microscopy figures.

GT-SPATIAL-01 retains every source row, exact-ID mappings and class/context gaps.
Marker fields are reviewed candidates, never automatically experimental truth.
The acquisition and audit run inside an annotated executed notebook.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import urllib.parse
import xml.etree.ElementTree as ET

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'scripts'))
from starplast import organisms as O  # noqa: E402
from starplast.provenance import SourceFile  # noqa: E402
from recover_pmc_sources import named_media  # noqa: E402
from recover_source_files import download  # noqa: E402

TABLES = {
    'lopit_tgon': (O.TOXOPLASMA, 'mmc4.xls', 'S3 - Feature metadata', 'Accession',
        [('markers', 'tagm.map.allocation.pred', 'extracellular_tachyzoites')]),
    'pf_spatial_proteome': (O.FALCIPARUM, None, 'PlasmoLOPIT_data_summary', 'Gene accession',
        [('markers (S1-S2)', 'Final location (S1-S2)', 'late_schizonts_S1_S2'),
         ('markers (S1-S2-S3)', 'Final location (S1-S2-S3)', 'schizonts_plus_merozoites_S1_S2_S3')]),
}
FIGURES = {'lopit_tgon': ('fig1', 'fig2'), 'pf_spatial_proteome': ('Fig3',)}


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _verify(file):
    path = Path(file['path'])
    if _sha(path) != file['sha256'] or path.stat().st_size != file['bytes']:
        raise ValueError('Pinned input changed: '+str(path))
    return path


def _media(record, prior, root, filename):
    candidates = []
    for path in sorted(prior.glob(record['pmcid']+'.*.json')):
        metadata = json.loads(path.read_text())
        candidates.extend((metadata, url) for url in named_media(metadata, record['pmcid'], record['pmid'], filename))
    if len(candidates) != 1:
        raise ValueError('Primary exact filename/version association is not unique')
    metadata, url = candidates[0]
    target = root/'spatial_truth_review_inputs'/record['source_id']/filename
    expected = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)['md5'][0].lower()
    if not target.exists():
        download(url, target)
    if hashlib.md5(target.read_bytes()).hexdigest() != expected:
        raise ValueError('Primary published MD5 mismatch')
    if filename.endswith('.jpg') and not target.read_bytes().startswith(b'\xff\xd8\xff'):
        raise ValueError('Expected a JPEG figure')
    return {'filename': filename, 'file': asdict(SourceFile.inspect(target, 'raw_input', url,
        'primary_XML_address_PMCID_PMID_published_version_MD5_verified; assay_interpretation_pending')),
        'published_md5': expected, 'license': metadata.get('license_code', 'unresolved')}


def audit(prior, root, output):
    """Freeze complete marker rows/mapping counts and continue after individual source gaps."""
    prior, root, output = Path(prior), Path(root), Path(output)
    if output.exists() or not root.is_dir():
        raise ValueError('Use a new audit snapshot and the existing source archive')
    output.mkdir(parents=True)
    inputs = [prior/'sources.json', Path(__file__), ROOT/'starplast/organisms.py',
        ROOT/'scripts/recover_pmc_sources.py', ROOT/'scripts/recover_source_files.py']
    reports, media = [], []
    for record in json.loads((prior/'sources.json').read_text()):
        key = record['source_id']
        report = {'source_id': key, 'primary_article': record['primary_article'], 'benchmark_admitted': False,
            'negative_semantics': 'unknown marker/missing value/failed tagging or IFA are not negative truth',
            'gaps': ['Mixed marker origins and profile-based selection; row-level independence unverified',
                'Same-paper targeted microscopy is selected validation, not an unbiased proteome-wide sample',
                'Compartment taxonomy, strain/stage and protein-group ambiguity require assay review']}
        try:
            organism, filename, sheet, id_column, fields = TABLES[key]
            binding = record['existing_file'] if filename is None else next(c['file'] for c in record['companions']
                if c['declared_primary_filename'] == filename and c['status'].startswith('retrieved'))
            path = _verify(binding)
            nodes_path = Path(O.nodes_path(organism))
            nodes = pd.read_parquet(nodes_path, columns=['gene_id'])
            reference = set(nodes.gene_id.astype(str))
            frame = pd.read_excel(path, sheet_name=sheet, dtype=str)
            inputs.extend((path, nodes_path))
            report['table_rows'] = len(frame)
            ids = frame[id_column]
            report['missing_accessions'] = int(ids.isna().sum())
            report['duplicate_accession_rows'] = int(ids.duplicated(keep=False).sum())
            rows = []
            for marker_column, prediction_column, context in fields:
                marker = frame[marker_column]
                prediction = frame[prediction_column]
                # Exact source labels remain unchanged. Sentinel detection never infers a class.
                known = marker.notna() & ~marker.str.strip().str.lower().isin(['', 'unknown', 'unassigned'])
                known_prediction = prediction.notna() & ~prediction.str.strip().str.lower().isin(['', 'unknown', 'unassigned'])
                mapped = ids.isin(reference)
                sub = pd.DataFrame({'source_row': range(2, len(frame)+2), 'source_id': key,
                    'accession': ids, 'context': context, 'marker_field': marker_column,
                    'marker_label': marker, 'marker_present': known, 'prediction_field': prediction_column,
                    'source_prediction': prediction, 'prediction_present': known_prediction,
                    'exact_installed_id': mapped, 'duplicate_accession': ids.duplicated(keep=False),
                    'truth_grade': 'unresolved', 'benchmark_admitted': False})
                sub.to_parquet(output/(key+'_'+context+'_rows.parquet'), index=False)
                rows.append({'context': context, 'marker_field': marker_column, 'marker_rows': int(known.sum()),
                    'marker_classes': int(marker[known].nunique()), 'marker_exact_ID_rows': int((known & mapped).sum()),
                    'marker_unmapped_rows': int((known & ~mapped).sum()), 'all_exact_ID_rows': int(mapped.sum()),
                    'classifier_assignment_rows': int(known_prediction.sum()),
                    'marker_vs_prediction_literal_disagreements': int((known & known_prediction & marker.ne(prediction)).sum()),
                    'interpretation': 'Literal mismatch may be coarse/refined taxonomy, not biological error',
                    'marker_class_counts': {str(k): int(v) for k, v in marker[known].value_counts().items()}})
            report['fields'] = rows
            report['status'] = 'complete_source_rows_and_exact_ID_mapping_reviewed; biology_withheld'
        except Exception as error:
            report['status'] = 'table_review_unavailable_continue'
            report['gaps'].append(type(error).__name__+': '+str(error))
        # Independent acquisition still runs after a table error.
        try:
            xml_path = _verify(record['fulltext'])
            inputs.append(xml_path)
            article = ET.parse(xml_path).getroot()
            requested = []
            for figure_id in FIGURES[key]:
                figure = next(e for e in article.findall('.//fig') if e.get('id') == figure_id)
                names = [e.get('{http://www.w3.org/1999/xlink}href') for e in figure.iter()
                    if (e.get('{http://www.w3.org/1999/xlink}href') or '').endswith('.jpg')]
                if len(names) != 1:
                    raise ValueError('Figure has no unique primary JPEG address')
                requested.append((figure_id, names[0]))
            if key == 'lopit_tgon':
                # Address was published in the primary Table S9 supplement caption.
                associations = [e for e in record['supplement_captions'] if 'mmc10.xls' in e['hrefs']]
                if len(associations) != 1:
                    raise ValueError('Table S9 primary address unavailable')
                requested.append(('Table_S9', 'mmc10.xls'))
            for figure_id, filename in requested:
                row = {'source_id': key, 'figure_or_table': figure_id, 'filename': filename}
                try:
                    row.update(_media(record, prior, root, filename), status='retrieved_verified_primary_media')
                except Exception as error:
                    row.update(status='unavailable_continue', error=type(error).__name__+': '+str(error))
                media.append(row)
        except Exception as error:
            report['gaps'].append('Figure discovery: '+type(error).__name__+': '+str(error))
        reports.append(report)
        print(key, report['status'], flush=True)
    (output/'table_review.json').write_text(json.dumps(reports, indent=2, allow_nan=False)+'\n')
    (output/'media_receipts.json').write_text(json.dumps(media, indent=2, allow_nan=False)+'\n')
    inputs.extend(sorted(prior.glob('PMC*.json')))
    summary = {'source_addresses': len(reports), 'source_row_reviews': sum('fields' in r for r in reports),
        'verified_additional_media': sum(r['status'] == 'retrieved_verified_primary_media' for r in media),
        'biological_benchmarks_admitted': 0, 'installed_data_or_runtime_changes': 'none'}
    (output/'manifest.json').write_text(json.dumps({'created_utc': datetime.now(timezone.utc).isoformat(),
        'summary': summary, 'input_sha256': {str(p.resolve()): _sha(p) for p in inputs},
        'outputs': {p.name: _sha(p) for p in output.iterdir()}}, indent=2)+'\n')
    return summary


def main():
    """Execute the bounded complete-table audit with an annotated notebook."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior', type=Path, required=True)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    nb = ExecutedNotebook('Spatial marker truth: complete counts, exact-ID mapping and microscopy inputs')
    nb.md('GT-SPATIAL-01, bounded extension: the same two registered spatial papers. Read SHA-verified original tables; retain every marker/prediction row and separate S1/S2 from S1/S2/S3. Match only exact installed gene IDs; no guessed aliases. Discover microscopy figures from verified primary XML and Table S9 from its verified caption. Verify original PMC media names, published versions and MD5; individual failures become gaps while the next task proceeds.')
    nb.code('from scripts.audit_spatial_marker_tables import audit',
        f'summary=audit({str(args.prior)!r},{str(args.root)!r},{str(args.out)!r})', 'summary')
    nb.md('All labels remain unresolved truth candidates. Final marker sets combine prior knowledge, inferred/GO-derived assignments, profile selection and same-paper microscopy. Literal marker/prediction disagreements may reflect different compartment granularity. Classifier assignments are predictions. Failed tagging/detection and unassigned markers are unknown. Targeted validation cannot estimate unbiased proteome-wide performance; admission requires per-gene assay outcomes, timing/dependence and a reviewed compartment taxonomy.')
    nb.write(str(args.out/'audit.ipynb'))
    print(nb.ns['summary'])


if __name__ == '__main__':
    main()
