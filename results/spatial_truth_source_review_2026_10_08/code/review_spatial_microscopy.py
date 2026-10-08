"""Freeze directly inspected Plasmodium microscopy outcomes without truth admission.

Figure shorthand identifiers are resolved only to unique full identifiers already
published in the same paper's primer table. Unknown outcomes remain unknown.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from starplast import organisms as O  # noqa: E402

# Manual transcription of Figure 3 panel labels, visually inspected 2026-10-08.
# These are shorthand labels, not invented full gene identifiers or fine classes.
PANELS = (
    ('0932000', 'IMC', 'MTIP', 'top_left'),
    ('0920500', 'IMC', 'MTIP', 'top_right'),
    ('1137400', 'mitochondria', 'MitoTracker', 'second_left'),
    ('1137300', 'ER', 'BiP', 'second_right'),
    ('1013800', 'rhoptries', 'RAP1', 'third_left'),
    ('1136200', 'rhoptries', 'RAP1', 'third_right'),
    ('1401600', 'rhoptries', 'RAP1', 'fourth_left'),
    ('0723300', 'host PM', 'KAHRP', 'fourth_right'),
    ('0419500', 'plasma membrane', 'MSP3', 'bottom_left'),
)


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def review(prior, table_audit, output):
    """Retain all attempted targets, verified Figure 3 positives and unresolved assay scope."""
    prior, table_audit, output = Path(prior), Path(table_audit), Path(output)
    if output.exists():
        raise ValueError('Use a new microscopy-review snapshot')
    sources = json.loads((prior/'sources.json').read_text())
    source = next(r for r in sources if r['source_id'] == 'pf_spatial_proteome')
    primer = next(c['file'] for c in source['companions'] if c['declared_primary_filename'].endswith('MOESM4_ESM.xlsx'))
    figure = next(r['file'] for r in json.loads((table_audit/'media_receipts.json').read_text())
        if r['source_id'] == source['source_id'] and r['figure_or_table'] == 'Fig3')
    for binding in (primer, figure, source['existing_file']):
        if _sha(binding['path']) != binding['sha256']:
            raise ValueError('Pinned primary microscopy input changed')
    attempts = pd.read_excel(primer['path'], sheet_name='Oligonucleotides')['Gene'].dropna().astype(str).unique().tolist()
    if len(attempts) != 15:
        raise ValueError('Attempted-target table differs from reviewed primary cohort')
    reference_path = Path(O.nodes_path(O.FALCIPARUM))
    reference = set(pd.read_parquet(reference_path, columns=['gene_id']).gene_id.astype(str))
    table = pd.read_excel(source['existing_file']['path'], sheet_name='PlasmoLOPIT_data_summary', dtype=str)
    if table['Gene accession'].duplicated().any():
        raise ValueError('Do not choose between duplicate source accessions')
    table = table.set_index('Gene accession')
    positives = {}
    for shorthand, compartment, stain, panel in PANELS:
        matches = [gene for gene in attempts if gene.rsplit('_', 1)[-1] == shorthand]
        if len(matches) != 1 or matches[0] in positives:
            raise ValueError('Figure shorthand has no unique primary-table association')
        positives[matches[0]] = (shorthand, compartment, stain, panel)
    rows = []
    for gene in attempts:
        panel = positives.get(gene)
        row = {'source_id': source['source_id'], 'organism': O.FALCIPARUM, 'gene_id': gene,
            'attempted_tagging': True, 'IFA_detected': True if panel else None,
            'experimental_compartment': panel[1] if panel else None,
            'figure_shorthand': panel[0] if panel else None, 'co_stain': panel[2] if panel else None,
            'figure_panel': panel[3] if panel else None,
            'mapping_basis': 'unique suffix match to same-paper full primer-table ID' if panel else 'primary primer-table full ID',
            'exact_installed_id': gene in reference, 'protein_in_spatial_table': gene in table.index,
            'assay_outcome_grade': 'direct_experiment_candidate' if panel else 'unresolved',
            'benchmark_admitted': False,
            'outcome_scope': 'selected endogenous C-terminal HA-tagged reporter line; single experiment; primary Figure 3',
            'negative_semantics': 'no observed IFA outcome is unknown, not failed localization or biological absence'}
        for column in ('markers (S1-S2)', 'markers (S1-S2-S3)', 'Final location (S1-S2)', 'Final location (S1-S2-S3)'):
            value = table.at[gene, column] if gene in table.index else None
            row[column] = value if pd.notna(value) else None
        rows.append(row)
    output.mkdir(parents=True)
    (output/'outcomes.json').write_text(json.dumps(rows, indent=2, allow_nan=False)+'\n')
    summary = {'primary_attempted_targets': len(rows), 'Figure_3_detected_IFA_targets': len(positives),
        'attempted_without_resolved_IFA_outcome': len(rows)-len(positives),
        'IFA_exact_installed_IDs': sum(r['IFA_detected'] is True and r['exact_installed_id'] for r in rows),
        'IFA_source_marker_rows_S1_S2': sum(r['IFA_detected'] is True and r['markers (S1-S2)'] not in (None, 'unknown') for r in rows),
        'biological_benchmarks_admitted': 0,
        'gaps': ['Figure labels are coarse compartments; no equivalence to refined classifier niches assumed',
            'Target selection used hyperLOPIT profiles; evaluate only this selected tagged cohort, not proteome-wide accuracy',
            'Reporter tagging/one experiment limits wild-type and replication interpretation',
            'Need explicit target taxonomy, feature/source exclusions, context and frozen truth before benchmark admission',
            'The six unresolved IFA outcomes are not negatives; PCR/Western status remains per-target unreviewed']}
    inputs = [Path(__file__), prior/'sources.json', table_audit/'media_receipts.json',
        reference_path, *(Path(b['path']) for b in (primer, figure, source['existing_file']))]
    (output/'manifest.json').write_text(json.dumps({'created_utc': datetime.now(timezone.utc).isoformat(),
        'summary': summary, 'input_sha256': {str(p.resolve()): _sha(p) for p in inputs},
        'output_sha256': {'outcomes.json': _sha(output/'outcomes.json')}}, indent=2)+'\n')
    return summary


def main():
    """Execute the Figure 3 transcription/mapping review in an annotated notebook."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--prior', type=Path, required=True)
    parser.add_argument('--table-audit', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    nb = ExecutedNotebook('Primary microscopy outcomes: selected Plasmodium Figure 3 cohort')
    nb.md('GT-SPATIAL-01: primary Figure 3 was visually inspected from the original MD5/SHA-verified JPEG. Nine panel shorthand IDs, compartment labels and co-stains are transcribed exactly in PANELS. Full gene IDs come from unique matching suffixes in the same paper\'s original primer table; the prefix is never invented. Preserve all 15 attempted targets, with six unknown IFA outcomes. No classifier label, microscopy failure or missing value supplies biological negatives.')
    nb.code('from scripts.review_spatial_microscopy import review',
        f'summary=review({str(args.prior)!r},{str(args.table_audit)!r},{str(args.out)!r})', 'summary')
    nb.md('The images establish measured reporter localization only for this selected tagged cohort. Compartment granularity, training-marker use, profiling-based target selection, source-feature dependence and assay context must be reviewed before truth admission. Single reporter-line experiment and >100 observed cells do not create >100 independent gene samples. No runtime data or calibration change; no proteome-wide accuracy claim.')
    nb.write(str(args.out/'review.ipynb'))
    print(nb.ns['summary'])


if __name__ == '__main__':
    main()
