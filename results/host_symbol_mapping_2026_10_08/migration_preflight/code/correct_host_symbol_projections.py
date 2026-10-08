"""Withdraw verified ambiguous rhoptry protein projections without losing gene evidence.

66.04 migration uses the frozen builder review, exact source scores and a
per-cell withdrawal ledger. Every unrelated installed cell and host row remains.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from starplast import organisms as O  # noqa: E402


def _sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def migrate(review, backup, output):
    """Apply only frozen, source-verified withdrawals and preserve original input backups."""
    review, backup, output = Path(review), Path(backup), Path(output)
    if output.exists() or backup.exists():
        raise ValueError('Use new migration and backup directories')
    data = ROOT/'starplast/data'
    target = data/O.HOST_TABLES[O.HUMAN]
    cache = data/'deposit_host_k562_rhoptry_screen.tsv'
    snapshot = json.loads((review/'manifest.json').read_text())
    expected = snapshot['input_sha256'][str(target)]
    if _sha(target) != expected:
        raise ValueError('Installed table differs from reviewed revision')
    for filename, digest in snapshot['outputs'].items():
        if _sha(review/filename) != digest:
            raise ValueError('Frozen mapping/gene review changed')
    evidence = pd.read_parquet(review/'rhoptry_gene_evidence.parquet')
    ambiguous = pd.read_parquet(review/'ambiguous_legacy_projection_review.parquet')
    corrected = pd.read_parquet(review/'rhoptry_corrected_projection_candidate.parquet')
    if len(evidence) != 20010 or len(ambiguous) != 39 or len(corrected) != 18700:
        raise ValueError('Frozen source/projection populations differ')
    before = pd.read_parquet(target)
    after = before.copy(deep=True)
    indexed = before.set_index('host_id')
    fields = ['rhoptry_discharge_score', 'rhoptry_discharge_beta', 'rhoptry_discharge_fdr']
    records = []
    for _, row in ambiguous.iterrows():
        accession = row.legacy_assigned_host_id
        alternatives = list(row.candidate_accessions)
        if len(alternatives) < 2 or accession not in alternatives or accession in set(corrected.host_id):
            raise ValueError('Withdrawal is not exclusively ambiguous')
        position = before.index[before.host_id == accession]
        if len(position) != 1:
            raise ValueError('Withdrawal has no unique installed endpoint')
        for field in fields:
            old = indexed.at[accession, field]
            if pd.isna(old) or old != row[field]:
                raise ValueError('Installed value differs from original ambiguous gene score')
            records.append({'organism': O.HUMAN, 'source_id': 'host_k562_rhoptry_screen',
                'host_id': accession, 'field': field, 'source_row': int(row.source_row),
                'gene_symbol': row.gene_symbol, 'old_value': float(old),
                'candidate_accessions': alternatives, 'reason': 'ambiguous_gene_to_reviewed_proteins',
                'new_semantics': 'withheld protein projection; original measured gene score retained'})
            after.at[position[0], field] = float('nan')
    # Unrelated columns, identities, order and unambiguous measurements are exact.
    pd.testing.assert_frame_equal(after.drop(columns=fields), before.drop(columns=fields))
    affected = set(ambiguous.legacy_assigned_host_id)
    pd.testing.assert_frame_equal(after.loc[~after.host_id.isin(affected)], before.loc[~before.host_id.isin(affected)])
    after_index = after.set_index('host_id')
    for field in fields:
        pd.testing.assert_series_equal(after_index.loc[corrected.host_id, field].sort_index(),
            corrected.set_index('host_id')[field].sort_index(), check_dtype=False)
    # Keep exact original TSV serialization for every retained cache row.
    lines = cache.read_text().splitlines(keepends=True)
    header = lines[0].rstrip('\n').split('\t')
    id_pos = header.index('host_id')
    retained = [line for line in lines[1:] if line.rstrip('\n').split('\t')[id_pos] not in affected]
    removed = [line for line in lines[1:] if line.rstrip('\n').split('\t')[id_pos] in affected]
    if len(removed) != 39 or len(retained) != 18700:
        raise ValueError('Derived cache does not match the reviewed withdrawal population')
    gene_target = data/'deposit_host_k562_rhoptry_gene_evidence.parquet'
    ledger_target = data/'host_projection_withdrawals.json'
    if gene_target.exists() or ledger_target.exists():
        raise ValueError('An existing source revision must not be overwritten')
    output.mkdir(parents=True)
    backup.mkdir(parents=True)
    originals = [target, cache, data/O.HOST_TABLES[O.MOUSE], data/'nodes.parquet', data/'pf_nodes.parquet',
        data/'graph.npz', data/'pf_graph.npz', data/'track_record.parquet', data/'claims.parquet']
    original_hashes = {str(path): _sha(path) for path in originals}
    shutil.copy2(target, backup/target.name)
    shutil.copy2(cache, backup/cache.name)
    shutil.copy2(review/'rhoptry_gene_evidence.parquet', gene_target)
    after.to_parquet(target, index=False)
    cache.write_text(lines[0]+''.join(retained))
    ledger = {'schema_version': 1, 'revision': '66.04_2026_10_08', 'source_id': 'host_k562_rhoptry_screen',
        'gene_evidence_file': gene_target.name, 'gene_evidence_sha256': _sha(gene_target),
        'original_gene_population': len(evidence), 'measured_unit': 'gene_perturbation',
        'target_file': target.name, 'before_sha256': original_hashes[str(target)], 'after_sha256': _sha(target),
        'original_source_sha256': snapshot['input_sha256'][next(p for p in snapshot['input_sha256'] if p.endswith('/v2_media-1.xlsx'))],
        'records': records}
    ledger_target.write_text(json.dumps(ledger, indent=2, allow_nan=False)+'\n')
    untouched = [path for path in originals if path not in (target, cache)]
    if any(_sha(path) != original_hashes[str(path)] for path in untouched):
        raise ValueError('An unrelated runtime input changed')
    pd.testing.assert_frame_equal(pd.read_parquet(target), after)
    summary = {'withdrawn_protein_projections': len(affected), 'withdrawn_cells': len(records),
        'preserved_gene_records': len(evidence), 'preserved_host_rows': len(after),
        'unambiguous_protein_values': len(corrected), 'unchanged_parasite_inference_inputs': True,
        'unchanged_mouse_table': True, 'calibration_scope': 'no strategy, parasite feature/graph, track-record or claim input changed'}
    manifest = {'summary': summary, 'original_runtime_sha256': original_hashes,
        'backup_sha256': {str(p): _sha(p) for p in backup.iterdir()},
        'changed_sha256': {str(p): _sha(p) for p in (target, cache, gene_target, ledger_target)},
        'input_sha256': {str(Path(__file__)): _sha(__file__), str(review/'manifest.json'): _sha(review/'manifest.json')}}
    (output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    return summary


def main():
    """Execute the explicit source revision in an annotated migration notebook."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--review', type=Path, required=True)
    parser.add_argument('--backup', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    nb = ExecutedNotebook('Correct installed ambiguous host projections and preserve original gene evidence')
    nb.md('66.04 second partition: withdraw only 39 verified ambiguous rhoptry protein projections, retain original gene scores and accession alternatives, and write a per-cell withdrawal ledger. Require exact original source scores and unique reviewed endpoints before changing anything. Preserve every host row and unrelated cell; retain exact TSV serialization of all unambiguous cache rows. Archive previous installed inputs externally. No inferred absence or numeric zero replaces a withheld value.')
    nb.code('from scripts.correct_host_symbol_projections import migrate',
        f'summary=migrate({str(args.review)!r},{str(args.backup)!r},{str(args.out)!r})', 'summary')
    nb.md('All 20,010 source gene records remain available. Reduced protein coverage corrects a false entity assignment; it does not erase a gene experiment. Historical non-loss counts must reconcile through the source-verified ledger. Parasite feature/graph, track-record/claim and mouse inputs are byte-identical; no existing calibration is refitted or newly claimed. Run affected host/deposit/provenance/leakage/layout/catalogue and unchanged-inference checks before item completion.')
    nb.write(str(args.out/'migration.ipynb'))
    print(nb.ns['summary'])


if __name__ == '__main__':
    main()
