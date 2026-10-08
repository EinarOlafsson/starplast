"""Independent receipts and exact denominator verification of the frozen Pfam pilot."""
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path('/media/carruthers/mnt3/claude/repo/starplast')
sys.path.insert(0, str(ROOT))
from scripts.notebook_runner import ExecutedNotebook
OUT = Path(__file__).resolve().parent
PILOT = ROOT / 'results/functional_profile_knn_2026_10_08_v4'


def sha(path):
    with Path(path).open('rb') as file:
        return hashlib.file_digest(file, 'sha256').hexdigest()


def verify():
    output_count = input_count = artifact_count = 0
    for line in (PILOT / 'SHA256SUMS.txt').read_text().splitlines():
        digest, name = line.split('  ', 1)
        assert sha(PILOT / name) == digest, name
        output_count += 1
    for name, digest in json.loads((PILOT / 'input_manifest.json').read_text()).items():
        assert sha(name) == digest, name
        input_count += 1
    manifest = json.loads((PILOT / 'held_out/manifest.json').read_text())
    for name, receipt in manifest['files'].items():
        path = PILOT / 'held_out' / name
        assert sha(path) == receipt['sha256'] and path.stat().st_size == receipt['bytes'], name
        artifact_count += 1
    summary = json.loads((PILOT / 'summary.json').read_text())
    rows = json.loads((PILOT / 'held_out/rows.json').read_text())
    assert len(rows) == len({row['entity'] for row in rows}) == 635
    assert [row['entity'] for row in rows] == manifest['spec']['entity_order']
    assert not set(manifest['spec']['fit_entities']) & {row['entity'] for row in rows}
    assert len(manifest['spec']['fit_entities']) == 2381
    answered = sum(row['prediction'] is not None for row in rows)
    correct = sum(row['prediction'] is not None and row['prediction'] == row['truth'] for row in rows)
    assert answered == 20 and correct == 12
    assert sum(not row['training_supported'] for row in rows) == 333
    assert all(row['abstained'] == (row['prediction'] is None) for row in rows)
    assert summary['counts'] == {'eligible':635, 'answered':20, 'abstained':615,
        'correct':12, 'wrong':8, 'unique_biological_entities':635}
    assert summary['metrics']['accuracy'] == 12 / 635
    assert summary['metrics']['coverage'] == 20 / 635
    assert summary['metrics']['precision_of_calls'] == 12 / 20
    assert summary['baseline_comparisons']['training_majority']['metrics']['accuracy'] == 19 / 635
    assert summary['baseline_comparisons']['training_prevalence']['metrics']['accuracy'] == 2 / 635
    assert all(all(value[key] for key in ('same_entity_order','same_truth','same_class_order'))
        for value in summary['baseline_comparisons'].values())
    assert len(summary['selected_feature_columns']) == len(summary['kept_feature_columns']) == 348
    assert summary['biological_admission'] is False and summary['biological_accuracy'] is None
    assert summary['artifact_identity'] == manifest['identity']
    return {'output_receipts':output_count, 'input_receipts':input_count,
        'artifact_receipts':artifact_count, 'counts':summary['counts'],
        'artifact_identity':manifest['identity'], 'exact_denominators':True}


nb = ExecutedNotebook('Independent acceptance of frozen complete-Pfam kNN recovery')
nb.ns.update(verify=verify)
nb.md('Verify original source/code and every output/artifact receipt. Independently recount exact test outcomes and controls. Annotation recovery remains separate from biological validation.')
nb.code('result = verify()', 'result')
(OUT / 'summary.json').write_text(json.dumps(nb.ns['result'], indent=2) + '\n')
nb.write(str(OUT / 'executed.ipynb'))
(OUT / 'SHA256SUMS.txt').write_text('\n'.join(f'{sha(path)}  {path.relative_to(OUT)}'
    for path in sorted(OUT.rglob('*')) if path.is_file() and path.name != 'SHA256SUMS.txt') + '\n')
print(json.dumps(nb.ns['result']))
