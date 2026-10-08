"""Merge validated frozen functional views without modifying either source artifact."""
import hashlib
import json
from pathlib import Path
import sys
ROOT = Path('/media/carruthers/mnt3/claude/repo/starplast')
sys.path.insert(0, str(ROOT))
from scripts.notebook_runner import ExecutedNotebook
from starplast import functional_results as F
OUT = Path(__file__).resolve().parent
OLD = ROOT / 'starplast/data/functional_results.json'
ADD = ROOT / 'results/functional_profile_ui_bundle_2026_10_08/functional_results.json'


def sha(path):
    with Path(path).open('rb') as file:
        return hashlib.file_digest(file, 'sha256').hexdigest()


def merge():
    old_sha = sha(OLD)
    assert old_sha == F.FUNCTIONAL_RESULTS_SHA256
    add_sha = json.loads((ADD.parent / 'receipt.json').read_text())['sha256']
    code = [Path(F.__file__), ROOT / 'scripts/notebook_runner.py']
    before = {str(path): sha(path) for path in (OLD, ADD, Path(__file__), *code)}
    old = F.load(OLD, expected_sha256=old_sha)
    added = F.load(ADD, expected_sha256=add_sha)
    assert len(old) == len(added) == 1
    assert old[0].namespace == 'ec_major' and added[0].namespace == 'pfam'
    entries = json.loads(OLD.read_text())['benchmarks'] + json.loads(ADD.read_text())['benchmarks']
    output = OUT / 'functional_results.json'
    output.write_text(json.dumps({'schema_version': 1, 'benchmarks': entries},
        sort_keys=True, ensure_ascii=False, allow_nan=False, separators=(',', ':')) + '\n')
    digest = sha(output)
    merged = F.load(output, expected_sha256=digest)
    assert [item.metadata for item in merged] == [item.metadata for item in old + added]
    assert [item.payloads for item in merged] == [item.payloads for item in old + added]
    assert before == {name: sha(name) for name in before}
    (OUT / 'original_ec_bundle.json').write_bytes(OLD.read_bytes())
    (OUT / 'code').mkdir()
    for path in code:
        (OUT / 'code' / path.name).write_bytes(path.read_bytes())
    (OUT / 'input_manifest.json').write_text(json.dumps(before, indent=2) + '\n')
    return {'sha256': digest, 'bytes': output.stat().st_size, 'benchmarks': 2,
        'original_ec_sha256': old_sha, 'addition_sha256': add_sha,
        'source_artifacts': [item.metadata['source_artifact_identity'] for item in merged],
        'rows': [len(item.rows) for item in merged], 'source_entries_exact': True,
        'biological_admission': False, 'fitting_performed': False}


nb = ExecutedNotebook('Merge exact enzyme and complete-Pfam frozen recovery views')
nb.ns.update(merge=merge)
nb.md('Verify both externally pinned source bundles, retain every native payload and scope, refuse changed input bytes, and validate the combined bundle before installation. No fitting or biological admission.')
nb.code('receipt = merge()', 'receipt')
(OUT / 'receipt.json').write_text(json.dumps(nb.ns['receipt'], indent=2) + '\n')
nb.write(str(OUT / 'executed.ipynb'))
(OUT / 'SHA256SUMS.txt').write_text('\n'.join(f'{sha(path)}  {path.relative_to(OUT)}'
    for path in sorted(OUT.rglob('*')) if path.is_file() and path.name != 'SHA256SUMS.txt') + '\n')
print(json.dumps(nb.ns['receipt']))
