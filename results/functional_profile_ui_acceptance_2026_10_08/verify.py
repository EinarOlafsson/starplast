"""Independent output, source and historical-code receipt verification."""
import hashlib, json, sys
from pathlib import Path
ROOT=Path('/media/carruthers/mnt3/claude/repo/starplast')
sys.path.insert(0,str(ROOT))
from scripts.notebook_runner import ExecutedNotebook
OUT=Path(__file__).resolve().parent


def sha(path):
    with Path(path).open('rb') as file:return hashlib.file_digest(file,'sha256').hexdigest()


def verify():
    counts={}
    for name in ('functional_profile_results_merged_2026_10_08','functional_profile_ui_2026_10_08_v3'):
        folder=ROOT/'results'/name
        lines=(folder/'SHA256SUMS.txt').read_text().splitlines()
        for line in lines:
            digest,path=line.split('  ',1)
            assert sha(folder/path)==digest,path
        inputs=json.loads((folder/'input_manifest.json').read_text())
        for path,digest in inputs.items():
            candidate=Path(path)
            if sha(candidate)!=digest:
                candidate=folder/('original_ec_bundle.json' if Path(path).name=='functional_results.json' else 'code/'+Path(path).name)
            assert sha(candidate)==digest,path
        counts[name]={'outputs':len(lines),'inputs':len(inputs)}
    folder=ROOT/'results/functional_profile_ui_bundle_2026_10_08'
    manifest=json.loads((folder/'manifest.json').read_text())
    for path,receipt in manifest['outputs'].items():
        assert sha(folder/path)==receipt['sha256'] and (folder/path).stat().st_size==receipt['bytes'],path
    for path,digest in manifest['input_sha256'].items():
        candidate=Path(path)
        if sha(candidate)!=digest:candidate=folder/'code'/candidate.name
        assert sha(candidate)==digest,path
    receipts=json.loads((folder/'source_receipts.json').read_text())
    for receipt in receipts:
        path=Path(receipt['path'])
        assert sha(path)==receipt['sha256'] and path.stat().st_size==receipt['bytes']
    counts[folder.name]={'outputs':len(manifest['outputs']),'inputs':len(manifest['input_sha256']),'original_artifact_receipts':len(receipts)}
    assert sha(ROOT/'starplast/data/functional_results.json')==sha(ROOT/'results/functional_profile_results_merged_2026_10_08/functional_results.json')
    return counts


nb=ExecutedNotebook('Independent domain UI integration receipt acceptance')
nb.ns.update(verify=verify)
nb.md('Check every packet output and original artifact receipt. Historical code/bundle inputs verify against archived exact bytes after intentional checksum-pin installation; actual UI inputs verify unchanged current bytes.')
nb.code('counts=verify()','counts')
(OUT/'summary.json').write_text(json.dumps(nb.ns['counts'],indent=2)+'\n')
nb.write(str(OUT/'executed.ipynb'))
(OUT/'SHA256SUMS.txt').write_text('\n'.join(f'{sha(path)}  {path.relative_to(OUT)}' for path in sorted(OUT.rglob('*')) if path.is_file() and path.name!='SHA256SUMS.txt')+'\n')
print(json.dumps(nb.ns['counts']))
