"""Package a checksum-pinned functional pilot view without altering its immutable model."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from starplast import functional_results as F


def package(pilot, output):
    """Verify every source payload byte and retain its scope in a compact offline view."""
    pilot=Path(pilot).resolve();output=Path(output)
    if output.exists():raise ValueError('Use a new immutable packaged-view output')
    manifest=json.loads((pilot/'held_out/manifest.json').read_text())
    for name,receipt in manifest['files'].items():
        data=(pilot/'held_out'/name).read_bytes()
        if hashlib.sha256(data).hexdigest()!=receipt['sha256'] or len(data)!=receipt['bytes']:
            raise ValueError('Source artifact payload changed: '+name)
    summary=json.loads((pilot/'summary.json').read_text())
    if summary['artifact_identity']!=manifest['identity']:
        raise ValueError('Pilot summary and source artifact differ')
    payloads={name:json.loads((pilot/'held_out'/name).read_text()) for name in F._PAYLOAD_NAMES}
    spec=manifest['spec'];evaluation=json.loads(spec['evaluation_scope_json'])
    metadata={'organism':summary['organism'],'target':summary['target'],'strategy':spec['strategy'],
        'benchmark_id':spec['benchmark_id'],'truth_grade':evaluation['truth_grade'],
        'source_targets':['ec_number'],'source_context':'Installed Toxoplasma gene EC annotations resolved through pinned ENZYME nomenclature; gene-level curation/prediction provenance unresolved.',
        'source_artifact_identity':manifest['identity'],'source_manifest':manifest,
        'source_directory':str(pilot.relative_to(ROOT)),
        'upstream_model_reference':{'path':'held_out/model_state.json',**manifest['files']['model_state.json']},
        'summary':summary,'biological_admission':False,'calibrated_confidence':None}
    document={'schema_version':1,'benchmarks':[{'metadata':metadata,'payloads':payloads}]}
    output.mkdir(parents=True)
    path=output/'functional_results.json'
    path.write_text(json.dumps(document,sort_keys=True,ensure_ascii=False,allow_nan=False,indent=2)+'\n')
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    loaded=F.load(path,expected_sha256=digest)
    assert len(loaded)==1 and loaded[0].rows==payloads['rows.json']
    (output/'receipt.json').write_text(json.dumps({'sha256':digest,'source_artifact_identity':manifest['identity'],
        'rows':len(loaded[0].rows),'source_payloads':'all source artifact payload bytes verified; displayed payloads retain exact source hashes',
        'biological_admission':False},indent=2)+'\n')
    (output/'code').mkdir()
    code=[Path(__file__),Path(F.__file__),ROOT/'scripts/notebook_runner.py']
    for source in code:(output/'code'/source.name).write_bytes(source.read_bytes())
    (output/'manifest.json').write_text(json.dumps({
        'input_sha256':{str(path):hashlib.sha256(path.read_bytes()).hexdigest()
            for path in [pilot/'held_out/manifest.json',pilot/'summary.json',*code]},
        'output_sha256':{str(path.relative_to(output)):hashlib.sha256(path.read_bytes()).hexdigest()
            for path in output.rglob('*') if path.is_file()}},indent=2)+'\n')
    return {'sha256':digest,'rows':len(loaded[0].rows),'source_artifact_identity':manifest['identity']}


def main():
    """Execute packaging and verification in an inspectable notebook."""
    from notebook_runner import ExecutedNotebook
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--pilot',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    notebook=ExecutedNotebook('Packaged offline functional annotation-recovery scorecards')
    notebook.md('Compact UI view of the immutable FN-EC-01 reference-recovery pilot. All original artifact payload hashes and exact displayed rows/cards are retained; model state remains upstream. No independent biological admission, calibrated confidence or new unknown-gene claims.')
    try:
        notebook.code('from scripts.package_functional_results import package',
            f'receipt=package({str(args.pilot)!r},{str(args.out)!r})','receipt')
    except Exception as exc:
        args.out.mkdir(parents=True,exist_ok=True)
        notebook.md('Diagnostic: '+type(exc).__name__+': '+str(exc))
        notebook.write(str(args.out/'diagnostic.ipynb'))
        raise
    notebook.write(str(args.out/'packaging.ipynb'))
    print(json.dumps(notebook.ns['receipt'],indent=2))


if __name__=='__main__':main()
