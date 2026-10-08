"""Package verified complete-Pfam recovery records with explicit source complement semantics.

This bounded reader-only operation verifies the original artifact bytes and
retains native rows, cards, profile classes and controls without fitting or
rewriting their source identities. Derived domain cards recover recorded source
presence only; they do not establish biological absence or gene confidence.
Large native profile-card payloads are streamed into compact JSON without model
state duplication. Every execution and diagnostic is retained in a notebook.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from starplast import functional_results as F,organisms as O

SOURCE=ROOT/'results/functional_profile_knn_2026_10_08_v4'
SOURCE_ID='38cfb6bf110cadb559eac41c798298d7d95882efe2524a3591ccf7ad63d0385f'
NATIVE=('rows.json','card.json','class_cards.json','baseline_cards.json')


def _sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def _write(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def package(output):
    """Verify source bytes and freeze a compact native-plus-derived offline view."""
    output=Path(output)
    if output.exists():raise ValueError('Use a new immutable profile UI bundle directory')
    output.mkdir(parents=True)
    manifest=json.loads((SOURCE/'held_out/manifest.json').read_text())
    if manifest['identity']!=SOURCE_ID or F._hash({k:v for k,v in manifest.items() if k!='identity'})!=SOURCE_ID:
        raise ValueError('Only the accepted verified Pfam V4 source is permitted')
    receipts=[]
    for name,receipt in manifest['files'].items():
        path=SOURCE/'held_out'/name
        if path.is_symlink() or _sha(path)!=receipt['sha256'] or path.stat().st_size!=receipt['bytes']:
            raise ValueError('Original source artifact payload changed: '+name)
        receipts.append({'name':name,'path':str(path),**receipt})
    spec=manifest['spec'];evaluation=json.loads(spec['evaluation_scope_json'])
    table=next(d for d in spec['dependencies'] if d['kind']=='table' and d['name']=='installed_nodes')
    if _sha(O.nodes_path(O.TOXOPLASMA))!=table['sha256']:raise ValueError('Installed source differs from original artifact')
    original=json.loads((SOURCE/'summary.json').read_text())
    if original['artifact_identity']!=SOURCE_ID or original['test_genes']!=635 or original['unsupported_test_genes']!=333:
        raise ValueError('Keep every accepted test gene and unsupported reference profile')
    mapping={'train':'train_genes','test':'test_genes','eligible_profiles':'eligible_genes'}
    summary={**original,**{name:original[source] for name,source in mapping.items()}}
    split_path=SOURCE/'split.json';split_document=json.loads(split_path.read_text())
    split=split_document['split'];split_identity=F._hash(split)
    split_body={k:v for k,v in split_document.items() if k!='content_sha256'}
    if (split_identity!=split_document['identity'] or not any(d['kind']=='split' and d['sha256']==split_identity for d in spec['dependencies'])
            or hashlib.sha256(json.dumps(split_body,sort_keys=True).encode()).hexdigest()!=split_document['content_sha256']):
        raise ValueError('Original protected split content/identity changed')
    role_mapping={role:role for role in ('train','tune','calibration','test')}
    summary.update({role:sum(row['role']==role for row in split['assignments']) for role in role_mapping})
    rows=json.loads((SOURCE/'held_out/rows.json').read_text())
    members=F.member_cards(rows,'pfam')
    metadata={'profile_namespace':'pfam','organism':O.TOXOPLASMA,'target':spec['target'],'strategy':spec['strategy'],
        'benchmark_id':spec['benchmark_id'],'truth_grade':evaluation['truth_grade'],'source_targets':['pfam_id'],
        'source_context':'Installed Toxoplasma complete Pfam annotation profiles; assignment evidence/release and source/homology independence unresolved.',
        'source_artifact_identity':SOURCE_ID,'source_manifest':manifest,'source_directory':str(SOURCE.relative_to(ROOT)),
        'upstream_model_reference':{'path':'held_out/model_state.json',**manifest['files']['model_state.json']},
        'source_summary':original,'source_summary_sha256':_sha(SOURCE/'summary.json'),'summary':summary,'summary_mapping':mapping,
        'summary_role_mapping':role_mapping,'split_manifest':split,
        'source_split_receipt':{'path':str(split_path.relative_to(ROOT)),'sha256':_sha(split_path),
            'bytes':split_path.stat().st_size,'identity':split_identity},
        'membership_derivation':{'recipe':'recorded_profile_member_cards_v1','namespace':'pfam',
            'source_rows_sha256':manifest['files']['rows.json']['sha256'],'payload_sha256':F._hash(members)},
        'biological_admission':False,'calibrated_confidence':None,
        'control_interpretation':'Original control names/scopes retained; same split, target, cohort and test population, distinct control benchmark/strategy address'}
    _write(output/'scope.json',{'source_artifact_identity':SOURCE_ID,'namespace':'pfam','test_genes':635,
        'unsupported_test_genes':333,'native_payloads':list(NATIVE),'derived_payload':'major_class_cards.json',
        'fitting_performed':False,'biological_admission':False,'memory_cap_bytes':400*1024*1024})
    path=output/'functional_results.json'
    with path.open('wb') as stream:
        stream.write(b'{"schema_version":1,"benchmarks":[{"metadata":')
        stream.write(F._canonical(metadata).encode());stream.write(b',"payloads":{')
        for index,name in enumerate(NATIVE):
            if index:stream.write(b',')
            stream.write(json.dumps(name).encode()+b':')
            with (SOURCE/'held_out'/name).open('rb') as source:shutil.copyfileobj(source,stream,1024*1024)
        stream.write(b',"major_class_cards.json":'+F._canonical(members).encode()+b'}}]}\n')
    digest=_sha(path)
    benchmark=F.load(path,expected_sha256=digest)[0]
    if benchmark.namespace!='pfam' or benchmark.rows!=rows or len(benchmark.rows)!=635:
        raise ValueError('Packaged profile population differs')
    if sum(row['training_supported'] is False for row in benchmark.rows)!=333:raise ValueError('Unsupported profiles were dropped')
    for receipt in receipts:
        if _sha(receipt['path'])!=receipt['sha256']:raise ValueError('Original source changed during packaging')
    receipt={'sha256':digest,'source_artifact_identity':SOURCE_ID,'rows':635,'unsupported_test_genes':333,
        'profile_classes':len(benchmark.profile_class_cards),'domain_member_classes':len(members),
        'source_payload_receipts':len(receipts),'native_displayed_payloads':list(NATIVE),
        'member_derivation':metadata['membership_derivation'],'biological_admission':False,'calibrated_confidence':None,
        'counts':benchmark.card['counts'],'bytes':path.stat().st_size}
    _write(output/'receipt.json',receipt);_write(output/'source_receipts.json',receipts)
    return receipt


def main():
    """Execute packaging and record actual runtime, source receipts and output hashes."""
    from notebook_runner import ExecutedNotebook
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();started=time.monotonic();nb=ExecutedNotebook('Verified Pfam profile recovery: native records and domain presence cards')
    nb.md('Reader-only packaging of accepted V4: retain all 635 test genes, including 333 training-unsupported profiles and all abstentions. Original native payload filenames and hashes are unchanged. Derived domain-membership metrics describe recorded profile presence/complement only. Original control benchmark scopes remain distinct, with matched split/target/cohort. No fitting, calibration, biological admission or deployment. External MemoryMax400M.')
    try:
        nb.code('from scripts.package_functional_profile_results import package',f'receipt=package({str(args.out.resolve())!r})','receipt')
    except Exception as error:
        args.out.mkdir(parents=True,exist_ok=True);nb.md('Diagnostic: '+type(error).__name__+': '+str(error))
        nb.write(str(args.out/'diagnostic.ipynb'));raise
    nb.write(str(args.out/'packaging.ipynb'))
    code=args.out/'code';code.mkdir()
    inputs=[Path(__file__),Path(F.__file__),ROOT/'scripts/notebook_runner.py',SOURCE/'held_out/manifest.json',SOURCE/'summary.json',SOURCE/'split.json']
    for source in inputs[:3]:(code/source.name).write_bytes(source.read_bytes())
    peak=next(int(line.split()[1])*1024 for line in Path('/proc/self/status').read_text().splitlines() if line.startswith('VmHWM:'))
    _write(args.out/'runtime.json',{'seconds':time.monotonic()-started,'process_peak_bytes':peak,'external_memory_cap_bytes':400*1024*1024})
    _write(args.out/'manifest.json',{'input_sha256':{str(path):_sha(path) for path in inputs},
        'outputs':{str(path.relative_to(args.out)):{'sha256':_sha(path),'bytes':path.stat().st_size}
            for path in args.out.rglob('*') if path.is_file()}})
    print(json.dumps(nb.ns['receipt'],indent=2))


if __name__=='__main__':main()
