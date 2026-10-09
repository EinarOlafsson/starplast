"""Install an externally pinned Pf reader candidate, retaining existing entries exactly."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

from starplast import organisms as _ORGANISMS
CANDIDATE=ROOT/'results/pf_functional_ec_knn_acceptance_2026_10_08_v2/functional_results.json'
CANDIDATE_SHA='c89b317a36e01c788b56985675a97f5494ba6c798ccb9ea06be1b83f009e7bf3'


def install(output, *, candidate_path=CANDIDATE, expected_sha256=CANDIDATE_SHA):
    from starplast import functional_results as F
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    shipped=ROOT/'starplast/data/functional_results.json';reader=Path(F.__file__)
    original=shipped.read_bytes();old_pin=F.FUNCTIONAL_RESULTS_SHA256
    previous=F.load(shipped,expected_sha256=old_pin)
    candidate_path=Path(candidate_path)
    candidate,=F.load(candidate_path,expected_sha256=expected_sha256)
    assert candidate.organism==_ORGANISMS.FALCIPARUM and candidate.metadata['target']=='ec_direct_complete_major_profile'
    assert len(candidate.rows)==152 and sum(not row['training_supported'] for row in candidate.rows)==4
    (output/'original_functional_results.json').write_bytes(original)
    (output/'original_functional_results.py').write_bytes(reader.read_bytes())
    (output/'installer.py').write_bytes(Path(__file__).read_bytes())
    document=json.loads(original);incoming=json.loads(candidate_path.read_bytes())
    document['benchmarks'].extend(incoming['benchmarks'])
    merged=(F._canonical(document)+'\n').encode();digest=hashlib.sha256(merged).hexdigest()
    path=output/'functional_results.json';path.write_bytes(merged)
    results=F.load(path,expected_sha256=digest)
    assert [(item.metadata,item.payloads) for item in results[:-1]]==[(item.metadata,item.payloads) for item in previous]
    assert results[-1].metadata==candidate.metadata and results[-1].payloads==candidate.payloads
    assert shipped.read_bytes()==original
    code=reader.read_text();assert code.count(old_pin)==1
    reader.write_text(code.replace(old_pin,digest));shipped.write_bytes(merged)
    result={'old_sha256':old_pin,'candidate_sha256':expected_sha256,'merged_sha256':digest,
        'organism_counts':{organism:sum(item.organism==organism for item in results) for organism in (_ORGANISMS.TOXOPLASMA,_ORGANISMS.FALCIPARUM)},'original_entries_retained_exactly':True,
        'native_candidate_retained_exactly':True,'fitting_performed':False,'biological_admission':False}
    (output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    import argparse
    from scripts.notebook_runner import ExecutedNotebook
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--candidate',type=Path,default=CANDIDATE)
    parser.add_argument('--expected-sha256',default=CANDIDATE_SHA)
    args=parser.parse_args();nb=ExecutedNotebook('Install verified Pf functional recovery')
    nb.ns.update(install=install,output=args.out,candidate_path=args.candidate,expected_sha256=args.expected_sha256)
    nb.md('Verified immutable candidate, unchanged existing entries, native controls and unsupported profiles. No fitting.')
    try:
        nb.code('result=install(output,candidate_path=candidate_path,expected_sha256=expected_sha256)','result');print(json.dumps(nb.ns['result'],indent=2))
    finally:
        args.out.mkdir(parents=True,exist_ok=True);nb.write(str(args.out/'executed.ipynb'))
        (args.out/'SHA256SUMS.txt').write_text(''.join(f'{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(args.out)}\n'
            for path in sorted(args.out.rglob('*')) if path.is_file() and path.name!='SHA256SUMS.txt'))
