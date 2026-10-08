"""Install the accepted Pf reader candidate, retaining both existing Tg entries exactly."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
CANDIDATE=ROOT/'results/pf_functional_ec_knn_acceptance_2026_10_08_v2/functional_results.json'
CANDIDATE_SHA='c89b317a36e01c788b56985675a97f5494ba6c798ccb9ea06be1b83f009e7bf3'


def install(output):
    from starplast import functional_results as F
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    shipped=ROOT/'starplast/data/functional_results.json';reader=Path(F.__file__)
    original=shipped.read_bytes();old_pin=F.FUNCTIONAL_RESULTS_SHA256
    previous=F.load(shipped,expected_sha256=old_pin)
    candidate,=F.load(CANDIDATE,expected_sha256=CANDIDATE_SHA)
    assert len(previous)==2 and all(item.organism=='Tg' for item in previous)
    assert candidate.organism=='Pf' and candidate.card['counts']['correct']==40
    assert len(candidate.rows)==152 and sum(not row['training_supported'] for row in candidate.rows)==4
    (output/'original_functional_results.json').write_bytes(original)
    (output/'original_functional_results.py').write_bytes(reader.read_bytes())
    (output/'installer.py').write_bytes(Path(__file__).read_bytes())
    document=json.loads(original);incoming=json.loads(CANDIDATE.read_bytes())
    document['benchmarks'].extend(incoming['benchmarks'])
    merged=(F._canonical(document)+'\n').encode();digest=hashlib.sha256(merged).hexdigest()
    path=output/'functional_results.json';path.write_bytes(merged)
    results=F.load(path,expected_sha256=digest)
    assert [(item.metadata,item.payloads) for item in results[:2]]==[(item.metadata,item.payloads) for item in previous]
    assert results[2].metadata==candidate.metadata and results[2].payloads==candidate.payloads
    assert shipped.read_bytes()==original
    code=reader.read_text();assert code.count(old_pin)==1
    reader.write_text(code.replace(old_pin,digest));shipped.write_bytes(merged)
    result={'old_sha256':old_pin,'candidate_sha256':CANDIDATE_SHA,'merged_sha256':digest,
        'organism_counts':{'Tg':2,'Pf':1},'original_entries_retained_exactly':True,
        'native_candidate_retained_exactly':True,'fitting_performed':False,'biological_admission':False}
    (output/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    import argparse
    from scripts.notebook_runner import ExecutedNotebook
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();nb=ExecutedNotebook('Install verified Pf functional recovery')
    nb.ns.update(install=install,output=args.out)
    nb.md('Verified immutable candidate, unchanged existing entries, native controls and unsupported profiles. No fitting.')
    try:
        nb.code('result=install(output)','result');print(json.dumps(nb.ns['result'],indent=2))
    finally:
        args.out.mkdir(parents=True,exist_ok=True);nb.write(str(args.out/'executed.ipynb'))
        (args.out/'SHA256SUMS.txt').write_text(''.join(f'{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(args.out)}\n'
            for path in sorted(args.out.rglob('*')) if path.is_file() and path.name!='SHA256SUMS.txt'))
