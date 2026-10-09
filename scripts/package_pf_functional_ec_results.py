"""Verify and package frozen Pf EC recovery without fitting or changing native cards."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
from pathlib import Path
import shutil
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))

from starplast import organisms as _ORGANISMS


def sha(path):
    with Path(path).open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def package(source,output,strategy='feature_knn'):
    import pandas as pd
    from scripts import freeze_pf_functional_ec_knn as K
    from starplast import artifacts as A,functional_results as F,strategies as S
    from starplast.splits import read_split
    if strategy=='feature_knn':method=K
    elif strategy=='random_forest':
        from scripts import freeze_pf_functional_ec_forest as method
    else:raise ValueError('Unsupported Pf recovery adapter')
    source=Path(source).resolve();output=Path(output)
    receipt_count=0
    for line in (source/'SHA256SUMS.txt').read_text().splitlines():
        digest,name=line.split('  ',1);path=source/name
        assert not path.is_symlink() and path.resolve().is_relative_to(source.resolve())
        assert sha(path)==digest,name
        receipt_count+=1
    inputs=json.loads((source/'input_manifest.json').read_text())
    for name,digest in inputs.items():assert sha(name)==digest,name
    packet=K.load_packet()
    split,exclusions=read_split(source/'split.json')
    assert split==packet['split']
    manifest=json.loads((source/'held_out/manifest.json').read_text())
    original=json.loads((source/'summary.json').read_text())
    deps=tuple(A.Dependency(**value) for value in manifest['spec']['dependencies'])
    feature_scope=json.loads((source/'feature_scope.json').read_text())
    assert not {'ortholog_number','has_paralog'}&set(feature_scope['provenance']['selected_columns'])
    assert manifest['spec']['strategy']==strategy
    expected=method.artifact_spec(packet,deps,gaps=feature_scope['provenance']['gaps'],
        model_available=original['artifact_status']=='ok')
    artifact=A.read_artifact(source/'held_out',expected=expected,split=split)
    assert artifact.identity==original['artifact_identity']
    assert feature_scope['exclusions']==json.loads(F._canonical(asdict(exclusions)))
    expected_inputs={('table','installed_nodes'):packet['source_hash'],
        ('source','prepared_EC_packet'):K.PREPARATION_SHA,
        ('truth','complete_direct_profiles'):sha(K.PREPARATION/'direct_profiles.parquet'),
        ('split','prepared_protected'):split.identity,
        ('exclusions','training_full_source_closure'):F._hash(asdict(exclusions)),
        ('model','native_training_state'):F._hash(artifact.payloads['model_state.json'])}
    for dependency in deps:
        if dependency.kind=='code':assert dependency.sha256==sha(ROOT/dependency.name)
        else:assert dependency.sha256==expected_inputs.pop((dependency.kind,dependency.name))
    assert not expected_inputs
    assert artifact.payloads['baseline_cards.json']==packet['controls']
    rows=artifact.payloads['rows.json'];controls=packet['control_rows']
    for key in ('entity','truth','group','training_supported'):
        assert [row[key] for row in rows]==controls[key].tolist()
    assert len(rows)==152 and sum(not row['training_supported'] for row in rows)==4
    assert original['train']==600 and original['tune']==151 and original['calibration']==147
    assert original['eligible_profiles']==1050
    parquet=pd.read_parquet(source/'rows.parquet')
    assert parquet.astype(object).where(parquet.notna(),None).to_dict('records')==rows
    assert sum(row['prediction'] is not None and row['prediction']==row['truth'] for row in rows)==original['counts']['correct']
    spec=manifest['spec'];evaluation=json.loads(spec['evaluation_scope_json'])
    metadata={'profile_namespace':'ec_major','control_format':'native_ec_controls_v1',
        'organism':_ORGANISMS.FALCIPARUM,'target':spec['target'],'strategy':spec['strategy'],'benchmark_id':spec['benchmark_id'],
        'truth_grade':'unresolved','source_targets':['ec_number'],
        'source_context':'Installed direct Plasmodium EC annotations; individual curation/prediction provenance unresolved. Orthology-derived annotations kept separate.',
        'source_artifact_identity':artifact.identity,'source_manifest':manifest,
        'source_directory':str(source.relative_to(ROOT)),'source_summary':original,'summary':original,
        'source_summary_sha256':sha(source/'summary.json'),'biological_admission':False,'calibrated_confidence':None,
        'evaluation_scope':evaluation}
    payloads={name:artifact.payloads[name] for name in F._PAYLOAD_NAMES}
    candidate=output/'functional_results.json'
    candidate.write_text(F._canonical({'schema_version':1,'benchmarks':[{'metadata':metadata,'payloads':payloads}]})+'\n')
    benchmark,=F.load(candidate,expected_sha256=sha(candidate))
    assert benchmark.payloads==payloads and benchmark.rows==rows
    nodes=pd.read_parquet(packet['source_path'])
    F.require_context(benchmark,S.Context(nodes,graph={},organism=_ORGANISMS.FALCIPARUM))
    for name,digest in inputs.items():assert sha(name)==digest,name
    result={'source_artifact_identity':artifact.identity,'source_outputs_verified':receipt_count,
        'source_inputs_verified':len(inputs),'artifact_receipts_verified':len(manifest['files']),
        'candidate_sha256':sha(candidate),'test_genes':152,'unsupported_test_genes':4,
        'counts':original['counts'],'metrics':original['metrics'],'controls_retained_exactly':True,
        'installed_context_verified':True,'biological_admission':False,'shipped_bundle_changed':False}
    (output/'acceptance.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def main():
    from scripts.notebook_runner import ExecutedNotebook
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--strategy',choices=('feature_knn','random_forest'),default='feature_knn')
    args=parser.parse_args();args.out.mkdir(parents=True,exist_ok=False)
    (args.out/'code').mkdir()
    for path in (Path(__file__),ROOT/'starplast/functional_results.py',ROOT/'tests/test_functional_native_ec_reader.py'):
        shutil.copyfile(path,args.out/'code'/path.name)
    nb=ExecutedNotebook('Independent Pf EC receipt, population and reader acceptance')
    nb.ns.update(package=package,source=args.source,output=args.out,strategy=args.strategy)
    nb.md('Native controls and unsupported genes retained exactly. No fitting or shipped bundle change.')
    try:
        nb.code('result=package(source,output,strategy)','result')
        print(json.dumps(nb.ns['result'],indent=2))
    except Exception as exc:
        (args.out/'failure.json').write_text(json.dumps({'error':type(exc).__name__,'detail':str(exc)})+'\n')
        raise
    finally:
        nb.write(str(args.out/'executed.ipynb'))
        (args.out/'SHA256SUMS.txt').write_text(''.join(f'{sha(path)}  {path.relative_to(args.out)}\n'
            for path in sorted(args.out.rglob('*')) if path.is_file() and path.name!='SHA256SUMS.txt'))


if __name__=='__main__':main()
