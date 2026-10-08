"""Execute bounded numeric operators from the verified Tg EC V4 source lineage.

Only train-fit numeric preprocessing runs here: no classifier, voting, biological
evaluation, deployment, source acquisition or hidden-label feature selection. The
frozen source split and historical numeric columns constrain this pilot; current
functional source closure may deliberately withhold additional coordinates.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict,replace
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
import pandas as pd
from starplast import artifacts as A,functional_exclusions as X,organisms as O,precompute_features as F,precompute_jobs as P
from starplast.splits import read_split,write_split

SOURCE=ROOT/'results/functional_ec_knn_pilot_2026_10_08_v4'
SOURCE_ID='349317afd8f24b860d94ae1e11280a6130bf451700b02b9d2145fab85fc99162'
BUDGET=P.Budget(2,120,32*1024*1024,32*1024*1024,400*1024*1024)
LIMITS=F.FeatureLimits(1226,400,1226*400,12*1024*1024)


def _sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def _write(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def _source():
    manifest=json.loads((SOURCE/'manifest.json').read_text())
    receipts=[]
    for relative,digest in manifest['output_sha256'].items():
        path=SOURCE/relative
        if path.is_symlink() or not path.is_file() or _sha(path)!=digest:
            raise ValueError('Verified V4 source output changed: '+relative)
        receipts.append({'path':str(path),'sha256':digest,'bytes':path.stat().st_size})
    split,old_exclusions=read_split(SOURCE/'split.json')
    spec=A.spec_from_dict(json.loads((SOURCE/'held_out/manifest.json').read_text())['spec'])
    artifact=A.read_artifact(SOURCE/'held_out',expected=spec,split=split)
    if artifact.identity!=SOURCE_ID or spec.organism!=O.TOXOPLASMA:
        raise ValueError('Only the verified Tg V4 artifact is permitted')
    if F.exclusions_identity(old_exclusions) not in [d.sha256 for d in spec.dependencies if d.kind=='exclusions']:
        raise ValueError('Frozen V4 exclusion identity differs')
    nodes=Path(O.nodes_path(O.TOXOPLASMA)).resolve()
    if manifest['input_sha256'].get(str(nodes))!=_sha(nodes):
        raise ValueError('Installed numeric source differs from verified V4 input')
    if len(split.assignments)!=1226 or len(split.entities('train'))!=679:
        raise ValueError('Retain the frozen 1226-entity/679-training population')
    return split,old_exclusions,spec,artifact.payloads['model_state.json'],nodes,receipts


def prepare(output):
    """Freeze exact source, current closure, recipe and byte dependencies before jobs."""
    output=Path(output)
    if output.exists():raise ValueError('Use a new immutable numeric-pilot directory')
    output.mkdir(parents=True)
    started=time.monotonic()
    split,old_exclusions,source_spec,model,nodes,receipts=_source()
    ids=[row.entity for row in split.assignments]
    train_ids=list(split.entities('train'))
    training=pd.read_parquet(nodes,filters=[('gene_id','in',train_ids)]).set_index('gene_id',drop=False).loc[train_ids].copy()
    target=source_spec.target
    # The helper needs only training target values for empirical source closure.
    training_truth=pd.read_parquet(SOURCE/'truth.parquet',filters=[('gene_id','in',train_ids)]).set_index('gene_id')[target]
    training[target]=training_truth.loc[train_ids]
    exclusions=X.make_exclusions(training,(target,),source_targets=('ec_number',),benchmark_id=split.benchmark_id,split=split)
    # Preserve every old exclusion; the current policy can only broaden closure.
    exclusions=replace(exclusions,columns=tuple(sorted(set(exclusions.columns)|set(old_exclusions.columns))),
        layers=tuple(sorted(set(exclusions.layers)|set(old_exclusions.layers))),
        targets=tuple(dict.fromkeys((*old_exclusions.targets,*exclusions.targets))))
    historical=list(model['input_columns'])
    columns=[column for column in historical if column not in exclusions.columns]
    dropped=[column for column in historical if column not in columns]
    del training,training_truth
    features=pd.read_parquet(nodes,columns=['gene_id',*columns],filters=[('gene_id','in',ids)]).set_index('gene_id').loc[ids,columns]
    source_path=output/'numeric_source.json'
    source_path.write_text(A.canonical_object(F.snapshot_payload(features,organism=split.organism))+'\n')
    table=P.Snapshot(A.Dependency('table','permitted_numeric_inputs',_sha(source_path)),str(source_path.resolve()))
    builder=F.NumericFeatureBuilder(features,table,split=split,exclusions=exclusions,limits=LIMITS)
    write_split(output/'split.json',split,exclusions=exclusions)
    code_dir=output/'code';code_dir.mkdir()
    code_paths=[Path(__file__),ROOT/'scripts/notebook_runner.py',*(ROOT/'starplast'/name for name in
        ('precompute_features.py','precompute_jobs.py','functional_exclusions.py','artifacts.py','splits.py',
         'strategies.py','search.py','datasets.py','slots.py','discovery_labels.py','functional_domains.py',
         'capabilities.py','query.py','organisms.py','paths.py'))]
    extra=[]
    for path in code_paths:
        frozen=code_dir/path.name;frozen.write_bytes(path.read_bytes())
        extra.append(P.Snapshot(A.Dependency('code','pilot:'+path.stem,_sha(frozen)),str(frozen.resolve())))
    receipts.extend({'path':str(path),'sha256':_sha(path),'bytes':path.stat().st_size}
        for path in (SOURCE/'manifest.json',nodes))
    for receipt in receipts:
        path=Path(receipt['path'])
        name='installed_nodes' if path==nodes else 'v4:'+str(path.relative_to(SOURCE))
        extra.append(P.Snapshot(A.Dependency('source',name,receipt['sha256']),str(path)))
    # The full old artifact content identity includes the native model state.
    snapshots=builder.snapshots+tuple(extra)
    deps=tuple(snapshot.dependency for snapshot in snapshots)+(A.Dependency('artifact','verified_ec_v4',SOURCE_ID),
        A.Dependency('split','nested',split.identity),A.Dependency('exclusions','functional_source_closure',F.exclusions_identity(exclusions)))
    available=builder.availability()
    evaluation=json.loads(source_spec.evaluation_scope_json)
    evaluation.update(eligible_population=len(ids),biological_admission=False,
        interpretation='Intermediate numeric operator; no biological accuracy, inference or calibrated confidence')
    version='numeric-operator-pilot-v1;python='+sys.version.split()[0]+';numpy='+np.__version__+';pandas='+pd.__version__
    spec=replace(source_spec,entity_order=tuple(ids),role='reusable_base',partition_id='numeric_ec_v4:operator',
        dependencies=deps,settings_json=A.canonical_object(F.settings(LIMITS)),evaluation_scope_json=A.canonical_object(evaluation),
        code_version=version,confidence_kind='none',evaluation_partition='none',status=available['status'],gaps=tuple(available['gaps']))
    child=replace(spec,partition_id='numeric_ec_v4:reference',dependencies=deps+(A.Dependency('artifact','numeric_operator','0'*64),))
    plan=P.Plan((P.Job('numeric_operator',spec,split=split),
        P.Job('operator_reference',child,(P.Upstream('numeric_operator','artifact','numeric_operator'),),split)),snapshots)
    scope={'source_artifact_identity':SOURCE_ID,'source_spec_identity':source_spec.identity,'organism':O.TOXOPLASMA,
        'target':target,'split_identity':split.identity,'ordered_entities':len(ids),'training_entities':len(train_ids),
        'historical_input_columns':len(historical),'permitted_input_columns':len(columns),'deliberate_extra_exclusions':dropped,
        'added_exclusion_columns':sorted(set(exclusions.columns)-set(old_exclusions.columns)),
        'added_exclusion_layers':sorted(set(exclusions.layers)-set(old_exclusions.layers)),
        'exclusion_policy':X.POLICY_VERSION,'exclusion_identity':F.exclusions_identity(exclusions),
        'historical_exclusion_identity':F.exclusions_identity(old_exclusions),'numeric_source_sha256':table.dependency.sha256,
        'raw_installed_nodes_sha256':_sha(nodes),'plan_identity':plan.identity,'job_fingerprints':dict(plan.fingerprints()),
        'budget':asdict(BUDGET),'feature_limits':asdict(LIMITS),'snapshots':[asdict(s) for s in snapshots],
        'phase':F.PHASE,'inference_outputs':False,'biological_admission':False,'biological_accuracy':None,
        'calibrated_confidence':None,'deployment':'unavailable','runtime_versions':version,
        'preparation_seconds':time.monotonic()-started,'preparation_process_peak_rss_bytes':P._peak_rss(),
        'measurement_scope':'Runner RSS at checkpoints and executed-process lifetime peak; not isolated per-job peak',
        'limitations':['Train-fit ranks and numeric operator only; no classifier, graph, fitting of biological inference or acquisition',
            'Historical source annotation evidence remains unresolved; this software parity does not establish biological accuracy',
            'External MemoryMax=400M required; 64.25 remains incomplete']}
    _write(output/'scope.json',scope);_write(output/'source_receipts.json',receipts)
    return {'output':output,'builder':builder,'plan':plan,'scope':scope,'model':model,'calls':[]}


def stage(session,name):
    """Run controlled interruption, dependent-only resume and verified no-builder replay."""
    session['calls'].clear()
    def numeric(context):
        session['calls'].append('numeric_operator')
        return session['builder'](context)
    def reference(context):
        session['calls'].append('operator_reference');context.checkpoint()
        parent=context.parents['numeric_operator'];parent.verify_contents()
        payload=parent.payloads['operator.json']
        return {'reference.json':{'upstream_identity':parent.identity,'entity_order':payload['entity_order'],
            'columns':payload['columns'],'phase':payload['phase'],'inference_outputs':False,
            'biological_accuracy':None,'calibrated_confidence':None}}
    budget=replace(BUDGET,max_jobs=1) if name=='stop_after_operator' else BUDGET
    result=P.run(session['plan'],{} if name=='verified_replay' else {'numeric_operator':numeric,'operator_reference':reference},
        session['output']/'run',budget)
    expected={'stop_after_operator':['numeric_operator'],'resume_reference':['operator_reference'],'verified_replay':[]}
    if session['calls']!=expected[name]:raise ValueError('Unexpected builder calls: '+repr(session['calls']))
    if name=='stop_after_operator':
        assert dict(result.states)=={'numeric_operator':'succeeded','operator_reference':'pending'}
        assert result.stop_reason=='Run job attempt budget exhausted'
    else:assert set(result.states.values())=={'succeeded'} and not result.stop_reason
    measured=json.loads(result.journal.read_text())['journal']['runs'][-1]
    record={'stage':name,'states':dict(result.states),'builder_calls':list(session['calls']),
        'artifact_identities':{name:artifact.identity for name,artifact in result.artifacts.items()},'measured_run':measured}
    _write(session['output']/(name+'.json'),record);session['result']=result
    return record


def verify(session):
    """Require exact native rank, all-row ECDF and historical training-state parity."""
    result=session['result'];artifact=result.artifacts['numeric_operator'];artifact.verify_contents()
    payload=artifact.payloads['operator.json'];frame=session['builder'].features;split=session['builder'].split
    columns=payload['columns'];matrix=np.array(payload['matrix'],dtype=float)
    row_at={entity:i for i,entity in enumerate(payload['entity_order'])}
    train_at=[row_at[entity] for entity in split.entities('train')]
    native=(frame.loc[list(split.entities('train')),columns].rank(pct=True)-.5).fillna(0.).to_numpy(dtype=float)
    np.testing.assert_array_equal(matrix[train_at],native)
    independent=np.zeros_like(matrix)
    for j,column in enumerate(columns):
        training=np.array(payload['training_distributions'][column],dtype=float)
        values=frame[column].to_numpy(dtype=float,na_value=np.nan)
        for i,value in enumerate(values):
            if np.isfinite(value):
                ties=np.flatnonzero(training==value)
                rank=(ties[0]+ties[-1]+2)/2 if len(ties) else int((training<=value).sum())
                independent[i,j]=rank/len(training)-.5
    np.testing.assert_array_equal(matrix,independent)
    model=session['model'];old_columns=model['kept_columns']
    if any(column not in old_columns for column in columns):raise ValueError('Pilot adds an unreviewed historical feature')
    old_at=[old_columns.index(column) for column in columns]
    np.testing.assert_array_equal(matrix[train_at],np.array(model['training_vectors'])[:,old_at])
    for column in columns:
        np.testing.assert_array_equal(payload['training_distributions'][column],model['training_distributions'][column])
    np.testing.assert_array_equal(np.array(payload['observed']),frame[columns].notna().to_numpy())
    reference=result.artifacts['operator_reference'];reference.verify_contents()
    assert reference.payloads['reference.json']['upstream_identity']==artifact.identity
    assert A.Dependency('artifact','numeric_operator',artifact.identity) in reference.spec.dependencies
    assert payload['entity_order']==[row.entity for row in split.assignments]
    assert payload['biological_accuracy'] is None and payload['calibrated_confidence'] is None
    summary={'source_artifact_identity':SOURCE_ID,'ordered_entities_retained':len(matrix),'training_entities':len(train_at),
        'permitted_input_columns':len(frame.columns),'retained_training_observed_columns':len(columns),
        'withheld_all_missing_training_columns':payload['withheld_all_missing_training_columns'],
        'deliberate_extra_exclusions':session['scope']['deliberate_extra_exclusions'],
        'native_training_ranks':'exact','independent_all_row_ecdf':'exact','historical_shared_training_vectors':'exact',
        'historical_shared_distributions':'exact','observed_masks':'exact','typed_child_binding':'verified',
        'artifact_identities':{name:artifact.identity for name,artifact in result.artifacts.items()},
        'instruction_64_25_complete':False,'inference_outputs':False,'biological_admission':False,
        'biological_accuracy':None,'calibrated_confidence':None,
        'remaining':'Graph, representation, strategy fitting, calibration and deployment builders; full organism/target coverage'}
    _write(session['output']/'summary.json',summary)
    return summary


def main():
    """Persist an actually executed notebook, operational logs and immutable checksums."""
    from notebook_runner import ExecutedNotebook
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();started=time.monotonic()
    nb=ExecutedNotebook('Tg EC V4 bounded numeric precomputation: train ranks, frozen ECDF and exact resume')
    nb.md('Freeze the existing verified FN-EC-01 V4 population and historical numeric input columns. '
        'Recompute current source exclusion closure on the 679 training rows only; preserve all historical exclusions. '
        'This intermediate feature operator fits numeric ranks only. No classifier, voting, new data or biological admission. '
        'External MemoryMax400M; two job attempts, 120 seconds per run, 32 MiB storage/package, bounded 1226x400 matrices.')
    try:
        nb.code('from scripts.freeze_precompute_features import prepare, stage, verify',
            f'session=prepare({str(args.out.resolve())!r})',"session['scope']")
        nb.md('Stop after the first reusable operator artifact. Preserve the pending child and the explicit stop reason.')
        nb.code("first=stage(session,'stop_after_operator')",'first')
        nb.md('Resume only the dependent reference. Its dependency must bind the complete verified operator identity.')
        nb.code("second=stage(session,'resume_reference')",'second')
        nb.md('Replay with no builders. Compare training ranks, all-row ECDF, shared old-model state and observed masks exactly.')
        nb.code("third=stage(session,'verified_replay')",'third')
        nb.code('verification=verify(session)','verification')
        nb.write(str(args.out/'pilot.ipynb'))
        _write(args.out/'runtime.json',{'total_seconds':time.monotonic()-started,'process_peak_rss_bytes':P._peak_rss(),
            'process_rss_bytes':P._rss(),'external_memory_cap_bytes':BUDGET.memory_limit_bytes,
            'measurement_scope':'Executed Python process lifetime; not individual isolated-job peaks'})
        files={str(path.relative_to(args.out)):{'sha256':_sha(path),'bytes':path.stat().st_size}
            for path in args.out.rglob('*') if path.is_file()}
        evidence_bytes=sum(record['bytes'] for record in files.values())
        _write(args.out/'manifest.json',{'files':files,'evidence_bytes_before_manifest':evidence_bytes,
            'package_budget_bytes':BUDGET.max_package_bytes,'source_artifact_identity':SOURCE_ID})
        total=sum(path.stat().st_size for path in args.out.rglob('*') if path.is_file())
        if total>BUDGET.max_package_bytes:raise ValueError('Entire evidence directory exceeds package budget')
    except Exception as error:
        nb.md('Diagnostic: '+type(error).__name__+': '+str(error))
        nb.write(str(args.out/'diagnostic.ipynb'))
        _write(args.out/'failure.json',{'type':type(error).__name__,'reason':str(error),'seconds':time.monotonic()-started})
        raise
    print(json.dumps(nb.ns['verification'],indent=2));print('Evidence bytes including manifest:',total)


if __name__=='__main__':main()
