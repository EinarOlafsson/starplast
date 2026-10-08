"""Freeze a native functional-profile recovery pilot with explicit biological truth gaps."""
from __future__ import annotations

import argparse
from dataclasses import asdict,replace
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import numpy as np
import pandas as pd
from starplast import artifacts as A, baselines as B, capabilities as C, functional_ontology as F, label_records as L, organisms as O, record_scorecards as R, scorecard as SC, strategies as S
from starplast.query import Query
from starplast.splits import make_exclusions,make_split,read_split,write_split

REVIEW=ROOT/'results/functional_source_review_2026_10_08'
TARGET='ec_major_classes'


def _sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _json(value):
    return json.loads(json.dumps(value,allow_nan=False,default=lambda v:v.item() if isinstance(v,np.generic) else _unsupported(v)))


def _unsupported(value):raise ValueError('Unsupported result scalar: '+type(value).__name__)


def _inputs():
    manifest=json.loads((REVIEW/'manifest.json').read_text())
    for path,digest in manifest['output_sha256'].items():assert _sha(REVIEW/path)==digest,path
    # Review code is replayable from its pinned snapshot. New benchmark/card
    # functions need not pretend to be the earlier source-review implementation.
    for path,digest in manifest['input_sha256'].items():
        original=Path(path)
        snapshot=REVIEW/'code'/original.name
        assert _sha(snapshot if snapshot.exists() else original)==digest,path
    receipts=json.loads((REVIEW/'receipts.json').read_text())
    enzyme=Path(next(r['path'] for r in receipts if r['name']=='enzyme.dat' and r['status']=='available'))
    nodes=Path(O.nodes_path(O.TOXOPLASMA))
    ctx=S.Context.shipped(O.TOXOPLASMA)
    profiles,resolutions=F.ec_profiles(ctx,'ec_number',F.parse_enzyme(enzyme.read_text()))
    frozen=pd.read_parquet(REVIEW/f'{O.TOXOPLASMA}_ec_number_profiles.parquet')
    assert profiles.gene_id.tolist()==frozen.gene_id.tolist()
    assert profiles.profile.fillna('UNRESOLVED').tolist()==frozen.profile.fillna('UNRESOLVED').tolist()
    assert profiles.eligible.tolist()==frozen.eligible.tolist()
    return ctx,profiles,enzyme,nodes


def _native(features,labels,split,state):
    ids=[a.entity for a in split.assignments]
    matrix=np.zeros((len(ids),len(state['kept_columns'])))
    for j,column in enumerate(state['kept_columns']):
        ordered=np.array(state['training_distributions'][column],dtype=float)
        values=features[column].to_numpy(dtype=float)
        for i,value in enumerate(values):
            if np.isfinite(value):
                equal=np.flatnonzero(ordered==value)
                rank=(equal[0]+equal[-1]+2)/2 if len(equal) else int((ordered<=value).sum())
                matrix[i,j]=rank/len(ordered)-.5
    at={entity:i for i,entity in enumerate(ids)}
    train=np.array([at[e] for e in split.entities('train')])
    expected=(features.loc[list(labels.index),state['kept_columns']].rank(pct=True)-.5).fillna(0.).to_numpy()
    np.testing.assert_array_equal(matrix[train],expected)
    np.testing.assert_array_equal(matrix[train],state['training_vectors'])
    visible=pd.Series([None]*len(ids),dtype=object)
    visible.iloc[train]=labels.to_numpy()
    query=np.array([at[e] for e in split.entities('test')])
    native,shares=S.knn_vote(matrix,visible,state['k'],query=query)
    calls=native.where(shares>=state['min_share']).iloc[query]
    scores=S.class_scores_of(native).iloc[query]
    return calls,shares.iloc[query],scores


def freeze(output):
    """Fit one predeclared EC-profile adapter without seeing hidden truth or changing runtime data."""
    output=Path(output)
    if output.exists():raise ValueError('Use a new immutable functional pilot directory')
    ctx,profiles,enzyme,nodes=_inputs()
    frame=ctx.nodes.set_index('gene_id',drop=False).copy()
    truth=profiles.set_index('gene_id').profile
    eligible=profiles.gene_id[profiles.eligible].tolist()
    groups=ctx.groups()[profiles.eligible.to_numpy()].tolist()
    benchmark='FN-EC-01:'+O.TOXOPLASMA+':major_class_profiles'
    split=make_split(eligible,groups,organism=O.TOXOPLASMA,benchmark_id=benchmark,seed=23,feature_access='inductive')
    frame[TARGET]=truth
    training=frame.loc[list(split.entities('train'))].reset_index(drop=True)
    exclusions=make_exclusions(training,('ec_number',TARGET),benchmark_id=benchmark,split=split)
    # The gene-level assignment methods are unresolved. Conservatively withhold
    # functional annotations and homology summaries that could have generated EC.
    withheld={'ec_number','ec_number_orthology','has_ec','interpro_id','interpro_ids','interpro_desc','pfam_id','pfam_ids',
        'n_interpro','has_domain','orthogroup','paralog_number','has_pf_ortholog','has_cp_ortholog','lineage_specific'}
    exclusions=replace(exclusions,columns=tuple(sorted(set(exclusions.columns)|withheld)))
    train_ctx=S.Context(training,graph={},organism=O.TOXOPLASMA)
    columns=train_ctx.numeric_columns(exclude=exclusions.columns)
    features=frame.loc[eligible,columns]
    labels=truth.loc[list(split.entities('train'))].astype(str)
    batch=L.feature_knn(features,labels,split=split,exclusions=exclusions,k=15,min_share=.3)
    calls,support,scores=_native(features,labels,split,batch.model_state)
    assert batch.rows.prediction.fillna('ABSTAIN').tolist()==calls.fillna('ABSTAIN').tolist()
    np.testing.assert_array_equal(batch.rows.support.to_numpy(dtype=float),support.to_numpy(dtype=float))
    np.testing.assert_array_equal(batch.class_scores.to_numpy(),scores.to_numpy())
    rows=batch.rows.copy()
    rows['truth']=truth.loc[list(rows.entity)].to_numpy()
    rows['group']=rows.entity.map(dict(zip(eligible,groups)))
    gaps=tuple(sorted(set((*batch.gaps,
        'Gene-level EC curation/experimental/predicted evidence grades remain unresolved',
        'Nomenclature resolution is not a new gene activity measurement',
        'Unannotated/unresolvable profiles remain unknown; no genome-wide verified negatives',
        'Exact profile recovery preserves multiple major classes; novel combinations may be unsupported',
        'Source release and complete homology/source independence remain unresolved',
        str(sum(str(g).startswith('__') for g in groups))+' eligible genes lack resolved homology groups'))))
    settings={'k':15,'min_share':.3,'target':'complete major-class profile','selection':'fixed before testing',
        'feature_transform':'native training ranks; frozen query ECDF','ontology':'pinned ENZYME 02-Sep-2026'}
    scope=R.RecordScope(O.TOXOPLASMA,'feature_knn',TARGET,SC.T_LABEL,A.canonical_object(settings),23,split.identity,
        'outer_test',benchmark,'unresolved','gene','complete recorded annotation profiles; biological absence unknown',gaps)
    card=R.aggregate(rows,scope,parameters={'class_scores':batch.class_scores.reset_index(drop=True)})
    baseline=B.label_baselines(labels,split,seed=23)
    baseline_cards={}
    for column in ('majority','prevalence_call'):
        base_rows=rows.drop(columns=['prediction','abstained']).assign(prediction=baseline[column].to_numpy())
        baseline_cards[column]=R.aggregate(base_rows,scope)
    code=[Path(__file__),*(ROOT/'starplast'/name for name in ('functional_ontology.py','discovery_labels.py','label_records.py',
        'strategies.py','splits.py','search.py','slots.py','embedding.py','datasets.py','artifacts.py','record_scorecards.py','scorecard.py','baselines.py','capabilities.py'))]
    inputs=[nodes,enzyme,REVIEW/'manifest.json',REVIEW/f'{O.TOXOPLASMA}_ec_number_profiles.parquet',*code]
    dependencies=[A.Dependency('code',p.stem,_sha(p)) for p in code]
    dependencies.extend((A.Dependency('table','installed_nodes',_sha(nodes)),A.Dependency('source','ENZYME',_sha(enzyme)),
        A.Dependency('truth','complete_major_class_profiles',_sha(inputs[3])),A.Dependency('split','nested',split.identity),
        A.Dependency('exclusions','training_source_closure',hashlib.sha256(A.canonical_object(asdict(exclusions)).encode()).hexdigest())))
    evaluation={'unit':'gene','eligible_population':len(rows),'truth_grade':'unresolved','negative_semantics':scope.negative_semantics,
        'limitations':list(gaps),'benchmark_status':'candidate','source_ids':['toxodb_ec_numbers'],
        'source_annotated_genes':int(profiles.status.ne('unannotated').sum()),'unknown_genes':int(profiles.status.eq('unannotated').sum()),
        'unresolved_annotated_genes':int(profiles.status.eq('unresolved').sum()),'biological_admission':False}
    version=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()+'+functional-ec-pilot'
    spec=A.ArtifactSpec(Query(O.TOXOPLASMA,'label',target=TARGET).to_json(),'feature_knn',SC.T_LABEL,C.get('feature_knn').outputs[0],
        TARGET,tuple(rows.entity),'held_out','outer:test:'+split.identity,tuple(dependencies),A.canonical_object(settings),
        A.canonical_object(evaluation),23,version,fit_entities=split.entities('train'),fit_role='train',benchmark_id=benchmark,
        confidence_kind='method_support',evaluation_partition='test',gaps=gaps)
    payload=_json({'rows.json':rows.astype(object).where(rows.notna(),None).to_dict('records'),
        'row_columns.json':rows.columns.tolist(),'class_scores.json':{'classes':batch.class_scores.columns.tolist(),'scores':batch.class_scores.to_numpy().tolist()},
        'model_state.json':batch.model_state,'card.json':card,'profile_class_cards.json':R.class_cards(rows,scope),
        'major_class_cards.json':F.member_class_cards(rows),'baseline_cards.json':baseline_cards})
    output.mkdir(parents=True);write_split(output/'split.json',split,exclusions=exclusions)
    identity=A.write_artifact(output/'held_out',spec,payload,split=split)
    assert A.read_artifact(output/'held_out',expected=spec,split=split).payloads==payload
    summary={'organism':O.TOXOPLASMA,'target':TARGET,'eligible_profiles':len(eligible),'train':len(labels),
        'tune':len(split.entities('tune')),'calibration':len(split.entities('calibration')),'test':len(rows),
        'feature_columns':len(columns),'profile_classes':len(batch.class_scores.columns),
        'multi_major_class_test_genes':int(rows.truth.map(lambda v:len(F.profile_members(v))).gt(1).sum()),
        'metrics':card['metrics'],'counts':card['counts'],'majority_baseline_metrics':baseline_cards['majority']['metrics'],
        'major_class_cards':payload['major_class_cards.json'],'unsupported_test_profiles':rows[~rows.truth.isin(batch.class_scores.columns)].truth.value_counts().to_dict(),
        'artifact_identity':identity,'biological_benchmark_admitted':False,'confidence':'uncalibrated native support',
        'interpretation':'Functional annotation-profile recovery; independent biological accuracy unknown'}
    (output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    (output/'code').mkdir()
    for path in code:(output/'code'/path.name).write_bytes(path.read_bytes())
    (output/'manifest.json').write_text(json.dumps({'input_sha256':{str(p):_sha(p) for p in inputs},
        'output_sha256':{str(p.relative_to(output)):_sha(p) for p in output.rglob('*') if p.is_file()}},indent=2)+'\n')
    return summary


def verify(output):
    """Replay every native prediction/score and recompute matched profile/member cards."""
    output=Path(output);manifest=json.loads((output/'manifest.json').read_text())
    for path,digest in manifest['input_sha256'].items():assert _sha(path)==digest,path
    for path,digest in manifest['output_sha256'].items():assert _sha(output/path)==digest,path
    split,exclusions=read_split(output/'split.json')
    spec=A.spec_from_dict(json.loads((output/'held_out/manifest.json').read_text())['spec'])
    artifact=A.read_artifact(output/'held_out',expected=spec,split=split)
    payload=artifact.payloads
    ctx,profiles,enzyme,nodes=_inputs()
    truth=profiles.set_index('gene_id').profile
    state=payload['model_state.json'];ids=[a.entity for a in split.assignments]
    features=ctx.nodes.set_index('gene_id').loc[ids,state['input_columns']]
    labels=truth.loc[list(split.entities('train'))].astype(str)
    batch=L.feature_knn(features,labels,split=split,exclusions=exclusions,k=15,min_share=.3)
    assert _json(batch.model_state)==state
    rows=pd.DataFrame(payload['rows.json'],columns=payload['row_columns.json'])
    calls,support,scores=_native(features,labels,split,state)
    assert rows.prediction.fillna('ABSTAIN').tolist()==calls.fillna('ABSTAIN').tolist()
    np.testing.assert_array_equal(rows.support.to_numpy(dtype=float),support.to_numpy(dtype=float))
    np.testing.assert_array_equal(payload['class_scores.json']['scores'],scores.to_numpy())
    assert payload['class_scores.json']['classes']==scores.columns.tolist()
    assert rows.truth.tolist()==truth.loc[list(rows.entity)].tolist()
    scope_data=dict(payload['card.json']['scope']);scope_data['gaps']=tuple(scope_data['gaps'])
    scope=R.RecordScope(**scope_data)
    assert _json(R.aggregate(rows,scope,parameters={'class_scores':batch.class_scores.reset_index(drop=True)}))==payload['card.json']
    assert _json(R.class_cards(rows,scope))==payload['profile_class_cards.json']
    assert F.member_class_cards(rows)==payload['major_class_cards.json']
    baseline=B.label_baselines(labels,split,seed=23)
    for column in ('majority','prevalence_call'):
        base_rows=rows.drop(columns=['prediction','abstained']).assign(prediction=baseline[column].to_numpy())
        assert _json(R.aggregate(base_rows,scope))==payload['baseline_cards.json'][column]
    return {'native_predictions_and_scores':'exact','profile_and_member_cards':'exact','all_test_genes_retained':len(rows),
        'source_code_hashes':'verified','biological_admission':False}


def main():
    """Execute fitting and native replay; preserve any failed wrapper as a diagnostic notebook."""
    from notebook_runner import ExecutedNotebook
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();nb=ExecutedNotebook('Functional EC profiles: frozen native kNN recovery and capacity')
    nb.md('FN-EC-01, fixed seed 23/k=15/min-share=.3, protected homology groups, complete major-class profiles. All assigned EC entries must uniquely resolve into active pinned nomenclature; overlapping major classes stay in the full profile. Unknown or unresolved genes never become negatives. Train-only feature selection/ranks/imputation; conservatively exclude enzyme/domain/homology annotations. Native profile calls and class scores remain uncalibrated method support. Candidate reference recovery is not independent measured gene activity.')
    try:
        nb.code('from scripts.freeze_functional_ec_pilot import freeze, verify',f'summary=freeze({str(args.out)!r})','summary')
        nb.md('Independent rank construction and native vote replay, exact float scores, retained abstentions/unseen profile combinations, profile/member class cards and matched train-only baselines. No original data/calibration changes and no new verified unknown-gene claims.')
        nb.code(f'verification=verify({str(args.out)!r})','verification')
    except Exception as exc:
        nb.md('Diagnostic: '+type(exc).__name__+': '+str(exc));nb.write(str(args.out/'diagnostic.ipynb'));raise
    nb.write(str(args.out/'pilot.ipynb'));print(json.dumps(nb.ns['summary'],indent=2));print(nb.ns['verification'])


if __name__=='__main__':main()
