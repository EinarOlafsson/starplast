#!/usr/bin/env python3
"""Export benchmark metrics, reliability bins and paired family-bootstrap intervals.

Public artifacts contain aggregate results and reproducible run contracts. They do
not turn unknown-gene hypotheses into an annotation table.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import accuracy_score, f1_score, mean_squared_error


def paired_interval(first, second, classification, repetitions=500):
    """Bootstrap whole groups jointly, conditional on the already fitted OOF models."""
    left=first[first.role.eq('held_out_evaluation')].set_index('gene_id').sort_index()
    right=second[second.role.eq('held_out_evaluation')].set_index('gene_id').sort_index()
    if not left.index.equals(right.index) or not left.group.equals(right.group) or not left.fold.equals(right.fold):
        raise ValueError('paired comparison requires identical evaluation genes, groups and folds')
    if not left.truth.equals(right.truth): raise ValueError('paired comparison requires identical truth')
    groups=[np.flatnonzero(left.group.to_numpy()==g) for g in sorted(left.group.unique())]
    truth=left.truth.astype(str).to_numpy() if classification else left.truth.to_numpy(dtype=float)
    a=left.prediction.fillna('__abstain__').astype(str).to_numpy() if classification else left.prediction.to_numpy(dtype=float)
    b=right.prediction.fillna('__abstain__').astype(str).to_numpy() if classification else right.prediction.to_numpy(dtype=float)
    labels=sorted(set(truth)) if classification else None
    def score(indices, calls):
        if classification:
            return f1_score(truth[indices],calls[indices],labels=labels,average='macro',zero_division=0)
        return np.sqrt(mean_squared_error(truth[indices],calls[indices]))
    rng=np.random.default_rng(42); differences=[]
    for _ in range(repetitions):
        indices=np.concatenate([groups[i] for i in rng.integers(0,len(groups),len(groups))])
        differences.append(score(indices,b)-score(indices,a))
    low,high=np.quantile(differences,[.025,.975])
    all_rows=np.arange(len(left))
    return dict(metric='macro_f1' if classification else 'rmse',
                difference=score(all_rows,b)-score(all_rows,a),low_95=float(low),high_95=float(high),
                repetitions=repetitions,unit='held-out orthogroup',
                limitation='Conditional on fixed OOF fits; does not include training or model-selection uncertainty')


def summarize(inputs, networks, output):
    """Write aggregate benchmark tables and the exact split/recipe contracts."""
    output.mkdir(parents=True,exist_ok=True)
    summary=[]; classes=[]; reliability=[]; strata=[]; contracts={}; predictions={}
    for path in sorted(inputs.glob('*/run.json')):
        name=path.parent.name; target,method,view=name.split('__')
        record=json.loads(path.read_text()); frame=pd.read_csv(path.parent/'predictions.csv')
        contracts[name]=record; predictions[name]=frame
        summary.append(dict(task=target,method=method,features=view,**record['metrics']))
        per_class=pd.read_csv(path.parent/'per_class.csv') if record['spec']['kind']=='classification' else pd.DataFrame()
        if len(per_class): classes.append(per_class.assign(task=target,method=method,features=view))
        known=frame[frame.role.eq('held_out_evaluation')]
        if record['spec']['kind']=='classification':
            labels=record['provenance']['classes']
            probabilities=known[[f'probability::{label}' for label in labels]].to_numpy()
            confidence=probabilities.max(axis=1)
            correct=np.asarray(labels)[probabilities.argmax(axis=1)]==known.truth.astype(str).to_numpy()
            for i in range(10):
                mask=(confidence>=i/10)&(confidence<((i+1)/10) if i<9 else confidence<=1)
                if mask.any(): reliability.append(dict(task=target,method=method,features=view,
                    bin_lower=i/10,count=int(mask.sum()),mean_confidence=float(confidence[mask].mean()),
                    accuracy=float(correct[mask].mean())))
    for name,frame in predictions.items():
        target,method,view=name.split('__')
        if target!='compartment' or view=='shuffled_labels': continue
        baseline=predictions[f'{target}__linear__published'].set_index('gene_id')
        known=frame[frame.role.eq('held_out_evaluation')].copy()
        coverage=baseline.loc[known.gene_id,'feature_coverage'].to_numpy()
        # Define sparsity using the same published-feature coverage across methods;
        # ESM's dense dimensions must not move a sparse gene into the dense group.
        low,high=np.quantile(coverage,[1/3,2/3])
        for label,mask in [('lower third',coverage<=low),('middle third',(coverage>low)&(coverage<=high)),
                           ('upper third',coverage>high)]:
            if mask.any():
                truth=known.truth.astype(str).to_numpy()[mask]
                calls=known.prediction.fillna('__abstain__').astype(str).to_numpy()[mask]
                strata.append(dict(task=target,method=method,features=view,evidence_coverage=label,
                                   genes=int(mask.sum()),accuracy=float(accuracy_score(truth,calls)),
                                   supported_fraction=float(known.supported.to_numpy()[mask].mean())))
    differences=[]
    for target in sorted({name.split('__')[0] for name in predictions}):
        baseline=f'{target}__boosted__published'
        for view in ('published_af3','published_af3_esm'):
            candidate=f'{target}__boosted__{view}'
            if baseline in predictions and candidate in predictions:
                differences.append(dict(task=target,baseline=baseline,comparison=candidate,
                                        **paired_interval(predictions[baseline],predictions[candidate],target=='compartment')))
    pd.DataFrame(summary).to_csv(output/'summary.csv',index=False)
    (pd.concat(classes,ignore_index=True) if classes else pd.DataFrame()).to_csv(output/'per_class.csv',index=False)
    pd.DataFrame(reliability).to_csv(output/'reliability_bins.csv',index=False)
    pd.DataFrame(strata).to_csv(output/'coverage_strata.csv',index=False)
    pd.DataFrame(differences).to_csv(output/'paired_group_bootstrap.csv',index=False)
    network_rows=[]
    for path in sorted(networks.glob('*/run.json')):
        record=json.loads(path.read_text()); name=path.parent.name
        contracts['network__'+name]=record
        network_rows.append(dict(layers=name,**record['metrics']))
    pd.DataFrame(network_rows).to_csv(output/'network_summary.csv',index=False)
    with (output/'run_contracts.json.gz').open('wb') as raw:
        with gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as archive:
            archive.write(json.dumps(contracts,sort_keys=True).encode())
    root=Path(__file__).resolve().parents[1]
    hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in
            [root/'starplast'/f for f in ('prediction.py','network_prediction.py','data/nodes.parquet',
                                         'data/esm_features.parquet','data/af3_features.parquet','data/graph.npz')]}
    manifest={'data_and_implementation_sha256':hashes,'run_count':len(contracts),
              'evaluation':'Development cross-validation; no prospective claim',
              'bootstrap':'500 paired orthogroup resamples, seed 42, percentile 95% intervals, fixed OOF fits',
              'artifacts':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in output.iterdir()
                           if p.is_file() and p.name!='manifest.json'}}
    (output/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(pd.DataFrame(differences).to_string(index=False),flush=True)


def main():
    """Summarize completed feature and network benchmarks without rerunning models."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--networks',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    summarize(args.input,args.networks,args.output)


if __name__=='__main__': main()
