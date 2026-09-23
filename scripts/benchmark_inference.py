#!/usr/bin/env python3
"""Reproducible family-held-out gene prediction comparisons on bundled data.

Outputs include predictions, full split memberships, feature exclusions, metrics,
source ablations and shuffled-label controls. This is a development benchmark;
it does not turn hypotheses into experimentally validated gene functions.
"""
from pathlib import Path
import argparse
from dataclasses import asdict
import hashlib
import json
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from starplast.prediction import TaskSpec, run


def main():
    """Compare all models using the same outcomes, family groups and random seed."""
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,default=ROOT/'results/inference_043')
    p.add_argument('--sequence',type=Path)
    p.add_argument('--tasks',nargs='+',default=['compartment','protein_ibaq_log2'])
    p.add_argument('--regression-targets',nargs='*',default=['protein_ibaq_log2'],
                   help='Which selected outcome columns are continuous measurements')
    p.add_argument('--methods',nargs='+',default=['prior','linear','boosted','neighbors','pca','umap','multiview'])
    p.add_argument('--controls',action='store_true')
    args=p.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    nodes=pd.read_parquet(ROOT/'starplast/data/nodes.parquet')
    if args.sequence:
        nodes=nodes.merge(pd.read_parquet(args.sequence),on='gene_id',how='left',validate='one_to_one')
    numeric=tuple(nodes.select_dtypes(include=np.number).columns)
    baseline=tuple(c for c in numeric if not c.startswith(('esm_','af3_')))
    structured=tuple(c for c in numeric if not c.startswith('esm_'))
    summaries=[]
    for target in args.tasks:
        kind='regression' if target in args.regression_targets else 'classification'
        configurations=[(method,'published',baseline,False,20) for method in args.methods]
        configurations += [(method,'published_af3',structured,False,20) for method in ('linear','boosted')]
        if args.sequence:
            configurations += [('linear','published_af3_esm',numeric,False,20),
                               ('boosted','published_af3_esm',numeric,False,20)]
        if target=='compartment' and 'umap' in args.methods:
            configurations += [('umap','published_3d',baseline,False,3)]
        if args.controls:
            configurations += [('linear','shuffled_labels',baseline,True,20)]
        for method,view,features,permute,dimensions in configurations:
            name=f'{target}__{method}__{view}';dest=args.output/name
            started=time.monotonic()
            data=nodes.copy()
            if permute:
                from starplast.prediction import known_labels
                known=data[target].notna() if kind=='regression' else known_labels(data[target])
                values=data.loc[known,target].to_numpy(copy=True)
                np.random.default_rng(42).shuffle(values)
                data.loc[known,target]=values
            spec=TaskSpec(target,kind=kind,method=method,features=features,dimensions=dimensions)
            if (dest/'run.json').exists():
                previous=json.loads((dest/'run.json').read_text())
                digest=hashlib.sha256(pd.util.hash_pandas_object(data,index=True).values.tobytes()).hexdigest()
                if previous['provenance']['data_sha256']!=digest or previous['spec']!=json.loads(json.dumps(asdict(spec))):
                    raise ValueError(f'{name}: cached run has different inputs/settings; use a fresh output directory')
                summaries.append({'task':target,'method':method,'features':view,**previous['metrics']})
                continue
            print(f'\n{name}',flush=True)
            try:
                result=run(data,spec,
                           log=lambda message:print(message,flush=True))
                result.provenance['benchmark_feature_set']=view
                result.provenance['permuted_labels']=permute
                result.save(dest)
                summary={'task':target,'method':method,'features':view,**result.metrics,
                         'seconds':round(time.monotonic()-started,2)}
            except Exception as exc:
                summary={'task':target,'method':method,'features':view,'error':f'{type(exc).__name__}: {exc}'}
                print(summary,flush=True)
                dest.mkdir(parents=True,exist_ok=True)
                (dest/'failure.json').write_text(json.dumps(summary,indent=2)+'\n')
            summaries.append(summary)
            pd.DataFrame(summaries).to_csv(args.output/'summary.csv',index=False)
            print(summary,flush=True)
    pd.DataFrame(summaries).to_csv(args.output/'summary.csv',index=False)


if __name__=='__main__':
    main()
