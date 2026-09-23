#!/usr/bin/env python3
"""Compare fixed, sourced network layers for held-out localization labels."""
import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def main():
    """Use identical family splits for single layers and a learned combination."""
    from starplast.prediction import TaskSpec
    from starplast.network_prediction import run_network,bundled_network
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/network_043')
    args=parser.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    nodes=pd.read_parquet(ROOT/'starplast/data/nodes.parquet')
    graph,sources=bundled_network(nodes)
    rows=[]
    for layers in (['domain'],['coexpression'],['struct'],list(sources)):
        name='+'.join(layers)
        result=run_network(nodes,TaskSpec('compartment'),graph,layers,sources,
                           log=lambda message:print(f'{name}: {message}',flush=True))
        result.save(args.output/name)
        rows.append({'layers':name,**result.metrics})
        print(rows[-1],flush=True)
        pd.DataFrame(rows).to_csv(args.output/'summary.csv',index=False)


if __name__=='__main__':
    main()
