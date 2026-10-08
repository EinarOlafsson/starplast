"""Execute a bounded offline census of functional labels and verified result coverage."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))


def _sha(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024*1024),b''):digest.update(chunk)
    return digest.hexdigest()


def freeze(output):
    """Verify installed sources, inventory their functional columns, and retain scoped coverage."""
    import pandas as pd
    import pyarrow.parquet as pq
    from starplast import capabilities as C,discovery_labels as D,functional_coverage as F,functional_results as R,organisms as O,strategies as S,track_record as T
    output=Path(output)
    if output.exists():raise ValueError('Use a new immutable functional census directory')
    previous=ROOT/'results/discoveries_label_coverage_2026_10_08_v4'
    manifest=json.loads((previous/'manifest.json').read_text())
    catalogues={};inputs=[]
    for organism in (O.TOXOPLASMA,O.FALCIPARUM):
        source=Path(O.nodes_path(organism));inputs.append(source)
        if _sha(source)!=manifest['input_sha256'][str(source)]:raise ValueError('Installed functional source changed after membership audit')
        labels_path=previous/organism/'labels.parquet';inputs.append(labels_path)
        relative=str(labels_path.relative_to(previous))
        if _sha(labels_path)!=manifest['output_sha256'][relative]:raise ValueError('Prior functional catalogue changed')
        old=pd.read_parquet(labels_path);old=old[old.family.eq('Function')]
        columns=['gene_id',*old.target.tolist()]
        nodes=pd.read_parquet(source,columns=columns)
        catalogue,_=D.catalogue(S.Context(nodes,graph={},organism=organism),pd.DataFrame(),pd.DataFrame())
        catalogue=catalogue[catalogue.family.eq('Function')]
        check=['target','genes','annotated_genes','unannotated_genes','classes']
        pd.testing.assert_frame_equal(catalogue[check].sort_values('target').reset_index(drop=True),
            old[check].sort_values('target').reset_index(drop=True),check_exact=True)
        catalogues[organism]=catalogue
    ledger=Path(T.shipped_path());inputs.append(ledger)
    if _sha(ledger)!=manifest['input_sha256'][str(ledger)]:raise ValueError('Historical held-out ledger changed')
    counts=Counter();scanned=0
    for batch in pq.ParquetFile(ledger).iter_batches(batch_size=65536,columns=['organism','target','strategy']):
        rows=batch.to_pydict();scanned+=len(rows['organism'])
        counts.update(zip(rows['organism'],rows['target'],rows['strategy']))
    legacy=[{'organism':organism,'target':target,'strategy':strategy,'task':'label calls','rows':number}
        for (organism,target,strategy),number in sorted(counts.items())]
    benchmarks=[];availability={}
    for organism in F.DEFAULT_ORGANISMS:
        loaded,reason=R.shipped(organism);benchmarks.extend(loaded);availability[organism]=reason
    document=F.build(catalogues,capabilities=C.catalog(),benchmarks=benchmarks,legacy_records=legacy)
    output.mkdir(parents=True)
    (output/'matrix.json').write_text(json.dumps(document,sort_keys=True,indent=2,allow_nan=False)+'\n')
    (output/'summary.json').write_text(json.dumps({**document['summary'],'legacy_rows_scanned':scanned,
        'legacy_scoped_counts':len(legacy),'packaged_result_availability':availability,
        'annotation_counts':'exact current source replay against previous full label audit',
        'biology_calibration_deployment':'all unavailable; no inferred admission'},indent=2,allow_nan=False)+'\n')
    code=[Path(__file__),Path(F.__file__),Path(R.__file__),Path(D.__file__),Path(C.__file__)]
    inputs.extend([ROOT/'starplast/data/functional_results.json',ROOT/'starplast/data/functional_domain_names.json',*code])
    (output/'code').mkdir()
    for source in code:(output/'code'/source.name).write_bytes(source.read_bytes())
    (output/'manifest.json').write_text(json.dumps({'input_sha256':{str(path):_sha(path) for path in inputs},
        'output_sha256':{str(path.relative_to(output)):_sha(path) for path in output.rglob('*') if path.is_file()}},indent=2)+'\n')
    return json.loads((output/'summary.json').read_text())


def verify(output):
    """Replay every coverage address exactly from retained source labels and current declarations."""
    from starplast import capabilities as C,functional_coverage as F,functional_results as R
    output=Path(output);manifest=json.loads((output/'manifest.json').read_text())
    for path,digest in manifest['input_sha256'].items():assert _sha(path)==digest,path
    for path,digest in manifest['output_sha256'].items():assert _sha(output/path)==digest,path
    document=json.loads((output/'matrix.json').read_text())
    # The small inventory is independent of the repeated matrix cells.
    catalogues={}
    for source in document['source_label_inventory']:
        catalogues.setdefault(source['organism'],[]).append({'organism':source['organism'],
            'target':source['source_label'],'family':'Function','genes':source['genes'],
            'annotated_genes':source['known_annotation_genes'],'unannotated_genes':source['unknown_annotation_genes'],
            'classes':source['known_classes'],'source_ids':source['source_ids']})
    benchmarks=[benchmark for organism in F.DEFAULT_ORGANISMS for benchmark in R.shipped(organism)[0]]
    legacy=[]
    assert not any(row['legacy_result']['status']!='unavailable' for row in document['rows'])
    replay=F.build(catalogues,capabilities=C.catalog(),benchmarks=benchmarks,legacy_records=legacy)
    assert replay==document
    assert all(row['independent_biological_test']['admitted'] is False and row['calibration']['calibrated_confidence'] is None
        and row['deployment']['unknown_gene_claims'] is None for row in document['rows'])
    return {'matrix_replay':'exact','addresses':len(document['rows']),'pooled_accuracy':None,
        'independent_biological_admission':False,'calibration_deployment':'unavailable'}


def main():
    """Write an executed, replayed notebook rather than an unexecuted coverage claim."""
    from notebook_runner import ExecutedNotebook
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();notebook=ExecutedNotebook('Functional source, target and strategy coverage census')
    notebook.md('Bounded offline census: exact current functional source memberships, current capability declarations, historical row-count scopes and checksum-pinned functional profile recovery. Hosts remain uninventoried unsupported gene adapters. All independent biology/calibration/deployment remain unavailable; no pooled accuracy, gene count or admission is inferred.')
    try:
        notebook.code('from scripts.freeze_functional_coverage import freeze, verify',f'summary=freeze({str(args.out)!r})','summary')
        notebook.code(f'replay=verify({str(args.out)!r})','replay')
    except Exception as exc:
        args.out.mkdir(parents=True,exist_ok=True);notebook.md('Diagnostic: '+type(exc).__name__+': '+str(exc))
        notebook.write(str(args.out/'diagnostic.ipynb'));raise
    notebook.write(str(args.out/'census.ipynb'));print(json.dumps(notebook.ns['summary'],indent=2));print(notebook.ns['replay'])


if __name__=='__main__':main()
