"""Acquire pinned functional nomenclature and audit installed source semantics."""
from __future__ import annotations

import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import urllib.request

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
import pandas as pd
from starplast import datasets as D, functional_ontology as F, organisms as O, strategies as S

DOWNLOADS = {
    'enzyme.dat':'https://ftp.expasy.org/databases/enzyme/enzyme.dat',
    'enzclass.txt':'https://ftp.expasy.org/databases/enzyme/enzclass.txt',
    'enzuser.txt':'https://ftp.expasy.org/databases/enzyme/enzuser.txt',
    'interpro_index.html':'https://ftp.ebi.ac.uk/pub/databases/interpro/current_release/',
    'Pfam-A.hmm.dat.gz':'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz',
    'Pfam.version.gz':'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam.version.gz',
    'pfam_relnotes.txt':'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/relnotes.txt',
    'pfam_userman.txt':'https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/userman.txt',
}


def _sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def fetch(name,url,archive):
    """Archive one bounded public metadata response with receipts; no access workaround."""
    path=Path(archive)/name
    if path.exists():raise ValueError('Do not overwrite an existing source acquisition')
    receipt={'requested_url':url,'retrieved_utc':datetime.now(timezone.utc).isoformat(),'name':name}
    try:
        request=urllib.request.Request(url,headers={'User-Agent':'Starplast-functional-source-review/0.54'})
        with urllib.request.urlopen(request,timeout=30) as response:
            payload=response.read(20_000_001)
            if len(payload)>20_000_000:raise ValueError('Metadata exceeds frozen 20 MB limit')
            if name!='interpro_index.html' and b'<html' in payload[:2000].lower():raise ValueError('HTML is not ontology data')
            path.write_bytes(payload)
            receipt.update(status='available',path=str(path.resolve()),sha256=_sha(path),bytes=len(payload),
                final_url=response.url,http_status=response.status,headers=dict(response.headers))
    except Exception as exc:receipt.update(status='unavailable',error=type(exc).__name__+': '+str(exc))
    return receipt


def review(output,archive,receipts):
    """Audit current EC resolution and explicit unresolved annotation/assay lineage."""
    output,archive=Path(output),Path(archive)
    entries=F.parse_enzyme((archive/'enzyme.dat').read_text())
    header=(archive/'enzyme.dat').read_text().split('ID   ',1)[0]
    summaries=[]
    for organism in (O.TOXOPLASMA,O.FALCIPARUM):
        ctx=S.Context.shipped(organism)
        for target in ('ec_number','ec_number_orthology'):
            if target not in ctx.nodes:continue
            profiles,resolutions=F.ec_profiles(ctx,target,entries)
            profiles.to_parquet(output/f'{organism}_{target}_profiles.parquet',index=False)
            resolutions.to_parquet(output/f'{organism}_{target}_resolution.parquet',index=False)
            owners=[d for d in D.REGISTRY if d.organism==organism and target in d.columns]
            summaries.append({'organism':organism,'target':target,'genes':ctx.n,
                'source_annotated_genes':int(ctx.truth(target).notna().sum()),'eligible_complete_profiles':int(profiles.eligible.sum()),
                'profile_status':profiles.status.value_counts().to_dict(),'resolution_status':resolutions.status.value_counts().to_dict(),
                'multi_major_class_genes':int(profiles.major_classes.map(len).gt(1).sum()),
                'major_class_profiles':profiles.profile.dropna().value_counts().to_dict(),
                'source_ids':[d.key for d in owners],
                'declared_source_grade':'orthology_transfer' if target=='ec_number_orthology' else 'registry annotation; individual evidence grades unresolved',
                'biological_benchmark_admitted':False,'negative_semantics':'unannotated_is_unknown; recovery of complete recorded profiles only',
                'gaps':['Gene-level curation/experimental/predicted evidence not independently verified',
                    'Original source annotation release and homology/source independence unresolved',
                    'ENZYME nomenclature defines reactions; it is not new gene activity evidence']})
    source_rows=[]
    for d in D.REGISTRY:
        if d.key in {'interpro','plasmodb_pf3d7_domains','toxodb_ec_numbers','pf_enzyme_classification'}:
            source_rows.append({'key':d.key,'organism':d.organism,'columns':d.columns,'url':d.url,'path':d.path,'note':d.note,
                'independent_gene_activity_truth':False,'source_grade':'orthology_transfer for explicitly derived EC; all other individual grades unresolved'})
    (output/'source_review.json').write_text(json.dumps(source_rows,indent=2)+'\n')
    (output/'summary.json').write_text(json.dumps({'ec_release_header':header,'ec_entry_counts':pd.Series([e.status for e in entries.values()]).value_counts().to_dict(),
        'sources':summaries,'domain_note':'Gene domain matches are annotations/predictions, not an independently measured activity assay',
        'ontology_selection':'Official current nomenclature metadata; no gene dataset replacement or literature-best claim'},indent=2)+'\n')
    inputs=[Path(O.nodes_path(o)) for o in (O.TOXOPLASMA,O.FALCIPARUM)]
    inputs += [Path(r['path']) for r in receipts if r['status']=='available']
    code=[Path(__file__),ROOT/'starplast/functional_ontology.py',ROOT/'starplast/discovery_labels.py',ROOT/'starplast/datasets.py']
    inputs+=code
    (output/'code').mkdir()
    for path in code:(output/'code'/path.name).write_bytes(path.read_bytes())
    (output/'manifest.json').write_text(json.dumps({'input_sha256':{str(p):_sha(p) for p in inputs},
        'output_sha256':{str(p.relative_to(output)):_sha(p) for p in output.rglob('*') if p.is_file()},
        'code_version':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()+'+functional-source-review'},indent=2)+'\n')
    return summaries


def main():
    """Execute source acquisition and profile review as an annotated notebook."""
    from notebook_runner import ExecutedNotebook
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--archive',type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists() or args.archive.exists():raise ValueError('Use new immutable review and archive directories')
    args.out.mkdir(parents=True);args.archive.mkdir(parents=True)
    nb=ExecutedNotebook('Functional sources: nomenclature, obsolete entries and truth gaps')
    nb.md('FN-EC-01 source partition. Acquire only public ontology metadata, not bulk alignments/models. Original gene data stays unchanged. Resolve only a unique transfer chain into an active EC; multiple alternatives do not become extra gene assignments. Preserve complete multi-valued major-class profiles; unresolved/unannotated genes stay outside annotation-recovery truth. Registry annotation is not independent measured activity. Sources/assays, ontology release and individual evidence grades remain explicit.')
    nb.code('from scripts.review_functional_sources import fetch, review','receipts=[]')
    try:
        for name,url in DOWNLOADS.items():
            nb.code(f'receipt=fetch({name!r},{url!r},{str(args.archive)!r})','receipts.append(receipt)','receipt')
            print(name,nb.ns['receipt']['status'],flush=True)
        receipts=nb.ns['receipts']
        (args.out/'receipts.json').write_text(json.dumps(receipts,indent=2)+'\n')
        nb.code(f'summary=review({str(args.out)!r},{str(args.archive)!r},receipts)','summary')
    except Exception as exc:
        nb.md('Diagnostic: '+type(exc).__name__+': '+str(exc));nb.write(str(args.out/'diagnostic.ipynb'));raise
    nb.write(str(args.out/'review.ipynb'))
    print(json.dumps(nb.ns['summary'],indent=2),flush=True)


if __name__=='__main__':main()
