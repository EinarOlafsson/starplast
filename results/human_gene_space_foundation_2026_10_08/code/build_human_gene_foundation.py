"""Execute the verified human source-gene builder and strict source-text replay."""
from __future__ import annotations

import argparse
import csv
from dataclasses import asdict
import gzip
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from starplast import organisms as O, provenance as P  # noqa: E402
from starplast.spaces import human as H, validate_nodes  # noqa: E402


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):digest.update(block)
    return digest.hexdigest()


def _inputs():
    review = ROOT/'results/host_expression_source_review_2026_10_08/comparison_v3/manifest.json'
    bindings = ROOT/'results/source_recovery_final_2026_10_07_v2/bindings.json'
    source = json.loads(bindings.read_text())['host_gtex_transcriptome']
    previous = json.loads(review.read_text())['input_sha256']
    inputs = {'gtex':Path(source['path'])}
    for key,name in (('idmap','HUMAN_9606_idmapping.dat.gz'),('reviewed','reviewed_human.tsv.gz')):
        inputs[key] = Path(next(path for path in previous if path.endswith('/'+name)))
    for key,path in inputs.items():
        expected = source['sha256'] if key=='gtex' else previous[str(path)]
        if _sha(path)!=expected:raise ValueError('Verified source snapshot is stale: '+str(path))
    return inputs,review,bindings,source


def build(output):
    """Retain canonical and qualified gene evidence without promoting an installed pack."""
    output = Path(output)
    if output.exists():raise ValueError('Use a new immutable human gene foundation directory')
    inputs,review,bindings,receipt = _inputs()
    records = H.read_gtex_gene_records(inputs['gtex'])
    links = H.reviewed_gene_protein_links(inputs['idmap'],inputs['reviewed'])
    batch = H.build_gene_foundation(records,links)
    assert len(records)==59033 and len(batch.nodes)==58988
    assert batch.mapping_summary['qualified_PAR_Y_records']==45
    code = [Path(__file__),ROOT/'starplast/spaces/human.py',ROOT/'starplast/spaces/__init__.py',
        ROOT/'starplast/provenance.py',ROOT/'starplast/organisms.py']
    source_files = tuple(P.SourceFile(str(path),_sha(path),path.stat().st_size,
        'processed_input' if key=='gtex' else 'mapping_reference') for key,path in inputs.items())
    transform = P.Transform('GTEx source-gene identity and summary preservation','starplast.spaces.human',_sha(code[1]),
        ('host_gtex_transcriptome:Name',*H.FIELDS),tuple(batch.nodes.columns),
        {'reference':'GENCODE39/GRCh38','qualified_PAR_Y':'retain separately, no ordinary-gene overwrite',
         'source_float_conversion':'Python binary64 round-trip','protein_measurement_projection':False})
    gaps = ('Independent biological benchmark admission remains unresolved',
        'GTEx v10 release-to-publication and sample-summary/context lineage review remains pending',
        'Complete pack redistribution/reference/graph/opening gates remain pending',
        'Cultured adult fibroblasts are not validated HFF measurements',
        'Gene annotation scope is this source, not genome completeness')
    traces = [P.MeasurementTrace('host_gtex_transcriptome',O.HUMAN,O.HUMAN,'gene',column,
        evidence_grade='direct_experiment',quantity_unit='median transcripts per million (TPM)',
        context=(source,),source_files=(source_files[0],),transforms=(transform,),
        publication_identity='PMID:32913098',publication_status='registered; release lineage unresolved',
        license='GTEx open-access acknowledgement terms; specific redistribution review unresolved',
        redistribution='unresolved',gaps=gaps,source_species=('Homo sapiens',),
        source_species_status='source_file_verified',source_species_reference=receipt['url'])
        for source,column in H.FIELDS.items()]
    mapping_audit = P.audit_mapping(batch.nodes.gene_id,links[['gene_id','protein_id']].itertuples(index=False,name=None),
        source_organism=O.HUMAN,target_organism=O.HUMAN,source_namespace='Ensembl gene',target_namespace='reviewed UniProt protein',
        source_version='GTEx v10 GENCODE39',target_version='UniProt 2026_03',reference=_sha(inputs['idmap']))
    # Mapping audits describe whether a hypothetical unique projection exists;
    # every explicit candidate edge remains retained and nothing is projected.
    column_metadata = {}
    for column in batch.nodes:
        numeric = column in H.FIELDS.values()
        column_metadata[column] = {'organism':O.HUMAN,'storage_unit':'gene',
            'quantity_unit':'median TPM' if numeric else 'identifier or source-provided name',
            'evidence_grade':'direct_experiment' if numeric else 'curation',
            'source_id':'host_gtex_transcriptome','source_sha256':_sha(inputs['gtex']),
            'source_field':next((k for k,v in H.FIELDS.items() if v==column),'Name' if column=='gene_id' else 'Description'),
            'context':next((k for k,v in H.FIELDS.items() if v==column),'GENCODE39 source-gene annotation'),
            'license':'GTEx public terms; pack redistribution review pending','redistribution':'unresolved',
            'zero_semantics':'observed zero median; not missing' if numeric else 'not a measurement',
            'missingness':'source missing retained; unmapped protein does not remove gene',
            'targets':'candidate numeric benchmark only, no independent accuracy admission' if numeric else 'identity only',
            'name_status':'source-provided; not asserted HGNC validation' if column=='gene_name' else 'not applicable'}
    output.mkdir(parents=True)
    batch.nodes.to_parquet(output/'hs_gene_nodes.parquet',index=False)
    batch.source_records.to_parquet(output/'gtex_source_gene_records.parquet',index=False)
    batch.protein_links.to_parquet(output/'gene_protein_crossrefs.parquet',index=False)
    batch.gene_mapping.to_parquet(output/'gene_mapping_status.parquet',index=False)
    P.write_traces(traces,output/'measurement_provenance.json')
    (output/'mapping_audit.json').write_text(json.dumps(asdict(mapping_audit),indent=2)+'\n')
    (output/'column_metadata.json').write_text(json.dumps(column_metadata,indent=2)+'\n')
    summary = dict(batch.mapping_summary,benchmark_admitted=False,installed_space_registered=False,
        distributable_pack_built=False,redistribution='unresolved',source_sha256=_sha(inputs['gtex']),gaps=list(gaps))
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    (output/'code').mkdir()
    for path in code:(output/'code'/path.name).write_bytes(path.read_bytes())
    installed = [ROOT/'starplast/data'/name for name in (*O.HOST_TABLES.values(),'nodes.parquet','pf_nodes.parquet',
        'graph.npz','pf_graph.npz','gene_track_record.parquet','pf_gene_track_record.parquet') if (ROOT/'starplast/data'/name).exists()]
    all_inputs = [*inputs.values(),review,bindings,*code,*installed]
    (output/'manifest.json').write_text(json.dumps({'input_sha256':{str(p):_sha(p) for p in all_inputs},
        'outputs':{str(p.relative_to(output)):_sha(p) for p in output.rglob('*') if p.is_file()},
        'source_inputs':{k:str(v) for k,v in inputs.items()},'unchanged_installed_files':[str(p) for p in installed]},indent=2)+'\n')
    return summary


def verify(output):
    """Replay exact source values, gene preservation, links and every artifact hash."""
    output = Path(output)
    manifest = json.loads((output/'manifest.json').read_text())
    for path,digest in manifest['input_sha256'].items():assert _sha(path)==digest,path
    for path,digest in manifest['outputs'].items():assert _sha(output/path)==digest,path
    inputs = {k:Path(v) for k,v in manifest['source_inputs'].items()}
    batch = H.build_gene_foundation(H.read_gtex_gene_records(inputs['gtex']),
        H.reviewed_gene_protein_links(inputs['idmap'],inputs['reviewed']))
    for name,frame in (('hs_gene_nodes.parquet',batch.nodes),('gtex_source_gene_records.parquet',batch.source_records),
        ('gene_protein_crossrefs.parquet',batch.protein_links),('gene_mapping_status.parquet',batch.gene_mapping)):
        pd.testing.assert_frame_equal(pd.read_parquet(output/name),frame,check_exact=True)
    validate_nodes(H.HUMAN_SPACE,batch.nodes)
    # Independent line-level source replay does not call the builder's reader.
    records = batch.source_records.set_index('source_row')
    checked = 0
    with gzip.open(inputs['gtex'],'rt',newline='') as stream:
        stream.readline();stream.readline()
        reader = csv.DictReader(stream,delimiter='\t')
        for line,row in enumerate(reader,4):
            retained = records.loc[line]
            assert row['Name']==retained.source_gene_id
            assert row['Description']==retained.source_gene_name
            for field,column in H.FIELDS.items():
                original = float(row[field])
                assert original==retained[column]
                checked+=1
    ordinary = batch.source_records[batch.source_records.source_qualifier.eq('ordinary')].set_index('gene_id')
    np.testing.assert_array_equal(batch.nodes.gene_id,ordinary.index)
    for column in H.FIELDS.values():np.testing.assert_array_equal(batch.nodes[column],ordinary[column])
    summary = json.loads((output/'summary.json').read_text())
    assert all(summary[k]==v for k,v in batch.mapping_summary.items())
    for trace in P.read_traces(output/'measurement_provenance.json'):
        assert trace.output_organism==O.HUMAN and trace.storage_unit=='gene'
        assert trace.quantity_unit=='median transcripts per million (TPM)'
        for source in trace.source_files:source.verify()
    # Preserve all reviewed mapping candidates; ambiguous joins are reported by
    # the independent common MappingAudit contract, not chosen by row order.
    audit = P.audit_mapping(batch.nodes.gene_id,batch.protein_links[['gene_id','protein_id']].itertuples(index=False,name=None),
        source_organism=O.HUMAN,target_organism=O.HUMAN,source_namespace='Ensembl gene',target_namespace='reviewed UniProt protein',
        source_version='GTEx v10 GENCODE39',target_version='UniProt 2026_03',reference=_sha(inputs['idmap']))
    assert asdict(audit)==json.loads((output/'mapping_audit.json').read_text())
    return {'source_measurement_cells_exact':checked,'canonical_genes_exact':len(batch.nodes),'qualified_Y_preserved':45,
        'crossrefs_and_reverse_cardinality':'exact replay','all_installed_inputs':'hash unchanged',
        'graph_and_pack_gates':'pending, not fabricated','biological_admission':False}


def main():
    """Write the actual executed source build and independent strict replay notebook."""
    from notebook_runner import ExecutedNotebook
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args = parser.parse_args()
    nb = ExecutedNotebook('Human canonical source-gene foundation and explicit protein mappings')
    nb.md('64.21 first partition: verified GTEx v10 processed medians, selected cultured-fibroblast and LCM-hepatocyte contexts, existing UniProt 2026_03 reviewed cross-references. All source gene rows and zero/missingness semantics remain intact; PAR_Y copies cannot overwrite ordinary-gene values. This is a source-scoped gene foundation, not a promoted full genome or distributable pack. Protein measurements remain in their own reference tables. Source/unit/context, license/redistribution, graph/opening and independent benchmark gates stay explicit.')
    nb.code('from scripts.build_human_gene_foundation import build, verify',f'summary=build({str(args.out)!r})','summary')
    nb.md('Replay every source text value and source/derived artifact hash; exact checks use binary64 identity, never default dataframe tolerance. Verify one-to-many/reverse mapping census and all qualified source records. Numeric source positions/versions and mapping counts stay outside the biological measurement table.')
    nb.code(f'verification=verify({str(args.out)!r})','verification')
    nb.write(str(args.out/'build.ipynb'))
    print(nb.ns['summary']);print(nb.ns['verification'])


if __name__=='__main__':main()
