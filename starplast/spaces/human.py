"""Build human gene records from measured source genes, keeping protein links separate.

This foundation does not relabel the installed host protein table or activate a
gene pack. GENCODE pseudoautosomal Y records retain their qualified source identity
and cannot overwrite the corresponding ordinary Ensembl gene's measurements.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
import gzip
import re

import numpy as np
import pandas as pd

from .. import organisms as O
from . import validate_nodes

GENE = re.compile(r'ENSG\d{11}')
SOURCE_GENE = re.compile(r'(ENSG\d{11})\.(\d+)(_PAR_Y)?')
PROTEIN = re.compile(r'(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})')
FIELDS = {'Cells_Cultured_fibroblasts':'gtex_fibroblast_median_tpm',
          'Liver_Hepatocyte':'gtex_lcm_hepatocyte_median_tpm'}
HUMAN_SPACE = O.Space(O.HUMAN,'Homo sapiens','GTEx v10 / GENCODE 39 / GRCh38',O.HOST,
    r'ENSG\d{11}',('ENSG',),'Ensembl','https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={id}',
    'hs_gene_nodes.parquet',graph='hs_gene_graph.npz',contexts=frozenset(FIELDS),
    numbers=tuple(FIELDS.values()),distribution='pack')


@dataclass(frozen=True)
class HumanGeneBuild:
    """Source records, canonical gene rows and explicit many-to-many protein cross-references."""

    nodes: pd.DataFrame
    source_records: pd.DataFrame
    protein_links: pd.DataFrame
    gene_mapping: pd.DataFrame
    mapping_summary: dict


def read_gtex_gene_records(path) -> pd.DataFrame:
    """Read both selected TPM summaries exactly, checking IDs and declared GCT dimensions.

    Python float conversion rounds the original text to binary64 without pandas'
    fast-parser rounding. Genuine zero is retained; explicit missing tokens remain
    NaN. Source rows are one-based physical line numbers, including GCT headers.
    """
    rows = []
    with gzip.open(path,'rt',newline='') as stream:
        if stream.readline().strip()!='#1.2':
            raise ValueError('Expected GTEx GCT format 1.2')
        dimensions = stream.readline().strip().split('\t')
        if len(dimensions)!=2 or not all(v.isdigit() for v in dimensions):
            raise ValueError('Invalid declared GCT dimensions')
        n,m = map(int,dimensions)
        reader = csv.reader(stream,delimiter='\t')
        header = next(reader,None)
        if (not header or header[:2]!=['Name','Description'] or len(header)!=m+2
                or len(header)!=len(set(header)) or not set(FIELDS)<=set(header)):
            raise ValueError('Source columns differ from the declared GTEx context contract')
        positions = {source:header.index(source) for source in FIELDS}
        for line,values in enumerate(reader,4):
            if len(values)!=len(header):
                raise ValueError('Incomplete or excessive GCT source fields at line '+str(line))
            match = SOURCE_GENE.fullmatch(values[0])
            if match is None:
                raise ValueError('Noncanonical human source gene at line '+str(line))
            row = {'source_row':line,'source_gene_id':values[0],'gene_id':match.group(1),
                'annotation_version':int(match.group(2)),'source_qualifier':'PAR_Y' if match.group(3) else 'ordinary',
                'source_gene_name':values[1] or None}
            for source,column in FIELDS.items():
                text = values[positions[source]]
                try:
                    value = np.nan if text in ('','NA','NaN') else float(text)
                except ValueError as error:
                    raise ValueError('Invalid numeric TPM at line '+str(line)) from error
                if np.isinf(value) or value<0:
                    raise ValueError('TPM must be nonnegative finite or explicitly missing')
                row[column] = value
            rows.append(row)
    if len(rows)!=n or not rows:
        raise ValueError('GCT row count differs from declared nonempty source population')
    records = pd.DataFrame(rows)
    if not records.source_gene_id.is_unique:
        raise ValueError('Repeated qualified source gene IDs are ambiguous')
    ordinary = records[records.source_qualifier.eq('ordinary')]
    if not ordinary.gene_id.is_unique:
        raise ValueError('Multiple annotation versions cannot share a canonical gene row')
    if not set(records.gene_id)<=set(ordinary.gene_id):
        raise ValueError('Qualified PAR_Y gene lacks its ordinary canonical source record')
    return records


def reviewed_gene_protein_links(mapping_path, reviewed_path) -> pd.DataFrame:
    """Retain all reviewed human Ensembl cross-references and their isoform/version evidence.

    Nothing is selected by symbol or source row order. Isoform accessions retain
    their original value alongside the reviewed parent protein accession. The
    caller verifies source hashes/species/release; IDs alone cannot prove species.
    """
    reviewed = set()
    with gzip.open(reviewed_path,'rt',newline='') as stream:
        reader = csv.reader(stream,delimiter='\t')
        if next(reader,None)!=['Entry','Gene Names (primary)']:
            raise ValueError('Expected the verified reviewed human entry-list schema')
        for row in reader:
            if len(row)!=2 or not PROTEIN.fullmatch(row[0]) or row[0] in reviewed:
                raise ValueError('Reviewed protein IDs must be canonical and unique')
            reviewed.add(row[0])
    if not reviewed:
        raise ValueError('Reviewed protein population is empty')
    records = set()
    with gzip.open(mapping_path,'rt',newline='') as stream:
        for line in stream:
            row = line.rstrip('\r\n').split('\t')
            if len(row)!=3:
                raise ValueError('Incomplete UniProt cross-reference line')
            if row[1]!='Ensembl':
                continue
            protein = row[0].split('-')[0]
            if protein not in reviewed:
                continue
            match = re.fullmatch(r'(ENSG\d{11})(?:\.(\d+))?',row[2])
            if match is None or not re.fullmatch(PROTEIN.pattern+r'(?:-\d+)?',row[0]):
                raise ValueError('Reviewed Ensembl reference contains a noncanonical human identifier')
            records.add((match.group(1),protein,row[0],row[2]))
    return pd.DataFrame(sorted(records),columns=['gene_id','protein_id','source_protein_id','source_gene_xref'])


def build_gene_foundation(source_records, protein_links) -> HumanGeneBuild:
    """Keep source gene measurements and quantify mapping ambiguity without projection."""
    required = {'source_row','source_gene_id','gene_id','annotation_version','source_qualifier','source_gene_name',*FIELDS.values()}
    if set(source_records)!=required or source_records.empty or not source_records.source_gene_id.is_unique:
        raise ValueError('Complete unique source gene records are required')
    records = source_records.copy()
    if records.source_row.tolist()!=list(range(4,len(records)+4)):
        raise ValueError('Source record positions/order must preserve the complete GCT population')
    for row in records.itertuples(index=False):
        match = SOURCE_GENE.fullmatch(row.source_gene_id)
        if (not match or match.group(1)!=row.gene_id or int(match.group(2))!=row.annotation_version
                or row.source_qualifier!=('PAR_Y' if match.group(3) else 'ordinary')):
            raise ValueError('Source identity and canonical gene projection disagree')
    values = records[list(FIELDS.values())].to_numpy(dtype=float)
    if np.isinf(values).any() or (values[np.isfinite(values)]<0).any():
        raise ValueError('TPM must be nonnegative finite or missing')
    nodes = records[records.source_qualifier.eq('ordinary')][['gene_id','source_gene_name',*FIELDS.values()]].rename(
        columns={'source_gene_name':'gene_name'}).reset_index(drop=True)
    validate_nodes(HUMAN_SPACE,nodes)
    if not set(records.gene_id)<=set(nodes.gene_id):
        raise ValueError('Qualified source record has no ordinary gene row')
    link_columns = ['gene_id','protein_id','source_protein_id','source_gene_xref']
    if list(protein_links)!=link_columns or protein_links.duplicated().any():
        raise ValueError('Unique explicit gene/protein cross-reference records are required')
    for row in protein_links.itertuples(index=False):
        if (not GENE.fullmatch(str(row.gene_id)) or not PROTEIN.fullmatch(str(row.protein_id))
                or not re.fullmatch(PROTEIN.pattern+r'(?:-\d+)?',str(row.source_protein_id))
                or row.source_protein_id.split('-')[0]!=row.protein_id
                or not re.fullmatch(re.escape(row.gene_id)+r'(?:\.\d+)?',str(row.source_gene_xref))):
            raise ValueError('Invalid explicit gene/protein cross-reference identity')
    links = protein_links.copy()
    links['in_gene_universe'] = links.gene_id.isin(nodes.gene_id)
    pairs = links[['gene_id','protein_id']].drop_duplicates()
    gene_counts = pairs.groupby('gene_id').protein_id.nunique()
    protein_counts = pairs.groupby('protein_id').gene_id.nunique()
    gene_mapping = nodes[['gene_id']].copy()
    gene_mapping['reviewed_protein_count'] = nodes.gene_id.map(gene_counts).fillna(0).astype(int)
    gene_mapping['protein_mapping_status'] = np.where(gene_mapping.reviewed_protein_count.eq(0),'unmapped',
        np.where(gene_mapping.reviewed_protein_count.eq(1),'one_reviewed_protein','multiple_reviewed_proteins'))
    shared = set(pairs[pairs.protein_id.map(protein_counts).gt(1)].gene_id)
    gene_mapping['has_shared_protein'] = gene_mapping.gene_id.isin(shared)
    summary = {'organism':O.HUMAN,'unit':'gene','source_records':len(records),'canonical_genes':len(nodes),
        'qualified_PAR_Y_records':int(records.source_qualifier.eq('PAR_Y').sum()),
        'unmapped_genes':int(gene_mapping.reviewed_protein_count.eq(0).sum()),
        'single_protein_genes':int(gene_mapping.reviewed_protein_count.eq(1).sum()),
        'multiple_protein_genes':int(gene_mapping.reviewed_protein_count.gt(1).sum()),
        'source_genes_with_shared_protein':int(gene_mapping.has_shared_protein.sum()),
        'reference_mapping_records':len(links),'reference_gene_protein_pairs':len(pairs),
        'mapping_genes_outside_source_universe':int((~gene_counts.index.isin(nodes.gene_id)).sum()),
        'mapped_proteins_in_source_universe':int(pairs[pairs.gene_id.isin(nodes.gene_id)].protein_id.nunique()),
        'reference_proteins_with_multiple_genes':int(protein_counts.gt(1).sum()),
        'source_measurement_cells':int(records[list(FIELDS.values())].notna().sum().sum()),
        'canonical_gene_measurement_cells':int(nodes[list(FIELDS.values())].notna().sum().sum()),
        'protein_measurements_projected':0,'source_records_lost':0,
        'universe_scope':'ordinary genes in this GTEx annotation; qualified Y records preserved separately'}
    return HumanGeneBuild(nodes,records,links,gene_mapping,summary)
