"""Source-level non-loss, pseudoautosomal qualifiers and many-to-many host identity."""
import gzip

import numpy as np
import pandas as pd
import pytest

from starplast import organisms as O, strategies as S
from starplast.spaces import human as H


def source(tmp_path,text=None):
    path = tmp_path/'source.gct.gz'
    content = text or ('#1.2\n4\t2\nName\tDescription\tCells_Cultured_fibroblasts\tLiver_Hepatocyte\n'
        'ENSG00000000001.2\tGeneA\t0.0\t0.12345678901234567\n'
        'ENSG00000000001.2_PAR_Y\tGeneA\t99.0\t88.0\n'
        'ENSG00000000002.3\tGeneB\t2.0\tNA\n'
        'ENSG00000000003.1\tGeneC\t3.0\t4.0\n')
    with gzip.open(path,'wt') as f:f.write(content)
    return path


def links():
    return pd.DataFrame([
        ['ENSG00000000001','P11111','P11111','ENSG00000000001.2'],
        ['ENSG00000000001','P11111','P11111-2','ENSG00000000001.2'],
        ['ENSG00000000001','Q11111','Q11111','ENSG00000000001.2'],
        ['ENSG00000000002','Q11111','Q11111','ENSG00000000002.3'],
        ['ENSG00000000009','P99999','P99999','ENSG00000000009.1']],
        columns=['gene_id','protein_id','source_protein_id','source_gene_xref'])


def test_source_text_values_zero_and_qualified_y_are_preserved_without_overwriting(tmp_path):
    records = H.read_gtex_gene_records(source(tmp_path))
    batch = H.build_gene_foundation(records,links())
    assert batch.nodes.gene_id.tolist()==['ENSG00000000001','ENSG00000000002','ENSG00000000003']
    assert batch.nodes.gtex_fibroblast_median_tpm.tolist()==[0.,2.,3.]
    assert batch.nodes.gtex_lcm_hepatocyte_median_tpm.iloc[0]==float('0.12345678901234567')
    assert np.isnan(batch.nodes.gtex_lcm_hepatocyte_median_tpm.iloc[1])
    pd.testing.assert_frame_equal(batch.source_records,records,check_exact=True)
    assert batch.source_records.source_row.tolist()==[4,5,6,7]
    assert batch.mapping_summary['qualified_PAR_Y_records']==1
    assert batch.mapping_summary['source_records_lost']==0
    assert batch.mapping_summary['canonical_gene_measurement_cells']==5
    assert batch.mapping_summary['source_measurement_cells']==7
    assert batch.mapping_summary['protein_measurements_projected']==0
    # Numeric source positions/versions and mapping degree must not become
    # biological features when these nodes later enter a host Context.
    assert set(batch.nodes)-{'gene_id','gene_name'}==set(H.FIELDS.values())
    context = S.Context(batch.nodes,graph={},organism=O.HUMAN)
    assert context.numeric_columns(min_values=2)==list(H.FIELDS.values())
    assert context.categorical_columns(min_labelled=2)==[]
    # Shape/content validation does not register or activate a host pack.
    assert H.HUMAN_SPACE.gene_regex==r'ENSG\d{11}'


def test_many_to_many_reverse_links_and_unmapped_genes_are_quantified(tmp_path):
    batch = H.build_gene_foundation(H.read_gtex_gene_records(source(tmp_path)),links())
    summary = batch.mapping_summary
    assert summary['unmapped_genes']==summary['single_protein_genes']==summary['multiple_protein_genes']==1
    assert summary['reference_mapping_records']==5 and summary['reference_gene_protein_pairs']==4
    assert summary['reference_proteins_with_multiple_genes']==1 and summary['source_genes_with_shared_protein']==2
    assert summary['mapping_genes_outside_source_universe']==1
    assert batch.protein_links.source_protein_id.tolist()==links().source_protein_id.tolist()
    assert batch.protein_links.in_gene_universe.tolist()==[True,True,True,True,False]
    altered = H.build_gene_foundation(batch.source_records,links().iloc[::-1].reset_index(drop=True))
    pd.testing.assert_frame_equal(batch.nodes,altered.nodes,check_exact=True)
    assert batch.mapping_summary==altered.mapping_summary


@pytest.mark.parametrize('mutation',[
    lambda s:s.replace('4\t2','5\t2'),
    lambda s:s.replace('GeneC\t3.0\t4.0','GeneC\t3.0'),
    lambda s:s.replace('ENSG00000000003.1','P11111'),
    lambda s:s.replace('ENSG00000000003.1','ENSG00000000001.3'),
    lambda s:s.replace('ENSG00000000001.2_PAR_Y','ENSG00000000004.2_PAR_Y'),
    lambda s:s.replace('GeneC\t3.0','GeneC\t-1.0'),
    lambda s:s.replace('GeneC\t3.0','GeneC\tinf')])
def test_dimension_identity_collision_or_invalid_measurement_is_refused(tmp_path,mutation):
    path = source(tmp_path)
    with gzip.open(path,'rt') as f:text=f.read()
    with pytest.raises(ValueError):H.read_gtex_gene_records(source(tmp_path,mutation(text)))


def test_tampered_source_or_cross_reference_cannot_change_gene_identity(tmp_path):
    records = H.read_gtex_gene_records(source(tmp_path))
    changed = records.copy();changed.loc[0,'gene_id']='ENSG00000000002'
    with pytest.raises(ValueError,match='projection disagree'):H.build_gene_foundation(changed,links())
    with pytest.raises(ValueError,match='positions/order'):H.build_gene_foundation(records.iloc[::-1],links())
    changed_links = links();changed_links.loc[0,'source_gene_xref']='ENSG00000000002.2'
    with pytest.raises(ValueError,match='cross-reference identity'):H.build_gene_foundation(records,changed_links)


def test_reviewed_reference_retains_versions_isoforms_and_discards_unreviewed(tmp_path):
    reviewed = tmp_path/'reviewed.tsv.gz';mapping = tmp_path/'mapping.dat.gz'
    with gzip.open(reviewed,'wt') as f:f.write('Entry\tGene Names (primary)\nP11111\tA\nQ11111\tB\n')
    with gzip.open(mapping,'wt') as f:f.write('P11111\tEnsembl\tENSG00000000001.2\nP11111-2\tEnsembl\tENSG00000000001.2\n'
        'Q11111\tEnsembl\tENSG00000000001.2\nP99999\tEnsembl\tENSG00000000003.1\nP11111\tGene_Name\tA\n')
    parsed = H.reviewed_gene_protein_links(mapping,reviewed)
    assert len(parsed)==3 and set(parsed.source_protein_id)=={'P11111','P11111-2','Q11111'}
    assert parsed.source_gene_xref.eq('ENSG00000000001.2').all()
    # A truncated gzip must fail rather than quietly demote missing proteins.
    reviewed.write_bytes(reviewed.read_bytes()[:-8])
    with pytest.raises((EOFError,OSError)):H.reviewed_gene_protein_links(mapping,reviewed)
