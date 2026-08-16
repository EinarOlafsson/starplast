#!/usr/bin/env python3
"""Strain variation, and the ordering that says the column is what it claims to be.

The real table is shipped, so the biology check runs against it rather than against a fixture: SRS
antigens above the virulence loci above the genome above the ribosomal core. Everything else here is
fixture-based, because the failure modes -- a renamed header, a missing file, an unresolved
accession -- are about the reader and not about Toxoplasma.
"""
from __future__ import annotations

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import variation as V  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HEADER = "\t".join(["Gene ID"] + list(V.COLUMNS))


def _table(tmp_path, body, header=HEADER):
    data = tmp_path / "starplast" / "data"
    data.mkdir(parents=True, exist_ok=True)
    (data / V.TABLE).write_text(header + "\n" + body)
    return tmp_path


def test_the_five_effect_classes_are_read_and_renamed(tmp_path):
    _table(tmp_path, "TGME49_200010\t119\t30\t20\t69\t0\n")
    out = V.strain_snps(str(tmp_path), log=lambda *_: None)
    assert list(out.columns) == list(V.COLUMNS.values())
    assert out.loc["TGME49_200010", "snp_nonsynonymous"] == 30


def test_zero_is_kept_as_a_measurement(tmp_path):
    """690 genes carry no SNP in any sequenced strain. That is a fact about the gene, and writing it
    as missing would move them from 'identical everywhere' to 'nobody looked'."""
    _table(tmp_path, "TGME49_200010\t0\t0\t0\t0\t0\n")
    out = V.strain_snps(str(tmp_path), log=lambda *_: None)
    assert out.loc["TGME49_200010", "snp_total_all_strains"] == 0
    assert out.notna().all().all()


def test_accessions_go_through_the_identity_layer(tmp_path):
    _table(tmp_path, "TGGT1_100010\t5\t1\t2\t2\t0\n")
    out = V.strain_snps(str(tmp_path),
                        resolve=lambda g: {"TGGT1_100010": "TGME49_200010"}.get(g),
                        log=lambda *_: None)
    assert list(out.index) == ["TGME49_200010"]


def test_an_unresolvable_accession_keeps_its_own_name(tmp_path):
    _table(tmp_path, "TGME49_999999\t5\t1\t2\t2\t0\n")
    out = V.strain_snps(str(tmp_path), resolve=lambda g: None, log=lambda *_: None)
    assert list(out.index) == ["TGME49_999999"]


def test_two_rows_for_one_gene_take_the_larger_count(tmp_path):
    """Two transcripts of one gene report the same locus twice; summing would double it."""
    _table(tmp_path, "TGME49_200010\t119\t30\t20\t69\t0\nTGME49_200010\t119\t30\t20\t69\t0\n")
    out = V.strain_snps(str(tmp_path), log=lambda *_: None)
    assert out.loc["TGME49_200010", "snp_total_all_strains"] == 119


def test_a_missing_table_yields_nothing(tmp_path):
    assert V.strain_snps(str(tmp_path), log=lambda *_: None).empty


def test_a_report_with_none_of_the_expected_columns_yields_nothing(tmp_path):
    """ToxoDB ships display names as the header, so a renamed attribute must not read as zero SNPs."""
    _table(tmp_path, "TGME49_200010\t1\n", header="Gene ID\tSomething Else")
    assert V.strain_snps(str(tmp_path), log=lambda *_: None).empty


def test_a_report_whose_first_column_is_not_the_gene_yields_nothing(tmp_path):
    _table(tmp_path, "chr1\t1\t1\t1\t1\t1\n",
           header="\t".join(["Sequence ID"] + list(V.COLUMNS)))
    assert V.strain_snps(str(tmp_path), log=lambda *_: None).empty


@pytest.mark.skipif(not os.path.exists(os.path.join(ROOT, "starplast", "data", V.TABLE)),
                    reason="shipped SNP table not present")
def test_polymorphism_orders_the_families_the_way_biology_does():
    """The check that no metadata field can give: does the column rank genes the way it must?

    SRS surface antigens are the most polymorphic thing in the genome because host immunity looks at
    them; ROP5/ROP18/GRA15 are the classic strain-typing virulence loci; ribosomal proteins are the
    conserved core. If a join went wrong this ordering collapses, and nothing else would say so.
    """
    snps = V.strain_snps(ROOT, log=lambda *_: None)
    n = pd.read_parquet(os.path.join(ROOT, "starplast", "data", "nodes.parquet")).set_index("gene_id")
    density = (snps["snp_nonsynonymous"] / n["length"].reindex(snps.index).clip(lower=1)) * 1000
    prod = n["product"].astype(str)

    def median_for(pattern):
        ids = [g for g in n.index[prod.str.contains(pattern, case=False, regex=True, na=False)]
               if g in density.index]
        return density.loc[ids].median()

    srs = median_for(r"SRS\d|SAG1|SAG2")
    virulence = median_for(r"ROP5|ROP18|GRA15")
    ribosome = median_for(r"ribosomal protein")
    genome = density.median()
    assert srs > virulence > genome > ribosome, (
        f"SRS {srs:.1f}, virulence {virulence:.1f}, genome {genome:.1f}, ribosome {ribosome:.1f}")
    assert ribosome * 5 < genome, "the conserved core is not conserved in this column"
