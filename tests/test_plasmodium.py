#!/usr/bin/env python3
"""The Plasmodium table: the second species, and the two ways it could quietly be wrong.

The report is served per TRANSCRIPT, so a gene with two transcripts is two rows and 71 genes would
be double-weighted by anyone who took the row count at face value. And the piggyBac screen's two
scores have a direction that is easy to invert -- a low mutagenesis index means the gene resists
disruption, which is to say it is essential. Inverting it would turn the essential genome into the
dispensable one and nothing would crash. Both are checked here against biology rather than against
the file: ribosomal proteins have to come out essential and the variant surface families have to
come out dispensable in culture and hypervariable between strains.
"""
from __future__ import annotations

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import plasmodium as P  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPORT = os.path.join(ROOT, "datasets", "reference", "plasmodb",
                      "plasmodb_pf3d7_gene_attributes.tsv")


def _report(tmp_path, rows=None, drop_id=False):
    """A PlasmoDB attributes report shaped like the real one."""
    rows = rows if rows is not None else [
        {"Gene ID": "PF3D7_0100100", "Gene Type": "protein coding gene",
         "Product Description": "erythrocyte membrane protein 1", "Protein Length": "2163",
         "Chromosome": "01", "Transcript Length": "6492", "# Exons in Transcript": "2",
         "Molecular Weight": "245807", "Isoelectric Point": "5.24", "# TM Domains": "0",
         "SignalP Peptide": "N/A", "Ortholog count": "2211", "Paralog count": "65",
         "Ortholog Group": "OG6_104345", "Total SNPs All Strains": "2687",
         "NonSynonymous SNPs All Strains": "1546", "Synonymous SNPs All Strains": "648",
         "Non-Coding SNPs All Strains": "474", "SNPs with Stop Codons All Strains": "19",
         "NonSyn/Syn SNP Ratio All Strains": "2.39",
         "P.falciparum 3D7 piggyBac insertion mutagenesis - mutant fitness score": "-1.706",
         "P.falciparum 3D7 piggyBac insertion mutagenesis - mutagenesis index score": "1",
         "Interpro ID": "IPR004258;IPR008602", "PFam ID": "PF03011"},
        {"Gene ID": "PF3D7_0100200", "Gene Type": "protein coding gene",
         "Product Description": "rifin", "Protein Length": "331", "Chromosome": "01",
         "Transcript Length": "996", "# Exons in Transcript": "2", "Molecular Weight": "37415",
         "Isoelectric Point": "8.86", "# TM Domains": "1", "SignalP Peptide": "Yes",
         "Ortholog count": "4601", "Paralog count": "0", "Ortholog Group": "OG6_100719",
         "Total SNPs All Strains": "10", "NonSynonymous SNPs All Strains": "5",
         "Synonymous SNPs All Strains": "3", "Non-Coding SNPs All Strains": "2",
         "SNPs with Stop Codons All Strains": "0", "NonSyn/Syn SNP Ratio All Strains": "1.67",
         "P.falciparum 3D7 piggyBac insertion mutagenesis - mutant fitness score": "null",
         "P.falciparum 3D7 piggyBac insertion mutagenesis - mutagenesis index score": "null",
         "Interpro ID": "", "PFam ID": ""}]
    frame = pd.DataFrame(rows)
    if drop_id:
        frame = frame.drop(columns=["Gene ID"])
    path = tmp_path / "report.tsv"
    frame.to_csv(path, sep="\t", index=False)
    return str(path)


# --------------------------------------------------------------------------- reading
def test_the_report_becomes_one_row_per_gene(tmp_path):
    d = P.build(_report(tmp_path))
    assert len(d) == 2 and d["gene_id"].is_unique


def test_numbers_are_numbers(tmp_path):
    d = P.build(_report(tmp_path)).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "length"] == 2163
    assert d.loc["PF3D7_0100100", "snp_nonsynonymous"] == 1546
    assert d.loc["PF3D7_0100100", "piggybac_mfs"] == pytest.approx(-1.706)


def test_the_three_ways_plasmodb_writes_a_missing_value_all_read_as_missing(tmp_path):
    d = P.build(_report(tmp_path)).set_index("gene_id")
    assert pd.isna(d.loc["PF3D7_0100200", "piggybac_mfs"])       # "null"
    assert not d.loc["PF3D7_0100100", "has_signal_peptide"]      # "N/A"
    assert d.loc["PF3D7_0100200", "n_interpro"] == 0             # ""


# --------------------------------------------------------------------------- the transcript trap
def test_a_gene_with_two_transcripts_yields_one_row(tmp_path):
    """The report is per transcript. 5,791 rows describe 5,720 genes."""
    base = pd.read_csv(_report(tmp_path), sep="\t", dtype=str).to_dict("records")
    second = dict(base[0])
    second["Transcript Length"] = "3000"
    second["Protein Length"] = "999"
    d = P.build(_report(tmp_path, rows=base + [second]))
    assert len(d) == 2
    row = d[d.gene_id == "PF3D7_0100100"]
    assert len(row) == 1


def test_the_longest_transcript_wins(tmp_path):
    """So protein length and exon count describe one transcript rather than a mix of two."""
    base = pd.read_csv(_report(tmp_path), sep="\t", dtype=str).to_dict("records")
    short = dict(base[0])
    short["Transcript Length"] = "10"
    short["Protein Length"] = "5"
    d = P.build(_report(tmp_path, rows=[short] + base)).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "length"] == 2163


# --------------------------------------------------------------------------- derived columns
def test_domain_content_is_counted_from_the_identifier_list(tmp_path):
    d = P.build(_report(tmp_path)).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "n_interpro"] == 2
    assert d.loc["PF3D7_0100100", "has_domain"]
    assert not d.loc["PF3D7_0100200", "has_domain"]


def test_an_uncalled_signal_peptide_is_not_read_as_no_signal_peptide(tmp_path):
    """Only 10.5% of genes carry a call. That is a fact about what was run."""
    d = P.build(_report(tmp_path)).set_index("gene_id")
    assert d.loc["PF3D7_0100200", "has_signal_peptide"]
    assert "signalp" not in d.columns


def test_membrane_and_paralogy_flags_follow_their_counts(tmp_path):
    d = P.build(_report(tmp_path)).set_index("gene_id")
    assert not d.loc["PF3D7_0100100", "is_tm"] and d.loc["PF3D7_0100200", "is_tm"]
    assert d.loc["PF3D7_0100100", "has_paralog"] and not d.loc["PF3D7_0100200", "has_paralog"]


# --------------------------------------------------------------------------- refusals
def test_a_missing_report_yields_nothing(tmp_path):
    assert P.build(str(tmp_path / "absent.tsv")).empty


def test_a_report_without_a_gene_column_is_refused(tmp_path):
    assert P.build(_report(tmp_path, drop_id=True)).empty


def test_collapse_refuses_a_frame_with_no_gene_column():
    assert P._collapse(pd.DataFrame({"length": [1]})).empty


def test_a_report_with_a_header_and_no_genes_yields_nothing(tmp_path):
    """PlasmoDB answers an empty search with column names and no rows."""
    path = _report(tmp_path)
    header = open(path).readline()
    empty = tmp_path / "empty.tsv"
    empty.write_text(header)
    d = P.build(str(empty))
    assert d.empty and "n_interpro" not in d.columns


def test_load_returns_nothing_when_the_table_is_not_built(tmp_path):
    assert P.load(str(tmp_path)).empty


# --------------------------------------------------------------------------- the real table
@pytest.mark.skipif(not os.path.exists(REPORT), reason="PlasmoDB report not fetched")
def test_the_real_report_gives_the_gene_count_the_registry_claims():
    d = P.build(REPORT)
    assert len(d) == 5720 and d["gene_id"].is_unique
    assert d["gene_id"].str.startswith("PF3D7_").all()


@pytest.mark.skipif(not os.path.exists(REPORT), reason="PlasmoDB report not fetched")
def test_the_piggybac_direction_is_not_inverted():
    """A low mutagenesis index means the gene resists disruption, which means it is essential.

    Checked against biology rather than against the file: ribosomal proteins are essential and the
    variant surface antigen families are dispensable in culture. If the sign were flipped this is
    the only thing that would notice.
    """
    d = P.build(REPORT)
    d = d[d["piggybac_mis"].notna()]
    ribosomal = d[d["product"].str.contains("ribosomal protein", case=False, na=False)]
    surface = d[d["product"].str.contains("erythrocyte membrane protein 1|rifin|stevor",
                                          case=False, na=False)]
    assert len(ribosomal) > 100 and len(surface) > 100
    assert ribosomal["piggybac_mis"].median() < 0.3
    assert surface["piggybac_mis"].median() > 0.8
    assert ribosomal["piggybac_mfs"].median() < surface["piggybac_mfs"].median()


@pytest.mark.skipif(not os.path.exists(REPORT), reason="PlasmoDB report not fetched")
def test_the_variant_surface_families_carry_the_strain_variation():
    """The other direction check: var, rifin and stevor are the hypervariable families."""
    d = P.build(REPORT)
    surface = d[d["product"].str.contains("erythrocyte membrane protein 1|rifin|stevor",
                                          case=False, na=False)]
    ribosomal = d[d["product"].str.contains("ribosomal protein", case=False, na=False)]
    assert surface["snp_nonsynonymous"].median() > 100
    assert ribosomal["snp_nonsynonymous"].median() < 10


@pytest.mark.skipif(not os.path.exists(os.path.join(ROOT, "starplast", "data", P.TABLE)),
                    reason="Plasmodium table not built")
def test_the_shipped_table_is_falciparum_and_not_gondii():
    """Nothing is merged. A TGME49 accession in here would mean two species in one table."""
    d = P.load(ROOT)
    assert len(d) > 5000
    assert not d["gene_id"].str.contains("TGME49_|TGGT1_").any()
