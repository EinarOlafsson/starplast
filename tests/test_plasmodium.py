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


# --------------------------------------------------------------------------- expression
EXPR = os.path.join(ROOT, "datasets", "reference", "plasmodb", P.EXPRESSION_TABLE)

#: Stage markers whose timing is textbook, used to check that the columns are labelled correctly.
#: Pfs25 is the interesting one: its transcript is stockpiled in mature female gametocytes and only
#: translated in the ookinete, so a TRANSCRIPT peak in gametocyte V is right and a peak in ookinete
#: would suggest the columns had been shifted by one.
MARKERS = {
    "PF3D7_0304600": ("expr_sporozoite", "circumsporozoite protein"),
    "PF3D7_0930300": ("expr_schizont", "merozoite surface protein 1"),
    "PF3D7_0406200": ("expr_gametocyte_ii", "early gametocyte marker Pfs16"),
    "PF3D7_1031000": ("expr_gametocyte_v", "Pfs25, transcript stockpiled before the ookinete"),
}
STAGES = ("expr_ring", "expr_early_trophozoite", "expr_late_trophozoite", "expr_schizont",
          "expr_gametocyte_ii", "expr_gametocyte_v", "expr_ookinete", "expr_oocyst",
          "expr_sporozoite")


def _expr_report(tmp_path, header=True):
    rows = {"Gene ID": ["PF3D7_0100100", "PF3D7_0100200"]}
    for study, sample, _column in P.EXPRESSION:
        rows[f"{study} - {sample} - unique" if sample else study] = ["1.5", "N/A"]
    path = tmp_path / "expr.tsv"
    frame = pd.DataFrame(rows)
    if not header:
        frame = frame.drop(columns=["Gene ID"])
    frame.to_csv(path, sep="\t", index=False)
    return str(path)


def test_the_expression_report_maps_every_sample_to_a_column(tmp_path):
    d = P.expression(_expr_report(tmp_path))
    for _study, _sample, column in P.EXPRESSION:
        assert column in d.columns, column


def test_expression_values_are_numeric_and_missing_stays_missing(tmp_path):
    d = P.expression(_expr_report(tmp_path)).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "expr_ring"] == pytest.approx(1.5)
    assert pd.isna(d.loc["PF3D7_0100200", "expr_ring"])


def test_a_missing_expression_report_yields_nothing(tmp_path):
    assert P.expression(str(tmp_path / "absent.tsv")).empty


def test_an_expression_report_without_a_gene_column_is_refused(tmp_path):
    assert P.expression(_expr_report(tmp_path, header=False)).empty


def test_an_expression_report_with_no_recognised_sample_is_refused(tmp_path):
    path = tmp_path / "e.tsv"
    pd.DataFrame({"Gene ID": ["PF3D7_0100100"], "something else": ["1"]}).to_csv(
        path, sep="\t", index=False)
    assert P.expression(str(path)).empty


@pytest.mark.skipif(not os.path.exists(EXPR), reason="expression report not fetched")
def test_the_life_stage_columns_are_not_mislabelled():
    """Marker genes have to peak where a hundred years of malaria biology says they do.

    The columns come out of PlasmoDB in one wide report and are matched by substring, so a mistake
    here would silently shift a stage. Nothing else in the pipeline would notice.
    """
    d = P.expression(EXPR).set_index("gene_id")
    for gene, (stage, why) in MARKERS.items():
        assert gene in d.index, gene
        peak = d.loc[gene, list(STAGES)].idxmax()
        assert peak == stage, f"{gene} ({why}) peaks in {peak}, expected {stage}"


@pytest.mark.skipif(not os.path.exists(EXPR), reason="expression report not fetched")
def test_polysomal_and_steady_state_are_kept_apart():
    """One is what is on ribosomes and the other is what is in the cell. Averaging them would
    destroy the only comparison in this table that separates transcription from translation."""
    d = P.expression(EXPR)
    poly = [c for c in d.columns if c.startswith("polysomal_")]
    steady = [c for c in d.columns if c.startswith("steady_state_")]
    assert len(poly) == 3 and len(steady) == 3
    assert not set(poly) & set(steady)
    # They must also disagree: identical columns would mean the substring match caught one twice.
    for a, b in zip(sorted(poly), sorted(steady)):
        assert not d[a].equals(d[b]), f"{a} and {b} are the same column"


@pytest.mark.skipif(not os.path.exists(EXPR), reason="expression report not fetched")
def test_the_proteome_columns_are_compositional_and_say_so_in_their_names():
    """PlasmoDB serves that TMT study row-normalised, and the name has to carry that.

    Called `protein_ring` it reads as an abundance, and the slot catalog would have claimed it for
    `protein abundance · asexual blood stage`. It is a share: the three values sum to a constant per
    gene and the columns are anti-correlated with one another by construction. The check is on the
    data rather than on the name so that a future PlasmoDB release serving true abundances fails
    here instead of silently keeping a wrong label.
    """
    d = P.expression(EXPR)
    share = [c for c in d.columns if c.startswith("protein_stage_share_")]
    assert len(share) == 3
    rows = d[share].dropna()
    assert rows.sum(axis=1).std() < 0.5, "no longer compositional: rename and re-check the slot"
    corr = rows.corr(method="spearman")
    off = [corr.iloc[i, j] for i in range(3) for j in range(3) if i != j]
    assert sum(v < 0 for v in off) >= 4, "compositional columns should mostly anti-correlate"
    assert not any(c.startswith("protein_ring") for c in d.columns), (
        "a bare `protein_` name would let the abundance slot claim a share")


def test_no_plasmodium_slot_claims_the_compositional_proteome():
    """The slot this was nearly given asks how much protein there is, which this cannot answer."""
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "gst", os.path.join(ROOT, "scripts", "generate_slot_table.py"))
    gst = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gst)
    for slot, patterns in gst.PF_PATTERNS.items():
        for pattern in patterns:
            assert not pattern.startswith("protein_stage_share"), f"{slot} claims a share"
        if "protein abundance" in slot:
            assert not patterns, f"{slot} must stay empty until a true abundance arrives"


# --------------------------------------------------------------------------- export prediction
EXPORT = os.path.join(ROOT, "datasets", "reference", "plasmodb", P.EXPORT_DIR)
NODES = os.path.join(ROOT, "starplast", "data", P.TABLE)


def _export_dir(tmp_path, tiers=None):
    """PlasmoDB answers ExportPred as a gene list per score threshold, so the tiers are nested."""
    tiers = tiers if tiers is not None else {1: ["A", "B", "C"], 5: ["A", "B"], 10: ["A"]}
    folder = tmp_path / P.EXPORT_DIR
    folder.mkdir()
    for threshold, genes in tiers.items():
        pd.DataFrame({"Gene ID": genes, "Product Description": ["x"] * len(genes)}).to_csv(
            folder / f"exportpred_score_ge_{threshold}.tsv", sep="\t", index=False)
    return str(folder)


def test_export_tiers_are_ordinal(tmp_path):
    d = P.export_prediction(_export_dir(tmp_path)).set_index("gene_id")
    assert d.loc["A", "export_pred_tier"] == 3
    assert d.loc["B", "export_pred_tier"] == 2
    assert d.loc["C", "export_pred_tier"] == 1


def test_only_the_top_tier_counts_as_exported_at_the_default(tmp_path):
    d = P.export_prediction(_export_dir(tmp_path)).set_index("gene_id")
    assert d.loc["A", "is_exported"]
    assert not d.loc["B", "is_exported"] and not d.loc["C", "is_exported"]


def test_a_missing_export_folder_yields_nothing(tmp_path):
    assert P.export_prediction(str(tmp_path / "absent")).empty


def test_an_export_folder_with_no_readable_answer_yields_nothing(tmp_path):
    folder = tmp_path / P.EXPORT_DIR
    folder.mkdir()
    pd.DataFrame({"wrong": ["A"]}).to_csv(folder / "exportpred_score_ge_1.tsv", sep="\t", index=False)
    assert P.export_prediction(str(folder)).empty


def test_a_partly_fetched_export_folder_still_gives_the_tiers_it_has(tmp_path):
    d = P.export_prediction(_export_dir(tmp_path, tiers={1: ["A", "B"]}))
    assert set(d["gene_id"]) == {"A", "B"}
    assert (d["export_pred_tier"] == 1).all() and not d["is_exported"].any()


@pytest.mark.skipif(not os.path.isdir(EXPORT), reason="ExportPred answers not fetched")
def test_the_export_prediction_finds_the_known_exportome():
    """Composition check: the exported set has to be the families every malaria textbook lists."""
    d = P.export_prediction(EXPORT)
    assert 150 < int(d["is_exported"].sum()) < 300
    top = d[d["is_exported"]]["gene_id"]
    assert "PF3D7_0202000" in set(top), "KAHRP is not called exported"


@pytest.mark.skipif(not os.path.isdir(EXPORT), reason="ExportPred answers not fetched")
def test_the_tiers_keep_exported_proteins_the_default_threshold_loses():
    """Why this is a tier and not a boolean. MESA and PfEMP3 are exported and score below 10."""
    d = P.export_prediction(EXPORT).set_index("gene_id")
    for gene, name in (("PF3D7_0500800", "MESA"), ("PF3D7_0201900", "PfEMP3")):
        assert gene in d.index, f"{name} is not called exported at any threshold"
        assert d.loc[gene, "export_pred_tier"] < 3, f"{name} was expected below the default"


@pytest.mark.skipif(not os.path.exists(NODES), reason="Plasmodium table not built")
def test_a_gene_absent_from_exportpred_is_predicted_not_exported_rather_than_unknown():
    """A sequence model was run on every protein, so its silence is an answer."""
    d = pd.read_parquet(NODES, columns=["gene_id", "export_pred_tier", "is_exported"])
    assert d["export_pred_tier"].notna().all()
    assert (d["export_pred_tier"] == 0).sum() > 5000


# --------------------------------------------------------------------------- whole-table assembly
def _dataset_root(tmp_path, with_stages=True, with_export=True, genes=True):
    """A dataset tree shaped the way `build_all` expects to find one."""
    base = tmp_path / "reference" / "plasmodb"
    base.mkdir(parents=True)
    if genes:
        src = _report(tmp_path)
        os.replace(src, base / "plasmodb_pf3d7_gene_attributes.tsv")
    if with_stages:
        rows = {"Gene ID": ["PF3D7_0100100", "PF3D7_0100200"]}
        for study, sample, _column in P.EXPRESSION:
            rows[f"{study} - {sample} - unique" if sample else study] = ["1.5", "2.5"]
        pd.DataFrame(rows).to_csv(base / P.EXPRESSION_TABLE, sep="\t", index=False)
    if with_export:
        folder = base / P.EXPORT_DIR
        folder.mkdir()
        for threshold, members in {1: ["PF3D7_0100100"], 5: ["PF3D7_0100100"],
                                   10: ["PF3D7_0100100"]}.items():
            pd.DataFrame({"Gene ID": members}).to_csv(
                folder / f"exportpred_score_ge_{threshold}.tsv", sep="\t", index=False)
    return str(tmp_path)


def test_build_all_assembles_the_three_reports(tmp_path):
    d = P.build_all(_dataset_root(tmp_path), log=lambda *a: None)
    assert {"length", "expr_ring", "export_pred_tier", "is_exported"} <= set(d.columns)
    assert len(d) == 2


def test_build_all_records_an_unpredicted_gene_as_not_exported(tmp_path):
    d = P.build_all(_dataset_root(tmp_path), log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "is_exported"]
    assert d.loc["PF3D7_0100200", "export_pred_tier"] == 0
    assert not d.loc["PF3D7_0100200", "is_exported"]


def test_build_all_survives_a_dataset_root_with_only_the_gene_report(tmp_path):
    d = P.build_all(_dataset_root(tmp_path, with_stages=False, with_export=False),
                    log=lambda *a: None)
    assert len(d) == 2 and "expr_ring" not in d.columns and "is_exported" not in d.columns


def test_build_all_without_a_gene_report_yields_nothing(tmp_path):
    assert P.build_all(_dataset_root(tmp_path, genes=False), log=lambda *a: None).empty
