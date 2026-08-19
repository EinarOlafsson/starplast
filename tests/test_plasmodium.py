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

import numpy as np
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


def test_the_compositional_proteome_is_a_share_and_never_an_abundance():
    """The slot this was nearly given asks how much protein there is, which this cannot answer.

    It used to assert that NOTHING claimed these columns, which was right while they had no slot and
    wrong afterwards: a column no slot claims can never be held out or audited, and the slot tree's
    orphan alarm found these three sitting unclaimed. They now have a slot that says `share` in its
    name, and the rule that matters is the narrower one -- the ABUNDANCE slot must still be empty,
    because row sums are constant and the stages anti-correlate, so these say what fraction of a
    protein's signal falls in each stage rather than how much of it there is.
    """
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "gst", os.path.join(ROOT, "scripts", "generate_slot_table.py"))
    gst = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gst)
    claimants = [slot for slot, patterns in gst.PF_PATTERNS.items()
                 for p in patterns if p.startswith("protein_stage_share")]
    assert claimants == ["protein stage share"], claimants
    for slot, patterns in gst.PF_PATTERNS.items():
        if slot.startswith("protein abundance"):
            assert not patterns, f"{slot} must stay empty until a true abundance arrives"
    from starplast import slots as S
    abundance = [s for s in S.all_slots("Pf") if s.name.startswith("protein abundance")]
    assert abundance, "the abundance slots vanished rather than staying empty"
    for slot in abundance:
        assert not slot.patterns, f"{slot.name} started claiming columns"


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


def test_build_all_folds_in_the_phosphosites_and_sets_the_flag_everywhere(tmp_path):
    root = _dataset_root(tmp_path)
    folder = os.path.join(root, "reference", "plasmodb", P.PHOSPHO_DIR)
    os.makedirs(folder)
    pd.DataFrame([("PF3D7_0100100.1-p1", "12", "Phospho", "1")],
                 columns=["Proteins", "Protein Modification Positions", "Modification",
                          P.PHOSPHO_THRESHOLD]).to_csv(
        os.path.join(folder, "x_Site_Peptidoform_centric.tsv"), sep="\t", index=False)
    d = P.build_all(root, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "n_phosphosites"] == 1
    assert d.loc["PF3D7_0100100", "has_phospho"]
    assert not d.loc["PF3D7_0100200", "has_phospho"]
    assert pd.isna(d.loc["PF3D7_0100200", "n_phosphosites"])


def test_build_all_folds_in_the_palmitome_and_flags_the_rest_false(tmp_path):
    root = _dataset_root(tmp_path)
    folder = os.path.join(root, "post_translation", P.PALMITOME[0], P.PALMITOME[1])
    os.makedirs(folder)
    with pd.ExcelWriter(os.path.join(folder, P.PALMITOME[2])) as writer:
        pd.DataFrame({P.PALMITOME_PREDICTED: ["PF3D7_0100200"],
                      P.PALMITOME_OBSERVED: ["PF3D7_0100100"]}).to_excel(
            writer, sheet_name=P.PALMITOME_SHEET, index=False)
    d = P.build_all(root, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "is_palmitoylated"]
    assert not d.loc["PF3D7_0100200", "is_palmitoylated"], "the predicted column leaked in"


def test_build_all_folds_in_the_chromatin_perturbation(tmp_path):
    root = _dataset_root(tmp_path)
    base = os.path.join(root, "reference", "plasmodb")
    rows = {"Gene ID": ["PF3D7_0100100", "PF3D7_0100200"]}
    for label, _column in P.SIR2:
        rows[f"{label} (Sir2 KO Marray)"] = ["2.7", "2.9"]
    pd.DataFrame(rows).to_csv(os.path.join(base, P.SIR2_TABLE), sep="\t", index=False)
    d = P.build_all(root, log=lambda *a: None)
    assert "sir2_wt_ring" in d.columns and "sir2a_ko_ring" in d.columns


def test_build_all_folds_in_the_model_confidence(tmp_path):
    root = _dataset_root(tmp_path)
    base = os.path.join(root, "reference", "plasmodb")
    rows = {"gene_id": ["PF3D7_0100100"], "uniprot_used": ["Q1"]}
    for source, _column in P.ALPHAFOLD:
        rows[source] = ["70.5"] if source == "globalMetricValue" else ["0.25"]
    pd.DataFrame(rows).to_csv(os.path.join(base, P.ALPHAFOLD_TABLE), sep="\t", index=False)
    d = P.build_all(root, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "mean_plddt"] == pytest.approx(70.5)
    assert pd.isna(d.loc["PF3D7_0100200", "mean_plddt"])


def test_build_all_leaves_unassayed_genes_missing_for_myristoylation(tmp_path):
    root = _dataset_root(tmp_path)
    folder = os.path.join(root, "post_translation", P.MYRISTOYLOME[0], P.MYRISTOYLOME[1])
    os.makedirs(folder)
    pd.DataFrame([("PF3D7_0100100.1-p1", "+", -2.4)],
                 columns=["Protein IDs", "Significant", "Difference"]).to_excel(
        os.path.join(folder, P.MYRISTOYLOME[2]), index=False)
    d = P.build_all(root, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "is_myristoylated"]
    assert pd.isna(d.loc["PF3D7_0100200", "is_myristoylated"]), "an unassayed gene got a call"


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


# --------------------------------------------------------------------------- phosphosites
PHOS = os.path.join(ROOT, "datasets", "reference", "plasmodb", P.PHOSPHO_DIR)


def _phos_dir(tmp_path, rows=None):
    rows = rows if rows is not None else [
        # the same site seen twice, in two peptidoforms -- one site, not two
        ("PF3D7_0100100.1-p1", "12", "Phospho", "1"),
        ("PF3D7_0100100.1-p1", "12", "Phospho", "1"),
        ("PF3D7_0100100.1-p1", "44", "Phospho", "1"),
        ("PF3D7_0100200.1-p1", "7", "Phospho", "1"),
        ("PF3D7_0100200.1-p1", "9", "Phospho", "0"),      # fails the q-value cut
        ("PF3D7_0100200.1-p1", "9", "Acetyl", "1"),       # not phosphorylation
        ("CONTAM_HUMAN", "3", "Phospho", "1"),            # not a Plasmodium gene
    ]
    folder = tmp_path / P.PHOSPHO_DIR
    folder.mkdir()
    pd.DataFrame(rows, columns=["Proteins", "Protein Modification Positions", "Modification",
                                P.PHOSPHO_THRESHOLD]).to_csv(
        folder / "PXD000001_Site_Peptidoform_centric.tsv", sep="\t", index=False)
    return str(folder)


def test_a_site_seen_twice_is_counted_once(tmp_path):
    """The tables are one row per peptidoform per run, so rows count looking, not sites."""
    d = P.phosphosites(_phos_dir(tmp_path), log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "n_phosphosites"] == 2


def test_sites_below_the_q_value_and_other_modifications_are_excluded(tmp_path):
    d = P.phosphosites(_phos_dir(tmp_path), log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100200", "n_phosphosites"] == 1


def test_non_plasmodium_accessions_are_dropped(tmp_path):
    d = P.phosphosites(_phos_dir(tmp_path), log=lambda *a: None)
    assert d["gene_id"].str.startswith("PF3D7_").all()


def test_the_transcript_and_product_suffix_is_stripped(tmp_path):
    """`PF3D7_0100100.1-p1` is one gene; two products would otherwise count their sites twice."""
    d = P.phosphosites(_phos_dir(tmp_path), log=lambda *a: None)
    assert set(d["gene_id"]) == {"PF3D7_0100100", "PF3D7_0100200"}


def test_a_missing_phospho_folder_yields_nothing(tmp_path):
    assert P.phosphosites(str(tmp_path / "absent"), log=lambda *a: None).empty


def test_a_phospho_table_without_a_protein_column_is_skipped(tmp_path):
    folder = tmp_path / P.PHOSPHO_DIR
    folder.mkdir()
    pd.DataFrame({"wrong": ["x"]}).to_csv(folder / "a_Site_Peptidoform_centric.tsv",
                                          sep="\t", index=False)
    assert P.phosphosites(str(folder), log=lambda *a: None).empty


@pytest.mark.skipif(not os.path.isdir(PHOS), reason="phosphosite tables not fetched")
def test_the_pooled_phosphoproteome_behaves_like_a_phosphoproteome():
    """Two orderings that have to hold, and would not if rows were being counted instead of sites.

    Kinases are phosphorylated more often than proteins at large -- activation loops and
    autophosphorylation -- and longer proteins carry more sites. Both survive the numbers changing.
    """
    import scipy.stats as st
    nodes = pd.read_parquet(NODES, columns=["gene_id", "product", "length", "n_phosphosites",
                                            "has_phospho"])
    assert 1000 < int(nodes["has_phospho"].sum()) < 4000
    kinase = nodes["product"].str.contains("kinase", case=False, na=False)
    rate_kinase = nodes.loc[kinase, "has_phospho"].mean()
    rate_all = nodes["has_phospho"].mean()
    assert rate_kinase > rate_all * 1.3, (rate_kinase, rate_all)
    sub = nodes.dropna(subset=["n_phosphosites", "length"])
    assert st.spearmanr(sub["n_phosphosites"], sub["length"]).statistic > 0.3


@pytest.mark.skipif(not os.path.exists(NODES), reason="Plasmodium table not built")
def test_the_count_stays_missing_where_nothing_was_detected_but_the_flag_does_not():
    """Mirrors the Toxoplasma convention: how many sites is unknown, whether any were seen is no."""
    d = pd.read_parquet(NODES, columns=["n_phosphosites", "has_phospho"])
    assert d["has_phospho"].notna().all()
    assert d["n_phosphosites"].isna().sum() > 2000
    assert (d["n_phosphosites"] == 0).sum() == 0


# --------------------------------------------------------------------------- palmitome
PALM = os.path.join(ROOT, "datasets", "post_translation", P.PALMITOME[0], P.PALMITOME[1],
                    P.PALMITOME[2])


def _palm_root(tmp_path, sheet=None, observed=None):
    folder = tmp_path / "post_translation" / P.PALMITOME[0] / P.PALMITOME[1]
    folder.mkdir(parents=True)
    frame = pd.DataFrame({
        P.PALMITOME_PREDICTED: ["PF3D7_0100100", "PF3D7_0100200", "PF3D7_0100300"],
        P.PALMITOME_OBSERVED: (observed if observed is not None
                               else ["PF3D7_0100200", None, None])})
    with pd.ExcelWriter(folder / P.PALMITOME[2]) as writer:
        frame.to_excel(writer, sheet_name=sheet or P.PALMITOME_SHEET, index=False)
    return str(tmp_path)


def test_only_the_observed_column_is_read_not_the_predicted_one(tmp_path):
    """The trap this loader exists for.

    The source workbook also has a sheet CALLED `nrPalmitoylatedProteins` whose contents are the
    union of observed and motif-predicted -- 3,105 of 5,720 genes. Reading by name would have
    called 54% of the proteome palmitoylated against published palmitomes of 400 to 500.
    """
    d = P.palmitome(_palm_root(tmp_path), log=lambda *a: None)
    assert list(d["gene_id"]) == ["PF3D7_0100200"]
    assert "PF3D7_0100100" not in set(d["gene_id"]), "the predicted column leaked in"


def test_a_missing_palmitome_yields_nothing(tmp_path):
    assert P.palmitome(str(tmp_path), log=lambda *a: None).empty


def test_a_palmitome_without_the_expected_sheet_is_refused(tmp_path):
    assert P.palmitome(_palm_root(tmp_path, sheet="Something Else"), log=lambda *a: None).empty


def test_a_palmitome_sheet_without_the_observed_column_is_refused(tmp_path):
    folder = tmp_path / "post_translation" / P.PALMITOME[0] / P.PALMITOME[1]
    folder.mkdir(parents=True)
    with pd.ExcelWriter(folder / P.PALMITOME[2]) as writer:
        pd.DataFrame({P.PALMITOME_PREDICTED: ["PF3D7_0100100"]}).to_excel(
            writer, sheet_name=P.PALMITOME_SHEET, index=False)
    assert P.palmitome(str(tmp_path), log=lambda *a: None).empty


def test_a_palmitome_with_no_plasmodium_accessions_is_refused(tmp_path):
    assert P.palmitome(_palm_root(tmp_path, observed=["CONTAM", None, None]),
                       log=lambda *a: None).empty


@pytest.mark.skipif(not os.path.exists(PALM), reason="palmitome not fetched")
def test_the_real_palmitome_is_the_size_a_palmitome_should_be():
    d = P.palmitome("datasets" if os.path.isdir("datasets") else
                    os.path.join(ROOT, "datasets"), log=lambda *a: None)
    assert 300 < len(d) < 800, f"{len(d)} proteins: published palmitomes are 400-500"


@pytest.mark.skipif(not os.path.exists(NODES), reason="Plasmodium table not built")
def test_palmitoylation_lands_on_the_proteins_it_should():
    """GAP45 and CDPK1 are the canonical Plasmodium palmitoylation substrates, and the modification
    anchors proteins to membranes -- so membrane proteins have to be enriched among them."""
    from scipy.stats import fisher_exact
    n = pd.read_parquet(NODES, columns=["gene_id", "n_tm", "is_palmitoylated"]).set_index("gene_id")
    assert n.loc["PF3D7_1222700", "is_palmitoylated"], "GAP45"
    assert n.loc["PF3D7_0217500", "is_palmitoylated"], "CDPK1"
    tm = n["n_tm"].fillna(0) > 0
    p = n["is_palmitoylated"]
    odds, pvalue = fisher_exact([[int((p & tm).sum()), int((p & ~tm).sum())],
                                 [int((~p & tm).sum()), int((~p & ~tm).sum())]])
    assert odds > 1.5 and pvalue < 1e-6, (odds, pvalue)


# --------------------------------------------------------------------------- chromatin perturbation
SIR2 = os.path.join(ROOT, "datasets", "reference", "plasmodb", P.SIR2_TABLE)


def _sir2_report(tmp_path, header=True):
    rows = {"Gene ID": ["PF3D7_0100100", "PF3D7_0100200"]}
    for label, _column in P.SIR2:
        rows[f"{label} (Sir2 KO Marray)"] = ["2.7", "N/A"]
    frame = pd.DataFrame(rows)
    if not header:
        frame = frame.drop(columns=["Gene ID"])
    path = tmp_path / P.SIR2_TABLE
    frame.to_csv(path, sep="\t", index=False)
    return str(path)


def test_every_sir2_condition_becomes_its_own_column(tmp_path):
    d = P.sir2_perturbation(_sir2_report(tmp_path))
    for _label, column in P.SIR2:
        assert column in d.columns, column


def test_the_knockout_and_wild_type_stay_separate_columns(tmp_path):
    """The contrast is not precomputed, so both arms have to survive into the table."""
    d = P.sir2_perturbation(_sir2_report(tmp_path))
    assert any(c.startswith("sir2_wt_") for c in d.columns)
    assert any(c.startswith("sir2a_ko_") for c in d.columns)
    assert any(c.startswith("sir2b_ko_") for c in d.columns)


def test_a_missing_sir2_report_yields_nothing(tmp_path):
    assert P.sir2_perturbation(str(tmp_path / "absent.tsv")).empty


def test_a_sir2_report_without_a_gene_column_is_refused(tmp_path):
    assert P.sir2_perturbation(_sir2_report(tmp_path, header=False)).empty


def test_a_sir2_report_with_no_recognised_condition_is_refused(tmp_path):
    path = tmp_path / P.SIR2_TABLE
    pd.DataFrame({"Gene ID": ["PF3D7_0100100"], "unrelated": ["1"]}).to_csv(
        path, sep="\t", index=False)
    assert P.sir2_perturbation(str(path)).empty


@pytest.mark.skipif(not os.path.exists(SIR2), reason="Sir2 report not fetched")
def test_the_sir2_arrays_are_on_a_comparable_scale():
    """Differencing them is only meaningful if they are, and the module leaves that to the caller.

    This asserts the precondition rather than the conclusion: the arrays are log intensities whose
    medians line up, which is what makes a knockout-minus-wild-type difference interpretable at all.
    """
    d = P.sir2_perturbation(SIR2)
    values = d[[c for c in d.columns if c != "gene_id"]]
    medians = values.median()
    assert medians.max() - medians.min() < 0.1, medians.to_dict()
    assert values.min().min() > 1 and values.max().max() < 20


# --------------------------------------------------------------------------- fold confidence
AF = os.path.join(ROOT, "datasets", "reference", "plasmodb", P.ALPHAFOLD_TABLE)


def _af_report(tmp_path, header=True):
    rows = {"gene_id": ["PF3D7_0100100", "PF3D7_0100200"], "uniprot_used": ["Q1", "Q2"]}
    for source, _column in P.ALPHAFOLD:
        rows[source] = ["70.5", "40.0"] if source == "globalMetricValue" else ["0.25", "0.25"]
    frame = pd.DataFrame(rows)
    if not header:
        frame = frame.drop(columns=["gene_id"])
    path = tmp_path / P.ALPHAFOLD_TABLE
    frame.to_csv(path, sep="\t", index=False)
    return str(path)


def test_the_alphafold_summary_becomes_a_mean_and_a_shape(tmp_path):
    d = P.alphafold(_af_report(tmp_path)).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "mean_plddt"] == pytest.approx(70.5)
    for _source, column in P.ALPHAFOLD:
        assert column in d.columns


def test_the_accession_a_model_came_from_is_recorded(tmp_path):
    """A gene can carry eight UniProt accessions; a number should name the structure it came from."""
    d = P.alphafold(_af_report(tmp_path)).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "alphafold_accession"] == "Q1"


def test_a_missing_alphafold_report_yields_nothing(tmp_path):
    assert P.alphafold(str(tmp_path / "absent.tsv")).empty


def test_an_alphafold_report_without_a_gene_column_is_refused(tmp_path):
    assert P.alphafold(_af_report(tmp_path, header=False)).empty


def test_an_alphafold_report_without_the_metric_is_refused(tmp_path):
    path = tmp_path / P.ALPHAFOLD_TABLE
    pd.DataFrame({"gene_id": ["PF3D7_0100100"], "uniprot_used": ["Q1"]}).to_csv(
        path, sep="\t", index=False)
    assert P.alphafold(str(path)).empty


@pytest.mark.skipif(not os.path.exists(AF), reason="AlphaFold summaries not fetched")
def test_the_confidence_fractions_are_internally_consistent():
    """Four fractions of one model have to account for all of it."""
    d = P.alphafold(AF)
    fractions = d[[c for c in d.columns if c.startswith("plddt_fraction_")]]
    assert len(fractions.columns) == 4
    total = fractions.sum(axis=1)
    assert (total - 1.0).abs().max() < 0.01


@pytest.mark.skipif(not os.path.exists(NODES), reason="Plasmodium table not built")
def test_proteins_with_a_recognised_domain_are_modelled_more_confidently():
    """A domain is a thing that folds, so this ordering has to hold whatever the numbers are."""
    from scipy.stats import mannwhitneyu
    n = pd.read_parquet(NODES, columns=["mean_plddt", "has_domain"]).dropna(subset=["mean_plddt"])
    with_domain = n.loc[n["has_domain"] == True, "mean_plddt"]      # noqa: E712
    without = n.loc[n["has_domain"] == False, "mean_plddt"]         # noqa: E712
    assert with_domain.median() > without.median() + 5
    assert mannwhitneyu(with_domain, without).pvalue < 1e-20


# --------------------------------------------------------------------------- myristoylome
MYR = os.path.join(ROOT, "datasets", "post_translation", P.MYRISTOYLOME[0], P.MYRISTOYLOME[1],
                   P.MYRISTOYLOME[2])


def _myr_root(tmp_path, rows=None):
    rows = rows if rows is not None else [
        ("PF3D7_0100100.1-p1", "+", -2.4),    # substrate: significant AND depleted
        ("PF3D7_0100200.1-p1", "+", 2.4),     # significant the wrong way
        ("PF3D7_0100300.1-p1", None, -3.0),   # depleted but not significant
        ("CONTAM;OTHER", "+", -2.0),          # not a Plasmodium gene
    ]
    folder = tmp_path / "post_translation" / P.MYRISTOYLOME[0] / P.MYRISTOYLOME[1]
    folder.mkdir(parents=True)
    pd.DataFrame(rows, columns=["Protein IDs", "Significant", "Difference"]).to_excel(
        folder / P.MYRISTOYLOME[2], index=False)
    return str(tmp_path)


def test_a_substrate_must_be_significant_and_depleted(tmp_path):
    """Depletion when the transferase is blocked is the evidence. Enrichment is not."""
    d = P.myristoylome(_myr_root(tmp_path), log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "is_myristoylated"]
    assert not d.loc["PF3D7_0100200", "is_myristoylated"], "wrong direction was accepted"
    assert not d.loc["PF3D7_0100300", "is_myristoylated"]


def test_assayed_and_unassayed_are_different_states(tmp_path):
    """A protein in the pulldown and not a substrate is a result; one absent from it is not."""
    d = P.myristoylome(_myr_root(tmp_path), log=lambda *a: None)
    assert set(d["gene_id"]) == {"PF3D7_0100100", "PF3D7_0100200", "PF3D7_0100300"}


def test_a_missing_myristoylome_yields_nothing(tmp_path):
    assert P.myristoylome(str(tmp_path), log=lambda *a: None).empty


def test_a_myristoylome_without_the_expected_columns_is_refused(tmp_path):
    folder = tmp_path / "post_translation" / P.MYRISTOYLOME[0] / P.MYRISTOYLOME[1]
    folder.mkdir(parents=True)
    pd.DataFrame({"wrong": [1]}).to_excel(folder / P.MYRISTOYLOME[2], index=False)
    assert P.myristoylome(str(tmp_path), log=lambda *a: None).empty


def test_a_myristoylome_with_no_plasmodium_accessions_is_refused(tmp_path):
    assert P.myristoylome(_myr_root(tmp_path, rows=[("CONTAM", "+", -2.0)]),
                          log=lambda *a: None).empty


@pytest.mark.skipif(not os.path.exists(MYR), reason="myristoylome not fetched")
def test_the_myristoylome_is_the_families_that_are_myristoylated():
    """ARF and Rab GTPases, GAP45, ARO, CDPK1 and the ISP family are the known substrates."""
    d = P.myristoylome(os.path.join(ROOT, "datasets"), log=lambda *a: None).set_index("gene_id")
    for gene, name in (("PF3D7_1222700", "GAP45"), ("PF3D7_0414900", "ARO"),
                       ("PF3D7_0217500", "CDPK1"), ("PF3D7_1020900", "ARF1")):
        assert d.loc[gene, "is_myristoylated"], name
    assert 5 < int(d["is_myristoylated"].sum()) < 40


@pytest.mark.skipif(not os.path.exists(NODES), reason="Plasmodium table not built")
def test_only_the_assayed_proteins_carry_a_myristoylation_call():
    n = pd.read_parquet(NODES, columns=["is_myristoylated"])["is_myristoylated"]
    assert n.notna().sum() < 1000, "the call spread beyond the 609 proteins in the pulldown"
    assert int((n == True).sum()) < 50                            # noqa: E712


# --------------------------------------------------------------------------- acetylome
ACET = os.path.join(ROOT, "datasets", "post_translation", P.ACETYLOME[0], P.ACETYLOME[1],
                    P.ACETYLOME[2])


def _acet_root(tmp_path, rows=None, sheet=None):
    rows = rows if rows is not None else [
        ("Acetyl", "AAA", "PF3D7_0100100", "K10", 1.0),
        ("Acetyl", "AAA", "PF3D7_0100100", "K20", 0.9),
        ("Acetyl", "AAA", "PF3D7_0100100", "K20", 0.9),   # the same site twice
        ("Acetyl", "AAA", "PF3D7_0100200", "K5", 0.30),   # identified, not localised
        ("Acetyl", "AAA", "CONTAM", "K1", 1.0),
    ]
    folder = tmp_path / "post_translation" / P.ACETYLOME[0] / P.ACETYLOME[1]
    folder.mkdir(parents=True)
    body = pd.DataFrame(rows, columns=["Modification", "Surrounding Sequence", "Accession Number",
                                       "Site", "Ascore Localization Probability"])
    with pd.ExcelWriter(folder / P.ACETYLOME[2].replace(".xls", ".xlsx")) as writer:
        title = pd.DataFrame([["Acetyl-lysine sites in Plasmodium"]])
        title.to_excel(writer, sheet_name=sheet or P.ACETYLOME_SHEET, index=False, header=False)
        body.to_excel(writer, sheet_name=sheet or P.ACETYLOME_SHEET, index=False, startrow=1)
    os.rename(folder / P.ACETYLOME[2].replace(".xls", ".xlsx"), folder / P.ACETYLOME[2])
    return str(tmp_path)


def test_a_gene_seen_acetylated_is_flagged_even_if_no_site_is_localised(tmp_path):
    """Identifying an acetylated peptide and localising the acetyl group are different claims."""
    d = P.acetylome(_acet_root(tmp_path), log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100200", "has_acetyl"]
    assert pd.isna(d.loc["PF3D7_0100200", "n_acetylsites"])


def test_the_site_count_uses_localised_sites_only_and_counts_each_once(tmp_path):
    d = P.acetylome(_acet_root(tmp_path), log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "n_acetylsites"] == 2


def test_non_plasmodium_accessions_are_dropped_from_the_acetylome(tmp_path):
    d = P.acetylome(_acet_root(tmp_path), log=lambda *a: None)
    assert d["gene_id"].str.startswith("PF3D7_").all()


def test_a_missing_acetylome_yields_nothing(tmp_path):
    assert P.acetylome(str(tmp_path), log=lambda *a: None).empty


def test_an_acetylome_without_the_expected_sheet_is_refused(tmp_path):
    assert P.acetylome(_acet_root(tmp_path, sheet="Other"), log=lambda *a: None).empty


def test_an_acetylome_with_no_plasmodium_accessions_is_refused(tmp_path):
    assert P.acetylome(_acet_root(tmp_path, rows=[("Acetyl", "A", "CONTAM", "K1", 1.0)]),
                       log=lambda *a: None).empty


@pytest.mark.skipif(not os.path.exists(ACET), reason="acetylome not fetched")
def test_the_acetylome_is_led_by_the_acetylation_machinery():
    """Chromatin readers and writers are the most heavily acetylated proteins in any acetylome."""
    d = P.acetylome(os.path.join(ROOT, "datasets"), log=lambda *a: None)
    n = pd.read_parquet(NODES, columns=["gene_id", "product"])
    m = d.merge(n, on="gene_id", how="left")
    assert 800 < len(m) < 2000
    top = " ".join(m.nlargest(10, "n_acetylsites")["product"].astype(str)).lower()
    assert "phd" in top or "acetyltransferase" in top or "histone" in top
    assert m["product"].str.contains("histone", case=False, na=False).sum() >= 5


def test_an_acetylome_sheet_without_the_expected_columns_is_refused(tmp_path):
    folder = tmp_path / "post_translation" / P.ACETYLOME[0] / P.ACETYLOME[1]
    folder.mkdir(parents=True)
    with pd.ExcelWriter(folder / "tmp.xlsx") as writer:
        pd.DataFrame([["title"]]).to_excel(writer, sheet_name=P.ACETYLOME_SHEET, index=False,
                                           header=False)
        pd.DataFrame({"wrong": ["x"]}).to_excel(writer, sheet_name=P.ACETYLOME_SHEET, index=False,
                                                startrow=1)
    os.rename(folder / "tmp.xlsx", folder / P.ACETYLOME[2])
    assert P.acetylome(str(tmp_path), log=lambda *a: None).empty


def test_build_all_splits_the_acetyl_flag_from_the_acetyl_count(tmp_path):
    root = _dataset_root(tmp_path)
    folder = os.path.join(root, "post_translation", P.ACETYLOME[0], P.ACETYLOME[1])
    os.makedirs(folder)
    body = pd.DataFrame([("Acetyl", "A", "PF3D7_0100100", "K10", 1.0)],
                        columns=["Modification", "Surrounding Sequence", "Accession Number",
                                 "Site", "Ascore Localization Probability"])
    tmp = os.path.join(folder, "tmp.xlsx")
    with pd.ExcelWriter(tmp) as writer:
        pd.DataFrame([["title"]]).to_excel(writer, sheet_name=P.ACETYLOME_SHEET, index=False,
                                           header=False)
        body.to_excel(writer, sheet_name=P.ACETYLOME_SHEET, index=False, startrow=1)
    os.rename(tmp, os.path.join(folder, P.ACETYLOME[2]))
    d = P.build_all(root, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "has_acetyl"] and d.loc["PF3D7_0100100", "n_acetylsites"] == 1
    assert not d.loc["PF3D7_0100200", "has_acetyl"]
    assert pd.isna(d.loc["PF3D7_0100200", "n_acetylsites"])


# --------------------------------------------------------------------------- strain identity
NF54 = os.path.join(ROOT, "datasets", "reference", "plasmodb", P.NF54_TABLE)


def _nf54_report(tmp_path, rows=None):
    rows = rows if rows is not None else [
        ("PfNF54_A", "OG6_1", "100"),       # one on each side, same length -> maps
        ("PfNF54_B", "OG6_2", "100"),       # group has two 3D7 genes -> dropped
        ("PfNF54_C", "OG6_3", "500"),       # length disagrees wildly -> dropped
        ("PfNF54_D", "OG6_4", None),        # no length -> kept, nothing to contradict
        ("PfNF54_E", "NOT_OG", "100"),      # not an orthogroup
    ]
    path = tmp_path / P.NF54_TABLE
    pd.DataFrame(rows, columns=["Gene ID", "Ortholog Group", "Protein Length"]).to_csv(
        path, sep="\t", index=False)
    return str(path)


def _three_d7():
    return pd.DataFrame({
        "gene_id": ["PF3D7_A", "PF3D7_B1", "PF3D7_B2", "PF3D7_C", "PF3D7_D"],
        "orthogroup": ["OG6_1", "OG6_2", "OG6_2", "OG6_3", "OG6_4"],
        "length": [100.0, 100.0, 100.0, 100.0, 250.0]})


def test_a_one_to_one_orthogroup_maps_and_a_family_does_not(tmp_path):
    m = P.strain_map(_nf54_report(tmp_path), _three_d7())
    assert m["PfNF54_A"] == "PF3D7_A"
    assert "PfNF54_B" not in m, "a group with two 3D7 genes was guessed at"


def test_a_pair_whose_proteins_are_different_sizes_is_dropped(tmp_path):
    """Orthology is a claim about ancestry. This needs a claim about identity."""
    m = P.strain_map(_nf54_report(tmp_path), _three_d7())
    assert "PfNF54_C" not in m


def test_a_pair_with_no_length_to_check_is_kept(tmp_path):
    """Unknown is not a contradiction; refusing it would drop genes for missing metadata."""
    m = P.strain_map(_nf54_report(tmp_path), _three_d7())
    assert m.get("PfNF54_D") == "PF3D7_D"


def test_a_missing_or_empty_strain_report_maps_nothing(tmp_path):
    assert P.strain_map(str(tmp_path / "absent.tsv"), _three_d7()) == {}
    assert P.strain_map(_nf54_report(tmp_path), pd.DataFrame()) == {}


def test_a_strain_report_with_no_shared_orthogroups_maps_nothing(tmp_path):
    other = _three_d7().assign(orthogroup="OG6_999")
    assert P.strain_map(_nf54_report(tmp_path, rows=[("PfNF54_A", "OG6_1", "100")]), other) == {}


@pytest.mark.skipif(not (os.path.exists(NF54) and os.path.exists(NODES)),
                    reason="strain report not fetched")
def test_the_strain_map_really_is_an_identity_map():
    """The check that licenses using orthology as identity: 3D7 was cloned from NF54, so the paired
    proteins have to be the same size. If a PlasmoDB release ever breaks that, this fails first."""
    nodes = pd.read_parquet(NODES)
    m = P.strain_map(NF54, nodes)
    assert len(m) > 3000
    nf = pd.read_csv(NF54, sep="\t", dtype=str).replace({"N/A": None})
    nf.columns = ["nf54", "orthogroup", "length"]
    nf["length"] = pd.to_numeric(nf["length"], errors="coerce")
    lengths = dict(zip(nodes["gene_id"], nodes["length"]))
    pairs = [(a, lengths.get(b)) for a, b in m.items()]
    nf_len = dict(zip(nf["nf54"], nf["length"]))
    same = [abs(nf_len[a] - b) < 1 for a, b in pairs
            if b == b and nf_len.get(a) == nf_len.get(a) and nf_len.get(a) is not None]
    assert sum(same) / len(same) > 0.9, f"only {100*sum(same)/len(same):.0f}% of pairs match in length"


# --------------------------------------------------------------------------- lactylome
def test_lactylation_needs_a_strain_map(tmp_path):
    assert P.lactylome(str(tmp_path), {}, log=lambda *a: None).empty


def test_a_missing_lactylome_yields_nothing(tmp_path):
    assert P.lactylome(str(tmp_path), {"PfNF54_A": "PF3D7_A"}, log=lambda *a: None).empty


@pytest.mark.skipif(not os.path.exists(NODES), reason="Plasmodium table not built")
def test_lactylation_landed_on_histones():
    """Histone lactylation is the paper's headline, so histones have to be in the result."""
    n = pd.read_parquet(NODES, columns=["product", "has_lactyl", "n_lactylsites"])
    assert 50 < int(n["has_lactyl"].sum()) < 400
    marked = n[n["has_lactyl"]]
    assert marked["product"].str.contains("histone", case=False, na=False).sum() >= 3


def _lac_root(tmp_path, sheet=None, cols=True):
    folder = tmp_path / "post_translation" / P.LACTYLOME[0] / P.LACTYLOME[1]
    folder.mkdir(parents=True)
    body = (pd.DataFrame({"Gene ID": ["PfNF54_A;PfNF54_X", "PfNF54_A", "PfNF54_B", "PfNF54_Z"],
                          "Localization prob": [1.0, 1.0, 0.10, 1.0],
                          "Position within protein ": [10, 20, 5, 1]})
            if cols else pd.DataFrame({"wrong": [1]}))
    body.to_excel(folder / P.LACTYLOME[2], sheet_name=sheet or P.LACTYLOME_SHEET, index=False)
    return str(tmp_path)


LMAP = {"PfNF54_A": "PF3D7_A", "PfNF54_B": "PF3D7_B"}


def test_lactylation_sites_are_counted_per_resolved_gene(tmp_path):
    d = P.lactylome(_lac_root(tmp_path), LMAP, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_A", "n_lactylsites"] == 2


def test_a_gene_seen_lactylated_without_a_localised_site_is_still_flagged(tmp_path):
    d = P.lactylome(_lac_root(tmp_path), LMAP, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_B", "has_lactyl"] and pd.isna(d.loc["PF3D7_B", "n_lactylsites"])


def test_an_unmappable_strain_accession_is_dropped(tmp_path):
    d = P.lactylome(_lac_root(tmp_path), LMAP, log=lambda *a: None)
    assert set(d["gene_id"]) == {"PF3D7_A", "PF3D7_B"}


def test_a_lactylome_without_the_expected_sheet_or_columns_is_refused(tmp_path):
    assert P.lactylome(_lac_root(tmp_path, sheet="Other"), LMAP, log=lambda *a: None).empty
    assert P.lactylome(_lac_root(tmp_path / "b", cols=False), LMAP, log=lambda *a: None).empty


def test_a_lactylome_whose_accessions_all_fail_to_map_is_refused(tmp_path):
    assert P.lactylome(_lac_root(tmp_path), {"PfNF54_Q": "PF3D7_Q"}, log=lambda *a: None).empty


def test_build_all_resolves_lactylation_through_the_strain_map(tmp_path):
    root = _dataset_root(tmp_path)
    base = os.path.join(root, "reference", "plasmodb")
    pd.DataFrame([("PfNF54_A", "OG6_104345", "2163")],
                 columns=["Gene ID", "Ortholog Group", "Protein Length"]).to_csv(
        os.path.join(base, P.NF54_TABLE), sep="\t", index=False)
    folder = os.path.join(root, "post_translation", P.LACTYLOME[0], P.LACTYLOME[1])
    os.makedirs(folder)
    pd.DataFrame({"Gene ID": ["PfNF54_A"], "Localization prob": [1.0],
                  "Position within protein ": [10]}).to_excel(
        os.path.join(folder, P.LACTYLOME[2]), sheet_name=P.LACTYLOME_SHEET, index=False)
    d = P.build_all(root, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "has_lactyl"], "the NF54 accession did not resolve"
    assert not d.loc["PF3D7_0100200", "has_lactyl"]


def test_a_strain_report_missing_a_column_maps_nothing(tmp_path):
    """The report needs an accession, an orthogroup and a length; two of three is not enough."""
    path = tmp_path / P.NF54_TABLE
    pd.DataFrame({"Gene ID": ["PfNF54_A"], "Ortholog Group": ["OG6_1"]}).to_csv(
        path, sep="\t", index=False)
    assert P.strain_map(str(path), _three_d7()) == {}


# --------------------------------------------------------------------------- isoforms
ISO = os.path.join(ROOT, "datasets", "transcription", P.ISOFORMS[0], P.ISOFORMS[1], P.ISOFORMS[2])


def _iso_root(tmp_path, rows=None, sheet=None, cols=True):
    rows = rows if rows is not None else [
        ("PF3D7_0100100", "full-splice_match"),
        ("PF3D7_0100100", "novel_in_catalog"),
        ("PF3D7_0100200", "full-splice_match"),
        ("PF3D7_0100300_novel_gene_1", "novel_not_in_catalog"),
        ("novelGene_123", "genic"),
    ]
    folder = tmp_path / "transcription" / P.ISOFORMS[0] / P.ISOFORMS[1]
    folder.mkdir(parents=True)
    body = (pd.DataFrame(rows, columns=["associated_gene", "structural_category"])
            if cols else pd.DataFrame({"wrong": [1]}))
    body.to_excel(folder / P.ISOFORMS[2], sheet_name=sheet or P.ISOFORM_SHEET, index=False)
    return str(tmp_path)


def test_transcript_models_are_counted_per_gene(tmp_path):
    d = P.isoforms(_iso_root(tmp_path), log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "n_transcript_models"] == 2
    assert d.loc["PF3D7_0100200", "n_transcript_models"] == 1


def test_only_models_the_annotation_lacks_count_as_novel(tmp_path):
    """`full-splice_match` is the reference transcript recovered, which is not a discovery."""
    d = P.isoforms(_iso_root(tmp_path), log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "novel_transcript_models"] == 1
    assert d.loc["PF3D7_0100200", "novel_transcript_models"] == 0


def test_a_novel_gene_suffix_is_stripped_back_to_the_gene(tmp_path):
    d = P.isoforms(_iso_root(tmp_path), log=lambda *a: None)
    assert "PF3D7_0100300" in set(d["gene_id"])


def test_models_not_anchored_on_a_gene_are_dropped(tmp_path):
    d = P.isoforms(_iso_root(tmp_path), log=lambda *a: None)
    assert d["gene_id"].str.startswith("PF3D7_").all()


def test_a_missing_isoform_table_yields_nothing(tmp_path):
    assert P.isoforms(str(tmp_path), log=lambda *a: None).empty


def test_an_isoform_table_without_the_sheet_or_columns_is_refused(tmp_path):
    assert P.isoforms(_iso_root(tmp_path, sheet="Other"), log=lambda *a: None).empty
    assert P.isoforms(_iso_root(tmp_path / "b", cols=False), log=lambda *a: None).empty


def test_an_isoform_table_with_no_plasmodium_genes_is_refused(tmp_path):
    assert P.isoforms(_iso_root(tmp_path, rows=[("novelGene_1", "genic")]),
                      log=lambda *a: None).empty


@pytest.mark.skipif(not os.path.exists(NODES), reason="Plasmodium table not built")
def test_a_gene_with_no_long_read_model_stays_missing_rather_than_reading_as_one():
    """Absence is sequencing depth, not a statement that the gene has a single transcript."""
    d = pd.read_parquet(NODES, columns=["n_transcript_models", "novel_transcript_models"])
    assert d["n_transcript_models"].isna().sum() > 3000
    assert d["n_transcript_models"].min() >= 1


def test_build_all_folds_in_the_transcript_models(tmp_path):
    root = _dataset_root(tmp_path)
    folder = os.path.join(root, "transcription", P.ISOFORMS[0], P.ISOFORMS[1])
    os.makedirs(folder)
    pd.DataFrame([("PF3D7_0100100", "novel_in_catalog")],
                 columns=["associated_gene", "structural_category"]).to_excel(
        os.path.join(folder, P.ISOFORMS[2]), sheet_name=P.ISOFORM_SHEET, index=False)
    d = P.build_all(root, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "novel_transcript_models"] == 1
    assert pd.isna(d.loc["PF3D7_0100200", "n_transcript_models"])


# --------------------------------------------------------------------------- derived labels
def _stage_frame(n=80, seed=0):
    rng = np.random.default_rng(seed)
    frame = pd.DataFrame({"gene_id": [f"PF3D7_{i:06d}" for i in range(n)]})
    for cols in P.PF_STAGE_COLUMNS.values():
        for c in cols:
            frame[c] = rng.lognormal(size=n)
    return frame


def test_the_maximum_across_stages_is_on_a_log_scale():
    frame = _stage_frame()
    d = P.derived_labels(frame, log=lambda *a: None)
    expected = np.log2(frame[[c for cols in P.PF_STAGE_COLUMNS.values() for c in cols]].max(axis=1) + 1)
    assert np.allclose(d["expr_max"], expected)


def test_the_stage_call_comes_with_its_margin():
    d = P.derived_labels(_stage_frame(), log=lambda *a: None)
    assert {"stage_enriched_derived", "stage_margin_derived"} <= set(d.columns)
    called = d["stage_enriched_derived"].notna()
    assert d.loc[called, "stage_margin_derived"].notna().all()
    assert d.loc[~called, "stage_margin_derived"].isna().all()


def test_a_gene_with_no_clear_winner_is_left_unlabelled():
    """A label that is really a coin toss looks like a measurement in every table it reaches."""
    flat = _stage_frame(n=60)
    for cols in P.PF_STAGE_COLUMNS.values():
        for c in cols:
            flat[c] = 1.0
    d = P.derived_labels(flat, log=lambda *a: None)
    assert d["stage_enriched_derived"].isna().all()


def test_derived_labels_need_at_least_two_stage_columns():
    frame = pd.DataFrame({"gene_id": ["PF3D7_000001"], "expr_ring": [1.0]})
    assert P.derived_labels(frame, log=lambda *a: None).empty or \
        "expr_max" not in P.derived_labels(frame, log=lambda *a: None).columns


def test_the_two_arms_label_stages_by_the_same_construction():
    """`cellcycle.stage_enrichment` is shared on purpose, so a difference means biology."""
    from starplast import cellcycle
    frame = _stage_frame()
    mine = cellcycle.stage_enrichment(frame, log=lambda *a: None, stages=P.PF_STAGE_COLUMNS)
    assert "stage_enriched_derived" in mine.columns
    # And the Toxoplasma default is untouched by passing a different map.
    assert set(cellcycle.STAGE_COLUMNS) == {"tachyzoite", "bradyzoite", "oocyst"}


@pytest.mark.skipif(not os.path.exists(NODES), reason="Plasmodium table not built")
def test_the_shipped_stage_calls_are_conservative():
    d = pd.read_parquet(NODES, columns=["expr_max", "stage_enriched_derived"])
    assert d["expr_max"].notna().all()
    called = d["stage_enriched_derived"].notna().sum()
    assert 50 < called < 2000, f"{called} calls: the margin rule is not doing its job"


# --------------------------------------------------------------------------- antibody epitopes
IEDB = os.path.join(ROOT, "datasets", "reference", "plasmodb", P.IEDB_TABLE)


def _iedb_root(tmp_path, epitopes=None, uniprot=None):
    base = tmp_path / "reference" / "plasmodb"
    base.mkdir(parents=True)
    pd.DataFrame(uniprot if uniprot is not None else
                 [("PF3D7_A", "Q1"), ("PF3D7_B", "Q2,Q3"), ("PF3D7_C", "Q3")],
                 columns=["Gene ID", "UniProt ID(s)"]).to_csv(
        base / P.UNIPROT_TABLE, sep="\t", index=False)
    pd.DataFrame(epitopes if epitopes is not None else
                 [("Q1", "AAA"), ("Q1", "AAA"), ("Q1", "BBB"), ("Q3", "CCC"), ("Q9", "DDD")],
                 columns=["uniprot", "epitope"]).to_csv(
        base / P.IEDB_TABLE, sep="\t", index=False)
    return str(tmp_path)


def test_an_accession_naming_two_genes_is_dropped_from_the_index(tmp_path):
    """An epitope belongs to a protein; picking whichever paralogue sorted first invents it."""
    idx = P.uniprot_index(os.path.join(_iedb_root(tmp_path), "reference", "plasmodb",
                                       P.UNIPROT_TABLE))
    assert idx["Q1"] == "PF3D7_A" and idx["Q2"] == "PF3D7_B"
    assert "Q3" not in idx, "an accession shared by two genes was assigned to one"


def test_distinct_epitopes_are_counted_not_assay_records(tmp_path):
    d = P.bcell_epitopes(_iedb_root(tmp_path), log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_A", P.IEDB_COLUMN] == 2


def test_an_unmappable_antigen_is_dropped(tmp_path):
    d = P.bcell_epitopes(_iedb_root(tmp_path), log=lambda *a: None)
    assert set(d["gene_id"]) == {"PF3D7_A"}


def test_epitopes_need_both_tables(tmp_path):
    base = tmp_path / "reference" / "plasmodb"
    base.mkdir(parents=True)
    pd.DataFrame([("Q1", "AAA")], columns=["uniprot", "epitope"]).to_csv(
        base / P.IEDB_TABLE, sep="\t", index=False)
    assert P.bcell_epitopes(str(tmp_path), log=lambda *a: None).empty
    assert P.uniprot_index(str(base / "absent.tsv")) == {}


def test_an_epitope_table_without_the_expected_columns_is_refused(tmp_path):
    root = _iedb_root(tmp_path)
    pd.DataFrame({"wrong": ["x"]}).to_csv(
        os.path.join(root, "reference", "plasmodb", P.IEDB_TABLE), sep="\t", index=False)
    assert P.bcell_epitopes(root, log=lambda *a: None).empty


def test_a_uniprot_table_with_one_column_maps_nothing(tmp_path):
    base = tmp_path / "reference" / "plasmodb"
    base.mkdir(parents=True)
    pd.DataFrame({"Gene ID": ["PF3D7_A"]}).to_csv(base / P.UNIPROT_TABLE, sep="\t", index=False)
    assert P.uniprot_index(str(base / P.UNIPROT_TABLE)) == {}


def test_an_epitope_table_whose_antigens_all_fail_to_map_is_refused(tmp_path):
    root = _iedb_root(tmp_path, epitopes=[("Q9", "AAA")])
    assert P.bcell_epitopes(root, log=lambda *a: None).empty


@pytest.mark.skipif(not os.path.exists(IEDB), reason="IEDB epitopes not fetched")
def test_the_antigens_are_the_ones_malaria_serology_has_studied():
    """MSP1 and CSP are the two most studied antigens in the organism's history."""
    d = P.bcell_epitopes(os.path.join(ROOT, "datasets"), log=lambda *a: None).set_index("gene_id")
    assert d[P.IEDB_COLUMN].idxmax() == "PF3D7_0930300", "MSP1 is not the top antigen"
    assert "PF3D7_0304600" in d.index, "CSP has no antibody epitope"
    assert 200 < len(d) < 1500


def test_build_all_folds_in_the_antibody_epitopes(tmp_path):
    root = _dataset_root(tmp_path)
    base = os.path.join(root, "reference", "plasmodb")
    pd.DataFrame([("PF3D7_0100100", "Q1")], columns=["Gene ID", "UniProt ID(s)"]).to_csv(
        os.path.join(base, P.UNIPROT_TABLE), sep="\t", index=False)
    pd.DataFrame([("Q1", "AAA")], columns=["uniprot", "epitope"]).to_csv(
        os.path.join(base, P.IEDB_TABLE), sep="\t", index=False)
    d = P.build_all(root, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", P.IEDB_COLUMN] == 1
    assert pd.isna(d.loc["PF3D7_0100200", P.IEDB_COLUMN]), "absent was read as zero"


def test_the_two_iedb_halves_are_separate_columns(tmp_path):
    """Antibodies and T cells are different questions and get different slots."""
    root = _iedb_root(tmp_path)
    pd.DataFrame([("Q1", "TTT"), ("Q1", "TTT")], columns=["uniprot", "epitope"]).to_csv(
        os.path.join(root, "reference", "plasmodb", P.IEDB_TABLES["n_tcell_epitopes"]),
        sep="\t", index=False)
    d = P.bcell_epitopes(root, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_A", "n_bcell_epitopes"] == 2
    assert d.loc["PF3D7_A", "n_tcell_epitopes"] == 1


def test_one_iedb_half_alone_is_enough(tmp_path):
    """The T-cell fetch failing must not cost the antibody column."""
    d = P.bcell_epitopes(_iedb_root(tmp_path), log=lambda *a: None)
    assert "n_bcell_epitopes" in d.columns and "n_tcell_epitopes" not in d.columns


@pytest.mark.skipif(not os.path.exists(os.path.join(
    ROOT, "datasets", "reference", "plasmodb", P.IEDB_TABLES["n_tcell_epitopes"])),
    reason="T-cell epitopes not fetched")
def test_the_two_halves_are_not_interchangeable():
    """434 antigens carry an antibody epitope and only 44 a T-cell one; pooling would blur that."""
    d = P.bcell_epitopes(os.path.join(ROOT, "datasets"), log=lambda *a: None)
    b = d["n_bcell_epitopes"].notna().sum()
    t = d["n_tcell_epitopes"].notna().sum()
    assert b > t * 3, (b, t)
    assert "PF3D7_0304600" in set(d.loc[d["n_tcell_epitopes"].notna(), "gene_id"]), "CSP"


# --------------------------------------------------------------------------- febrile stress
FEB = os.path.join(ROOT, "datasets", "reference", "plasmodb", P.FEBRILE_TABLE)


def _feb_report(tmp_path, header=True):
    rows = {"Gene ID": ["PF3D7_0100100", "PF3D7_0100200"]}
    for label, _column in P.FEBRILE:
        rows[f"sense - {label} - unique only"] = ["10.0", "N/A"]
    frame = pd.DataFrame(rows)
    if not header:
        frame = frame.drop(columns=["Gene ID"])
    path = tmp_path / P.FEBRILE_TABLE
    frame.to_csv(path, sep="\t", index=False)
    return str(path)


def test_every_febrile_condition_becomes_its_own_column(tmp_path):
    d = P.febrile(_feb_report(tmp_path))
    for _label, column in P.FEBRILE:
        assert column in d.columns, column


def test_both_temperatures_and_all_three_lines_survive(tmp_path):
    """The 41-versus-37 contrast is the caller's to make, so both arms have to be there."""
    d = P.febrile(_feb_report(tmp_path))
    assert {"febrile_wt_37c", "febrile_wt_41c"} <= set(d.columns)
    assert sum(c.endswith("_41c") for c in d.columns) == 3


def test_a_missing_febrile_report_yields_nothing(tmp_path):
    assert P.febrile(str(tmp_path / "absent.tsv")).empty


def test_a_febrile_report_without_a_gene_column_is_refused(tmp_path):
    assert P.febrile(_feb_report(tmp_path, header=False)).empty


def test_a_febrile_report_with_no_recognised_condition_is_refused(tmp_path):
    path = tmp_path / P.FEBRILE_TABLE
    pd.DataFrame({"Gene ID": ["PF3D7_0100100"], "unrelated": ["1"]}).to_csv(
        path, sep="\t", index=False)
    assert P.febrile(str(path)).empty


@pytest.mark.skipif(not os.path.exists(FEB), reason="febrile report not fetched")
def test_the_febrile_arms_are_on_a_comparable_scale():
    """Asserting the precondition rather than the conclusion, as with the Sir2 arrays."""
    d = P.febrile(FEB)
    values = d[[c for c in d.columns if c != "gene_id"]]
    medians = values.median()
    assert medians.min() > 5 and medians.max() < 200, medians.to_dict()
    assert values.min().min() >= 0


def test_build_all_folds_in_the_febrile_conditions(tmp_path):
    root = _dataset_root(tmp_path)
    base = os.path.join(root, "reference", "plasmodb")
    rows = {"Gene ID": ["PF3D7_0100100", "PF3D7_0100200"]}
    for label, _column in P.FEBRILE:
        rows[f"sense - {label} - unique only"] = ["10.0", "12.0"]
    pd.DataFrame(rows).to_csv(os.path.join(base, P.FEBRILE_TABLE), sep="\t", index=False)
    d = P.build_all(root, log=lambda *a: None)
    assert "febrile_wt_41c" in d.columns


# --------------------------------------------------------------------------- complexes
CPLX = os.path.join(ROOT, "datasets", "reference", "plasmodb", P.COMPLEXES[0], P.COMPLEXES[1],
                    P.COMPLEXES[2])


def _cplx_root(tmp_path, rows=None, sheet=None):
    rows = rows if rows is not None else [
        ("Q1", 1, "Plasmodium falciparum (isolate 3D7)", "Pathogen-only", 2),
        ("Q2", 1, "Plasmodium falciparum (isolate 3D7)", "Pathogen-only", 2),
        ("Q4", 2, "Plasmodium falciparum (isolate 3D7)", "Host-Pathogen", 3),
        ("P99", 2, "Homo sapiens", "Host-Pathogen", 3),
    ]
    base = tmp_path / "reference" / "plasmodb"
    (base / P.COMPLEXES[0] / P.COMPLEXES[1]).mkdir(parents=True)
    pd.DataFrame([("PF3D7_A", "Q1"), ("PF3D7_B", "Q2"), ("PF3D7_D", "Q4")],
                 columns=["Gene ID", "UniProt ID(s)"]).to_csv(
        base / P.UNIPROT_TABLE, sep="\t", index=False)
    pd.DataFrame(rows, columns=["uniprot", "complex nr. in supp. fig 5", "organism_name",
                                "interaction_type", "csize"]).to_excel(
        base / P.COMPLEXES[0] / P.COMPLEXES[1] / P.COMPLEXES[2],
        sheet_name=sheet or P.COMPLEXES_SHEET, index=False)
    return str(tmp_path)


def test_only_parasite_members_get_a_row(tmp_path):
    d = P.complexes(_cplx_root(tmp_path), log=lambda *a: None)
    assert set(d["gene_id"]) == {"PF3D7_A", "PF3D7_B", "PF3D7_D"}


def test_a_complex_reaching_the_host_is_flagged_rather_than_dropped(tmp_path):
    """Seven of the 47 real clusters contain human proteins. That is a finding, not contamination --
    but a parasite gene in one has partners this table cannot name."""
    d = P.complexes(_cplx_root(tmp_path), log=lambda *a: None).set_index("gene_id")
    assert not d.loc["PF3D7_A", "complex_spans_host"]
    assert d.loc["PF3D7_D", "complex_spans_host"]


def test_complex_size_comes_from_the_source(tmp_path):
    d = P.complexes(_cplx_root(tmp_path), log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_A", "complex_size"] == 2 and d.loc["PF3D7_D", "complex_size"] == 3


def test_a_missing_or_unreadable_complex_table_yields_nothing(tmp_path):
    assert P.complexes(str(tmp_path), log=lambda *a: None).empty
    assert P.complexes(_cplx_root(tmp_path, sheet="Other"), log=lambda *a: None).empty


def test_a_complex_table_with_no_mappable_parasite_member_is_refused(tmp_path):
    root = _cplx_root(tmp_path, rows=[("P99", 1, "Homo sapiens", "Host-only", 2)])
    assert P.complexes(root, log=lambda *a: None).empty


@pytest.mark.skipif(not os.path.exists(CPLX), reason="complexes not fetched")
def test_the_real_complexes_are_mostly_parasite_only():
    d = P.complexes(os.path.join(ROOT, "datasets"), log=lambda *a: None)
    assert 80 < len(d) < 300
    assert d["complex_spans_host"].mean() < 0.5, "most complexes should not reach the host"
    assert d["complex_size"].min() >= 2, "a complex of one is not a complex"


def test_a_complex_sheet_missing_its_columns_is_refused(tmp_path):
    base = tmp_path / "reference" / "plasmodb"
    (base / P.COMPLEXES[0] / P.COMPLEXES[1]).mkdir(parents=True)
    pd.DataFrame({"wrong": [1]}).to_excel(
        base / P.COMPLEXES[0] / P.COMPLEXES[1] / P.COMPLEXES[2],
        sheet_name=P.COMPLEXES_SHEET, index=False)
    assert P.complexes(str(tmp_path), log=lambda *a: None).empty


def test_complexes_need_the_uniprot_index(tmp_path):
    root = _cplx_root(tmp_path)
    os.remove(os.path.join(root, "reference", "plasmodb", P.UNIPROT_TABLE))
    assert P.complexes(root, log=lambda *a: None).empty


def test_build_all_folds_in_the_complexes(tmp_path):
    root = _dataset_root(tmp_path)
    base = os.path.join(root, "reference", "plasmodb")
    pd.DataFrame([("PF3D7_0100100", "Q1")], columns=["Gene ID", "UniProt ID(s)"]).to_csv(
        os.path.join(base, P.UNIPROT_TABLE), sep="\t", index=False)
    folder = os.path.join(base, P.COMPLEXES[0], P.COMPLEXES[1])
    os.makedirs(folder)
    pd.DataFrame([("Q1", 1, "Plasmodium falciparum (isolate 3D7)", "Pathogen-only", 2)],
                 columns=["uniprot", "complex nr. in supp. fig 5", "organism_name",
                          "interaction_type", "csize"]).to_excel(
        os.path.join(folder, P.COMPLEXES[2]), sheet_name=P.COMPLEXES_SHEET, index=False)
    d = P.build_all(root, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "complex_id"] == 1
    assert pd.isna(d.loc["PF3D7_0100200", "complex_id"])


# --------------------------------------------------------------------------- the substring collision
def test_a_sample_name_does_not_match_a_longer_name_containing_it():
    """The bug that deleted three columns by adding three.

    `sense - asexual blood stages` is a SUBSTRING of `antisense - asexual blood stages`. Plain
    containment made each sense entry match two headers, fail the one-match test, and vanish -- so
    fetching antisense removed sense, silently, and the slot count went DOWN by two while a slot was
    being added. The label must start the header or be preceded by a non-alphanumeric character.
    """
    assert P._matches("sense - asexual blood stages - unique only", "sense - asexual blood stages")
    assert not P._matches("antisense - asexual blood stages - unique only",
                          "sense - asexual blood stages")
    assert P._matches("antisense - asexual blood stages - unique only",
                      "antisense - asexual blood stages")
    # And a label that begins mid-header after a separator is still found.
    assert P._matches("Study X - Ring Ave (proteome)", "Ring Ave")
    assert not P._matches("SpringRing Ave", "Ring Ave")


@pytest.mark.skipif(not os.path.exists(EXPR), reason="expression report not fetched")
def test_every_declared_expression_column_survives_the_real_report():
    """The guard that would have caught it: no declared sample may go missing from the built table."""
    d = P.expression(EXPR)
    missing = [c for _study, _sample, c in P.EXPRESSION if c not in d.columns]
    assert not missing, f"declared samples absent from the table: {missing}"


@pytest.mark.skipif(not os.path.exists(NODES), reason="Plasmodium table not built")
def test_antisense_is_the_minority_strand_and_sits_beside_its_sense_partner():
    """Antisense far below sense is the check; both present is the point of not shipping a ratio."""
    d = pd.read_parquet(NODES, columns=["expr_asexual_blood", "antisense_asexual_blood"]).dropna()
    assert d["antisense_asexual_blood"].median() < d["expr_asexual_blood"].median() / 3
    assert len(d) > 4000


# --------------------------------------------------------------------------- codon usage
CDS = os.path.join(ROOT, "starplast", "data", P.CDS_TABLE)


@pytest.mark.skipif(not os.path.exists(CDS), reason="Plasmodium CDS not fetched")
def test_codon_usage_comes_from_the_shared_construction():
    """Not reimplemented. ENC and GC3 are definitions and the CAI reference set is the ribosomal
    proteins in both arms, so two implementations could only differ by being wrong in one."""
    d = P.codon_usage(ROOT, log=lambda *a: None)
    assert set(d.columns) == {"codon_enc", "codon_gc3", "codon_cai_ribosomal"}
    assert len(d) > 4000


@pytest.mark.skipif(not (os.path.exists(NODES) and os.path.exists(
    os.path.join(ROOT, "starplast", "data", "nodes.parquet"))), reason="both caches needed")
def test_the_two_arms_differ_in_codon_bias_the_way_their_genomes_do():
    """The first measurement the two arms can be compared on, and the reason for sharing the code.

    P. falciparum has the most AT-rich genome of any eukaryote, so its GC3 has to be far below
    Toxoplasma's and its ENC far lower -- lower ENC being MORE biased. If these ever converge, either
    a genome was mixed up or the two arms stopped computing the same quantity.
    """
    pf = pd.read_parquet(NODES, columns=["codon_gc3", "codon_enc"])
    tg = pd.read_parquet(os.path.join(ROOT, "starplast", "data", "nodes.parquet"),
                         columns=["codon_gc3", "codon_enc"])
    assert pf["codon_gc3"].median() < 0.30, "Plasmodium GC3 is not AT-rich"
    assert tg["codon_gc3"].median() > 0.45, "Toxoplasma GC3 has moved"
    assert pf["codon_enc"].median() < tg["codon_enc"].median() - 10, (
        "the two arms no longer differ in codon bias")


@pytest.mark.skipif(not os.path.exists(CDS), reason="Plasmodium CDS not fetched")
def test_cai_is_measured_against_plasmodium_ribosomal_proteins():
    """The reference set must come from the Pf node table, not the Toxoplasma one."""
    d = P.codon_usage(ROOT, log=lambda *a: None)
    cai = d["codon_cai_ribosomal"].dropna()
    assert len(cai) > 4000
    assert 0 < cai.min() and cai.max() <= 1.0


# --------------------------------------------------------------------------- enzyme classification
EC = os.path.join(ROOT, "datasets", "reference", "plasmodb", P.EC_TABLE)


def _ec_root(tmp_path, rows=None, drop_id=False):
    rows = rows if rows is not None else [
        ("PF3D7_0100100", "2.7.11.1", None),        # curated
        ("PF3D7_0100200", None, "3.1.3.2"),         # orthology only
        ("PF3D7_0100300", None, None),              # neither
    ]
    base = tmp_path / "reference" / "plasmodb"
    base.mkdir(parents=True)
    frame = pd.DataFrame(rows, columns=["Gene ID", "EC numbers", "EC numbers from OrthoMCL"])
    if drop_id:
        frame = frame.drop(columns=["Gene ID"])
    frame.to_csv(base / P.EC_TABLE, sep="\t", index=False)
    return str(tmp_path)


def test_curated_and_derived_ec_stay_in_separate_columns(tmp_path):
    """Different KINDS of evidence. Merging them would put inference where annotation is, and the
    merged column would give no way to tell which genes were never annotated in this organism."""
    d = P.enzyme_classification(_ec_root(tmp_path), log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "has_ec"]
    assert not d.loc["PF3D7_0100200", "has_ec"], "an orthology-derived EC counted as curated"
    assert d.loc["PF3D7_0100200", "ec_number_orthology"] == "3.1.3.2"


def test_a_gene_with_no_ec_at_all_is_flagged_false(tmp_path):
    d = P.enzyme_classification(_ec_root(tmp_path), log=lambda *a: None).set_index("gene_id")
    assert not d.loc["PF3D7_0100300", "has_ec"]
    assert pd.isna(d.loc["PF3D7_0100300", "ec_number"])


def test_a_missing_or_headerless_ec_table_yields_nothing(tmp_path):
    assert P.enzyme_classification(str(tmp_path), log=lambda *a: None).empty
    assert P.enzyme_classification(_ec_root(tmp_path, drop_id=True), log=lambda *a: None).empty


@pytest.mark.skipif(not os.path.exists(EC), reason="EC table not fetched")
def test_the_enzyme_count_is_the_size_a_proteome_of_this_kind_gives():
    d = P.enzyme_classification(os.path.join(ROOT, "datasets"), log=lambda *a: None)
    assert 800 < int(d["has_ec"].sum()) < 2500
    # Orthology adds genes the curation does not have; if it added none the second column is pointless.
    extra = (~d["has_ec"] & d["ec_number_orthology"].notna()).sum()
    assert extra > 50, f"only {extra} genes gain an EC from orthology"


def test_build_all_folds_in_the_enzyme_classification(tmp_path):
    root = _dataset_root(tmp_path)
    base = os.path.join(root, "reference", "plasmodb")
    pd.DataFrame([("PF3D7_0100100", "2.7.11.1", None)],
                 columns=["Gene ID", "EC numbers", "EC numbers from OrthoMCL"]).to_csv(
        os.path.join(base, P.EC_TABLE), sep="\t", index=False)
    d = P.build_all(root, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "has_ec"] and not d.loc["PF3D7_0100200", "has_ec"]


# --------------------------------------------------------------------------- gene identity
def _identity(tmp_path, rows=None):
    """A PlasmoDB identity table, in the shape `fetch_names` writes."""
    rows = rows if rows is not None else [
        ("PF3D7_0102200", "X", "Previous IDs: PFA0110w;MAL1P4.09", "a protein"),
        ("PF3D7_0102300", "Y", "Previous IDs: PF13_0222", "another protein"),
    ]
    path = tmp_path / "plasmodb_identity.tsv"
    pd.DataFrame(rows, columns=["gene_id", "gene_name", "previous_ids", "product"]).to_csv(
        path, sep="\t", index=False)
    return str(path)


def test_a_previous_accession_resolves_to_the_current_gene(tmp_path):
    index = P.previous_id_index(_identity(tmp_path))
    assert index["pfa0110w"] == "PF3D7_0102200"
    assert index["mal1p4.09"] == "PF3D7_0102200"
    assert index["pf13_0222"] == "PF3D7_0102300"


def test_an_old_id_claimed_by_two_genes_resolves_to_neither(tmp_path):
    """66 real ones, and every one is a gene model that was SPLIT, seen from the other side.

    A measurement made against the old model belongs to both halves and to neither in particular.
    Assigning it to whichever row came first would put real numbers on an arbitrary gene, which is
    the failure the Toxoplasma layer withdraws 153 strings to avoid.
    """
    index = P.previous_id_index(_identity(tmp_path, rows=[
        ("PF3D7_0109850", "A", "Previous IDs: PFA0485w", "half one"),
        ("PF3D7_0109950", "B", "Previous IDs: PFA0485w", "half two"),
        ("PF3D7_0102300", "C", "Previous IDs: PF13_0222", "unambiguous"),
    ]))
    assert "pfa0485w" not in index
    assert index["pf13_0222"] == "PF3D7_0102300"


def test_a_current_accession_outranks_another_genes_history(tmp_path):
    """Files mix both forms, and a live id must not be captured by some other gene's past."""
    index = P.previous_id_index(_identity(tmp_path, rows=[
        ("PF3D7_0102200", "X", "Previous IDs: PF3D7_0102300", "renamed onto a live id"),
        ("PF3D7_0102300", "Y", "N/A", "the live gene itself"),
    ]))
    assert index["pf3d7_0102300"] == "PF3D7_0102300"


def test_a_missing_or_headerless_identity_table_is_an_empty_index(tmp_path):
    assert P.previous_id_index(str(tmp_path / "absent.tsv")) == {}
    path = tmp_path / "wrong.tsv"
    pd.DataFrame({"something": ["else"]}).to_csv(path, sep="\t", index=False)
    assert P.previous_id_index(str(path)) == {}


# --------------------------------------------------------------------------- ribosome profiling
def _riboseq_dir(tmp_path, tables=None):
    """The deposit's per-sample RPKM files, keyed the way 2014 keyed them."""
    folder = tmp_path / "translation" / P.RIBOSEQ[0] / P.RIBOSEQ[1]
    folder.mkdir(parents=True)
    tables = tables if tables is not None else {
        "GSM1410291_mRNA_1_rpkm.txt.gz": {"pfa0110w": 100.0, "pf13_0222": 10.0},
        "GSM1410296_ribosome_footprints_1_rpkm.txt.gz": {"pfa0110w": 400.0, "pf13_0222": 5.0},
    }
    for name, values in tables.items():
        pd.DataFrame(list(values.items())).to_csv(
            folder / name, sep="\t", header=False, index=False, compression="gzip")
    return folder


def test_riboseq_resolves_pre_2012_accessions_and_names_its_arms(tmp_path, monkeypatch):
    """The join that found zero rows before the identity layer existed."""
    _riboseq_dir(tmp_path)
    monkeypatch.setattr(__import__("starplast.paths", fromlist=["x"]),
                        "cache_file", lambda name: _identity(tmp_path))
    d = P.riboseq(str(tmp_path), log=lambda *a: None).set_index("gene_id")
    assert list(d.columns) == ["riboseq_mrna_ring", "riboseq_rpf_ring"]
    assert d.loc["PF3D7_0102200", "riboseq_rpf_ring"] == 400.0
    assert d.loc["PF3D7_0102300", "riboseq_mrna_ring"] == 10.0


def test_riboseq_withdraws_two_source_ids_landing_on_one_gene(tmp_path, monkeypatch):
    """The deposit splits some genes into `-a` and `-b` segments, and merges are the same problem.

    RPKM is already length-normalised, so neither summing two segments nor averaging them means
    anything. Both rows are dropped and the gene keeps no value for that stage.
    """
    _riboseq_dir(tmp_path, tables={"GSM1410291_mRNA_1_rpkm.txt.gz": {
        "pfa0110w": 100.0, "mal1p4.09": 60.0, "pf13_0222": 10.0}})
    monkeypatch.setattr(__import__("starplast.paths", fromlist=["x"]),
                        "cache_file", lambda name: _identity(tmp_path))
    d = P.riboseq(str(tmp_path), log=lambda *a: None).set_index("gene_id")
    assert "PF3D7_0102200" not in d.index, "two ids for one gene were averaged rather than withdrawn"
    assert d.loc["PF3D7_0102300", "riboseq_mrna_ring"] == 10.0


def test_riboseq_without_its_directory_or_its_index_is_empty(tmp_path, monkeypatch):
    assert P.riboseq(str(tmp_path), log=lambda *a: None).empty
    _riboseq_dir(tmp_path)
    monkeypatch.setattr(__import__("starplast.paths", fromlist=["x"]),
                        "cache_file", lambda name: str(tmp_path / "absent.tsv"))
    said = []
    assert P.riboseq(str(tmp_path), log=said.append).empty
    assert any("fetch_names" in m for m in said), "a missing identity table must say how to get one"


def test_riboseq_ignores_files_that_are_not_the_per_gene_tables(tmp_path, monkeypatch):
    """The deposit also carries coverage tracks and UTR beds; only the RPKM tables are read."""
    folder = _riboseq_dir(tmp_path)
    (folder / "GSM1410291_mRNA_1_minus.wig.gz").write_bytes(b"not a table")
    (folder / "filelist.txt").write_text("Archive/File\tName\n")
    monkeypatch.setattr(__import__("starplast.paths", fromlist=["x"]),
                        "cache_file", lambda name: _identity(tmp_path))
    d = P.riboseq(str(tmp_path), log=lambda *a: None)
    assert list(d.columns) == ["gene_id", "riboseq_mrna_ring", "riboseq_rpf_ring"]


def test_build_all_leaves_an_unmeasured_gene_missing_for_ribosome_profiling(tmp_path, monkeypatch):
    """61% coverage, uneven across the cycle. A gap is silence, not an absence of ribosomes."""
    root = _dataset_root(tmp_path)
    _riboseq_dir(tmp_path, tables={"GSM1410296_ribosome_footprints_1_rpkm.txt.gz": {
        "pfa0110w": 400.0}})
    ident = _identity(tmp_path, rows=[("PF3D7_0100100", "X", "Previous IDs: PFA0110w", "p")])
    monkeypatch.setattr(__import__("starplast.paths", fromlist=["x"]),
                        "cache_file", lambda name: ident)
    d = P.build_all(root, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "riboseq_rpf_ring"] == 400.0
    assert pd.isna(d.loc["PF3D7_0100200", "riboseq_rpf_ring"])


# --------------------------------------------------------------------------- what the data says
RIBO_STAGES = ("ring", "early_trophozoite", "late_trophozoite", "schizont", "merozoite")


def _ribosomal(nodes):
    """Ribosomal proteins by product annotation, whole-word, in whichever arm's table."""
    return nodes["product"].fillna("").str.contains(r"\bribosomal protein\b", case=False, regex=True)


@pytest.mark.skipif(not os.path.exists(NODES), reason="Plasmodium table not built")
def test_ribosome_footprints_pile_up_on_ribosomal_proteins_at_every_stage():
    """The check that says the measurement behaves, and it is an ordering rather than a total.

    Ribosomal proteins are the most heavily translated things in a growing cell. If this came out
    the other way the footprints would be measuring something else, and no amount of coverage would
    make the column worth shipping.
    """
    d = pd.read_parquet(NODES)
    ribo = _ribosomal(d)
    for stage in RIBO_STAGES:
        v = np.log1p(d[f"riboseq_rpf_{stage}"])
        seen = v.notna()
        assert v[seen & ribo].median() > v[seen & ~ribo].median() + 1.0, (
            f"{stage}: ribosomal proteins do not carry more footprint than the rest")


@pytest.mark.skipif(not os.path.exists(NODES), reason="Plasmodium table not built")
def test_the_riboseq_stage_labels_agree_with_the_independent_stage_series():
    """Labels are the deposit's own, so they are checked against a study that shares no sample.

    Four of the five stages have a counterpart in the PlasmoDB seven-stage series, and each has to
    correlate highest with its own. The merozoite has no counterpart there and is deliberately not
    asserted: its best match is the ring, which is the neighbouring point of the cycle -- a free
    merozoite is a ring that has not invaded yet -- and calling that a failure would be reading a
    missing column as a wrong label.
    """
    d = pd.read_parquet(NODES)
    for stage in RIBO_STAGES[:4]:
        mine = np.log1p(d[f"riboseq_mrna_{stage}"])
        theirs = {c: np.log1p(d[c]) for c in d.columns
                  if c.startswith("expr_") and c not in ("expr_max", "expr_asexual_blood")}
        scores = {c: mine.corr(v, method="spearman") for c, v in theirs.items()}
        best = max(scores, key=scores.get)
        assert best == f"expr_{stage}", f"{stage} looks most like {best} (rho {scores[best]:.2f})"


@pytest.mark.skipif(not os.path.exists(NODES), reason="Plasmodium table not built")
def test_translation_efficiency_is_not_a_column_on_either_arm():
    """Computed twice, from two instruments, and refused twice. The conditions ship; the ratio does not.

    The polysome arms gave a translation efficiency that correlated NEGATIVELY with codon adaptation
    where the textbook expects positive, and ribosome profiling -- a different instrument, a
    different strain and a decade earlier -- reproduces it (rho -0.01 to -0.12), while also putting
    ribosomal proteins BELOW the rest at the schizont. Two checks disagreeing is the contradictory
    case, whose answer is to ship the conditions rather than the contrast.
    """
    d = pd.read_parquet(NODES)
    assert not [c for c in d.columns if "translation_efficiency" in c or c.startswith("riboseq_te")]


@pytest.mark.skipif(not os.path.exists(NODES), reason="Plasmodium table not built")
def test_the_codon_index_separates_ribosomal_proteins_in_one_arm_and_not_the_other():
    """Where the contradiction above actually lives, found because the construction is SHARED.

    `codons.codon_usage` builds CAI against each arm's own ribosomal proteins, so the reference set
    should score at the top of its own index. In Toxoplasma it does. In *P. falciparum* it does not
    separate at all -- which says the quantity failing to behave in the translation-efficiency check
    is the codon index, not the footprints, and it is a property of an AT-rich genome rather than a
    bug in either arm. One implementation makes that comparable; two would have made it noise.
    """
    pf = pd.read_parquet(NODES)
    tg_path = os.path.join(ROOT, "starplast", "data", "nodes.parquet")
    if not os.path.exists(tg_path):
        pytest.skip("Toxoplasma table not built")
    tg = pd.read_parquet(tg_path)
    tg_gap = (tg.loc[_ribosomal(tg), "codon_cai_ribosomal"].median()
              - tg.loc[~_ribosomal(tg), "codon_cai_ribosomal"].median())
    pf_gap = (pf.loc[_ribosomal(pf), "codon_cai_ribosomal"].median()
              - pf.loc[~_ribosomal(pf), "codon_cai_ribosomal"].median())
    assert tg_gap > 0.03, "the Toxoplasma index stopped ranking its own reference set at the top"
    assert abs(pf_gap) < 0.01, "the Plasmodium index now separates: re-read the TE refusal"


# --------------------------------------------------------------------------- extracellular vesicles
def _secretome_book(tmp_path, rows=None):
    """The paper's compilation sheet: a union of two EV preparations, with its reference list."""
    folder = tmp_path / "post_translation" / P.SECRETOME[0] / P.SECRETOME[1]
    folder.mkdir(parents=True)
    rows = rows if rows is not None else [
        ("PFA0110w", 1, 1), ("PF13_0222", 1, 0), ("PF3D7_0102500", 0, 1),
        ("Mantel PY et al: Malaria-infected erythrocyte-derived microvesicles. "
         "Cell Host Microbe 2013, 13(5):521-534.", None, None)]
    pd.DataFrame(rows, columns=["geneid", *P.SECRETOME_STUDIES]).to_excel(
        folder / P.SECRETOME[2], sheet_name=P.SECRETOME_SHEET, index=False)
    return folder


def test_the_vesicle_sheet_counts_preparations_and_drops_its_own_references(tmp_path, monkeypatch):
    """A citation carries accessions, so resolution goes through the index rather than a regex."""
    _secretome_book(tmp_path)
    monkeypatch.setattr(__import__("starplast.paths", fromlist=["x"]), "cache_file",
                        lambda name: _identity(tmp_path, rows=[
                            ("PF3D7_0102200", "X", "Previous IDs: PFA0110w", "p"),
                            ("PF3D7_0102300", "Y", "Previous IDs: PF13_0222", "q"),
                            ("PF3D7_0102500", "Z", "N/A", "r")]))
    d = P.secretome(str(tmp_path), log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0102200", "ev_studies"] == 2
    assert d.loc["PF3D7_0102300", "ev_studies"] == 1
    assert len(d) == 3, "the reference line was read as a gene"


def test_a_missing_vesicle_file_sheet_or_column_is_an_empty_frame(tmp_path, monkeypatch):
    monkeypatch.setattr(__import__("starplast.paths", fromlist=["x"]), "cache_file",
                        lambda name: _identity(tmp_path))
    assert P.secretome(str(tmp_path), log=lambda *a: None).empty
    folder = _secretome_book(tmp_path, rows=[("PFA0110w", 1, 1)])
    pd.DataFrame({"geneid": ["PFA0110w"]}).to_excel(
        folder / P.SECRETOME[2], sheet_name="something else", index=False)
    assert P.secretome(str(tmp_path), log=lambda *a: None).empty


def test_a_gene_no_vesicle_preparation_reported_stays_unknown(tmp_path, monkeypatch):
    """The column that is deliberately NOT here is a False flag for the other 5,536 genes.

    The other mass-spectrometry columns on this arm do complete a boolean, because they pool
    proteome-wide assays where "never observed" is an answer. This is one preparation from one
    isolate, and the genes it did not report are six times less expressed than the ones it did, so a
    False would be a detection limit written down as a negative result -- and the atlas would have
    graded the slot A at 100% coverage for an experiment that identified 184 proteins.
    """
    root = _dataset_root(tmp_path)
    _secretome_book(tmp_path, rows=[("PFA0110w", 1, 1)])
    monkeypatch.setattr(__import__("starplast.paths", fromlist=["x"]), "cache_file",
                        lambda name: _identity(tmp_path, rows=[
                            ("PF3D7_0100100", "X", "Previous IDs: PFA0110w", "p")]))
    d = P.build_all(root, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "ev_studies"] == 2
    assert pd.isna(d.loc["PF3D7_0100200", "ev_studies"])
    assert not [c for c in d.columns if c.startswith("in_extracellular")]


@pytest.mark.skipif(not os.path.exists(NODES), reason="Plasmodium table not built")
def test_the_vesicle_proteome_leans_secretory_and_leans_abundant():
    """Both directions are asserted, because only one of them is a reason to trust the column.

    A secretome should carry signal peptides, and this one does. It is ALSO six times more expressed
    than the rest of the proteome, which is what mass spectrometry on a vesicle preparation returns
    and is why the column is named for detection rather than for secretion. Pinning the confound
    keeps it from quietly disappearing out of the prose.
    """
    d = pd.read_parquet(NODES)
    ev = d["ev_studies"].notna()
    assert d.loc[ev, "has_signal_peptide"].mean() > d.loc[~ev, "has_signal_peptide"].mean() * 1.5
    assert d.loc[ev, "expr_asexual_blood"].median() > d.loc[~ev, "expr_asexual_blood"].median() * 2


def test_a_riboseq_directory_holding_only_tracks_yields_nothing(tmp_path, monkeypatch):
    """The deposit's wiggle files are 40 MB and are not per-gene tables. Neither is filelist.txt."""
    folder = tmp_path / "translation" / P.RIBOSEQ[0] / P.RIBOSEQ[1]
    folder.mkdir(parents=True)
    (folder / "GSM1410291_mRNA_1_minus.wig.gz").write_bytes(b"not a table")
    monkeypatch.setattr(__import__("starplast.paths", fromlist=["x"]),
                        "cache_file", lambda name: _identity(tmp_path))
    assert P.riboseq(str(tmp_path), log=lambda *a: None).empty


def test_a_vesicle_sheet_without_its_study_columns_is_refused(tmp_path):
    """The membership columns ARE the measurement here; a sheet without them is a different table."""
    folder = tmp_path / "post_translation" / P.SECRETOME[0] / P.SECRETOME[1]
    folder.mkdir(parents=True)
    pd.DataFrame({"geneid": ["PF3D7_0102200"], "something else": [1]}).to_excel(
        folder / P.SECRETOME[2], sheet_name=P.SECRETOME_SHEET, index=False)
    assert P.secretome(str(tmp_path), log=lambda *a: None).empty


def test_the_vesicle_sheet_needs_the_identity_index_and_says_which_one(tmp_path, monkeypatch):
    _secretome_book(tmp_path)
    monkeypatch.setattr(__import__("starplast.paths", fromlist=["x"]),
                        "cache_file", lambda name: str(tmp_path / "absent.tsv"))
    said = []
    assert P.secretome(str(tmp_path), log=said.append).empty
    assert any("fetch_names" in m for m in said)


def test_a_vesicle_sheet_whose_ids_resolve_to_nothing_is_empty_rather_than_wrong(tmp_path,
                                                                                monkeypatch):
    """A join that matches nothing looks exactly like a dataset with no coverage, so it says so."""
    _secretome_book(tmp_path, rows=[("not_an_accession", 1, 1)])
    monkeypatch.setattr(__import__("starplast.paths", fromlist=["x"]),
                        "cache_file", lambda name: _identity(tmp_path))
    assert P.secretome(str(tmp_path), log=lambda *a: None).empty


# --------------------------------------------------------------------------- kinase substrates
def _kinase_book(tmp_path, windows=None):
    """The supplement's hypophosphorylation sheet, title row and all."""
    folder = tmp_path / "post_translation" / P.KINASE_SUBSTRATE[0] / P.KINASE_SUBSTRATE[1]
    folder.mkdir(parents=True)
    windows = windows if windows is not None else ["MAFKKMAF", "GGWWGG"]
    frame = pd.DataFrame({"Protein Description": ["x"] * len(windows),
                          P.KINASE_SUBSTRATE_WINDOW: windows})
    with pd.ExcelWriter(folder / P.KINASE_SUBSTRATE[2]) as writer:
        # Two title rows above the header, which is where the file keeps them.
        pd.DataFrame([["Supplementary Data 2a"], [None]]).to_excel(
            writer, sheet_name=P.KINASE_SUBSTRATE_SHEET, header=False, index=False)
        frame.to_excel(writer, sheet_name=P.KINASE_SUBSTRATE_SHEET, index=False,
                       startrow=P.KINASE_SUBSTRATE_HEADER)
    return folder


def _cds(tmp_path, sequences):
    path = tmp_path / P.CDS_TABLE
    pd.DataFrame({"gene_id": list(sequences), "cds": list(sequences.values())}).to_csv(
        path, sep="\t", index=False)
    return str(path)


def test_a_site_window_that_names_one_protein_is_counted(tmp_path, monkeypatch):
    """The only usable key in the file: its accessions are from an annotation nobody carries."""
    _kinase_book(tmp_path, windows=["MAFKK", "MAFKK", "GGWW"])
    monkeypatch.setattr(__import__("starplast.paths", fromlist=["x"]), "cache_file",
                        lambda name: _cds(tmp_path, {
                            "PF3D7_0100100": "ATGGCTTTTAAAAAA",     # MAFKK
                            "PF3D7_0100200": "GGAGGATGGTGGAAA"}))   # GGWWK
    d = P.kinase_substrates(str(tmp_path), log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "cdpk1_dependent_sites"] == 2
    assert d.loc["PF3D7_0100200", "cdpk1_dependent_sites"] == 1


def test_a_window_in_two_proteins_is_dropped_rather_than_assigned(tmp_path, monkeypatch):
    """Same rule as the accession index: a key that names two genes names neither."""
    _kinase_book(tmp_path, windows=["MAFKK", "GGWW"])
    monkeypatch.setattr(__import__("starplast.paths", fromlist=["x"]), "cache_file",
                        lambda name: _cds(tmp_path, {
                            "PF3D7_0100100": "ATGGCTTTTAAAAAA",
                            "PF3D7_0100200": "ATGGCTTTTAAAAAA",     # the same window, second gene
                            "PF3D7_0100300": "GGAGGATGGTGGAAA"}))
    said = []
    d = P.kinase_substrates(str(tmp_path), log=said.append).set_index("gene_id")
    assert list(d.index) == ["PF3D7_0100300"]
    assert any("1 matched more than one" in m for m in said), "the cost has to be reported"


def test_kinase_substrates_needs_the_file_the_sheet_the_column_and_the_proteome(tmp_path,
                                                                                monkeypatch):
    monkeypatch.setattr(__import__("starplast.paths", fromlist=["x"]), "cache_file",
                        lambda name: _cds(tmp_path, {"PF3D7_0100100": "ATGGCTTTTAAAAAA"}))
    assert P.kinase_substrates(str(tmp_path), log=lambda *a: None).empty
    folder = _kinase_book(tmp_path, windows=["NOTINANYPROTEIN"])
    assert P.kinase_substrates(str(tmp_path), log=lambda *a: None).empty, "no window matched"
    pd.DataFrame({"other": [1]}).to_excel(folder / P.KINASE_SUBSTRATE[2],
                                          sheet_name=P.KINASE_SUBSTRATE_SHEET, index=False)
    said = []
    assert P.kinase_substrates(str(tmp_path), log=said.append).empty
    assert any("not laid out" in m for m in said), (
        "a re-issued supplement with a different layout must report, not raise")
    pd.DataFrame({"other": [1] * 5}).to_excel(folder / P.KINASE_SUBSTRATE[2],
                                          sheet_name="another sheet", index=False)
    assert P.kinase_substrates(str(tmp_path), log=lambda *a: None).empty
    _kinase_book(tmp_path.parent / "no_cds", windows=["MAFKK"])
    monkeypatch.setattr(__import__("starplast.paths", fromlist=["x"]), "cache_file",
                        lambda name: str(tmp_path / "absent.tsv.gz"))
    said = []
    assert P.kinase_substrates(str(tmp_path.parent / "no_cds"), log=said.append).empty
    assert any("CDS table" in m for m in said)


@pytest.mark.skipif(not os.path.exists(NODES), reason="Plasmodium table not built")
def test_the_cdpk1_set_recovers_the_papers_own_conclusion():
    """The mapping is by sequence, so it is checked against something the sequence did not decide.

    The paper's finding is that CDPK1 signalling runs through the invasion motor. The genes this
    loader recovers are nine times more likely to be motor, IMC or invasion machinery than the
    proteome at large -- reached from the other end, through a 15-residue window and a translation.
    """
    d = pd.read_parquet(NODES)
    machinery = d["product"].fillna("").str.contains(
        r"myosin|glideosome|inner membrane complex|\bimc\b|actin|gap45|gap50|mtip", case=False,
        regex=True)
    hit = d["cdpk1_dependent_sites"].notna()
    assert hit.sum() > 50
    assert machinery[hit].mean() > machinery.mean() * 4


def test_a_kinase_sheet_without_the_window_column_is_refused(tmp_path):
    """The window IS the key here, so a sheet that lost it is not a version of this table."""
    folder = tmp_path / "post_translation" / P.KINASE_SUBSTRATE[0] / P.KINASE_SUBSTRATE[1]
    folder.mkdir(parents=True)
    pd.DataFrame({"a": range(6), "b": range(6)}).to_excel(
        folder / P.KINASE_SUBSTRATE[2], sheet_name=P.KINASE_SUBSTRATE_SHEET, index=False)
    assert P.kinase_substrates(str(tmp_path), log=lambda *a: None).empty


def test_build_all_folds_in_the_kinase_dependent_sites(tmp_path, monkeypatch):
    root = _dataset_root(tmp_path)
    _kinase_book(tmp_path, windows=["MAFKK"])
    monkeypatch.setattr(__import__("starplast.paths", fromlist=["x"]), "cache_file",
                        lambda name: (_cds(tmp_path, {"PF3D7_0100100": "ATGGCTTTTAAAAAA"})
                                      if str(name).endswith(".gz") else _identity(tmp_path)))
    d = P.build_all(root, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "cdpk1_dependent_sites"] == 1
    assert pd.isna(d.loc["PF3D7_0100200", "cdpk1_dependent_sites"])


def test_build_all_records_host_degree_where_the_crosslink_data_reaches(tmp_path):
    """The one line of the assembly that needed a crosslink file to run, and now has one.

    Three states, not two: a gene seen in the crosslink experiment gets a count -- zero included,
    since being crosslinked only to parasite proteins is an observation -- and a gene the experiment
    never saw stays missing rather than being called zero.
    """
    from starplast import pf_graph as G
    root = _dataset_root(tmp_path)
    base = tmp_path / "reference" / "plasmodb"
    (base / G.CROSSLINK[0] / G.CROSSLINK[1]).mkdir(parents=True)
    pd.DataFrame([("PF3D7_0100100", "Q1")], columns=["Gene ID", "UniProt ID(s)"]).to_csv(
        base / P.UNIPROT_TABLE, sep="\t", index=False)
    pd.DataFrame([("sp|Q1|Parasite", "sp|P02768|HumanAlbumin", 3)],
                 columns=["Protein1", "Protein2", "Num_Crosslinks"]).to_excel(
        base / G.CROSSLINK[0] / G.CROSSLINK[1] / G.CROSSLINK[2],
        sheet_name=G.CROSSLINK_SHEET, index=False)
    d = P.build_all(root, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "n_host_targets"] == 1
    assert pd.isna(d.loc["PF3D7_0100200", "n_host_targets"])


# --------------------------------------------------------------------------- the IP-MS interactome
class _Page:
    """One page of a supplementary PDF, as the loader consumes it."""

    def __init__(self, text):
        self._text = text

    def extract_text(self):
        return self._text


def _pdf(tmp_path, pages, monkeypatch):
    """Stand in for the published PDF: the loader only ever asks a page for its text."""
    folder = tmp_path / "reference" / "plasmodb" / P.IP_MS[0] / P.IP_MS[1]
    folder.mkdir(parents=True)
    (folder / P.IP_MS[2]).write_bytes(b"%PDF-1.4 stub")
    import pypdf
    monkeypatch.setattr(pypdf, "PdfReader",
                        lambda path: type("R", (), {"pages": [_Page(t) for t in pages]})())
    return str(tmp_path)


_PV1_PAGE = ("Gene ID Annotation Number of significant ms/ms spectra Experiment 1 Experiment 2 "
             "Avg. fold change PV1-HA WT PV1-HA WT 1 2 1 2 1 2 1 2 "
             "PF3D7_1129100 parasitophorous vacuolar protein 1 (PV1) 97 120 0 0 94 105 0 0 >104 "
             "PF3D7_1024800 exported protein 3 (EXP3) 90 107 0 0 68 77 0 0 >85.5")


def test_the_bait_is_read_from_the_column_header_and_must_top_its_own_table(tmp_path, monkeypatch):
    """The source's captions are wrong -- `parasite interacting proteins` heads a table of human
    ones -- so the bait comes from the header and the assignment is checked against the data: a
    pulldown's own bait is its first and strongest row."""
    root = _pdf(tmp_path, [_PV1_PAGE], monkeypatch)
    d = P.ip_ms(root, log=lambda *a: None)
    assert list(d["bait"]) == ["PF3D7_1129100"]
    assert list(d["prey"]) == ["PF3D7_1024800"], "the bait's own row is not an edge"


def test_a_table_that_does_not_start_with_its_bait_is_refused(tmp_path, monkeypatch):
    """Attributing one protein's partners to another is the failure this check exists to prevent."""
    swapped = _PV1_PAGE.replace("PF3D7_1129100 parasitophorous vacuolar protein 1 (PV1)",
                                "PF3D7_0000001 something else entirely")
    root = _pdf(tmp_path, [swapped], monkeypatch)
    said = []
    assert P.ip_ms(root, log=said.append).empty
    assert any("does not start with" in m for m in said)


def test_the_eight_counts_are_two_experiments_of_four_not_two_halves(tmp_path, monkeypatch):
    """The trap that shipped for ten minutes: bait, bait, control, control -- twice.

    Split down the middle it compares experiment 1 with experiment 2, and reports EXP3 as an
    enriched partner carrying 145 spectra in the untagged control. Read correctly it is 342 against
    zero, which is what an enrichment table looks like.
    """
    root = _pdf(tmp_path, [_PV1_PAGE], monkeypatch)
    row = P.ip_ms(root, log=lambda *a: None).iloc[0]
    assert row["spectra_bait"] == 90 + 107 + 68 + 77
    assert row["spectra_control"] == 0


def test_a_row_with_seven_counts_is_dropped_and_counted(tmp_path, monkeypatch):
    """The missing number could be either arm, so the row is refused rather than shifted.

    Its annotation ends in `(PIESP2)` in the file, which is what makes the anchored match refuse it:
    the eight counts have to be the LAST eight numbers in the record. An earlier version also
    required the annotation to end in a non-digit, and that refused `heat shock protein DnaJ
    homologue, Pfj2` -- a complete row whose name simply ends in a 2.
    """
    # Copied from the file: seven counts, and an annotation ending in a digit-bearing name, which
    # is the pair of properties that let this row look complete.
    short = (_PV1_PAGE + " PF3D7_0501200 parasite-infected erythrocyte surface protein 2 (PIESP2) "
             "0 0 0 3 2 0 0 >1.67")
    root = _pdf(tmp_path, [short], monkeypatch)
    said = []
    d = P.ip_ms(root, log=said.append)
    assert "PF3D7_0501200" not in set(d["prey"])
    assert any("wrong number of counts" in m for m in said)


def test_the_transgene_bait_contributes_no_pairs_and_says_how_many(tmp_path, monkeypatch):
    """PfEMP1B is built from a var gene and has no accession in the table; 38 real rows cannot be
    paired, and counting them is the difference between a limitation and a silence."""
    page = ("Gene ID Annotation Experiment 1 Experiment 2 PfEMP1B PfEMP1F PfEMP1B PfEMP1F "
            "PF3D7_1024800 exported protein 3 (EXP3) 5 6 0 0 7 8 0 0 >6 "
            "PF3D7_0917900 heat shock protein 70 (HSP70-2) 3 4 0 0 5 6 0 0 >4")
    root = _pdf(tmp_path, [page], monkeypatch)
    said = []
    assert P.ip_ms(root, log=said.append).empty
    assert any("2 rows from the transgene bait" in m for m in said)


def test_a_missing_interactome_pdf_is_an_empty_frame(tmp_path):
    assert P.ip_ms(str(tmp_path), log=lambda *a: None).empty


@pytest.mark.skipif(not os.path.exists(NODES), reason="Plasmodium table not built")
def test_no_shipped_ip_ms_pair_is_richer_in_its_control_than_in_its_bait():
    """The check that caught the column-order fault, kept against the shipped file."""
    import starplast.plasmodium as PL
    root = os.path.join(ROOT, "datasets")
    pairs = PL.ip_ms(root, log=lambda *a: None)
    if pairs.empty:
        pytest.skip("the interactome PDF is not on this machine")
    assert (pairs["spectra_bait"] > pairs["spectra_control"]).all()
    assert pairs["spectra_control"].median() == 0


def test_an_annotation_ending_in_a_digit_is_still_a_complete_row(tmp_path, monkeypatch):
    """`Pfj2` is a protein name, not a count. The regression the strict version would have caused."""
    page = ("Gene ID PV1-HA WT PV1-HA WT "
            "PF3D7_1129100 parasitophorous vacuolar protein 1 (PV1) 97 120 0 0 94 105 0 0 >104 "
            "PF3D7_1108700 heat shock protein DnaJ homologue, Pfj2 0 0 0 0 3 6 0 0 >2.25")
    root = _pdf(tmp_path, [page], monkeypatch)
    d = P.ip_ms(root, log=lambda *a: None)
    assert "PF3D7_1108700" in set(d["prey"])
    assert int(d.set_index("prey").loc["PF3D7_1108700", "spectra_bait"]) == 3 + 6


def test_build_all_counts_ip_ms_partners_both_ways_and_leaves_the_rest_missing(tmp_path,
                                                                               monkeypatch):
    """Degree over undirected pairs: a bait counts its partners, a partner counts its baits.

    Missing everywhere else, because four pulldowns are not a survey of the proteome and a zero
    would say `nothing binds this` about a protein nobody tested.
    """
    root = _dataset_root(tmp_path)
    monkeypatch.setattr(P, "ip_ms", lambda *a, **k: pd.DataFrame(
        [("PF3D7_0100100", "PF3D7_0100200", "x", 40, 0)],
        columns=["bait", "prey", "product", "spectra_bait", "spectra_control"]))
    d = P.build_all(root, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "n_ip_ms_partners"] == 1
    assert d.loc["PF3D7_0100200", "n_ip_ms_partners"] == 1


def test_pages_before_a_table_and_pages_with_no_readable_row_are_skipped(tmp_path, monkeypatch):
    """A PDF is mostly not tables: accessions appear in figure legends and in prose.

    Neither case may inherit a bait -- one has none yet, and the other would attach its unreadable
    rows to whatever table came before it.
    """
    prose = "Figure 3. Localisation of PF3D7_1129100 and PF3D7_1024800 in infected erythrocytes."
    unreadable = ("Gene ID PV1-HA WT PV1-HA WT "
                  "PF3D7_1129100 parasitophorous vacuolar protein 1 (PV1) not measured")
    root = _pdf(tmp_path, [prose, unreadable], monkeypatch)
    said = []
    assert P.ip_ms(root, log=said.append).empty
    assert any("1 rows carried the wrong number of counts" in m for m in said)


# --------------------------------------------------------------------------- the cell cycle
def _idc(tmp_path, samples, genes=("PF3D7_0100100", "PF3D7_0100200")):
    """A time course as the deposit ships it: one CSV per sample, plus its own metadata file."""
    folder = tmp_path / "transcription" / P.IDC[0] / P.IDC[1]
    folder.mkdir(parents=True)
    blocks = []
    for gsm, (hour, genotype, strain, values) in samples.items():
        pd.DataFrame({"ORF": list(genes), gsm: values}).to_csv(
            folder / f"{gsm}_sample.csv.gz", index=False, compression="gzip")
        blocks.append(f"^SAMPLE = {gsm}\n"
                      f"!Sample_characteristics_ch1 = hemoglobin genotype: {genotype}\n"
                      f"!Sample_characteristics_ch1 = strain: {strain}\n"
                      f"!Sample_characteristics_ch1 = hours post-invasion: {hour}\n")
    (folder / P.IDC_METADATA).write_text("".join(blocks))
    return str(tmp_path)


def _cycle(peak, hours, high=100.0, low=1.0):
    return [low + (high - low) * (0.5 + 0.5 * np.cos(2 * np.pi * (h - peak) / 48)) for h in hours]


def test_only_the_reference_arm_of_the_experiment_is_read(tmp_path):
    """The deposit crosses two haemoglobin genotypes with two parasite lines, and only one cell of
    that square is unperturbed: the laboratory strain in normal red cells. The rest is the study's
    variable, and reading it would time a gene under sickle-trait stress."""
    samples = {"GSM1": (3, "HbAA", "3D7", [1, 1]), "GSM2": (6, "HbAS", "3D7", [1, 1]),
               "GSM3": (9, "HbAA", "FUP", [1, 1]), "GSM4": (12, "HbAA", "3D7", [1, 1])}
    root = _idc(tmp_path, samples)
    hours = P.idc_samples(os.path.join(root, "transcription", P.IDC[0], P.IDC[1], P.IDC_METADATA))
    assert hours == {"GSM1": 3.0, "GSM4": 12.0}


def test_a_gene_is_timed_by_its_phase_and_a_flat_one_is_not(tmp_path):
    hours = list(range(3, 49, 3))
    samples = {f"GSM{i}": (h, "HbAA", "3D7", [_cycle(45, [h])[0], 50.0])
               for i, h in enumerate(hours)}
    d = P.idc_peak(_idc(tmp_path, samples), log=lambda *a: None).set_index("gene_id")
    assert abs(d.loc["PF3D7_0100100", "idc_peak_hour"] - 45) < 3
    assert "PF3D7_0100200" not in d.index, "a flat profile was given an hour"


def test_a_series_too_short_to_be_a_cycle_is_refused(tmp_path):
    samples = {f"GSM{i}": (h, "HbAA", "3D7", [10.0, 5.0]) for i, h in enumerate((3, 6, 9))}
    said = []
    assert P.idc_peak(_idc(tmp_path, samples), log=said.append).empty
    assert any("not a cycle" in m for m in said)


def test_no_time_course_is_an_empty_frame(tmp_path):
    assert P.idc_peak(str(tmp_path), log=lambda *a: None).empty


@pytest.mark.skipif(not os.path.exists(NODES), reason="Plasmodium table not built")
def test_the_shipped_phases_put_the_invasion_genes_together_across_the_wrap():
    """Read as a line these numbers look wrong, which is the point of asserting them.

    MSP1 lands at 45.7 h and AMA1 at 2.5 h -- five hours apart around a 48-hour circle, both in the
    invasion window, and forty-three hours apart to anyone reading the column as a line.
    """
    d = pd.read_parquet(NODES).set_index("gene_id")
    msp1 = d.loc["PF3D7_0930300", "idc_peak_hour"]
    ama1 = d.loc["PF3D7_1133400", "idc_peak_hour"]
    gap = abs(msp1 - ama1)
    assert min(gap, 48 - gap) < 8, "the two invasion transcripts are not adjacent on the circle"
    assert d.loc["PF3D7_0202000", "idc_peak_hour"] > 18, "KAHRP should peak after the ring"


@pytest.mark.skipif(not os.path.exists(NODES), reason="Plasmodium table not built")
def test_early_phase_genes_are_the_ring_enriched_ones():
    """Checked against a study that shares no sample with this one."""
    d = pd.read_parquet(NODES)
    early = d[d["idc_peak_hour"].between(6, 18)]
    late = d[d["idc_peak_hour"].between(33, 45)]
    assert early["expr_ring"].median() > late["expr_ring"].median() * 2
    assert late["expr_schizont"].median() > early["expr_schizont"].median()


def test_a_one_column_sample_file_is_skipped_rather_than_read_as_expression(tmp_path):
    """A CSV with no values is not a sample, and reading its index as numbers would be worse."""
    hours = list(range(3, 49, 3))
    samples = {f"GSM{i}": (h, "HbAA", "3D7", [_cycle(45, [h])[0], 50.0]) for i, h in enumerate(hours)}
    root = _idc(tmp_path, samples)
    folder = os.path.join(root, "transcription", P.IDC[0], P.IDC[1])
    pd.DataFrame({"ORF": ["PF3D7_0100100"]}).to_csv(
        os.path.join(folder, "GSM0_sample.csv.gz"), index=False, compression="gzip")
    d = P.idc_peak(root, log=lambda *a: None)
    assert len(d) == 1, "the empty file changed the answer"


def test_build_all_folds_in_the_cycle_timing(tmp_path, monkeypatch):
    root = _dataset_root(tmp_path)
    monkeypatch.setattr(P, "idc_peak", lambda *a, **k: pd.DataFrame(
        {"gene_id": ["PF3D7_0100100"], "idc_peak_hour": [45.7], "idc_cycling_amplitude": [0.86]}))
    d = P.build_all(root, log=lambda *a: None).set_index("gene_id")
    assert d.loc["PF3D7_0100100", "idc_peak_hour"] == 45.7
    assert pd.isna(d.loc["PF3D7_0100200", "idc_peak_hour"])
