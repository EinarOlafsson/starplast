#!/usr/bin/env python3
"""The published screens, built against synthetic supplements shaped like the real ones.

Fixtures rather than the real files, deliberately: the real supplements are not committed, so a test
that needs them runs on one machine. What is reproduced here is their *shape* -- the header on the
second row, the accession column that is not called gene_id, the strain accession, the gene listed
twice -- because the shape is what every historical bug in this module was about.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import screens as SC  # noqa: E402


@pytest.fixture(autouse=True)
def _no_leaked_resolver():
    """_RESOLVE is module-global, so a resolver left behind by one test would silently change the next."""
    yield
    SC._RESOLVE = None


# --------------------------------------------------------------------------- accessions
def test_an_accession_is_pulled_out_of_a_messy_cell():
    SC._RESOLVE = None
    s = pd.Series(["TGME49_200010", "gene TGME49_200020 (hypothetical)", "nothing here"])
    out = SC._acc(s)
    assert out.tolist()[:2] == ["TGME49_200010", "TGME49_200020"]
    assert pd.isna(out.iloc[2])


def test_superseded_accessions_go_through_the_identity_layer():
    """The 2019 in vivo screen uses pre-2012 TGME49_0xxxxx ids throughout. They match the regex, exist
    in no current node table, and silently drop all 181 of its genes unless resolved."""
    SC._RESOLVE = lambda a: {"TGME49_012345": "TGME49_249530"}.get(a)
    out = SC._acc(pd.Series(["TGME49_012345"]))
    assert out.tolist() == ["TGME49_249530"]


def test_an_unresolvable_accession_becomes_missing_not_a_guess():
    SC._RESOLVE = lambda a: None
    assert pd.isna(SC._acc(pd.Series(["TGME49_999999"])).iloc[0])


def test_without_a_resolver_the_raw_accession_is_kept():
    SC._RESOLVE = None
    assert SC._acc(pd.Series(["TGME49_200010"])).tolist() == ["TGME49_200010"]


# --------------------------------------------------------------------------- locating files
def test_a_supplement_is_found_flat_in_the_root(tmp_path):
    (tmp_path / "x.xlsx").write_text("")
    assert SC._find(str(tmp_path), "x.xlsx") == str(tmp_path / "x.xlsx")


def test_a_supplement_is_found_one_pmid_directory_down(tmp_path):
    """The tree was reorganised into datasets/<level>/<type>/<PMID>/ and both layouts still exist."""
    d = tmp_path / "10409377"
    d.mkdir()
    (d / "x.xlsx").write_text("")
    assert SC._find(str(tmp_path), "x.xlsx") == str(d / "x.xlsx")


def test_a_missing_supplement_returns_a_path_that_does_not_exist(tmp_path):
    """Callers test os.path.exists on the result, so returning None here would raise there instead."""
    p = SC._find(str(tmp_path), "absent.xlsx")
    assert isinstance(p, str) and not os.path.exists(p)


def test_find_survives_a_root_that_is_not_a_directory(tmp_path):
    f = tmp_path / "not_a_dir"
    f.write_text("")
    assert not os.path.exists(SC._find(str(f), "x.xlsx"))


def test_read_returns_none_for_a_missing_file(tmp_path):
    assert SC._read(str(tmp_path / "nope.xlsx")) is None


# --------------------------------------------------------------------------- the screens themselves
def _tree(tmp_path):
    d = tmp_path / "datasets" / "DNA" / "CRISPR_screen"
    d.mkdir(parents=True)
    return d


def test_no_screen_files_at_all_yields_an_empty_frame_not_an_error(tmp_path):
    _tree(tmp_path)
    out = SC.crispr_screens(str(tmp_path), log=lambda *_: None)
    assert out.empty or out.shape[1] == 0


def test_the_gra17_header_sits_on_the_second_row(tmp_path):
    """Read from row zero it collapses into one unnamed column and the screen contributes nothing."""
    d = _tree(tmp_path)
    frame = pd.DataFrame({
        "ME49ID": ["TGME49_200010", "TGME49_200020"],
        "RHΔgra17 AVG_P3/P4": [1.5, -2.0],
        "Δgra17-RH AVG_P3/P4": [0.5, -1.0],
        "Synthetic lethality/viability candidate": ["yes", None],
    })
    with pd.ExcelWriter(d / "gra17_synthlethal_PMC10409377_S1_phenotypes.xlsx") as w:
        # a title row above the real header, exactly as the published file has
        pd.DataFrame([["Table S1"]]).to_excel(w, sheet_name="TableS1_ALL_DATA",
                                              index=False, header=False)
        frame.to_excel(w, sheet_name="TableS1_ALL_DATA", index=False, startrow=1)

    out = SC.crispr_screens(str(tmp_path), log=lambda *_: None)
    assert "crispr_gra17ko_phenotype" in out.columns
    assert out.loc["TGME49_200010", "crispr_gra17ko_phenotype"] == pytest.approx(1.5)
    assert out.loc["TGME49_200010", "crispr_gra17_candidate"] == 1
    assert out.loc["TGME49_200020", "crispr_gra17_candidate"] == 0


def test_the_two_gra12_screens_are_kept_apart(tmp_path):
    """They are not replicates -- in-vivo L2FC correlates at r = 0.41 -- so merging them would average
    away a real difference and present it as one measurement."""
    d = _tree(tmp_path)
    for fn, sheet, v in (("gra12_PMC12003902_D3_gene_L2FC_screen1.xlsx", "2D.Gene L2FCs", 1.0),
                         ("gra12_PMC12003902_D4_gene_L2FC_screen2.xlsx", "3D.Gene L2FCs", -1.0)):
        pd.DataFrame({"GENE": ["TGME49_200010"], "MEDIAN_L2FC_IN_VITRO": [v],
                      "MEDIAN_L2FC_IN_VIVO": [v * 2], "DISCO_SCORE": [v * 3]}
                     ).to_excel(d / fn, sheet_name=sheet, index=False)
    out = SC.crispr_screens(str(tmp_path), log=lambda *_: None)
    assert out.loc["TGME49_200010", "crispr_gra12s1_l2fc_invitro"] == pytest.approx(1.0)
    assert out.loc["TGME49_200010", "crispr_gra12s2_l2fc_invitro"] == pytest.approx(-1.0)


def test_a_gene_listed_twice_is_averaged_not_taken_first(tmp_path):
    d = _tree(tmp_path)
    pd.DataFrame({"GENE": ["TGME49_200010", "TGME49_200010"],
                  "MEDIAN_L2FC_IN_VITRO": [1.0, 3.0],
                  "MEDIAN_L2FC_IN_VIVO": [0.0, 0.0], "DISCO_SCORE": [0.0, 0.0]}
                 ).to_excel(d / "gra12_PMC12003902_D3_gene_L2FC_screen1.xlsx",
                            sheet_name="2D.Gene L2FCs", index=False)
    out = SC.crispr_screens(str(tmp_path), log=lambda *_: None)
    assert out.loc["TGME49_200010", "crispr_gra12s1_l2fc_invitro"] == pytest.approx(2.0)


def test_the_older_flat_layout_is_still_searched(tmp_path):
    """A checkout from before the reorganisation has datasets/crispr_screens/, and silently finding
    nothing is the failure that move already caused once."""
    d = tmp_path / "datasets" / "crispr_screens"
    d.mkdir(parents=True)
    pd.DataFrame({"GENE": ["TGME49_200010"], "MEDIAN_L2FC_IN_VITRO": [1.0],
                  "MEDIAN_L2FC_IN_VIVO": [1.0], "DISCO_SCORE": [1.0]}
                 ).to_excel(d / "gra12_PMC12003902_D3_gene_L2FC_screen1.xlsx",
                            sheet_name="2D.Gene L2FCs", index=False)
    out = SC.crispr_screens(str(tmp_path), log=lambda *_: None)
    assert "crispr_gra12s1_l2fc_invitro" in out.columns


def test_rows_without_a_resolvable_accession_are_dropped(tmp_path):
    d = _tree(tmp_path)
    pd.DataFrame({"GENE": ["TGME49_200010", "not an accession"],
                  "MEDIAN_L2FC_IN_VITRO": [1.0, 2.0],
                  "MEDIAN_L2FC_IN_VIVO": [1.0, 2.0], "DISCO_SCORE": [1.0, 2.0]}
                 ).to_excel(d / "gra12_PMC12003902_D3_gene_L2FC_screen1.xlsx",
                            sheet_name="2D.Gene L2FCs", index=False)
    out = SC.crispr_screens(str(tmp_path), log=lambda *_: None)
    assert list(out.index) == ["TGME49_200010"]


def test_a_screen_missing_its_expected_column_contributes_what_it_has(tmp_path):
    d = _tree(tmp_path)
    frame = pd.DataFrame({"ME49ID": ["TGME49_200010"], "RHΔgra17 AVG_P3/P4": [1.5]})
    with pd.ExcelWriter(d / "gra17_synthlethal_PMC10409377_S1_phenotypes.xlsx") as w:
        pd.DataFrame([["Table S1"]]).to_excel(w, sheet_name="TableS1_ALL_DATA",
                                              index=False, header=False)
        frame.to_excel(w, sheet_name="TableS1_ALL_DATA", index=False, startrow=1)
    out = SC.crispr_screens(str(tmp_path), log=lambda *_: None)
    assert "crispr_gra17ko_phenotype" in out.columns
    assert "crispr_gra17_synthlethal_delta" not in out.columns


# --------------------------------------------------------------------------- the in vivo platform
def test_repeated_mean_lfc_columns_are_averaged_positionally(tmp_path):
    """These sheets repeat "mean LFC across replicates" once per condition, so selecting by NAME
    returns a frame rather than a column and the arithmetic silently changes meaning."""
    d = _tree(tmp_path)
    frame = pd.DataFrame({"Gene": ["TGME49_200010"],
                          "mean LFC across replicates": [2.0],
                          "mean LFC across replicates ": [4.0]})
    frame.to_excel(d / "invivo_platform_PMC6722137_D3_phenotype_scores.xlsx",
                   sheet_name="Phenotype scores", index=False)
    out = SC.crispr_screens(str(tmp_path), log=lambda *_: None)
    assert out.loc["TGME49_200010", "crispr_invivo_platform_lfc"] == pytest.approx(3.0)


def test_the_in_vivo_sheet_is_read_with_whichever_header_row_works(tmp_path):
    """Two of these workbooks put the header on row 0 and one on row 1; both layouts are tried."""
    d = _tree(tmp_path)
    with pd.ExcelWriter(d / "invivo_platform_PMC6722137_D3_phenotype_scores.xlsx") as w:
        pd.DataFrame([["title"]]).to_excel(w, sheet_name="Phenotype scores",
                                          index=False, header=False)
        pd.DataFrame({"Gene": ["TGME49_200010"],
                      "mean LFC across replicates": [1.0]}).to_excel(
            w, sheet_name="Phenotype scores", index=False, startrow=1)
    out = SC.crispr_screens(str(tmp_path), log=lambda *_: None)
    assert out.loc["TGME49_200010", "crispr_invivo_platform_lfc"] == pytest.approx(1.0)


def test_a_sheet_with_no_lfc_column_is_skipped(tmp_path):
    d = _tree(tmp_path)
    pd.DataFrame({"Gene": ["TGME49_200010"], "something_else": [1.0]}).to_excel(
        d / "invivo_platform_PMC6722137_D3_phenotype_scores.xlsx",
        sheet_name="Phenotype scores", index=False)
    out = SC.crispr_screens(str(tmp_path), log=lambda *_: None)
    assert out.empty or "crispr_invivo_platform_lfc" not in out.columns


# --------------------------------------------------------------------------- host-transcription
def test_a_target_named_by_symbol_is_mapped_through_the_papers_own_sgrna_table(tmp_path):
    """The statistic tables name targets by gene NAME, not accession. Without the map the whole screen
    resolves to nothing."""
    d = _tree(tmp_path)
    pd.DataFrame({"Gene_Name": ["GRA16"], "Gene_ID_Updated": ["TGME49_208830"]}).to_excel(
        d / "hosttx_effectors_PMC12033024_D3_sgRNA_gene_map.xlsx", index=False)
    pd.DataFrame({"Target": ["GRA16"], "T2_Statistic": [5.0],
                  "P_Value_Adjusted": [0.01]}).to_excel(
        d / "hosttx_effectors_PMC12033024_D4A_T2_statistic.xlsx", index=False)
    out = SC.crispr_screens(str(tmp_path), log=lambda *_: None)
    assert out.loc["TGME49_208830", "hosttx_T2"] == pytest.approx(5.0)
    assert out.loc["TGME49_208830", "hosttx_padj"] == pytest.approx(0.01)


def test_a_target_given_as_an_accession_does_not_need_the_symbol_map(tmp_path):
    d = _tree(tmp_path)
    pd.DataFrame({"Gene_Name": ["OTHER"], "Gene_ID_Updated": ["TGME49_111111"]}).to_excel(
        d / "hosttx_effectors_PMC12033024_D3_sgRNA_gene_map.xlsx", index=False)
    pd.DataFrame({"Target": ["TGME49_208830"], "T2_Statistic": [1.0],
                  "P_Value_Adjusted": [0.5]}).to_excel(
        d / "hosttx_effectors_PMC12033024_D5_T2_statistic.xlsx", index=False)
    out = SC.crispr_screens(str(tmp_path), log=lambda *_: None)
    assert "TGME49_208830" in out.index


def test_the_older_gene_id_column_is_used_when_the_updated_one_is_absent(tmp_path):
    d = _tree(tmp_path)
    pd.DataFrame({"Gene_Name": ["GRA16"], "Gene_ID_Old": ["TGME49_208830"]}).to_excel(
        d / "hosttx_effectors_PMC12033024_D3_sgRNA_gene_map.xlsx", index=False)
    pd.DataFrame({"Target": ["GRA16"], "T2_Statistic": [2.0],
                  "P_Value_Adjusted": [0.02]}).to_excel(
        d / "hosttx_effectors_PMC12033024_D4A_T2_statistic.xlsx", index=False)
    out = SC.crispr_screens(str(tmp_path), log=lambda *_: None)
    assert "TGME49_208830" in out.index


def test_a_statistic_table_without_a_target_column_is_skipped(tmp_path):
    d = _tree(tmp_path)
    pd.DataFrame({"Gene_Name": ["GRA16"], "Gene_ID_Updated": ["TGME49_208830"]}).to_excel(
        d / "hosttx_effectors_PMC12033024_D3_sgRNA_gene_map.xlsx", index=False)
    pd.DataFrame({"NotTarget": ["GRA16"]}).to_excel(
        d / "hosttx_effectors_PMC12033024_D4A_T2_statistic.xlsx", index=False)
    out = SC.crispr_screens(str(tmp_path), log=lambda *_: None)
    assert out.empty or "hosttx_T2" not in out.columns


def test_a_gene_map_without_any_id_column_yields_no_lookup(tmp_path):
    d = _tree(tmp_path)
    pd.DataFrame({"Gene_Name": ["GRA16"]}).to_excel(
        d / "hosttx_effectors_PMC12033024_D3_sgRNA_gene_map.xlsx", index=False)
    pd.DataFrame({"Target": ["GRA16"], "T2_Statistic": [1.0],
                  "P_Value_Adjusted": [0.1]}).to_excel(
        d / "hosttx_effectors_PMC12033024_D4A_T2_statistic.xlsx", index=False)
    out = SC.crispr_screens(str(tmp_path), log=lambda *_: None)
    assert out.empty or "hosttx_T2" not in out.columns


# --------------------------------------------------------------------------- proteomics
def _proteomics_tree(tmp_path):
    p = tmp_path / "toxo_stage_atlas" / "data" / "proteomics"
    p.mkdir(parents=True)
    return p


def test_ibaq_is_logged_and_the_median_taken_across_replicates(tmp_path):
    p = _proteomics_tree(tmp_path)
    with pd.ExcelWriter(p / "PXD065585_supp_S1_Pru_proteome_and_IP_iBAQ.xlsx") as w:
        for sh, v in (("Pru_Rep1", 2.0), ("Pru_Rep2", 8.0), ("Pru_Rep3", 32.0)):
            pd.DataFrame({"Protein_ID": ["TGME49_200010"], "iBAQ": [v]}).to_excel(
                w, sheet_name=sh, index=False)
    out = SC.proteomics(str(tmp_path), log=lambda *_: None)
    # log2 of 2, 8, 32 is 1, 3, 5 -- median 3
    assert out.loc["TGME49_200010", "protein_ibaq_log2"] == pytest.approx(3.0)


def test_a_zero_ibaq_is_missing_rather_than_negative_infinity(tmp_path):
    """log2(0) is -inf, which would then dominate every distance in the embedding."""
    p = _proteomics_tree(tmp_path)
    pd.DataFrame({"Protein_ID": ["TGME49_200010", "TGME49_200020"], "iBAQ": [0.0, 4.0]}).to_excel(
        p / "PXD065585_supp_S1_Pru_proteome_and_IP_iBAQ.xlsx", sheet_name="Pru_Rep1", index=False)
    out = SC.proteomics(str(tmp_path), log=lambda *_: None)
    assert "TGME49_200010" not in out.index
    assert np.isfinite(out.loc["TGME49_200020", "protein_ibaq_log2"])


def test_no_ibaq_tables_reports_and_returns_empty(tmp_path):
    _proteomics_tree(tmp_path)
    msgs = []
    assert SC.proteomics(str(tmp_path), log=msgs.append).empty
    assert any("no iBAQ tables" in m for m in msgs)


def test_a_sheet_without_an_ibaq_column_is_skipped(tmp_path):
    p = _proteomics_tree(tmp_path)
    pd.DataFrame({"Protein_ID": ["TGME49_200010"], "other": [1.0]}).to_excel(
        p / "PXD065585_supp_S1_Pru_proteome_and_IP_iBAQ.xlsx", sheet_name="Pru_Rep1", index=False)
    assert SC.proteomics(str(tmp_path), log=lambda *_: None).empty
