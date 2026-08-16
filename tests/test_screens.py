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


def test_hyperlopit_unassigned_libraries_remain_separate_and_missing_is_unknown(tmp_path):
    acquired = tmp_path / "datasets" / "toxoplasma_acquisition_2026_08_14"
    acquired.mkdir(parents=True)
    pd.DataFrame({"ME49_ID": ["TGME49_200010", "TGME49_200020"],
                  "In vivo fitness score": [1.0, 2.0]}).to_excel(
        acquired / "GSE253884_unassigned_1_summary.xlsx", sheet_name="Summary", index=False)
    pd.DataFrame({"ME49_ID": ["TGME49_200010"],
                  "In vivo fitness score": [-1.0]}).to_excel(
        acquired / "GSE253885_unassigned_2_summary.xlsx", sheet_name="Summary", index=False)
    out = SC.crispr_screens(str(tmp_path), log=lambda *_: None)
    assert out.loc["TGME49_200010", "fit_hyperlopit_unassigned_invivo_lib1"] == 1.0
    assert out.loc["TGME49_200010", "fit_hyperlopit_unassigned_invivo_lib2"] == -1.0
    assert pd.isna(out.loc["TGME49_200020", "fit_hyperlopit_unassigned_invivo_lib2"])


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


def test_full_host_response_signatures_are_reduced_per_parasite_effector(tmp_path):
    d = _tree(tmp_path) / "37827122"
    d.mkdir()
    pd.DataFrame({"Gene_Name": ["GRA16", "GRA24", "OTHER"],
                  "Gene_ID_Updated": ["TGME49_208830", "TGME49_230180", np.nan]}).to_excel(
        d / "hosttx_effectors_PMC12033024_D3_sgRNA_gene_map.xlsx", index=False)
    rows = []
    for target, shift in (("GRA16", 1.0), ("GRA24", -1.0), ("TGME49_200010", 0.25)):
        for i in range(5):
            rows.append({"target": target, "gene": f"HOST{i}",
                         "avg_log2FC": shift * (i + 1),
                         "p_val_bh": 0.01 if i < 3 else 0.5})
    pd.DataFrame(rows).to_csv(d / "hosttx_effectors_DE_host_genes.csv", index=False)
    resolve = lambda accession: accession if str(accession).startswith("TGME49_") else None
    out = SC.host_transcription_signatures(str(tmp_path), resolve=resolve, components=2,
                                           log=lambda *_: None)
    assert out.shape == (3, 4)
    assert set(out.index) == {"TGME49_208830", "TGME49_230180", "TGME49_200010"}
    assert out.loc["TGME49_208830", "hosttx_signature_n_de"] == 3
    assert np.isfinite(out.filter(like="_pc").to_numpy()).all()


def test_full_host_response_signature_absence_and_bad_schema_are_explicit(tmp_path):
    msgs = []
    assert SC.host_transcription_signatures(str(tmp_path), log=msgs.append).empty
    assert "absent" in msgs[-1]
    d = _tree(tmp_path) / "37827122"
    d.mkdir()
    pd.DataFrame({"Gene_Name": ["GRA16"], "Gene_ID_Updated": ["TGME49_208830"]}).to_excel(
        d / "hosttx_effectors_PMC12033024_D3_sgRNA_gene_map.xlsx", index=False)
    pd.DataFrame({"wrong": [1]}).to_csv(d / "hosttx_effectors_DE_host_genes.csv", index=False)
    msgs = []
    assert SC.host_transcription_signatures(str(tmp_path), log=msgs.append).empty
    assert "unexpected columns" in msgs[-1]


def test_full_host_response_drops_targets_that_cannot_be_resolved(tmp_path):
    d = _tree(tmp_path) / "37827122"
    d.mkdir()
    pd.DataFrame({"Gene_Name": ["UNKNOWN"], "Gene_ID_Updated": [np.nan]}).to_excel(
        d / "hosttx_effectors_PMC12033024_D3_sgRNA_gene_map.xlsx", index=False)
    pd.DataFrame({"target": ["UNKNOWN"], "gene": ["HOST1"], "avg_log2FC": [1.0],
                  "p_val_bh": [0.01]}).to_csv(
        d / "hosttx_effectors_DE_host_genes.csv", index=False)
    assert SC.host_transcription_signatures(
        str(tmp_path), resolve=lambda _a: None, log=lambda *_: None).empty


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


# ------------------------------------------------------- differentiation reporter screen
def _diff_archive(tmp_path, samples, *, matrix=True, matrix_body=None, plain=False):
    """A GSE132237-shaped deposit: guide count files in a tar, plus the series matrix.

    `samples` maps sample TITLE -> {gene: [per-guide counts]}. The archive names its members by GSM
    accession and the submitter's own suffix, exactly as GEO does, so a test that passes by matching
    the title against the filename would be testing the wrong thing.
    """
    import gzip as _gz
    import io
    import tarfile

    def _directory(name):
        info = tarfile.TarInfo(name)
        info.type = tarfile.DIRTYPE
        return info

    root = tmp_path / "deposit"
    root.mkdir(exist_ok=True)
    accessions = {title: f"GSM{3854790 + i}" for i, title in enumerate(samples)}
    with tarfile.open(root / SC.DIFFERENTIATION_SCREEN, "w") as archive:
        # A directory member: `extractfile` returns None for it and the loop must survive that.
        archive.addfile(_directory("subdir"))
        for title, genes in samples.items():
            rows = "".join(f"{gene}-{n}\t{count}\n"
                           for gene, counts in genes.items()
                           for n, count in enumerate(counts))
            rows += ("no-underscore-here\t5\n"          # a header or a control guide
                     "malformed_row\n"                  # one column
                     "TGME49_000001-9\tnotanumber\n")  # a count that is not a number
            blob = rows.encode() if plain else _gz.compress(rows.encode())
            name = f"{accessions[title]}_{title.replace(' ', '-')}_Counted.txt"
            info = tarfile.TarInfo(name if plain else name + ".gz")
            info.size = len(blob)
            archive.addfile(info, io.BytesIO(blob))
    if matrix:
        body = matrix_body if matrix_body is not None else (
            '!Sample_title\t' + "\t".join(f'"{t}"' for t in samples) + "\n"
            '!Sample_geo_accession\t' + "\t".join(f'"{accessions[t]}"' for t in samples) + "\n")
        with _gz.open(root / SC.DIFFERENTIATION_SCREEN.replace("_RAW.tar", "_series_matrix.txt.gz"),
                      "wt") as fh:
            fh.write(body)
    return root


def _two_arms(**over):
    """Both arms of both replicates, with per-gene guide counts. Defaults give a ratio of 0."""
    base = {"TGME49_000001": [10.0, 30.0], "TGME49_000002": [20.0, 20.0]}
    samples = {}
    for title in ("L1 mNG+ 10d", "L1 bulk brady 10d", "L2 mNG+ 10d", "L2 bulk brady 10d"):
        samples[title] = over.get(title, base)
    return samples


def test_guides_are_summed_per_gene_before_the_ratio(tmp_path):
    """Two guides at 10 and 30 is one gene at 40, not a gene whose ratio was averaged per guide."""
    root = _diff_archive(tmp_path, _two_arms(**{
        "L1 mNG+ 10d": {"TGME49_000001": [30.0, 30.0], "TGME49_000002": [20.0, 20.0]}}))
    out = SC.differentiation_screen(str(root), log=lambda *_: None, resolve=lambda s: s)
    # Arm 1 is 60/40 of a library of 100 vs 40/40 of a library of 80: 600000 vs 500000 CPM.
    per_gene = out["diff_reporter_log2_mNG_over_bulk"]
    assert per_gene["TGME49_000001"] > 0 and per_gene["TGME49_000002"] < 0


def test_the_two_replicate_arms_are_averaged(tmp_path):
    """L1 up and L2 down by the same amount is a gene with no phenotype, not one with L1's."""
    root = _diff_archive(tmp_path, _two_arms(**{
        "L1 mNG+ 10d": {"TGME49_000001": [40.0], "TGME49_000002": [10.0]},
        "L1 bulk brady 10d": {"TGME49_000001": [10.0], "TGME49_000002": [40.0]},
        "L2 mNG+ 10d": {"TGME49_000001": [10.0], "TGME49_000002": [40.0]},
        "L2 bulk brady 10d": {"TGME49_000001": [40.0], "TGME49_000002": [10.0]}}))
    out = SC.differentiation_screen(str(root), log=lambda *_: None, resolve=lambda s: s)
    assert out["diff_reporter_log2_mNG_over_bulk"].abs().max() < 1e-9


def test_counts_are_scaled_to_a_common_library_size(tmp_path):
    """One arm sequenced twice as deep is not every gene enriched twofold."""
    deep = {"TGME49_000001": [20.0], "TGME49_000002": [40.0]}
    shallow = {"TGME49_000001": [10.0], "TGME49_000002": [20.0]}
    out = SC.differentiation_screen(
        str(_diff_archive(tmp_path, _two_arms(**{"L1 mNG+ 10d": deep, "L1 bulk brady 10d": shallow,
                                                 "L2 mNG+ 10d": deep, "L2 bulk brady 10d": shallow}))),
        log=lambda *_: None, resolve=lambda s: s)
    assert out["diff_reporter_log2_mNG_over_bulk"].abs().max() < 0.02


def test_the_sample_mapping_comes_from_the_series_matrix_not_the_file_name(tmp_path):
    """Swap which GSM the titles point at and the ratio inverts -- proof the matrix is what is read."""
    samples = _two_arms(**{
        "L1 mNG+ 10d": {"TGME49_000001": [90.0], "TGME49_000002": [10.0]},
        "L1 bulk brady 10d": {"TGME49_000001": [10.0], "TGME49_000002": [90.0]}})
    root = _diff_archive(tmp_path, samples)
    straight = SC.differentiation_screen(str(root), log=lambda *_: None, resolve=lambda s: s)
    titles = list(samples)
    swapped = dict(zip(titles, [f"GSM{3854790 + i}" for i in (1, 0, 2, 3)]))
    import gzip as _gz
    with _gz.open(root / SC.DIFFERENTIATION_SCREEN.replace("_RAW.tar", "_series_matrix.txt.gz"),
                  "wt") as fh:
        fh.write('!Sample_title\t' + "\t".join(f'"{t}"' for t in titles) + "\n"
                 '!Sample_geo_accession\t' + "\t".join(f'"{swapped[t]}"' for t in titles) + "\n")
    inverted = SC.differentiation_screen(str(root), log=lambda *_: None, resolve=lambda s: s)
    c = "diff_reporter_log2_mNG_over_bulk"
    assert straight[c]["TGME49_000001"] == pytest.approx(-inverted[c]["TGME49_000001"])


def test_an_uncompressed_member_is_read_too(tmp_path):
    """GEO does not promise every supplementary file is gzipped."""
    out = SC.differentiation_screen(str(_diff_archive(tmp_path, _two_arms(), plain=True)),
                                    log=lambda *_: None, resolve=lambda s: s)
    assert len(out) == 2


def test_no_archive_returns_empty(tmp_path):
    assert SC.differentiation_screen(str(tmp_path), log=lambda *_: None).empty


def test_without_the_series_matrix_no_arm_can_be_identified(tmp_path):
    """Refusing is right: guessing the mapping from file names is what this design exists to avoid."""
    msgs = []
    assert SC.differentiation_screen(str(_diff_archive(tmp_path, _two_arms(), matrix=False)),
                                     log=msgs.append).empty
    assert any("no matching arms" in m for m in msgs)


def test_a_series_matrix_without_the_title_lines_is_not_used(tmp_path):
    root = _diff_archive(tmp_path, _two_arms(), matrix_body="!Series_title\tsomething else\n")
    assert SC.differentiation_screen(str(root), log=lambda *_: None).empty


def test_an_empty_arm_is_not_divided_by_its_own_zero_total(tmp_path):
    samples = _two_arms()
    samples["L1 mNG+ 10d"] = {"TGME49_000001": [0.0], "TGME49_000002": [0.0]}
    out = SC.differentiation_screen(str(_diff_archive(tmp_path, samples)),
                                    log=lambda *_: None, resolve=lambda s: s)
    assert np.isfinite(out["diff_reporter_log2_mNG_over_bulk"]).all()


def test_arms_sharing_no_genes_are_skipped(tmp_path):
    """L2 alone still yields a column; an arm pair with no gene in common contributes nothing."""
    samples = _two_arms()
    samples["L1 mNG+ 10d"] = {"TGME49_999999": [10.0]}
    samples["L1 bulk brady 10d"] = {"TGME49_888888": [10.0]}
    out = SC.differentiation_screen(str(_diff_archive(tmp_path, samples)),
                                    log=lambda *_: None, resolve=lambda s: s)
    assert set(out.index) == {"TGME49_000001", "TGME49_000002"}


def test_accessions_are_resolved_through_the_identity_layer_by_default(tmp_path):
    """resolve=None routes through _acc, which is how a retired id reaches its current one."""
    SC._RESOLVE = lambda s: "TGME49_111111" if s == "TGME49_000001" else s
    out = SC.differentiation_screen(str(_diff_archive(tmp_path, _two_arms())), log=lambda *_: None)
    assert "TGME49_111111" in out.index


def test_two_ids_resolving_to_one_gene_keep_a_single_row(tmp_path):
    out = SC.differentiation_screen(str(_diff_archive(tmp_path, _two_arms())),
                                    log=lambda *_: None, resolve=lambda s: "TGME49_111111")
    assert list(out.index) == ["TGME49_111111"]


def test_a_resolver_returning_nothing_leaves_the_original_id(tmp_path):
    out = SC.differentiation_screen(str(_diff_archive(tmp_path, _two_arms())),
                                    log=lambda *_: None, resolve=lambda s: None)
    assert set(out.index) == {"TGME49_000001", "TGME49_000002"}


# ------------------------------------------------------------- oxidative stress and cyst wall
def test_the_oxidative_screen_reads_the_authors_own_score(tmp_path):
    """Their score, not a recomputation from the guide counts beside it in the same workbook."""
    (tmp_path / SC.OXIDATIVE_SCREEN).write_text(
        "gene_id\toxidative_stress_screen_score\nTGGT1_100010\t-6.15\n")
    out = SC.oxidative_stress_screen(str(tmp_path), log=lambda *_: None,
                                     resolve=lambda g: "TGME49_200010")
    assert out.loc["TGME49_200010", "oxidative_stress_screen_score"] == pytest.approx(-6.15)


def test_a_gene_listed_twice_in_the_screen_is_averaged(tmp_path):
    (tmp_path / SC.OXIDATIVE_SCREEN).write_text(
        "gene_id\tscore\nTGME49_200010\t-2\nTGME49_200010\t-4\n")
    out = SC.oxidative_stress_screen(str(tmp_path), log=lambda *_: None, resolve=lambda g: g)
    assert out.loc["TGME49_200010", "oxidative_stress_screen_score"] == pytest.approx(-3.0)


def test_no_oxidative_screen_file_yields_nothing(tmp_path):
    assert SC.oxidative_stress_screen(str(tmp_path), log=lambda *_: None).empty


def test_a_one_column_oxidative_file_is_refused(tmp_path):
    (tmp_path / SC.OXIDATIVE_SCREEN).write_text("gene_id\nTGME49_200010\n")
    assert SC.oxidative_stress_screen(str(tmp_path), log=lambda *_: None).empty


def _cyst_wall(tmp_path, rows, baits=("CST1", "MAG1")):
    header = list(SC.CYST_WALL_META) + list(baits)
    lines = ["\t".join(header)]
    for acc, values in rows:
        blank = [""] * (len(SC.CYST_WALL_META) - 5)
        lines.append("\t".join(["1", "True", "", "protein", acc] + blank
                               + [str(v) for v in values]))
    (tmp_path / SC.CYST_WALL).write_text("\n".join(lines) + "\n")
    return tmp_path


def test_the_cyst_wall_keeps_the_strongest_bait_and_counts_them(tmp_path):
    _cyst_wall(tmp_path, [("TGME49_270240", (0.07, 0.02))])
    out = SC.cyst_wall_interactome(str(tmp_path), log=lambda *_: None, resolve=lambda g: g)
    assert out.loc["TGME49_270240", "cyst_wall_max_spectral"] == pytest.approx(0.07)
    assert out.loc["TGME49_270240", "cyst_wall_n_baits"] == 2


def test_a_bait_that_saw_nothing_is_not_counted(tmp_path):
    _cyst_wall(tmp_path, [("TGME49_270240", (0.07, 0.0))])
    out = SC.cyst_wall_interactome(str(tmp_path), log=lambda *_: None, resolve=lambda g: g)
    assert out.loc["TGME49_270240", "cyst_wall_n_baits"] == 1


def test_host_proteins_in_the_pulldown_are_dropped(tmp_path):
    """Most of the table is human -- the pulldowns were done on infected cultures."""
    _cyst_wall(tmp_path, [("ACACA_HUMAN", (0.5, 0.5)), ("TGME49_270240", (0.07, 0.02))])
    out = SC.cyst_wall_interactome(str(tmp_path), log=lambda *_: None, resolve=lambda g: g)
    assert list(out.index) == ["TGME49_270240"]


def test_the_type_i_namespace_is_recognised_too(tmp_path):
    _cyst_wall(tmp_path, [("TGGT1_100010", (0.07, 0.02))])
    out = SC.cyst_wall_interactome(str(tmp_path), log=lambda *_: None,
                                   resolve=lambda g: "TGME49_270240")
    assert list(out.index) == ["TGME49_270240"]


def test_no_cyst_wall_file_yields_nothing(tmp_path):
    assert SC.cyst_wall_interactome(str(tmp_path), log=lambda *_: None).empty


def test_a_cyst_wall_table_without_accessions_is_refused(tmp_path):
    (tmp_path / SC.CYST_WALL).write_text("something\telse\n1\t2\n")
    assert SC.cyst_wall_interactome(str(tmp_path), log=lambda *_: None).empty


def test_a_cyst_wall_table_with_no_baits_is_refused(tmp_path):
    (tmp_path / SC.CYST_WALL).write_text("\t".join(SC.CYST_WALL_META) + "\n")
    assert SC.cyst_wall_interactome(str(tmp_path), log=lambda *_: None).empty


def test_the_thermal_shift_score_is_read_per_protein(tmp_path):
    (tmp_path / SC.THERMAL_SHIFT).write_text("id\tEDscore\nTGGT1_100010\t0.411\n")
    out = SC.thermal_shift(str(tmp_path), log=lambda *_: None,
                           resolve=lambda g: "TGME49_200010")
    assert out.loc["TGME49_200010", "cetsa_calcium_ed_score"] == pytest.approx(0.411)


def test_a_protein_listed_twice_keeps_the_larger_shift(tmp_path):
    (tmp_path / SC.THERMAL_SHIFT).write_text(
        "id\tEDscore\nTGME49_200010\t0.1\nTGME49_200010\t0.4\n")
    out = SC.thermal_shift(str(tmp_path), log=lambda *_: None, resolve=lambda g: g)
    assert out.loc["TGME49_200010", "cetsa_calcium_ed_score"] == pytest.approx(0.4)


def test_the_thermal_shift_falls_back_to_the_identity_layer(tmp_path):
    """resolve=None routes through _acc, as every other loader here does."""
    SC._RESOLVE = None
    (tmp_path / SC.THERMAL_SHIFT).write_text("id\tEDscore\nTGME49_200010\t0.4\n")
    out = SC.thermal_shift(str(tmp_path), log=lambda *_: None)
    assert list(out.index) == ["TGME49_200010"]


def test_no_thermal_shift_file_yields_nothing(tmp_path):
    assert SC.thermal_shift(str(tmp_path), log=lambda *_: None).empty


def test_a_one_column_thermal_shift_file_is_refused(tmp_path):
    (tmp_path / SC.THERMAL_SHIFT).write_text("id\nTGME49_200010\n")
    assert SC.thermal_shift(str(tmp_path), log=lambda *_: None).empty
