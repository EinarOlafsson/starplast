#!/usr/bin/env python3
"""The six expression datasets, and the file shapes that made each of them wrong the first time.

Not one of these six supplements had its header on the first row, and two produced silently wrong
output rather than an error:

  * the in vivo table keys on StringTie identifiers (MSTRG.1) in its first column, which resolve for
    3,723 rows, while `Gene_ref` carries real accessions for 8,717. Taking column zero lost half the
    data and nothing said so.
  * the total proteome has a two-row merged header, so reading row zero alone collapsed a
    3,000-protein abundance block into one unnamed column.

Both are regression cases below. The fixtures reproduce the shapes, not the contents.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import expression as EX  # noqa: E402


# --------------------------------------------------------------------------- helpers
def test_an_accession_is_normalized_out_of_a_messy_string():
    assert EX._norm_acc("TGME49_200010") == "TGME49_200010"
    assert EX._norm_acc("tgme49-200010") == "TGME49_200010"
    assert EX._norm_acc("TGME49.200010") == "TGME49_200010"
    assert EX._norm_acc("gene TGME49_200010 here") == "TGME49_200010"


def test_a_missing_underscore_is_reinserted():
    """Typesetting drops it, and the accession then matches nothing."""
    assert EX._norm_acc("TGME49200010") == "TGME49_200010"


def test_a_string_with_no_accession_is_none():
    assert EX._norm_acc("no identifier here") is None
    assert EX._norm_acc(np.nan) is None


def test_resolution_is_applied_after_normalization():
    out = EX._resolve(pd.Series(["tgme49-200010"]), resolve=lambda a: a.replace("ME49", "GT1"))
    assert out.tolist() == ["TGGT1_200010"]


def test_unresolvable_values_survive_as_missing():
    out = EX._resolve(pd.Series(["TGME49_200010", "junk"]), resolve=lambda a: None)
    assert out.isna().all()


def test_a_gene_listed_once_per_peptide_is_collapsed_by_mean():
    """Supplements routinely list a gene several times, one row per peptide or transcript, and keeping
    the first would pick whichever the file happened to sort first."""
    df = pd.DataFrame({"v": [1.0, 3.0, 10.0]})
    out = EX._collapse(df, pd.Series(["g1", "g1", "g2"]))
    assert out.loc["g1", "v"] == pytest.approx(2.0)
    assert out.loc["g2", "v"] == pytest.approx(10.0)


def test_rows_with_no_gene_are_dropped_before_collapsing():
    df = pd.DataFrame({"v": [1.0, 2.0]})
    out = EX._collapse(df, pd.Series(["g1", None]))
    assert list(out.index) == ["g1"]


# --------------------------------------------------------------------------- in vivo brain
def _invivo(tmp_path, cols):
    d = tmp_path / "datasets" / "translation" / "proteomics" / "31726967"
    d.mkdir(parents=True)
    pd.DataFrame(cols).to_csv(d / "12864_2019_6213_MOESM4_ESM.csv", index=False)
    return d


def test_the_in_vivo_table_keys_on_gene_ref_not_the_first_column(tmp_path):
    """The regression case: column zero is StringTie ids that resolve for 3,723 rows, while Gene_ref
    carries real accessions for 8,717. Taking column zero silently lost half the dataset."""
    _invivo(tmp_path, {"gene_ID": ["MSTRG.1", "MSTRG.2"],
                       "Gene_ref": ["TGME49_200010", "TGME49_200020"],
                       "TZ_1": [1.0, 2.0], "WholeBrain_acute_1": [3.0, 4.0]})
    out = EX.invivo_brain(str(tmp_path), log=lambda *_: None)
    assert set(out.index) == {"TGME49_200010", "TGME49_200020"}
    assert out.shape[1] == 2


def test_the_in_vivo_loader_falls_back_to_the_first_column_if_gene_ref_is_absent(tmp_path):
    _invivo(tmp_path, {"gene_ID": ["TGME49_200010"], "TZ_1": [1.0]})
    out = EX.invivo_brain(str(tmp_path), log=lambda *_: None)
    assert list(out.index) == ["TGME49_200010"]


def test_only_the_named_condition_columns_are_taken(tmp_path):
    _invivo(tmp_path, {"Gene_ref": ["TGME49_200010"], "TZ_1": [1.0],
                       "some_annotation": ["x"], "length": [100]})
    out = EX.invivo_brain(str(tmp_path), log=lambda *_: None)
    assert list(out.columns) == ["invivo_TZ_1"]


def test_a_missing_in_vivo_file_yields_nothing(tmp_path):
    assert EX.invivo_brain(str(tmp_path), log=lambda *_: None).empty


# --------------------------------------------------------------------------- stress induction
def _atlas(tmp_path, kind):
    d = tmp_path / "toxo_stage_atlas" / "data" / kind
    d.mkdir(parents=True, exist_ok=True)
    return d


def test_the_gsm_prefix_is_stripped_because_it_is_a_sample_id_not_a_condition(tmp_path):
    d = _atlas(tmp_path, "transcriptomics")
    pd.DataFrame({"gene": ["TGME49_200010"], "GSM3856123_unstressed_24h": [100],
                  "GSM3856124_alkaline_48h": [200]}).to_csv(
        d / "GSE132248_STAR_counts_matrix.tsv", sep="\t", index=False)
    out = EX.stress_induction(str(tmp_path), log=lambda *_: None)
    assert sorted(out.columns) == ["stress_alkaline_48h", "stress_unstressed_24h"]


def test_non_numeric_columns_are_not_treated_as_samples(tmp_path):
    d = _atlas(tmp_path, "transcriptomics")
    pd.DataFrame({"gene": ["TGME49_200010"], "description": ["hypothetical"],
                  "GSM1_a": [10]}).to_csv(d / "GSE132248_STAR_counts_matrix.tsv",
                                          sep="\t", index=False)
    out = EX.stress_induction(str(tmp_path), log=lambda *_: None)
    assert list(out.columns) == ["stress_a"]


def test_a_missing_stress_file_yields_nothing(tmp_path):
    assert EX.stress_induction(str(tmp_path), log=lambda *_: None).empty


# --------------------------------------------------------------------------- MORC
def test_the_tpm_sheet_is_preferred_and_the_transform_follows_it(tmp_path):
    """The data sits on a second sheet behind a Legend, and which sheet was read decides whether the
    values are counts or TPM -- different transforms."""
    d = _atlas(tmp_path, "proteomics")
    with pd.ExcelWriter(d / "PXD058095_supp_DatasetEV1_MORC_RNAseq_counts_TPM.xlsx") as w:
        pd.DataFrame({"note": ["legend"]}).to_excel(w, sheet_name="Legend", index=False)
        pd.DataFrame({"gene": ["TGME49_200010"], "MORC_KD": [10.0]}).to_excel(
            w, sheet_name="TPM", index=False)
    out = EX.morc_depletion(str(tmp_path), log=lambda *_: None)
    assert list(out.columns) == ["morc_MORC_KD"]


def test_the_raw_data_sheet_is_used_when_there_is_no_tpm(tmp_path):
    d = _atlas(tmp_path, "proteomics")
    pd.DataFrame({"gene": ["TGME49_200010"], "MORC_KD": [10]}).to_excel(
        d / "PXD058095_supp_DatasetEV1_MORC_RNAseq_counts_TPM.xlsx",
        sheet_name="Raw Data", index=False)
    out = EX.morc_depletion(str(tmp_path), log=lambda *_: None)
    assert list(out.columns) == ["morc_MORC_KD"]


def test_a_missing_morc_file_yields_nothing(tmp_path):
    assert EX.morc_depletion(str(tmp_path), log=lambda *_: None).empty


# --------------------------------------------------------------------------- total proteome
def _proteome(tmp_path, header_top, header_sub, rows):
    d = _atlas(tmp_path, "proteomics")
    frame = pd.DataFrame([header_sub] + rows)
    with pd.ExcelWriter(d / "PXD039400_PXD042658_supp_SupplTable3_total_proteome.xlsx") as w:
        pd.DataFrame([header_top]).to_excel(w, sheet_name="MS DATA", index=False, header=False)
        frame.to_excel(w, sheet_name="MS DATA", index=False, header=False, startrow=1)
    return d


def test_the_two_row_merged_header_is_read_as_two_rows(tmp_path):
    """The regression case: reading row zero alone collapsed a 3,000-protein abundance block into one
    unnamed column, and the dataset appeared to contribute a single column."""
    _proteome(tmp_path,
              ["", "log2(normalized abundance)", "", "UT Vs T-24h"],
              ["Accession", "UT R1", "UT R2", "log2(fold change)"],
              [["TGME49_200010", 5.0, 6.0, 1.2], ["TGME49_200020", 7.0, 8.0, -0.4]])
    out = EX.total_proteome(str(tmp_path), log=lambda *_: None)
    assert len(out) == 2
    assert sum(c.startswith("proteome_lfc_") for c in out.columns) == 1
    assert sum(c.startswith("proteome_") and not c.startswith("proteome_lfc_")
               for c in out.columns) == 2


def test_abundances_and_fold_changes_are_normalized_differently(tmp_path):
    """One is an intensity already on a log scale, the other is a ratio: centring the ratio would move
    the zero, which is the reference condition."""
    _proteome(tmp_path,
              ["", "log2(normalized abundance)", "", "UT Vs T-24h"],
              ["Accession", "UT R1", "UT R2", "log2(fold change)"],
              [["TGME49_200010", 1.0, 1.0, 2.0], ["TGME49_200020", 3.0, 3.0, 4.0]])
    out = EX.total_proteome(str(tmp_path), log=lambda *_: None)
    lfc = [c for c in out.columns if c.startswith("proteome_lfc_")][0]
    assert out[lfc].tolist() == pytest.approx([2.0, 4.0])          # untouched
    ab = [c for c in out.columns if not c.startswith("proteome_lfc_")][0]
    assert out[ab].median() == pytest.approx(0.0)                  # centred


def test_a_proteome_with_no_recognisable_value_columns_yields_nothing(tmp_path):
    _proteome(tmp_path, ["", ""], ["Accession", "Description"],
              [["TGME49_200010", "hypothetical"]])
    assert EX.total_proteome(str(tmp_path), log=lambda *_: None).empty


def test_a_missing_proteome_file_yields_nothing(tmp_path):
    assert EX.total_proteome(str(tmp_path), log=lambda *_: None).empty


# --------------------------------------------------------------------------- phosphosites
def test_sites_are_counted_per_gene_and_the_ratio_summarised(tmp_path):
    """Site level collapsed to gene level: a per-gene table is the wrong shape for site positions, so
    what it keeps is how many were measured and how strongly they moved."""
    d = _atlas(tmp_path, "proteomics")
    for fn, ratios in (("PXD017032_supp_TableS1_upregulated_phosphosites_oocyst_vs_tachy.xlsx",
                        [2.0, 4.0]),
                       ("PXD017032_supp_TableS2_downregulated_phosphosites_oocyst_vs_tachy.xlsx",
                        [0.5])):
        pd.DataFrame({"Protein ID": ["TGME49_200010"] * len(ratios),
                      "Ratio": ratios}).to_excel(d / fn, index=False)
    out = EX.phosphosites(str(tmp_path), log=lambda *_: None)
    assert out.loc["TGME49_200010", "phospho_up_sites"] == 2
    assert out.loc["TGME49_200010", "phospho_up_ratio"] == pytest.approx(3.0)
    assert out.loc["TGME49_200010", "phospho_sites_measured"] == 3


def test_a_phospho_table_without_a_ratio_column_still_counts_sites(tmp_path):
    d = _atlas(tmp_path, "proteomics")
    pd.DataFrame({"Protein ID": ["TGME49_200010", "TGME49_200010"]}).to_excel(
        d / "PXD017032_supp_TableS1_upregulated_phosphosites_oocyst_vs_tachy.xlsx", index=False)
    out = EX.phosphosites(str(tmp_path), log=lambda *_: None)
    assert out.loc["TGME49_200010", "phospho_up_sites"] == 2


def test_missing_phospho_files_yield_nothing(tmp_path):
    assert EX.phosphosites(str(tmp_path), log=lambda *_: None).empty


# --------------------------------------------------------------------------- oocyst iTRAQ
def test_a_legacy_xls_that_cannot_be_converted_yields_nothing(tmp_path, monkeypatch):
    """libreoffice is not always installed, and its absence must cost this dataset rather than the
    build."""
    d = _atlas(tmp_path, "proteomics")
    (d / "PXD003765_supp_mmc2_iTRAQ_ratios_2095proteins.xls").write_text("not really an xls")
    monkeypatch.setattr(EX, "_libreoffice_xlsx", lambda p, log=print: None)
    assert EX.oocyst_itraq(str(tmp_path), log=lambda *_: None).empty


def test_the_converted_workbook_is_parsed_and_ratio_columns_found(tmp_path, monkeypatch):
    d = _atlas(tmp_path, "proteomics")
    (d / "PXD003765_supp_mmc2_iTRAQ_ratios_2095proteins.xls").write_text("stub")
    conv = tmp_path / "converted.xlsx"
    with pd.ExcelWriter(conv) as w:
        pd.DataFrame([["title"]]).to_excel(w, sheet_name="Sheet1", index=False, header=False)
        pd.DataFrame({"Accession": ["TGME49_200010"], "115:113": [1.5],
                      "116:113": [2.5]}).to_excel(w, sheet_name="Sheet1", index=False, startrow=1)
    monkeypatch.setattr(EX, "_libreoffice_xlsx", lambda p, log=print: str(conv))
    out = EX.oocyst_itraq(str(tmp_path), log=lambda *_: None)
    assert len(out.columns) == 2
    assert out.iloc[0].tolist() == pytest.approx([1.5, 2.5])       # ratios left untouched


def test_a_missing_itraq_file_yields_nothing(tmp_path):
    assert EX.oocyst_itraq(str(tmp_path), log=lambda *_: None).empty


def test_libreoffice_absence_is_reported_not_raised(tmp_path):
    import subprocess
    msgs = []

    def boom(*a, **k):
        raise FileNotFoundError("libreoffice")

    real = subprocess.run
    subprocess.run = boom
    try:
        assert EX._libreoffice_xlsx(str(tmp_path / "x.xls"), log=msgs.append) is None
    finally:
        subprocess.run = real
    assert any("libreoffice unavailable" in m for m in msgs)


def test_a_conversion_that_produces_nothing_returns_none(tmp_path, monkeypatch):
    import subprocess
    monkeypatch.setattr(subprocess, "run", lambda *a, **k: None)
    assert EX._libreoffice_xlsx(str(tmp_path / "x.xls"), log=lambda *_: None) is None


# --------------------------------------------------------------------------- composition
def test_load_all_joins_what_it_finds_and_survives_what_it_does_not(tmp_path):
    _invivo(tmp_path, {"Gene_ref": ["TGME49_200010"], "TZ_1": [1.0]})
    out = EX.load_all(str(tmp_path), log=lambda *_: None)
    assert "invivo_TZ_1" in out.columns


def test_one_loader_failing_does_not_lose_the_others(tmp_path, monkeypatch):
    """Six datasets, six file formats; one publisher reissuing a file should not cost the other five."""
    _invivo(tmp_path, {"Gene_ref": ["TGME49_200010"], "TZ_1": [1.0]})

    def boom(base, resolve=None, log=print):
        raise ValueError("this supplement changed shape")

    monkeypatch.setattr(EX, "LOADERS", (boom, EX.invivo_brain))
    msgs = []
    out = EX.load_all(str(tmp_path), log=msgs.append)
    assert "invivo_TZ_1" in out.columns
    assert any("failed" in m for m in msgs)


def test_no_datasets_at_all_returns_an_empty_frame(tmp_path):
    assert EX.load_all(str(tmp_path), log=lambda *_: None).empty


def test_a_space_separated_accession_is_deliberately_not_matched():
    """In running prose it would join two adjacent tokens into an accession that was never written."""
    assert EX._norm_acc("TGME49 200010") is None


def test_the_itraq_loader_falls_back_to_the_numeric_block_when_no_column_names_a_ratio(tmp_path,
                                                                                       monkeypatch):
    """Some releases of this file label the channels only in a merged header row, so the ratio columns
    arrive unnamed and have to be found by position instead."""
    d = _atlas(tmp_path, "proteomics")
    (d / "PXD003765_supp_mmc2_iTRAQ_ratios_2095proteins.xls").write_text("stub")
    conv = tmp_path / "converted2.xlsx"
    cols = {"Accession": ["TGME49_200010"]}
    for i in range(8):
        cols[f"col{i}"] = [float(i)]
    with pd.ExcelWriter(conv) as w:
        pd.DataFrame([["title"]]).to_excel(w, sheet_name="Sheet1", index=False, header=False)
        pd.DataFrame(cols).to_excel(w, sheet_name="Sheet1", index=False, startrow=1)
    monkeypatch.setattr(EX, "_libreoffice_xlsx", lambda p, log=print: str(conv))
    out = EX.oocyst_itraq(str(tmp_path), log=lambda *_: None)
    assert len(out.columns) == 4, "the first four numeric columns are annotation, not measurements"
