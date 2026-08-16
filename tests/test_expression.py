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
import gzip
import io
import tarfile

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


def test_gse22258_uses_accession_keyed_stage_values(tmp_path):
    d = tmp_path / "datasets"
    d.mkdir()
    pd.DataFrame({"ID_REF": ["TGME49_200010", "TGME49_200020"],
                  "GSM554066": [10.0, 11.0], "GSM554067": [12.0, 13.0]}).to_csv(
        d / "stagetranscriptome_GSE22258_series_matrix.txt.gz", sep="\t", index=False,
        compression="gzip")
    out = EX.gse22258_stage(str(tmp_path), log=lambda *_: None)
    assert list(out) == ["rna22258_tachyzoite", "rna22258_bradyzoite"]
    assert set(out.index) == {"TGME49_200010", "TGME49_200020"}


def test_gse168465_keeps_measurements_but_not_p_values(tmp_path):
    d = tmp_path / "datasets"
    d.mkdir()
    with pd.ExcelWriter(d / "stagetranscriptome_GSE168465_DESeq2-Toxo-all-time-points.xlsx") as w:
        pd.DataFrame({"gene": ["TGME49_200010"], "baseMean": [100.0],
                      "log2FoldChange": [2.0], "pvalue": [0.001], "padj": [0.01]}).to_excel(
            w, sheet_name="1d", index=False)
    out = EX.neuronal_differentiation(str(tmp_path), log=lambda *_: None)
    assert list(out) == ["brain168465_1d_base_mean", "brain168465_1d_lfc"]
    assert not any("pvalue" in c or "padj" in c for c in out)


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


# --------------------------------------------------------------------------- 2026 GEO acquisition
def _acquired(tmp_path):
    path = tmp_path / "datasets" / "toxoplasma_acquisition_2026_08_14"
    path.mkdir(parents=True)
    return path


def test_gse99395_keeps_rna_rpf_and_relative_translation_separate(tmp_path):
    path = _acquired(tmp_path)
    frame = pd.DataFrame({
        "gene": ["TGME49_200010", "TGME49_200020"],
        "Extracellular_Digested_1": [8, 16], "Extracellular_Digested_2": [4, 8],
        "Intracellular_Digested_1": [16, 8], "Intracellular_Digested_2": [8, 4],
        "Extracellular_Total_1": [4, 8], "Extracellular_Total_2": [2, 4],
        "Intracellular_Total_1": [8, 4], "Intracellular_Total_2": [4, 2],
    })
    frame.to_csv(path / "GSE99395_Raw_counts.txt.gz", sep="\t", index=False,
                 compression="gzip")
    out = EX.gse99395_ribosome_profiling(str(tmp_path), log=lambda *_: None)
    assert out.shape == (2, 12)
    assert {"rpf99395_extracellular_r1", "rna99395_extracellular_r1",
            "te99395_extracellular_r1"} <= set(out)
    assert np.allclose(out.te99395_extracellular_r1,
                       out.rpf99395_extracellular_r1 - out.rna99395_extracellular_r1)


def test_gse129869_reads_each_archived_assay_without_extracting_it(tmp_path):
    path = _acquired(tmp_path)
    archive_path = path / "GSE129869_RAW.tar"
    with tarfile.open(archive_path, "w") as archive:
        for condition in ("c", "s"):
            for assay in ("RFP", "RNA"):
                for replicate in (1, 2, 3):
                    payload = gzip.compress(b"TGME49_200010\t8\nTGME49_200020\t4\n")
                    info = tarfile.TarInfo(f"sample_{condition}{assay}.RH{replicate}_count.tab.gz")
                    info.size = len(payload)
                    archive.addfile(info, io.BytesIO(payload))
    out = EX.gse129869_host_context_ribosome_profiling(
        str(tmp_path), log=lambda *_: None)
    assert out.shape == (2, 18)
    assert {"rpf129869_confluent_r1", "rna129869_subconfluent_r3",
            "te129869_confluent_r2"} <= set(out)


def _write_geo_matrix(path, filename, titles):
    matrix = ["!Sample_title\t" + "\t".join(f'\"{title}\"' for title in titles),
              "!series_matrix_table_begin",
              '"ID_REF"\t' + "\t".join(f'\"GSM{i}\"' for i in range(len(titles))),
              '"probe1"\t' + "\t".join(str(i + 1) for i in range(len(titles))),
              "!series_matrix_table_end"]
    with gzip.open(path / filename, "wt") as handle:
        handle.write("\n".join(matrix) + "\n")
    platform = ["!platform_table_begin", "ID\tToxoDB", "probe1\t1.m00014",
                "!platform_table_end"]
    with gzip.open(path / "GPL7186_family.soft.gz", "wt") as handle:
        handle.write("\n".join(platform) + "\n")


def test_legacy_geo_arrays_resolve_probes_and_preserve_conditions(tmp_path):
    path = _acquired(tmp_path)
    resolve = lambda value: {"1.m00014": "TGME49_200010"}.get(value)
    _write_geo_matrix(path, "GSE19092_series_matrix.txt.gz",
                      ["asynchronous - 1", "1 hour release - 2"])
    cellcycle = EX.gse19092_cell_cycle(str(tmp_path), resolve=resolve,
                                      log=lambda *_: None)
    assert list(cellcycle) == ["cellcycle19092_async_r1", "cellcycle19092_1h_r2"]
    _write_geo_matrix(path, "GSE51780_series_matrix.txt.gz",
                      ["tachyzoite", "merozoite"])
    merozoite = EX.gse51780_merozoite(str(tmp_path), resolve=resolve,
                                     log=lambda *_: None)
    assert list(merozoite) == ["rna51780_tachy_r1", "rna51780_mero_r2"]


def test_gse168155_is_named_as_a_perturbation_not_direct_m6a(tmp_path):
    path = _acquired(tmp_path)
    frame = pd.DataFrame({"Name": ["TGME49_200010"],
                          "UT-1 - linear total RPKM": [4.0],
                          "IAA_7h-1 - linear total RPKM": [8.0],
                          "a statistic": [0.01]})
    frame.to_excel(path / "GSE168155_Matrix_table_processed_data.xlsx",
                   sheet_name="RPKM", index=False)
    out = EX.gse168155_rna_processing_perturbation(str(tmp_path), log=lambda *_: None)
    assert list(out) == ["cpsf4rna168155_ut_1", "cpsf4rna168155_iaa_7h_1"]
    assert not any("m6a" in column for column in out)


def test_gse200962_keeps_every_counted_condition(tmp_path):
    path = _acquired(tmp_path)
    pd.DataFrame({"gene": ["TGME49_200010", "TGME49_200020"],
                  "P2 with pH 8": [2, 4], "P5 without pH 8": [8, 16]}).to_csv(
        path / "GSE200962_gene_count_matrix_geo.csv.gz", index=False,
        compression="gzip")
    out = EX.gse200962_restriction_checkpoint(str(tmp_path), log=lambda *_: None)
    assert list(out) == ["restriction200962_p2_with_ph_8",
                         "restriction200962_p5_without_ph_8"]


# --------------------------------------------------------------------------- absent sources
def test_every_loader_returns_nothing_when_its_file_is_absent(tmp_path):
    """A dataset nobody has downloaded must not take the build down with it.

    Each of these reads one file from the archive; the cache is assembled from whatever is present,
    so a missing source is the ordinary case rather than the exceptional one. Returning an empty
    frame is what lets `load_all` compose the rest -- and the slot it would have filled then reads
    as empty, which is true.
    """
    import starplast.expression as E
    base = str(tmp_path)
    os.makedirs(os.path.join(base, "datasets"), exist_ok=True)
    loaders = [E.gse22258_stage, E.neuronal_differentiation,
               E.gse129869_host_context_ribosome_profiling, E.gse19092_cell_cycle,
               E.gse51780_merozoite]
    for loader in loaders:
        out = loader(base, log=lambda *a, **k: None)
        assert isinstance(out, pd.DataFrame) and out.empty, f"{loader.__name__} did not come back empty"


def test_a_series_matrix_with_a_header_and_no_table_yields_no_rows(tmp_path):
    """GEO series matrices carry their sample titles above the table. A file that has the titles and
    then ends -- a truncated download, or a series with no matrix published -- must give back the
    titles it did read and no data, rather than raising out of the csv parser."""
    import gzip
    import starplast.expression as E
    path = str(tmp_path / "truncated_series_matrix.txt.gz")
    with gzip.open(path, "wt") as fh:
        fh.write('!Sample_title\t"tachyzoite"\t"bradyzoite"\n')
        fh.write("!series_matrix_table_begin\n")
        fh.write("!series_matrix_table_end\n")
    frame, titles = E._geo_series_matrix(path)
    assert frame.empty
    assert titles == ["tachyzoite", "bradyzoite"], "the titles were lost with the table"


def test_a_source_that_is_present_but_unusable_yields_nothing_rather_than_half_a_dataset(tmp_path):
    """The other half of the ingestion contract, and the half that matters more.

    A missing file is obvious. A file that downloaded but is truncated, or whose table never
    appeared, is the one that quietly produces a column of nothing and fills a slot with it. Each of
    these loaders checks the shape it actually needs and returns empty when it is not there -- so
    the slot reads as empty, which is true, rather than as filled with a column of NaN.
    """
    import gzip
    import starplast.expression as E
    base = str(tmp_path)
    # The dated acquisition batch, named by `_acquired` -- not a directory called "acquired".
    acq = os.path.dirname(E._acquired(base, "x"))
    os.makedirs(acq, exist_ok=True)
    os.makedirs(os.path.join(base, "datasets"), exist_ok=True)

    # Two columns where three are needed: a series matrix with one sample instead of two.
    with gzip.open(os.path.join(base, "datasets",
                                "stagetranscriptome_GSE22258_series_matrix.txt.gz"), "wt") as fh:
        fh.write("ID_REF\tGSM1\n")
        fh.write("TGME49_000001\t5.0\n")
    assert E.gse22258_stage(base, log=lambda *a, **k: None).empty

    # A series matrix whose table never begins: the header is there and the data is not.
    empty_matrix = '!Sample_title\t"asynchronous"\n!series_matrix_table_begin\n' \
                   "!series_matrix_table_end\n"
    for name in ("GSE19092_series_matrix.txt.gz", "GSE51780_series_matrix.txt.gz"):
        with gzip.open(os.path.join(acq, name), "wt") as fh:
            fh.write(empty_matrix)
    with gzip.open(os.path.join(acq, "GPL7186_family.soft.gz"), "wt") as fh:
        fh.write("^PLATFORM = GPL7186\n!platform_table_begin\nID\tORF\n!platform_table_end\n")
    assert E.gse19092_cell_cycle(base, log=lambda *a, **k: None).empty
    assert E.gse51780_merozoite(base, log=lambda *a, **k: None).empty


def test_an_archive_with_no_count_tables_in_it_yields_nothing(tmp_path):
    """GSE129869 arrives as a tar of per-sample count tables, and the loader picks members by a
    name pattern. An archive that downloaded but holds none of them -- a partial upload, or a
    supplementary tar of something else entirely -- must produce no columns rather than an empty
    concat."""
    import tarfile
    import starplast.expression as E
    base = str(tmp_path)
    where = E._acquired(base, "GSE129869_RAW.tar")
    os.makedirs(os.path.dirname(where), exist_ok=True)
    readme = tmp_path / "README.txt"
    readme.write_text("nothing here matches the count-table pattern\n")
    with tarfile.open(where, "w") as archive:
        archive.add(str(readme), arcname="README.txt")
    assert E.gse129869_host_context_ribosome_profiling(base, log=lambda *a, **k: None).empty


def test_the_platform_map_skips_everything_before_its_table(tmp_path):
    """A GEO platform file opens with metadata lines before the probe table begins. They are not
    rows, and treating them as rows would map a probe id onto a licence string."""
    import gzip
    import starplast.expression as E
    path = str(tmp_path / "GPL7186_family.soft.gz")
    with gzip.open(path, "wt") as fh:
        fh.write("^PLATFORM = GPL7186\n")
        fh.write("!Platform_organism = Toxoplasma gondii\n")
        fh.write("!platform_table_begin\n")
        # The column is ToxoDB rather than ORF: the map reads a probe's PREVIOUS ToxoDB id and
        # resolves it forward, which is what makes a 2009 array usable against today's identifiers.
        fh.write("ID\tToxoDB\n")
        fh.write("probe1\tTGME49_000001\n")
        fh.write("!platform_table_end\n")
        fh.write("!Platform_trailing = ignored\n")
    out = E._gpl7186_gene_map(path)
    assert out.get("probe1") == "TGME49_000001"
    assert "^PLATFORM = GPL7186" not in out and len(out) == 1


def test_a_directory_inside_the_archive_is_not_read_as_a_count_table(tmp_path):
    """`extractfile` returns None for anything that is not a regular file, and a tar written by a
    submitter can carry a directory whose name matches the pattern. Reading None would raise from
    inside gzip rather than skipping the entry."""
    import tarfile
    import starplast.expression as E
    base = str(tmp_path)
    where = E._acquired(base, "GSE129869_RAW.tar")
    os.makedirs(os.path.dirname(where), exist_ok=True)
    folder = tmp_path / "GSM1_cRFP.RH1_count.tab.gz"
    folder.mkdir()
    with tarfile.open(where, "w") as archive:
        archive.add(str(folder), arcname="GSM1_cRFP.RH1_count.tab.gz", recursive=False)
    assert E.gse129869_host_context_ribosome_profiling(base, log=lambda *a, **k: None).empty


def test_a_workbook_sheet_without_fold_changes_is_skipped_and_an_empty_one_yields_nothing(tmp_path):
    """GSE168465 is one workbook with a sheet per time point, and submitters put other things in
    workbooks: a legend, a methods sheet, a sheet of counts with no comparison in it. A sheet with
    no log2FoldChange is not a time point and must be skipped rather than half-read -- and a
    workbook made entirely of those yields nothing rather than an empty concat."""
    import starplast.expression as E
    base = str(tmp_path)
    os.makedirs(os.path.join(base, "datasets"), exist_ok=True)
    path = os.path.join(base, "datasets",
                        "stagetranscriptome_GSE168465_DESeq2-Toxo-all-time-points.xlsx")
    with pd.ExcelWriter(path) as writer:
        pd.DataFrame({"gene": ["TGME49_000001"], "baseMean": [10.0]}).to_excel(
            writer, sheet_name="legend", index=False)
        pd.DataFrame({"gene": [], "note": []}).to_excel(writer, sheet_name="blank", index=False)
    assert E.neuronal_differentiation(base, log=lambda *a, **k: None).empty


# ------------------------------------------------------ GSE245775 differentiation ribosome profiling
def _gse245775(tmp_path, *, matrix=True, titles=None, counts=None, gsm_prefix=True,
               scramble_names=False):
    """A GSE245775-shaped deposit: per-gene count tables in a tar, plus the series matrix."""
    folder = (tmp_path / "datasets" / "quarantine" / "2026_08_16_unverified" / "Tg"
              / "stage_conversion_phenotype")
    folder.mkdir(parents=True, exist_ok=True)
    if titles is None:
        titles = [f"{strain}_{state}_{assay}_rep{rep}"
                  for strain in ("parental", "eIF1-2-KO")
                  for state in ("tachyzoite", "stress")
                  for assay in ("RIBOseq", "RNAseq")
                  for rep in (1, 2, 3)]
    accessions = [f"GSM{7848805 + i}" for i in range(len(titles))]
    body = counts or "TGME49_200010\t8\nTGME49_200020\t4\n"
    with tarfile.open(folder / "GSE245775_RAW.tar", "w") as archive:
        for i, (title, gsm) in enumerate(zip(titles, accessions)):
            payload = gzip.compress(body.encode())
            # Deliberately uninformative names: the mapping must come from the matrix.
            stem = f"sample_{i}" if scramble_names else title.replace("-", "_")
            name = (f"{gsm}_{stem}_counts.tabular.txt.gz" if gsm_prefix
                    else f"{stem}_counts.tabular.txt.gz")
            info = tarfile.TarInfo(name)
            info.size = len(payload)
            archive.addfile(info, io.BytesIO(payload))
    if matrix:
        with gzip.open(folder / "GSE245775_series_matrix.txt.gz", "wt") as fh:
            fh.write("!Sample_title\t" + "\t".join(f'"{t}"' for t in titles) + "\n")
            fh.write("!Sample_geo_accession\t" + "\t".join(f'"{a}"' for a in accessions) + "\n")
    return folder


def test_gse245775_yields_every_arm_and_its_efficiency(tmp_path):
    _gse245775(tmp_path)
    out = EX.gse245775_differentiation_ribosome_profiling(str(tmp_path), log=lambda *_: None)
    assert out.shape == (2, 36), "four arms x three assays x three replicates"
    assert {"rpf245775_parent_tachy_r1", "rna245775_parent_prebrady_r3",
            "te245775_eif12ko_tachy_r2"} <= set(out)


def test_the_stress_arm_is_named_pre_bradyzoite_not_bradyzoite(tmp_path):
    """The submitters call these `cell type: pre-bradyzoites` -- 48h at pH 8.3, not a tissue cyst.

    The column name is the only place a reader meets that distinction, so a rename to `brady` would
    quietly let this answer questions about chronic infection that 48 hours of alkaline stress cannot.
    """
    _gse245775(tmp_path)
    out = EX.gse245775_differentiation_ribosome_profiling(str(tmp_path), log=lambda *_: None)
    assert any("prebrady" in c for c in out)
    assert not any(c.split("245775_")[1].startswith("parent_brady") for c in out)


def test_samples_are_identified_by_accession_not_by_file_name(tmp_path):
    """File names here happen to be descriptive, which is exactly why they are not trusted."""
    _gse245775(tmp_path, scramble_names=True)
    out = EX.gse245775_differentiation_ribosome_profiling(str(tmp_path), log=lambda *_: None)
    assert out.shape == (2, 36)


def test_a_member_without_a_gsm_prefix_is_skipped(tmp_path):
    _gse245775(tmp_path, gsm_prefix=False)
    assert EX.gse245775_differentiation_ribosome_profiling(
        str(tmp_path), log=lambda *_: None).empty


def test_a_sample_title_of_an_unknown_shape_is_skipped(tmp_path):
    """A fifth arm nobody planned for is left out rather than guessed into an existing one."""
    titles = ["parental_tachyzoite_RIBOseq_rep1", "parental_tachyzoite_RNAseq_rep1",
              "some_other_condition_entirely"]
    _gse245775(tmp_path, titles=titles)
    out = EX.gse245775_differentiation_ribosome_profiling(str(tmp_path), log=lambda *_: None)
    assert set(out) == {"rpf245775_parent_tachy_r1", "rna245775_parent_tachy_r1",
                        "te245775_parent_tachy_r1"}


def test_efficiency_is_the_difference_of_the_logged_pair(tmp_path):
    """TE is RPF minus RNA in log space, which is the ratio -- not a second normalisation."""
    _gse245775(tmp_path, titles=["parental_tachyzoite_RIBOseq_rep1",
                                 "parental_tachyzoite_RNAseq_rep1"])
    out = EX.gse245775_differentiation_ribosome_profiling(str(tmp_path), log=lambda *_: None)
    expected = out["rpf245775_parent_tachy_r1"] - out["rna245775_parent_tachy_r1"]
    assert np.allclose(out["te245775_parent_tachy_r1"], expected)


def test_a_gsm_absent_from_the_matrix_contributes_nothing(tmp_path):
    folder = _gse245775(tmp_path)
    with gzip.open(folder / "GSE245775_series_matrix.txt.gz", "wt") as fh:
        fh.write('!Sample_title\t"parental_tachyzoite_RIBOseq_rep1"\n')
        fh.write('!Sample_geo_accession\t"GSM0000000"\n')
    assert EX.gse245775_differentiation_ribosome_profiling(
        str(tmp_path), log=lambda *_: None).empty


def test_no_archive_yields_nothing(tmp_path):
    assert EX.gse245775_differentiation_ribosome_profiling(
        str(tmp_path), log=lambda *_: None).empty


def test_no_series_matrix_yields_nothing(tmp_path):
    """Without the mapping there is no honest way to say which sample is which."""
    _gse245775(tmp_path, matrix=False)
    assert EX.gse245775_differentiation_ribosome_profiling(
        str(tmp_path), log=lambda *_: None).empty


def test_a_series_matrix_missing_the_accession_line_yields_nothing(tmp_path):
    folder = _gse245775(tmp_path)
    with gzip.open(folder / "GSE245775_series_matrix.txt.gz", "wt") as fh:
        fh.write('!Sample_title\t"parental_tachyzoite_RIBOseq_rep1"\n')
    assert EX.gse245775_differentiation_ribosome_profiling(
        str(tmp_path), log=lambda *_: None).empty


def test_accessions_are_resolved_before_joining(tmp_path):
    _gse245775(tmp_path, titles=["parental_tachyzoite_RIBOseq_rep1"],
               counts="TGGT1_100010\t8\nTGGT1_100020\t4\n")
    out = EX.gse245775_differentiation_ribosome_profiling(
        str(tmp_path), resolve=lambda a: {"TGGT1_100010": "TGME49_200010"}.get(a),
        log=lambda *_: None)
    assert list(out.index) == ["TGME49_200010"]


def test_a_directory_member_is_survived(tmp_path):
    folder = _gse245775(tmp_path)
    with tarfile.open(folder / "GSE245775_RAW.tar", "a") as archive:
        info = tarfile.TarInfo("GSM7848805_subdir")
        info.type = tarfile.DIRTYPE
        archive.addfile(info)
    out = EX.gse245775_differentiation_ribosome_profiling(str(tmp_path), log=lambda *_: None)
    assert out.shape == (2, 36)


# ------------------------------------------------------------------ GSE223620 BFD2 RIP-seq
def _bfd2(tmp_path, rows):
    folder = (tmp_path / "datasets" / "quarantine" / "2026_08_16_unverified" / "Tg"
              / "rna_binding_protein_targets")
    folder.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["GeneID", "HA_S_Input.ReadCount", "HA_S_IP.ReadCount"]).to_excel(
        folder / "GSE223620_ProcessedDataFile_BFD2.RIPseq.xls", index=False)
    return tmp_path


def test_bfd2_rip_is_the_ratio_of_the_two_libraries(tmp_path):
    """Equal share of each library is no enrichment, whatever the raw depths are."""
    _bfd2(tmp_path, [("TGME49_200010", 100, 400), ("TGME49_200020", 300, 1200)])
    out = EX.gse223620_bfd2_rip(str(tmp_path), log=lambda *_: None)
    assert abs(out["bfd2_rip_log2_ip_over_input"]).max() < 0.05


def test_an_enriched_transcript_is_positive(tmp_path):
    _bfd2(tmp_path, [("TGME49_200010", 10, 990), ("TGME49_200020", 990, 10)])
    out = EX.gse223620_bfd2_rip(str(tmp_path), log=lambda *_: None)
    assert out.loc["TGME49_200010", "bfd2_rip_log2_ip_over_input"] > 3
    assert out.loc["TGME49_200020", "bfd2_rip_log2_ip_over_input"] < -3


def test_genes_below_the_read_floor_are_dropped(tmp_path):
    """A ratio of two small numbers is noise, and a RIP is read for its tail."""
    _bfd2(tmp_path, [("TGME49_200010", 1, 2), ("TGME49_200020", 500, 500)])
    out = EX.gse223620_bfd2_rip(str(tmp_path), log=lambda *_: None)
    assert list(out.index) == ["TGME49_200020"]


def test_bfd2_accessions_are_resolved(tmp_path):
    _bfd2(tmp_path, [("TGGT1_100010", 100, 400), ("TGME49_200020", 100, 400)])
    out = EX.gse223620_bfd2_rip(str(tmp_path), resolve=lambda a: (
        "TGME49_200010" if a == "TGGT1_100010" else a), log=lambda *_: None)
    assert "TGME49_200010" in out.index


def test_no_bfd2_file_yields_nothing(tmp_path):
    assert EX.gse223620_bfd2_rip(str(tmp_path), log=lambda *_: None).empty


def test_a_bfd2_table_without_both_libraries_is_refused(tmp_path):
    folder = (tmp_path / "datasets" / "quarantine" / "2026_08_16_unverified" / "Tg"
              / "rna_binding_protein_targets")
    folder.mkdir(parents=True)
    pd.DataFrame({"GeneID": ["TGME49_200010"], "only": [1]}).to_excel(
        folder / "GSE223620_ProcessedDataFile_BFD2.RIPseq.xls", index=False)
    assert EX.gse223620_bfd2_rip(str(tmp_path), log=lambda *_: None).empty


def test_an_empty_library_does_not_divide_by_zero(tmp_path):
    _bfd2(tmp_path, [("TGME49_200010", 0, 40), ("TGME49_200020", 0, 40)])
    out = EX.gse223620_bfd2_rip(str(tmp_path), log=lambda *_: None)
    assert np.isfinite(out["bfd2_rip_log2_ip_over_input"]).all()
