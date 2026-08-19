#!/usr/bin/env python3
"""The host table and the first bridge, and the two filters that decide whether an IP is readable.

Both were arrived at by getting it wrong. A contaminant hides inside a protein group rather than at
its head, and a single peptide in a single run is an identification rather than an interaction.
"""
from __future__ import annotations

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import host as H  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _ip(tmp_path, rows):
    """A MaxQuant proteinGroups table shaped like the deposit's."""
    folder = tmp_path / H.MYR1_IP["folder"]
    folder.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(rows, columns=[
        "Protein IDs", "LFQ intensity M1", "LFQ intensity M2",
        "LFQ intensity R1", "LFQ intensity R2", "Unique peptides M1", "Unique peptides M2"])
    frame.to_csv(folder / H.MYR1_IP["file"], sep="\t", index=False)
    return tmp_path


def test_a_host_protein_enriched_over_the_control_is_kept(tmp_path):
    _ip(tmp_path, [("sp|O75340|PDCD6_HUMAN", 1000.0, 1000.0, 1.0, 1.0, 9, 3)])
    out = H.read_ip(str(tmp_path))
    assert list(out["host_id"]) == ["O75340"], "the row key is the accession, not the entry name"
    assert out["host_name"].iloc[0] == "PDCD6_HUMAN"
    assert out["host_ip_enrichment_log2"].iloc[0] > 8


def test_a_contaminant_inside_the_group_is_dropped(tmp_path):
    """MaxQuant prefixes a group only when the LEADING entry is a contaminant, so keratin arrives
    mid-group behind an `sp|` header and survives the obvious filter. It is otherwise the four most
    enriched host proteins in this deposit."""
    _ip(tmp_path, [("sp|P02533|K1C14_HUMAN;CON__P02533", 1000.0, 1000.0, 1.0, 1.0, 20, 20),
                   ("sp|O75340|PDCD6_HUMAN", 900.0, 900.0, 1.0, 1.0, 9, 3)])
    out = H.read_ip(str(tmp_path))
    assert list(out["host_name"]) == ["PDCD6_HUMAN"]


def test_one_peptide_in_one_run_is_not_an_interaction(tmp_path):
    _ip(tmp_path, [("sp|O75340|PDCD6_HUMAN", 1000.0, 1000.0, 1.0, 1.0, 9, 1),
                   ("sp|Q14315|FLNC_HUMAN", 900.0, 900.0, 1.0, 1.0, 5, 4)])
    out = H.read_ip(str(tmp_path))
    assert list(out["host_name"]) == ["FLNC_HUMAN"]


def test_the_parasite_side_of_the_experiment_is_not_a_host_row(tmp_path):
    """The bait's own partners belong to the parasite table, and a group naming both is not a host
    protein."""
    _ip(tmp_path, [("TGGT1_254470", 1000.0, 1000.0, 1.0, 1.0, 9, 9),
                   ("sp|O75340|PDCD6_HUMAN", 900.0, 900.0, 1.0, 1.0, 9, 3)])
    out = H.read_ip(str(tmp_path))
    assert list(out["host_name"]) == ["PDCD6_HUMAN"]


def test_a_bridge_names_both_ends_and_its_evidence(tmp_path):
    _ip(tmp_path, [("sp|O75340|PDCD6_HUMAN", 1000.0, 1000.0, 1.0, 1.0, 9, 3)])
    b = H.bridges(str(tmp_path))
    assert set(b.columns) == {"gene_id", "host_id", "host_name",
                              "host_ip_enrichment_log2", "evidence"}
    assert b["gene_id"].iloc[0] == H.MYR1_IP["bait"]
    assert "PXD016383" in b["evidence"].iloc[0]


def test_no_deposit_yields_nothing(tmp_path):
    assert H.read_ip(str(tmp_path)).empty and H.bridges(str(tmp_path)).empty


def test_a_table_without_protein_ids_is_refused(tmp_path):
    folder = tmp_path / H.MYR1_IP["folder"]
    folder.mkdir(parents=True)
    pd.DataFrame({"something": [1]}).to_csv(folder / H.MYR1_IP["file"], sep="\t", index=False)
    assert H.read_ip(str(tmp_path)).empty


def test_a_table_missing_the_quantitative_columns_is_refused(tmp_path):
    folder = tmp_path / H.MYR1_IP["folder"]
    folder.mkdir(parents=True)
    pd.DataFrame({"Protein IDs": ["sp|O75340|PDCD6_HUMAN"]}).to_csv(
        folder / H.MYR1_IP["file"], sep="\t", index=False)
    assert H.read_ip(str(tmp_path)).empty


def test_a_deposit_with_no_host_rows_yields_nothing(tmp_path):
    _ip(tmp_path, [("TGGT1_254470", 1000.0, 1000.0, 1.0, 1.0, 9, 9)])
    assert H.read_ip(str(tmp_path)).empty and H.bridges(str(tmp_path)).empty


def test_load_returns_nothing_when_the_table_is_not_built(tmp_path):
    assert H.load(str(tmp_path)).empty


# --------------------------------------------------------------------------- against the deposit
@pytest.mark.skipif(
    not os.path.exists(os.path.join(ROOT, "starplast", "data", H.BRIDGE_TABLE)),
    reason="host bridge not built")
def test_the_escrt_machinery_leads_the_host_side():
    """The check that decided this slot, and the one whose absence had it written off.

    A 2026 paper reports Toxoplasma GRA8 engaging host ALG-2 at the vacuole. ALG-2 is PDCD6, and it
    is rank 1 here; its partner ALIX and the ESCRT-I subunit VPS28 follow. Without an independent
    known-positive, an enriched host list cannot be told apart from abundant proteins that stick.
    """
    from scipy.stats import mannwhitneyu
    # Only the rows scored as a fold change: the SAINT rows carry a probability in their own
    # column, and ranking the two together would be comparing a 1.0 with a 4.11.
    b = H.load(ROOT, H.BRIDGE_TABLE).dropna(subset=["host_ip_enrichment_log2"])
    assert len(b) > 100
    top = b.sort_values("host_ip_enrichment_log2", ascending=False)["host_name"].tolist()
    assert "PDCD6" in str(top[0]), top[:3]
    escrt = b["host_name"].astype(str).str.match(
        r"^(PDCD6|PDC6I|PDCD6IP|VP37A|VPS28|CHMP\d|VPS4[AB]|TS101|TSG101|PEF1)")
    assert escrt.sum() >= 3, "the ESCRT machinery is not in the bridge"
    p = mannwhitneyu(b.loc[escrt, "host_ip_enrichment_log2"],
                     b.loc[~escrt, "host_ip_enrichment_log2"], alternative="greater").pvalue
    assert p < 0.05, f"ESCRT is not enriched over the rest, p = {p:.3f}"


def _dia(tmp_path, rows):
    """A differential-abundance table shaped like PXD080696's combined_results.csv."""
    import os as _os
    folder = tmp_path / _os.path.dirname(H.DIA_FILE)
    folder.mkdir(parents=True, exist_ok=True)
    spec = H.DIA_IPS[0]
    frame = pd.DataFrame(rows, columns=["gene_symbol", "uniprot_id", "organism",
                                        spec["logfc"], spec["padj"]])
    frame.to_csv(tmp_path / H.DIA_FILE, index=False)
    return tmp_path


def test_a_dia_bait_keeps_what_clears_both_thresholds(tmp_path):
    _dia(tmp_path, [("=\"PDCD6\"", "O75340", "Homo sapiens", 4.11, 0.001),
                    ("=\"NOISE\"", "P00001", "Homo sapiens", 4.11, 0.9),
                    ("=\"SMALL\"", "P00002", "Homo sapiens", 0.2, 0.001)])
    out = H.read_dia(str(tmp_path), H.DIA_IPS[0])
    assert list(out["host_name"]) == ["PDCD6"]


def test_the_spreadsheet_quoting_is_unwrapped(tmp_path):
    """Symbols arrive as ="PDCD6" -- a spreadsheet stopping Excel reading them as formulas."""
    _dia(tmp_path, [("=\"TSG101\"", "Q99816", "Homo sapiens", 3.09, 0.001)])
    out = H.read_dia(str(tmp_path), H.DIA_IPS[0])
    assert list(out["host_name"]) == ["TSG101"]


def test_parasite_rows_are_not_host_partners_in_the_dia_table(tmp_path):
    _dia(tmp_path, [("=\"GRA8\"", "A0A125", "Toxoplasma gondii", 4.78, 0.001),
                    ("=\"PDCD6\"", "O75340", "Homo sapiens", 4.11, 0.001)])
    out = H.read_dia(str(tmp_path), H.DIA_IPS[0])
    assert list(out["host_name"]) == ["PDCD6"]


def test_no_dia_file_yields_nothing(tmp_path):
    assert H.read_dia(str(tmp_path), H.DIA_IPS[0]).empty


def test_a_dia_table_without_the_contrast_columns_is_refused(tmp_path):
    import os as _os
    folder = tmp_path / _os.path.dirname(H.DIA_FILE)
    folder.mkdir(parents=True)
    pd.DataFrame({"gene_symbol": ["PDCD6"], "organism": ["Homo sapiens"]}).to_csv(
        tmp_path / H.DIA_FILE, index=False)
    assert H.read_dia(str(tmp_path), H.DIA_IPS[0]).empty


def test_a_dia_bait_with_nothing_significant_yields_nothing(tmp_path):
    _dia(tmp_path, [("=\"NOISE\"", "P00001", "Homo sapiens", 0.1, 0.9)])
    assert H.read_dia(str(tmp_path), H.DIA_IPS[0]).empty


def test_all_bridges_names_the_bait_on_every_row(tmp_path):
    _dia(tmp_path, [("=\"PDCD6\"", "O75340", "Homo sapiens", 4.11, 0.001)])
    b = H.all_bridges(str(tmp_path))
    assert set(b["gene_id"]) <= {s["bait"] for s in H.DIA_IPS}
    assert b["evidence"].str.contains("PXD080696").all()


def test_all_bridges_is_empty_when_nothing_is_downloaded(tmp_path):
    assert H.all_bridges(str(tmp_path)).empty


@pytest.mark.skipif(
    not os.path.exists(os.path.join(ROOT, "starplast", "data", H.BRIDGE_TABLE)),
    reason="host bridge not built")
def test_independent_baits_converge_on_the_escrt_machinery():
    """The claim a multi-bait bridge can make and a single IP cannot.

    PDCD6 -- ALG-2 -- is reached by MYR1, EAF1 and GRA35 independently, and a 2026 paper reports
    Toxoplasma GRA8 engaging it at the vacuole. Convergence across baits is the evidence; any one of
    them alone is a list.
    """
    from scipy.stats import fisher_exact
    b = H.load(ROOT, H.BRIDGE_TABLE)
    assert b["gene_id"].nunique() >= 4, "the bridge lost baits"
    reach = b.groupby("host_id")["gene_id"].nunique()
    names = b.drop_duplicates("host_id").set_index("host_id")["host_name"].astype(str).str.upper()
    escrt = names.str.match(r"^(PDCD6|PDC6I|PDCD6IP|ALIX|VP37A|VPS37|TSG101|TS101|VPS28|CHMP\d"
                            r"|VPS4[AB]|VP4[AB]|SNF8|VPS25|VPS36|IST1|VTA1|PEF1)")
    multi = reach >= 3
    odds, p = fisher_exact([[int((multi & escrt).sum()), int((multi & ~escrt).sum())],
                            [int((~multi & escrt).sum()), int((~multi & ~escrt).sum())]])
    assert p < 1e-6, f"odds {odds:.1f}, p {p:.1e}"
    assert int((multi & escrt).sum()) >= 5


@pytest.mark.skipif(
    not os.path.exists(os.path.join(ROOT, "starplast", "data", H.BRIDGE_TABLE)),
    reason="host bridge not built")
def test_the_eaf1_bait_pulls_the_escrt_pathway():
    """Five partners, all ESCRT, including the TSG101 the imaging screen measured."""
    b = H.load(ROOT, H.BRIDGE_TABLE)
    eaf1 = b[b["gene_id"] == "TGME49_225160"]
    if eaf1.empty:
        pytest.skip("the EAF1 bait is not in the bridge")
    got = set(eaf1["host_name"].astype(str).str.upper())
    assert {"PDCD6", "TSG101", "PDCD6IP", "CHMP4B"} <= got, got


def test_the_bridge_is_keyed_on_an_identifier_every_deposit_carries():
    """Three deposits, three naming conventions. Keyed on the readable name, one protein would count
    as three and the convergence claim would be an artefact of formatting."""
    b = H.load(ROOT, H.BRIDGE_TABLE)
    if b.empty:
        pytest.skip("host bridge not built")
    # UniProt accessions are 6 or 10 characters, letter-led, and never carry an underscore or a
    # lowercase letter -- which is exactly what distinguishes them from the entry names and symbols
    # the three deposits otherwise use.
    # UniProt's own accession grammar. An entry name or a gene symbol matches none of it, which is
    # what makes this a check on the key rather than on the formatting.
    ok = b["host_id"].str.fullmatch(
        r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")
    assert ok.all(), b.loc[~ok, "host_id"].head().tolist()


def _replicated(tmp_path, rows, spec=None, member=None):
    """A two-replicate-block sheet shaped like the GRA64 pulldowns."""
    import io as _io
    import zipfile as _zip
    spec = spec or H.REPLICATE_IPS[0]
    frame = pd.DataFrame(rows, columns=["Protein Accessions", "Protein Fold Change",
                                        "Protein Accessions ", "Protein Fold Change.1"])
    buf = _io.BytesIO()
    with pd.ExcelWriter(buf) as w:
        pd.DataFrame([["Table"]]).to_excel(w, sheet_name=spec["sheet"], index=False, header=False)
        frame.to_excel(w, sheet_name=spec["sheet"], index=False, startrow=1)
    folder = tmp_path / H.ESCRT_ROOT / os.path.dirname(spec["archive"])
    folder.mkdir(parents=True, exist_ok=True)
    with _zip.ZipFile(tmp_path / H.ESCRT_ROOT / spec["archive"], "w") as z:
        z.writestr(member or spec["member"], buf.getvalue())
    return tmp_path


def test_a_protein_enriched_in_both_replicates_is_kept(tmp_path):
    _replicated(tmp_path, [("O75340", 4.1, "O75340", 3.8)])
    out = H.read_replicated(str(tmp_path), H.REPLICATE_IPS[0])
    assert list(out["host_id"]) == ["O75340"]
    assert out["host_ip_enrichment_log2"].iloc[0] == pytest.approx(3.95)


def test_a_protein_enriched_in_only_one_replicate_is_dropped(tmp_path):
    _replicated(tmp_path, [("O75340", 4.1, "O75340", 0.2), ("Q99816", 3.0, "Q99816", 2.8)])
    out = H.read_replicated(str(tmp_path), H.REPLICATE_IPS[0])
    assert list(out["host_id"]) == ["Q99816"]


def test_the_two_replicate_blocks_are_intersected_not_read_row_wise(tmp_path):
    """The blocks are sorted differently, so row 1 of one is not row 1 of the other. Reading them
    row-wise would pair replicate 1 of one protein with replicate 2 of another."""
    _replicated(tmp_path, [("O75340", 4.1, "Q99816", 3.0), ("Q99816", 3.2, "O75340", 3.8)])
    out = H.read_replicated(str(tmp_path), H.REPLICATE_IPS[0])
    assert set(out["host_id"]) == {"O75340", "Q99816"}
    assert out.set_index("host_id").loc["O75340", "host_ip_enrichment_log2"] == pytest.approx(3.95)


def test_an_ambiguous_protein_group_keeps_its_leading_accession(tmp_path):
    """`P08134; P61586` is one peptide matching several proteins, not a protein called that."""
    _replicated(tmp_path, [("P08134; P61586", 4.1, "P08134; P61586", 3.8)])
    out = H.read_replicated(str(tmp_path), H.REPLICATE_IPS[0])
    assert list(out["host_id"]) == ["P08134"]


def test_parasite_accessions_are_not_host_rows_in_a_replicated_sheet(tmp_path):
    _replicated(tmp_path, [("TGME49_264660-t26_1-p1", 4.1, "TGME49_264660-t26_1-p1", 3.9),
                           ("O75340", 4.1, "O75340", 3.8)])
    out = H.read_replicated(str(tmp_path), H.REPLICATE_IPS[0])
    assert list(out["host_id"]) == ["O75340"]


def test_a_replicated_sheet_of_only_parasite_rows_yields_nothing(tmp_path):
    _replicated(tmp_path, [("TGME49_264660-t26_1-p1", 4.1, "TGME49_264660-t26_1-p1", 3.9)])
    assert H.read_replicated(str(tmp_path), H.REPLICATE_IPS[0]).empty


def test_a_replicated_sheet_with_nothing_reproducible_yields_nothing(tmp_path):
    _replicated(tmp_path, [("O75340", 0.1, "O75340", 0.2)])
    assert H.read_replicated(str(tmp_path), H.REPLICATE_IPS[0]).empty


def test_a_sheet_without_two_replicate_blocks_is_refused(tmp_path):
    import io as _io
    import zipfile as _zip
    spec = H.REPLICATE_IPS[0]
    buf = _io.BytesIO()
    with pd.ExcelWriter(buf) as w:
        pd.DataFrame([["Table"]]).to_excel(w, sheet_name=spec["sheet"], index=False, header=False)
        pd.DataFrame({"Protein Accessions": ["O75340"], "Protein Fold Change": [4.0]}).to_excel(
            w, sheet_name=spec["sheet"], index=False, startrow=1)
    folder = tmp_path / H.ESCRT_ROOT / os.path.dirname(spec["archive"])
    folder.mkdir(parents=True, exist_ok=True)
    with _zip.ZipFile(tmp_path / H.ESCRT_ROOT / spec["archive"], "w") as z:
        z.writestr(spec["member"], buf.getvalue())
    assert H.read_replicated(str(tmp_path), spec).empty


def test_a_missing_member_in_the_archive_is_refused(tmp_path):
    _replicated(tmp_path, [("O75340", 4.1, "O75340", 3.8)], member="something_else.xlsx")
    assert H.read_replicated(str(tmp_path), H.REPLICATE_IPS[0]).empty


def test_no_replicated_archive_yields_nothing(tmp_path):
    assert H.read_replicated(str(tmp_path), H.REPLICATE_IPS[0]).empty


def test_all_bridges_includes_the_replicated_pulldowns(tmp_path):
    """Four baits reach the bridge by three different readers; a bait whose reader returns nothing
    is skipped rather than contributing an empty block."""
    _replicated(tmp_path, [("O75340", 4.1, "O75340", 3.8)])
    b = H.all_bridges(str(tmp_path))
    assert not b.empty
    assert b["evidence"].str.contains("PMC9426488").any()
    assert set(b["host_id"]) == {"O75340"}


def _pv(tmp_path, rows):
    """The vacuole-proximity summary sheet, which lists parasite and host proteins together."""
    import io as _io
    import zipfile as _zip
    frame = pd.DataFrame(rows, columns=["Mammalian Protein", "Uniprot ID", "Tz-HFF", "Bz-HFF",
                                        "Neuron", "Average"])
    buf = _io.BytesIO()
    with pd.ExcelWriter(buf) as w:
        pd.DataFrame([["Log2 Fold Change vs Control"]]).to_excel(
            w, sheet_name=H.PV_UPTAKE["sheet"], index=False, header=False)
        frame.to_excel(w, sheet_name=H.PV_UPTAKE["sheet"], index=False, startrow=1)
    folder = tmp_path / H.ESCRT_ROOT / os.path.dirname(H.PV_UPTAKE["archive"])
    folder.mkdir(parents=True, exist_ok=True)
    with _zip.ZipFile(tmp_path / H.ESCRT_ROOT / H.PV_UPTAKE["archive"], "w") as z:
        z.writestr(H.PV_UPTAKE["member"], buf.getvalue())
    return tmp_path


def test_vacuole_enrichment_keeps_host_rows_only(tmp_path):
    """The parasite rows are the authors' positive control -- the dense granule proteins top the
    sheet -- and belong to no host table."""
    _pv(tmp_path, [("GRA1", "TGME49_270250", 6.0, 6.0, 6.0, 6.43),
                   ("PDCD6 / ALG-2", "O75340", 6.0, 6.0, 5.0, 5.90)])
    out = H.pv_enrichment(str(tmp_path))
    assert list(out["host_id"]) == ["O75340"]
    assert out["pv_enrichment_log2"].iloc[0] == pytest.approx(5.90)


def test_a_row_whose_identifier_is_not_an_accession_is_dropped(tmp_path):
    _pv(tmp_path, [("something", "not an id", 1.0, 1.0, 1.0, 1.0),
                   ("PDCD6", "O75340", 6.0, 6.0, 5.0, 5.90)])
    assert list(H.pv_enrichment(str(tmp_path))["host_id"]) == ["O75340"]


def test_no_vacuole_sheet_yields_nothing(tmp_path):
    assert H.pv_enrichment(str(tmp_path)).empty


def test_a_vacuole_archive_without_the_member_yields_nothing(tmp_path):
    import zipfile as _zip
    folder = tmp_path / H.ESCRT_ROOT / os.path.dirname(H.PV_UPTAKE["archive"])
    folder.mkdir(parents=True)
    with _zip.ZipFile(tmp_path / H.ESCRT_ROOT / H.PV_UPTAKE["archive"], "w") as z:
        z.writestr("other.xlsx", b"x")
    assert H.pv_enrichment(str(tmp_path)).empty


def test_a_vacuole_sheet_with_too_few_columns_is_refused(tmp_path):
    import io as _io
    import zipfile as _zip
    buf = _io.BytesIO()
    with pd.ExcelWriter(buf) as w:
        pd.DataFrame([["header"]]).to_excel(w, sheet_name=H.PV_UPTAKE["sheet"], index=False,
                                            header=False)
        pd.DataFrame({"a": ["PDCD6"], "b": ["O75340"]}).to_excel(
            w, sheet_name=H.PV_UPTAKE["sheet"], index=False, startrow=1)
    folder = tmp_path / H.ESCRT_ROOT / os.path.dirname(H.PV_UPTAKE["archive"])
    folder.mkdir(parents=True, exist_ok=True)
    with _zip.ZipFile(tmp_path / H.ESCRT_ROOT / H.PV_UPTAKE["archive"], "w") as z:
        z.writestr(H.PV_UPTAKE["member"], buf.getvalue())
    assert H.pv_enrichment(str(tmp_path)).empty


def test_the_host_table_joins_identity_and_properties(tmp_path):
    """A protein reached by a bait but never measured at the vacuole keeps a NaN, which is the
    difference between 'not at the vacuole' and 'this experiment did not see it'."""
    _dia(tmp_path, [("=\"PDCD6\"", "O75340", "Homo sapiens", 4.11, 0.001),
                    ("=\"TSG101\"", "Q99816", "Homo sapiens", 3.09, 0.001)])
    _pv(tmp_path, [("PDCD6", "O75340", 6.0, 6.0, 5.0, 5.90)])
    h = H.host_table(str(tmp_path)).set_index("host_id")
    assert h.loc["O75340", "pv_enrichment_log2"] == pytest.approx(5.90)
    assert pd.isna(h.loc["Q99816", "pv_enrichment_log2"])


def test_no_bridges_means_no_host_table(tmp_path):
    assert H.host_table(str(tmp_path)).empty


def _orthogonal(tmp_path, rows, spec=None):
    """A SAINT/MiST/CompPASS score table as the manuscript supplements publish them."""
    spec = spec or H.ORTHOGONAL_IPS[0]
    folder = tmp_path / H.ESCRT_ROOT / os.path.dirname(spec["file"])
    folder.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows, columns=["genesymbol", "SAINT_AvgP"]).to_csv(
        tmp_path / H.ESCRT_ROOT / spec["file"], index=False)
    return tmp_path


def test_saint_keeps_what_the_method_itself_calls(tmp_path):
    _dia(tmp_path, [("=\"PDCD6\"", "O75340", "Homo sapiens", 4.11, 0.001),
                    ("=\"WEAK\"", "P00001", "Homo sapiens", 4.11, 0.001)])
    _orthogonal(tmp_path, [("PDCD6", 1.0), ("WEAK", 0.4)])
    out = H.read_orthogonal(str(tmp_path), H.ORTHOGONAL_IPS[0])
    assert list(out["host_id"]) == ["O75340"]
    assert out["saint_avgp"].iloc[0] == 1.0


def test_a_saint_probability_is_not_put_in_the_log2_column(tmp_path):
    """A 1.0 beside a 4.11 in one column would read as the same quantity twice."""
    _dia(tmp_path, [("=\"PDCD6\"", "O75340", "Homo sapiens", 4.11, 0.001)])
    _orthogonal(tmp_path, [("PDCD6", 1.0)])
    out = H.read_orthogonal(str(tmp_path), H.ORTHOGONAL_IPS[0])
    assert "host_ip_enrichment_log2" not in out.columns
    b = H.all_bridges(str(tmp_path))
    saint = b[b["evidence"].str.contains("SAINT")]
    assert saint["host_ip_enrichment_log2"].isna().all()
    assert saint["saint_avgp"].notna().all()


def test_a_prey_with_no_accession_is_the_parasite_side(tmp_path):
    """The orthogonal tables score parasite preys too; those belong to no host table."""
    _dia(tmp_path, [("=\"PDCD6\"", "O75340", "Homo sapiens", 4.11, 0.001)])
    _orthogonal(tmp_path, [("PDCD6", 1.0), ("TGRH88_044880", 1.0)])
    out = H.read_orthogonal(str(tmp_path), H.ORTHOGONAL_IPS[0])
    assert list(out["host_id"]) == ["O75340"]


def test_the_symbol_map_comes_from_the_experiments_own_report(tmp_path):
    _dia(tmp_path, [("=\"PDCD6\"", "O75340", "Homo sapiens", 4.11, 0.001)])
    assert H._symbol_to_accession(str(tmp_path))["PDCD6"] == "O75340"


def test_no_symbol_map_without_the_report(tmp_path):
    assert H._symbol_to_accession(str(tmp_path)) == {}


def test_a_report_without_the_mapping_columns_gives_no_map(tmp_path):
    import os as _os
    folder = tmp_path / _os.path.dirname(H.DIA_FILE)
    folder.mkdir(parents=True)
    pd.DataFrame({"something": [1]}).to_csv(tmp_path / H.DIA_FILE, index=False)
    assert H._symbol_to_accession(str(tmp_path)) == {}


def test_no_orthogonal_file_yields_nothing(tmp_path):
    assert H.read_orthogonal(str(tmp_path), H.ORTHOGONAL_IPS[0]).empty


def test_an_orthogonal_table_without_saint_is_refused(tmp_path):
    spec = H.ORTHOGONAL_IPS[0]
    folder = tmp_path / H.ESCRT_ROOT / os.path.dirname(spec["file"])
    folder.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"genesymbol": ["PDCD6"], "logFC": [4.1]}).to_csv(
        tmp_path / H.ESCRT_ROOT / spec["file"], index=False)
    assert H.read_orthogonal(str(tmp_path), spec).empty


@pytest.mark.skipif(
    not os.path.exists(os.path.join(ROOT, "starplast", "data", H.BRIDGE_TABLE)),
    reason="host bridge not built")
def test_alg2_is_reached_by_every_bait():
    """Five parasite proteins, five experiments, one host protein in common.

    ALG-2 is the thing this map keeps arriving at from unrelated directions, and a bridge exists to
    make that visible rather than to hold any one pulldown.
    """
    b = H.load(ROOT, H.BRIDGE_TABLE)
    reach = b.groupby("host_id")["gene_id"].nunique()
    assert reach.get("O75340", 0) == b["gene_id"].nunique(), reach.sort_values().tail().to_dict()


# --------------------------------------------------------------------------- host tissue proteomes
def _rbc_book(tmp_path, rows=None, fractions=None):
    """The red-cell proteome as its paper ships it: one sheet per fraction."""
    folder = tmp_path / "host" / H.ERYTHROCYTE[0] / H.ERYTHROCYTE[1]
    folder.mkdir(parents=True)
    rows = rows if rows is not None else {
        "Membrane extract": [("P11277", "SPTB", 500), ("P02724", "GYPA", 120)],
        "Cytoplasmic extract": [("P69905", "HBA1; HBA2", 900), ("P11277", "SPTB", 40)]}
    with pd.ExcelWriter(folder / H.ERYTHROCYTE[2]) as writer:
        for sheet, values in rows.items():
            pd.DataFrame(values, columns=["Accession", "Gene", "# PSMs"]).to_excel(
                writer, sheet_name=sheet, index=False)
    return str(tmp_path)


def test_the_two_fractions_stay_two_measurements(tmp_path):
    """A protein at the surface the merozoite invades through and one in the haemoglobin around it
    are not the same observation, so they are not the same column."""
    d = H.erythrocyte_proteome(_rbc_book(tmp_path), log=lambda *a: None).set_index("host_id")
    assert d.loc["P11277", "rbc_membrane_psms"] == 500
    assert d.loc["P11277", "rbc_cytoplasm_psms"] == 40
    assert pd.isna(d.loc["P02724", "rbc_cytoplasm_psms"]), "absent from a fraction is not zero in it"


def test_a_row_naming_several_genes_keeps_the_string_it_was_given(tmp_path):
    """`HBA1; HBA2` is what the source says. Choosing one would be inventing a fact about which."""
    d = H.erythrocyte_proteome(_rbc_book(tmp_path), log=lambda *a: None).set_index("host_id")
    assert d.loc["P69905", "host_name"] == "HBA1; HBA2"


def test_a_missing_or_wrong_shaped_erythrocyte_file_is_empty(tmp_path):
    assert H.erythrocyte_proteome(str(tmp_path), log=lambda *a: None).empty
    root = _rbc_book(tmp_path, rows={"Membrane extract": []})
    folder = os.path.join(root, "host", H.ERYTHROCYTE[0], H.ERYTHROCYTE[1])
    pd.DataFrame({"something": [1]}).to_excel(os.path.join(folder, H.ERYTHROCYTE[2]),
                                              sheet_name="Membrane extract", index=False)
    assert H.erythrocyte_proteome(root, log=lambda *a: None).empty
