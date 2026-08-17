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
    b = H.load(ROOT, H.BRIDGE_TABLE)
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
