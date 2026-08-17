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
    assert list(out["host_id"]) == ["PDCD6_HUMAN"]
    assert out["host_accession"].iloc[0] == "O75340"
    assert out["host_ip_enrichment_log2"].iloc[0] > 8


def test_a_contaminant_inside_the_group_is_dropped(tmp_path):
    """MaxQuant prefixes a group only when the LEADING entry is a contaminant, so keratin arrives
    mid-group behind an `sp|` header and survives the obvious filter. It is otherwise the four most
    enriched host proteins in this deposit."""
    _ip(tmp_path, [("sp|P02533|K1C14_HUMAN;CON__P02533", 1000.0, 1000.0, 1.0, 1.0, 20, 20),
                   ("sp|O75340|PDCD6_HUMAN", 900.0, 900.0, 1.0, 1.0, 9, 3)])
    out = H.read_ip(str(tmp_path))
    assert list(out["host_id"]) == ["PDCD6_HUMAN"]


def test_one_peptide_in_one_run_is_not_an_interaction(tmp_path):
    _ip(tmp_path, [("sp|O75340|PDCD6_HUMAN", 1000.0, 1000.0, 1.0, 1.0, 9, 1),
                   ("sp|Q14315|FLNC_HUMAN", 900.0, 900.0, 1.0, 1.0, 5, 4)])
    out = H.read_ip(str(tmp_path))
    assert list(out["host_id"]) == ["FLNC_HUMAN"]


def test_the_parasite_side_of_the_experiment_is_not_a_host_row(tmp_path):
    """The bait's own partners belong to the parasite table, and a group naming both is not a host
    protein."""
    _ip(tmp_path, [("TGGT1_254470", 1000.0, 1000.0, 1.0, 1.0, 9, 9),
                   ("sp|O75340|PDCD6_HUMAN", 900.0, 900.0, 1.0, 1.0, 9, 3)])
    out = H.read_ip(str(tmp_path))
    assert list(out["host_id"]) == ["PDCD6_HUMAN"]


def test_a_bridge_names_both_ends_and_its_evidence(tmp_path):
    _ip(tmp_path, [("sp|O75340|PDCD6_HUMAN", 1000.0, 1000.0, 1.0, 1.0, 9, 3)])
    b = H.bridges(str(tmp_path))
    assert set(b.columns) == {"gene_id", "host_id", "host_ip_enrichment_log2", "evidence"}
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
    top = b.sort_values("host_ip_enrichment_log2", ascending=False)["host_id"].tolist()
    assert top[0] == "PDCD6_HUMAN", top[:3]
    escrt = b["host_id"].str.match(r"^(PDCD6|PDC6I|VPS28|CHMP\d|VPS4[AB]|TS101)")
    assert escrt.sum() >= 3, "the ESCRT machinery is not in the bridge"
    p = mannwhitneyu(b.loc[escrt, "host_ip_enrichment_log2"],
                     b.loc[~escrt, "host_ip_enrichment_log2"], alternative="greater").pvalue
    assert p < 0.05, f"ESCRT is not enriched over the rest, p = {p:.3f}"
