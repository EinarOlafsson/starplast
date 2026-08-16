#!/usr/bin/env python3
"""Antibody epitopes, and why they are not the epitope count already in the map."""
from __future__ import annotations

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import iedb as I  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _table(tmp_path, rows, header="gene_id\tn_bcell_epitopes"):
    data = tmp_path / "starplast" / "data"
    data.mkdir(parents=True, exist_ok=True)
    (data / I.TABLE).write_text(header + "\n" + "\n".join(f"{g}\t{v}" for g, v in rows) + "\n")
    return tmp_path


def test_counts_are_read_per_gene(tmp_path):
    _table(tmp_path, [("TGME49_233460", 47)])
    out = I.bcell_epitopes(str(tmp_path), log=lambda *_: None)
    assert out.loc["TGME49_233460", I.COLUMN] == 47


def test_accessions_go_through_the_identity_layer(tmp_path):
    _table(tmp_path, [("TGGT1_100010", 5)])
    out = I.bcell_epitopes(str(tmp_path), log=lambda *_: None,
                           resolve=lambda g: {"TGGT1_100010": "TGME49_200010"}.get(g))
    assert list(out.index) == ["TGME49_200010"]


def test_an_unresolvable_accession_keeps_its_own_name(tmp_path):
    _table(tmp_path, [("TGME49_999999", 5)])
    out = I.bcell_epitopes(str(tmp_path), resolve=lambda g: None, log=lambda *_: None)
    assert list(out.index) == ["TGME49_999999"]


def test_two_rows_for_one_gene_keep_the_larger(tmp_path):
    _table(tmp_path, [("TGME49_233460", 5), ("TGME49_233460", 47)])
    out = I.bcell_epitopes(str(tmp_path), log=lambda *_: None)
    assert out.loc["TGME49_233460", I.COLUMN] == 47


def test_a_one_column_table_is_refused(tmp_path):
    data = tmp_path / "starplast" / "data"
    data.mkdir(parents=True)
    (data / I.TABLE).write_text("gene_id\nTGME49_233460\n")
    assert I.bcell_epitopes(str(tmp_path), log=lambda *_: None).empty


def test_nothing_downloaded_yields_nothing(tmp_path):
    assert I.bcell_epitopes(str(tmp_path), log=lambda *_: None).empty


@pytest.mark.skipif(
    not os.path.exists(os.path.join(ROOT, "starplast", "data", "nodes.parquet")),
    reason="node table not present")
def test_the_serodiagnostic_antigens_are_the_top_of_the_column():
    """SAG1 -- ToxoDB calls it SRS29B -- GRA6, GRA7 and GRA1 are what Toxoplasma serology kits use."""
    n = pd.read_parquet(os.path.join(ROOT, "starplast", "data", "nodes.parquet")).set_index("gene_id")
    if I.COLUMN not in n.columns:
        pytest.skip("antibody epitope column not merged")
    top = n[I.COLUMN].dropna().sort_values(ascending=False).head(5).index
    products = " ".join(n.loc[top, "product"].astype(str)).upper()
    assert "SRS29B" in products, products
    assert sum(p in products for p in ("GRA6", "GRA7", "GRA1")) >= 2, products


@pytest.mark.skipif(
    not os.path.exists(os.path.join(ROOT, "starplast", "data", "nodes.parquet")),
    reason="node table not present")
def test_antibody_epitopes_are_not_the_all_types_count():
    """The two fill different slots, so they had better be different measurements.

    Most genes with an IEDB epitope have no ANTIBODY epitope: the ToxoDB column is dominated by
    T-cell assays. If that ever stopped being true, the two slots would be answering one question.
    """
    n = pd.read_parquet(os.path.join(ROOT, "starplast", "data", "nodes.parquet"))
    if I.COLUMN not in n.columns or "iedb_epitope_count" not in n.columns:
        pytest.skip("one of the epitope columns is not merged")
    any_type = n["iedb_epitope_count"].notna()
    antibody = n[I.COLUMN].notna()
    assert antibody.sum() < any_type.sum() / 3, "the antibody set is not the smaller one"
    assert (any_type & ~antibody).sum() > 100, "too few T-cell-only antigens to tell them apart"
