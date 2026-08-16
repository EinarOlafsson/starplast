#!/usr/bin/env python3
"""The palmitome, and the sign convention that would have silently halved it.

ToxoDB reports a signed fold DIFFERENCE, not a ratio: -3.12 means three-fold down. Passing that to
log2 gives NaN, which would have dropped every depleted protein while the column still looked like a
working measurement. That is the whole reason `signed_log2` exists and is tested first.
"""
from __future__ import annotations

import math
import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import palmitome as PM  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _tables(tmp_path, rows, filename="toxodb_palmitome_hydroxylamine.tsv"):
    data = tmp_path / "starplast" / "data"
    data.mkdir(parents=True, exist_ok=True)
    body = "\n".join(f"{g}\t{v}" for g, v in rows)
    (data / filename).write_text("Gene ID\tFold Difference\n" + body + "\n")
    return tmp_path


# --------------------------------------------------------------------------- the sign convention
def test_a_positive_fold_is_its_logarithm():
    assert PM.signed_log2(4) == pytest.approx(2.0)


def test_a_negative_fold_is_a_depletion_and_not_a_negative_quantity():
    """-4 means four-fold DOWN. log2(-4) is NaN, and that NaN would be half the table."""
    assert PM.signed_log2(-4) == pytest.approx(-2.0)


def test_no_change_is_zero():
    assert PM.signed_log2(1) == 0.0


def test_zero_and_nonsense_are_missing():
    assert math.isnan(PM.signed_log2(0))
    assert math.isnan(PM.signed_log2("n/a"))
    assert math.isnan(PM.signed_log2(None))


# --------------------------------------------------------------------------- the loader
def test_each_comparison_becomes_its_own_column(tmp_path):
    """They answer different questions -- one is thioester-specific and one is not -- so they are
    never averaged into a single 'palmitoylation' number."""
    _tables(tmp_path, [("TGGT1_200260", 4.0)])
    _tables(tmp_path, [("TGGT1_200260", 2.0)], filename="toxodb_palmitome_palmitate.tsv")
    out = PM.palmitome(str(tmp_path), log=lambda *_: None)
    assert sorted(out.columns) == sorted(PM.COMPARISONS.values())
    assert out.loc["TGGT1_200260", "palmitome_odya_vs_hydroxylamine_log2"] == pytest.approx(2.0)
    assert out.loc["TGGT1_200260", "palmitome_odya_vs_palmitate_log2"] == pytest.approx(1.0)


def test_accessions_go_through_the_identity_layer(tmp_path):
    """The whole table is TGGT1_; matched as strings it would join nothing."""
    _tables(tmp_path, [("TGGT1_200260", 4.0)])
    out = PM.palmitome(str(tmp_path),
                       resolve=lambda g: {"TGGT1_200260": "TGME49_200010"}.get(g),
                       log=lambda *_: None)
    assert list(out.index) == ["TGME49_200010"]


def test_an_unresolvable_accession_keeps_its_own_name(tmp_path):
    _tables(tmp_path, [("TGGT1_999999", 4.0)])
    out = PM.palmitome(str(tmp_path), resolve=lambda g: None, log=lambda *_: None)
    assert list(out.index) == ["TGGT1_999999"]


def test_two_rows_for_one_gene_keep_the_stronger_enrichment(tmp_path):
    """Averaging a hit with a non-hit erases the hit, which is the claim the study makes."""
    _tables(tmp_path, [("TGGT1_200260", 8.0), ("TGGT1_200260", 1.0)])
    out = PM.palmitome(str(tmp_path), log=lambda *_: None)
    assert out.loc["TGGT1_200260", "palmitome_odya_vs_hydroxylamine_log2"] == pytest.approx(3.0)


def test_genes_in_only_one_comparison_are_kept(tmp_path):
    _tables(tmp_path, [("A", 4.0)])
    _tables(tmp_path, [("B", 4.0)], filename="toxodb_palmitome_palmitate.tsv")
    out = PM.palmitome(str(tmp_path), log=lambda *_: None)
    assert set(out.index) == {"A", "B"}
    assert pd.isna(out.loc["B", "palmitome_odya_vs_hydroxylamine_log2"])


def test_a_missing_table_is_skipped(tmp_path):
    _tables(tmp_path, [("A", 4.0)])
    out = PM.palmitome(str(tmp_path), log=lambda *_: None)
    assert list(out.columns) == ["palmitome_odya_vs_hydroxylamine_log2"]


def test_a_one_column_table_is_refused(tmp_path):
    """The search can be asked for gene ids alone; that report has no measurement in it."""
    data = tmp_path / "starplast" / "data"
    data.mkdir(parents=True)
    (data / "toxodb_palmitome_hydroxylamine.tsv").write_text("Gene ID\nTGGT1_200260\n")
    assert PM.palmitome(str(tmp_path), log=lambda *_: None).empty


def test_nothing_downloaded_yields_nothing(tmp_path):
    assert PM.palmitome(str(tmp_path), log=lambda *_: None).empty


# --------------------------------------------------------------------------- against the real map
@pytest.mark.skipif(
    not os.path.exists(os.path.join(ROOT, "starplast", "data",
                                    "toxodb_palmitome_hydroxylamine.tsv")),
    reason="palmitome table not shipped")
def test_the_known_substrates_come_out_enriched():
    """The check a metadata field cannot give.

    GAP45 anchors the glideosome to the inner membrane complex through its palmitoyl group, MLC1's
    palmitoylation is the subject of a later paper by the same group, and ROP5 and AMA1 are both
    established. All four must sit above the median of the measured genes, or the column is not
    measuring palmitoylation.
    """
    n = pd.read_parquet(os.path.join(ROOT, "starplast", "data", "nodes.parquet")).set_index("gene_id")
    column = "palmitome_odya_vs_hydroxylamine_log2"
    if column not in n.columns:
        pytest.skip("palmitome column not merged")
    measured = n[column].dropna()
    prod = n["product"].astype(str)
    for name, pattern in (("GAP45", r"GAP45"), ("MLC1", r"myosin light chain MLC1\b"),
                          ("AMA1", r"apical membrane antigen"), ("ROP5", r"ROP5")):
        ids = [g for g in n.index[prod.str.contains(pattern, case=False, regex=True, na=False)]
               if g in measured.index]
        assert ids, f"{name} not measured in the palmitome"
        assert measured.loc[ids].max() > measured.median(), name
