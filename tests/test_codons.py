#!/usr/bin/env python3
"""Codon usage bias: the arithmetic against hand-checkable sequences, the biology against the map.

ENC and CAI are easy to get subtly wrong -- a stop codon counted as a synonymous choice, a
single-codon amino acid included, a reference set that quietly contains the answer -- and none of
those show up as an error. They show up as a column that looks fine and ranks genes wrongly.
"""
from __future__ import annotations

import gzip
import math
import os
import sys

import pandas as pd
import pytest
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import codons as C  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _cds(tmp_path, rows):
    data = tmp_path / "starplast" / "data"
    data.mkdir(parents=True, exist_ok=True)
    with gzip.open(data / C.TABLE, "wt") as fh:
        fh.write("Gene ID\tCoding Sequence\n")
        for gene, seq in rows:
            fh.write(f"{gene}\t{seq}\n")
    return tmp_path


# --------------------------------------------------------------------------- the genetic code
def test_the_code_has_sixty_four_codons_and_three_stops():
    assert len(C.CODE) == 64
    assert sum(1 for aa in C.CODE.values() if aa == "*") == 3
    assert C.CODE["ATG"] == "M" and C.CODE["TGG"] == "W" and C.CODE["TAA"] == "*"


def test_single_codon_amino_acids_are_not_degenerate():
    """Methionine and tryptophan have no synonymous choice, so they carry no information about bias."""
    assert "M" not in C.DEGENERATE and "W" not in C.DEGENERATE
    assert len(C.DEGENERATE["L"]) == 6 and len(C.DEGENERATE["F"]) == 2


# --------------------------------------------------------------------------- counting
def test_codons_are_read_in_frame_and_junk_is_skipped():
    counts = C.codon_counts("ATGNNNGCA")
    assert counts == {"ATG": 1, "GCA": 1}


def test_a_trailing_partial_codon_is_ignored():
    assert C.codon_counts("ATGGC") == {"ATG": 1}


def test_uracil_is_read_as_thymine():
    assert C.codon_counts("AUG") == {"ATG": 1}


# --------------------------------------------------------------------------- GC3
def test_gc3_counts_only_degenerate_third_positions():
    """TGG is tryptophan -- its G is not a synonymous choice and must not count as GC3."""
    assert C.gc3({"GCG": 1, "GCA": 1}) == pytest.approx(0.5)
    assert C.gc3({"TGG": 10}) != C.gc3({"TGG": 10, "GCG": 1})
    assert math.isnan(C.gc3({"TGG": 10, "ATG": 3}))


def test_stops_do_not_contribute_to_gc3():
    assert math.isnan(C.gc3({"TAA": 5, "TAG": 5}))


# --------------------------------------------------------------------------- ENC
def test_one_codon_per_amino_acid_is_maximum_bias():
    """A gene using a single codon for every family has an effective number near the floor of 20."""
    counts = {codons[0]: 40 for codons in C.DEGENERATE.values()}
    assert C.enc(counts) < 25


def test_even_use_of_every_codon_is_no_bias():
    counts = {codon: 40 for codons in C.DEGENERATE.values() for codon in codons}
    assert C.enc(counts) > 55


def test_enc_is_capped_at_the_number_of_sense_codons():
    """Wright's estimator can exceed 61 on small samples; more codons than exist is not a measurement."""
    counts = {codon: 2 for codons in C.DEGENERATE.values() for codon in codons}
    assert C.enc(counts) <= 61.0


def test_a_family_seen_once_cannot_estimate_homozygosity(tmp_path):
    """One observation gives a homozygosity of 1 and a spurious claim of total bias."""
    assert math.isnan(C.enc({"TTT": 1}))


def test_a_gene_with_no_degenerate_codons_at_all_is_refused():
    """NaN rather than 61: 'we could not tell' is not 'unbiased'."""
    assert math.isnan(C.enc({"ATG": 10, "TGG": 10}))


def test_a_missing_class_is_filled_from_another_rather_than_dropped():
    """A gene using only two-fold families still gets a number -- Wright's own substitution."""
    counts = {codon: 20 for aa, codons in C.DEGENERATE.items() if len(codons) == 2
              for codon in codons}
    assert not math.isnan(C.enc(counts))


# --------------------------------------------------------------------------- CAI
def test_the_reference_set_defines_the_optimal_codon():
    weights = C.reference_weights({"ref": {"TTT": 100, "TTC": 25}})
    assert weights["TTT"] == 1.0 and weights["TTC"] == pytest.approx(0.25)


def test_a_codon_the_reference_never_uses_gets_a_floor_not_a_zero():
    """log(0) would make CAI undefined for any gene that uses it once."""
    weights = C.reference_weights({"ref": {"TTT": 100}})
    assert weights["TTC"] == 0.01


def test_a_gene_matching_the_reference_scores_one():
    weights = C.reference_weights({"ref": {"TTT": 100, "TTC": 25}})
    assert C.cai({"TTT": 30}, weights) == pytest.approx(1.0)


def test_cai_is_a_geometric_and_not_arithmetic_mean():
    """One badly-adapted codon should pull the score down more than an average would allow."""
    weights = {"TTT": 1.0, "TTC": 0.01}
    assert C.cai({"TTT": 1, "TTC": 1}, weights) == pytest.approx(0.1)


def test_a_gene_of_only_non_degenerate_codons_has_no_cai():
    assert math.isnan(C.cai({"ATG": 5}, C.reference_weights({"r": {"TTT": 1}})))


# --------------------------------------------------------------------------- the loader
def test_the_table_yields_a_column_per_statistic(tmp_path):
    _cds(tmp_path, [("TGME49_200010", "ATG" + "GCG" * 40 + "TTT" * 40)])
    out = C.codon_usage(str(tmp_path), reference=["TGME49_200010"], log=lambda *_: None)
    assert list(out.columns) == ["codon_enc", "codon_gc3", "codon_cai_ribosomal"]
    assert out.loc["TGME49_200010", "codon_cai_ribosomal"] == pytest.approx(1.0)


def test_accessions_are_resolved_before_anything_is_computed(tmp_path):
    _cds(tmp_path, [("TGGT1_100010", "GCG" * 40)])
    out = C.codon_usage(str(tmp_path), resolve=lambda g: {"TGGT1_100010": "TGME49_200010"}.get(g),
                        reference=[], log=lambda *_: None)
    assert list(out.index) == ["TGME49_200010"]


def test_an_unresolvable_accession_keeps_its_own_name(tmp_path):
    _cds(tmp_path, [("TGME49_999999", "GCG" * 40)])
    out = C.codon_usage(str(tmp_path), resolve=lambda g: None, reference=[], log=lambda *_: None)
    assert list(out.index) == ["TGME49_999999"]


def test_a_row_with_no_sequence_is_skipped(tmp_path):
    _cds(tmp_path, [("TGME49_200010", ""), ("TGME49_200020", "GCG" * 40)])
    out = C.codon_usage(str(tmp_path), reference=[], log=lambda *_: None)
    assert list(out.index) == ["TGME49_200020"]


def test_a_row_of_only_junk_is_skipped(tmp_path):
    _cds(tmp_path, [("TGME49_200010", "NNNNNN"), ("TGME49_200020", "GCG" * 40)])
    out = C.codon_usage(str(tmp_path), reference=[], log=lambda *_: None)
    assert list(out.index) == ["TGME49_200020"]


def test_a_file_of_nothing_usable_yields_nothing(tmp_path):
    _cds(tmp_path, [("TGME49_200010", "NNN")])
    assert C.codon_usage(str(tmp_path), reference=[], log=lambda *_: None).empty


def test_no_cds_table_yields_nothing(tmp_path):
    assert C.codon_usage(str(tmp_path), log=lambda *_: None).empty


def test_without_a_reference_set_the_cai_column_is_left_out(tmp_path):
    """An arbitrary reference set would produce a number that means nothing."""
    _cds(tmp_path, [("TGME49_200010", "GCG" * 40)])
    said = []
    out = C.codon_usage(str(tmp_path), reference=[], log=said.append)
    assert "codon_cai_ribosomal" not in out
    assert any("no CAI reference" in m for m in said)


def test_the_reference_set_falls_back_to_the_node_table(tmp_path):
    """reference=None means 'the ribosomal proteins', read from the shipped table."""
    _cds(tmp_path, [("TGME49_200010", "GCG" * 40)])
    (tmp_path / "starplast" / "data" / "nodes.parquet").parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"gene_id": ["TGME49_200010"],
                  "product": ["40S ribosomal protein S3"]}).to_parquet(
        tmp_path / "starplast" / "data" / "nodes.parquet")
    out = C.codon_usage(str(tmp_path), log=lambda *_: None)
    assert "codon_cai_ribosomal" in out


def test_no_node_table_means_no_reference_set(tmp_path):
    _cds(tmp_path, [("TGME49_200010", "GCG" * 40)])
    assert "codon_cai_ribosomal" not in C.codon_usage(str(tmp_path), log=lambda *_: None)


# --------------------------------------------------------------------------- against the real map
@pytest.mark.skipif(not os.path.exists(os.path.join(ROOT, "starplast", "data", C.TABLE)),
                    reason="shipped CDS table not present")
def test_bias_tracks_expression_the_way_translational_selection_says_it_must():
    """The check no fixture can give, and the one that would catch a wrong reference set.

    Heavily translated genes use a narrower set of codons. So CAI must rise with expression and ENC
    must fall, and ribosomal proteins must be more biased than the genome. CAI's reference set is
    chosen by annotation rather than by expression, which is what makes this a finding instead of a
    restatement -- and what makes it able to fail.
    """
    X = C.codon_usage(ROOT, log=lambda *_: None)
    n = pd.read_parquet(os.path.join(ROOT, "starplast", "data", "nodes.parquet")).set_index("gene_id")
    ribo = [g for g in n.index[n["product"].astype(str).str.contains(
        "ribosomal protein", case=False, na=False)] if g in X.index]
    rest = [g for g in X.index if g not in set(ribo)]
    assert X.loc[ribo, "codon_enc"].median() < X.loc[rest, "codon_enc"].median() - 3

    j = X.join(n["expr_tachy"].dropna(), how="inner")
    assert spearmanr(j["codon_cai_ribosomal"], j["expr_tachy"]).statistic > 0.2
    assert spearmanr(j["codon_enc"], j["expr_tachy"]).statistic < -0.1
