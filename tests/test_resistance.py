#!/usr/bin/env python3
"""The curated resistance table, and the checks that make hand-typed data trustworthy.

This is the only source in the map that is read out of papers rather than out of a file, so it has
no parser to fail and no archive to disagree with. That removes the usual safety net, and these
tests are what replaces it: every allele has to name a gene that still exists, has to describe that
gene the way the current annotation describes it, and has to say how causality was established. The
product check is the important one -- it costs nothing and it catches a mistyped accession, which is
the failure a curated table is most exposed to and the one no amount of re-reading would find.
"""
from __future__ import annotations

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import resistance as R  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NODES = os.path.join(ROOT, "starplast", "data", "nodes.parquet")


# --------------------------------------------------------------------------- the rows themselves
def test_every_allele_is_fully_attributed():
    """A row with no paper behind it is an assertion, not a measurement."""
    for allele in R.ALLELES:
        assert allele.gene.startswith("TGME49_"), allele
        assert allele.pmid.isdigit(), allele
        assert allele.substitution and allele.compound and allele.product, allele
        assert len(allele.validation) > 20, f"{allele.gene} {allele.substitution} says how weakly"


def test_every_allele_records_how_causality_was_shown():
    """The bar is that the mutation was put back. The words that would mean otherwise are absent."""
    for allele in R.ALLELES:
        text = allele.validation.lower()
        assert not any(word in text for word in ("candidate", "associated with", "correlat")), (
            f"{allele.gene} {allele.substitution} reads as a candidate rather than a demonstration")


def test_substitutions_look_like_substitutions():
    import re
    for allele in R.ALLELES:
        assert re.fullmatch(r"[A-Z]\d{1,4}[A-Z]", allele.substitution), allele


# --------------------------------------------------------------------------- the table
def test_one_row_per_gene_with_its_alleles_counted():
    d = R.table()
    assert len(d) == d.index.nunique()
    row = d.loc["TGME49_312570"]
    assert row["resistance_allele_count"] == 3           # L162Q, I171N, S191Y
    assert row["resistance_compound_count"] == 3         # 1NM-PP1, 3BrB-PP1, 3MB-PP1
    assert "L162Q" in row["resistance_substitutions"]


def test_a_gene_with_no_allele_is_absent_rather_than_false():
    """Nobody selected resistance in most genes. Absence is ignorance, not a negative result."""
    d = R.resistance(log=lambda *a: None)
    assert "TGME49_200010" not in d.index
    assert not d.isna().any().any(), "a row that exists must be complete"


def test_accessions_are_resolved_when_a_resolver_is_given():
    d = R.resistance(log=lambda *a: None, resolve=lambda g: g.replace("TGME49_", "TGLATEST_"))
    assert all(g.startswith("TGLATEST_") for g in d.index)


def test_two_accessions_collapsing_to_one_gene_yield_one_row():
    d = R.resistance(log=lambda *a: None, resolve=lambda g: "TGME49_1")
    assert len(d) == 1


def test_the_loader_reports_what_it_carries():
    said = []
    R.resistance(log=said.append)
    assert said and "curated" in said[0]


def test_an_empty_curation_yields_an_empty_frame(monkeypatch):
    monkeypatch.setattr(R, "ALLELES", ())
    assert R.table().empty and R.resistance(log=lambda *a: None).empty


# --------------------------------------------------------------------------- against the map
@pytest.mark.skipif(not os.path.exists(NODES), reason="node table not built")
def test_every_curated_gene_still_exists_and_is_what_the_row_says_it_is():
    """The check that makes a hand-typed accession safe: a wrong one describes a different protein.

    A digit transposed in TGME49_312570 either names nothing or names some other gene, and either
    way the product stops matching. This is the same check that caught the artemisinin paper writing
    its own kinase accession two different ways.
    """
    nodes = pd.read_parquet(NODES, columns=["gene_id", "product"])
    product = dict(zip(nodes["gene_id"], nodes["product"]))
    for allele in R.ALLELES:
        assert allele.gene in product, f"{allele.gene} is not in the current annotation"
        assert product[allele.gene] == allele.product, (
            f"{allele.gene} is '{product[allele.gene]}', not '{allele.product}'")


@pytest.mark.skipif(not os.path.exists(NODES), reason="node table not built")
def test_the_alleles_deliberately_left_out_are_still_left_out():
    """DHFR-TS and DHODH are real and unreachable. If one is added it must bring a source with it."""
    for gene in R.DOCUMENTED_ELSEWHERE:
        assert gene not in {a.gene for a in R.ALLELES}, (
            f"{gene} was added without removing it from DOCUMENTED_ELSEWHERE")


@pytest.mark.skipif(not os.path.exists(NODES), reason="node table not built")
def test_the_absent_alleles_name_genes_that_actually_exist():
    """So the note stays actionable rather than becoming folklore about accessions nobody can find."""
    nodes = pd.read_parquet(NODES, columns=["gene_id"])
    known = set(nodes["gene_id"])
    for gene in R.DOCUMENTED_ELSEWHERE:
        assert gene in known, gene


@pytest.mark.skipif(not os.path.exists(NODES), reason="node table not built")
def test_the_shipped_column_did_not_spread_beyond_the_curated_genes():
    nodes = pd.read_parquet(NODES, columns=["gene_id", "resistance_allele_count"])
    assert nodes["resistance_allele_count"].notna().sum() == len({a.gene for a in R.ALLELES})
