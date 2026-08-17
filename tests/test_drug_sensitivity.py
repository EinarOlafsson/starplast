#!/usr/bin/env python3
"""The curated drug-sensitivity table, and the bar that keeps it honest.

Second curated source in the map, and it exists because six catalogue sweeps looked for a genome-wide
chemogenomic screen and correctly found none -- while the sentence they were testing, "no screen
exists", is not the same statement as "the question cannot be answered". A knockout with a measured
shift under a named compound answers it one gene at a time.

Same safety net as `resistance`: every row records the gene's product as the CURRENT annotation gives
it, so a mistyped accession stops matching and fails a test rather than describing a different
protein.
"""
from __future__ import annotations

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import drug_sensitivity as D  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NODES = os.path.join(ROOT, "starplast", "data", "nodes.parquet")


def test_every_effect_is_fully_attributed():
    for effect in D.EFFECTS:
        assert effect.gene.startswith("TGME49_"), effect
        assert effect.pmid.isdigit() and effect.compound and effect.product, effect
        assert effect.direction in D.DIRECTIONS, effect
        assert len(effect.evidence) > 20, f"{effect.gene} {effect.compound} says how weakly"


def test_target_engagement_is_not_admitted_as_drug_sensitivity():
    """"An inhibitor of this protein kills the parasite" is a different slot. Every row here has to be
    about a DISRUPTED GENE changing how the parasite responds to a compound."""
    for effect in D.EFFECTS:
        text = effect.evidence.lower()
        assert not any(word in text for word in ("inhibitor of", "inhibits the target",
                                                 "target engagement")), effect


def test_a_tested_and_unchanged_gene_is_kept():
    """A transporter deleted with no effect on analog sensitivity is a result. Dropping those rows
    would leave the column looking like a list of hits."""
    unchanged = [e for e in D.EFFECTS if e.direction == "unchanged"]
    assert unchanged, "no negative results are recorded"
    d = D.table()
    gene = unchanged[0].gene
    assert d.loc[gene, "drug_compounds_tested"] > 0
    assert d.loc[gene, "drug_sensitivity_shifts"] == 0


def test_a_double_knockout_effect_says_so():
    """TgENT3 was deleted in a ΔTgAT1 background, which makes it a genetic-interaction result.
    Reading it as a single-gene effect is the same error as reading a double mutant off one gene."""
    backgrounds = {e.background for e in D.EFFECTS}
    assert backgrounds != {"parental"}, "no row records a non-parental background"
    for effect in D.EFFECTS:
        assert effect.background, effect


def test_the_table_counts_compounds_and_shifts_separately():
    d = D.table()
    assert d.loc["TGME49_244440", "drug_sensitivity_shifts"] == 2
    assert "Ara-A:resistant" in d.loc["TGME49_244440", "drug_sensitivity_directions"]


def test_accessions_are_resolved_when_a_resolver_is_given():
    d = D.sensitivity(log=lambda *a: None, resolve=lambda g: g.replace("TGME49_", "TGLATEST_"))
    assert all(g.startswith("TGLATEST_") for g in d.index)


def test_two_accessions_collapsing_to_one_gene_yield_one_row():
    assert len(D.sensitivity(log=lambda *a: None, resolve=lambda g: "TGME49_1")) == 1


def test_the_loader_reports_what_it_carries():
    said = []
    D.sensitivity(log=said.append)
    assert said and "curated" in said[0]


def test_an_empty_curation_yields_an_empty_frame(monkeypatch):
    monkeypatch.setattr(D, "EFFECTS", ())
    assert D.table().empty and D.sensitivity(log=lambda *a: None).empty


@pytest.mark.skipif(not os.path.exists(NODES), reason="node table not built")
def test_every_curated_gene_is_what_the_row_says_it_is():
    """The check that makes a hand-typed accession safe, and it earned its keep here: `adenosine
    transporter AT1` in the annotation is independent confirmation that TGME49_244440 is TgAT1."""
    nodes = pd.read_parquet(NODES, columns=["gene_id", "product"])
    product = dict(zip(nodes["gene_id"], nodes["product"]))
    for effect in D.EFFECTS:
        assert effect.gene in product, f"{effect.gene} is not in the current annotation"
        assert product[effect.gene] == effect.product, (
            f"{effect.gene} is '{product[effect.gene]}', not '{effect.product}'")


@pytest.mark.skipif(not os.path.exists(NODES), reason="node table not built")
def test_the_shipped_column_did_not_spread_beyond_the_curated_genes():
    nodes = pd.read_parquet(NODES, columns=["drug_compounds_tested"])
    assert nodes["drug_compounds_tested"].notna().sum() == len({e.gene for e in D.EFFECTS})
