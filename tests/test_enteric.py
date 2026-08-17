#!/usr/bin/env python3
"""The curated enteric-fitness table, and the bar that keeps it honest.

Third curated source in the map. It exists because the sweeps that closed this slot looked for a
POOLED SCREEN through the enteroepithelial stages and correctly found none -- nobody has put a
barcoded library through a cat -- while the slot asks whether disrupting a gene costs the parasite
oocysts, which four labs have answered one knockout at a time.

Same safety net as `resistance` and `drug_sensitivity`: every row records the gene's product as the
CURRENT annotation gives it, so a mistyped accession stops matching and fails a test rather than
quietly describing a different protein.
"""
from __future__ import annotations

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import enteric as E  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NODES = os.path.join(ROOT, "starplast", "data", "nodes.parquet")


def test_every_measurement_is_fully_attributed():
    for m in E.MEASUREMENTS:
        assert m.gene.startswith("TGME49_"), m
        assert m.outcome in E.OUTCOMES, m
        assert m.sporulation in E.SPORULATION, m
        assert m.product and m.strain and m.background, m
        # A PMID where one has been assigned, a PMC id where the paper is too new to have one.
        assert m.pmid.isdigit() or m.pmid.startswith("PMC"), m
        assert len(m.evidence) > 25, f"{m.gene} does not say how it was measured"


def test_every_row_was_actually_carried_through_the_feline_stage():
    """The bar is a MEASURED enteric outcome. Two papers were found and refused for failing exactly
    this: one says the cat experiment "should be carried out", one says oocysts were seen but "not
    quantified". Neither may reappear as a row."""
    for m in E.MEASUREMENTS:
        text = m.evidence.lower()
        assert not any(p in text for p in ("should be carried out", "not quantified",
                                           "were not quantitated", "unclear whether")), m
        assert any(p in text for p in ("oocyst", "sporulation", "ploidy")), m


def test_an_in_vitro_only_phenotype_is_not_admitted():
    """Reduced cyst formation in culture is a different slot. Every row has to name an enteric
    readout, not a tachyzoite or bradyzoite one."""
    for m in E.MEASUREMENTS:
        assert "plaque assay" not in m.evidence.lower() or "oocyst" in m.evidence.lower(), m


def test_tested_and_unchanged_genes_are_kept():
    """A gene deleted with no effect on oocyst yield is a result. Dropping those rows would leave the
    column looking like a list of transmission-blocking hits."""
    unchanged = [m for m in E.MEASUREMENTS if m.outcome == "unchanged"]
    assert len(unchanged) >= 4, "the LEA cluster negatives are missing"
    d = E.table()
    for m in unchanged:
        assert "unchanged" in d.loc[m.gene, "enteric_oocyst_yield"]


def test_the_cluster_deletion_says_it_deleted_the_cluster():
    """All four LEA genes were removed together. That makes the null STRONGER -- redundancy cannot be
    hiding the phenotype -- but the row still has to disclose that no gene was deleted alone."""
    lea = [m for m in E.MEASUREMENTS if m.product.startswith("LEA domain")]
    assert len(lea) == 4
    for m in lea:
        assert "cluster" in m.background.lower() and "four" in m.background.lower(), m


def test_a_genetic_interaction_row_says_so():
    """AAH2 in a delta-aah1 background is not a clean single-gene result, and the double is no worse
    than delta-aah1 alone. Reading a severe defect onto AAH2 from it would be the same error as
    reading a double mutant's phenotype off one of its genes."""
    backgrounds = {m.background for m in E.MEASUREMENTS}
    assert any(b.startswith("d") and "cluster" not in b for b in backgrounds), backgrounds
    interaction = [m for m in E.MEASUREMENTS if m.background == "dAAH1"]
    assert interaction, "the double-mutant row has lost its background"
    for m in interaction:
        assert "indistinguishable" in m.evidence or "no evidence" in m.evidence, m


def test_the_hap2_caveat_is_carried_not_dropped():
    """That line has a second mutation the authors report and argue against. It stays in the record;
    a caveat the source discloses is not ours to silently drop."""
    hap2 = [m for m in E.MEASUREMENTS if "HAP2" in m.product]
    assert hap2
    assert any("second mutation" in m.evidence for m in hap2)


def test_sporulation_and_yield_are_separate_outcomes():
    """HAP2 sheds a few oocysts that never sporulate; Grx5 sheds fewer that sporulate poorly. Those
    are different events and collapsing them into one column would lose the distinction."""
    d = E.table()
    hap2 = d.loc["TGME49_285940"]
    assert "reduced" in hap2["enteric_oocyst_yield"]
    assert "abolished" in hap2["enteric_sporulation"]
    grx5 = d.loc["TGME49_227100"]
    assert "reduced" in grx5["enteric_oocyst_yield"] and "reduced" in grx5["enteric_sporulation"]


def test_unmeasured_sporulation_is_not_reported_as_a_result():
    d = E.table()
    assert d.loc["TGME49_287510", "enteric_sporulation"] == "not measured"


def test_measurement_counts_match_the_rows():
    d = E.table()
    for gene in {m.gene for m in E.MEASUREMENTS}:
        assert d.loc[gene, "enteric_measurements"] == sum(
            1 for m in E.MEASUREMENTS if m.gene == gene)


def test_table_is_one_row_per_gene():
    d = E.table()
    assert d.index.is_unique and d.index.name == "gene_id"
    assert len(d) == len({m.gene for m in E.MEASUREMENTS})


@pytest.mark.skipif(not os.path.exists(NODES), reason="built node table not present")
def test_every_curated_gene_is_what_the_row_says_it_is():
    """The check that earns its keep: the annotation names AAH1 and AAH2 and the four LEA proteins
    outright, so a transposed accession stops matching and this fails."""
    nodes = pd.read_parquet(NODES).set_index("gene_id")
    for m in E.MEASUREMENTS:
        assert m.gene in nodes.index, f"{m.gene} is not in the node table"
        assert str(nodes.loc[m.gene, "product"]).strip() == m.product, (
            f"{m.gene}: annotation says {nodes.loc[m.gene, 'product']!r}, row says {m.product!r}")


@pytest.mark.skipif(not os.path.exists(NODES), reason="built node table not present")
def test_the_curated_columns_reach_the_built_table():
    nodes = pd.read_parquet(NODES)
    for column in ("enteric_oocyst_yield", "enteric_sporulation", "enteric_measurements"):
        assert column in nodes.columns, f"{column} never reached the node table"
    filled = nodes["enteric_oocyst_yield"].notna().sum()
    assert filled == len({m.gene for m in E.MEASUREMENTS}), filled


def test_fitness_resolves_through_the_identity_index():
    out = E.fitness(resolve=lambda g: {"TGME49_227100": "TGME49_999999"}.get(g, g), log=lambda *_: None)
    assert "TGME49_999999" in out.index and "TGME49_227100" not in out.index


def test_fitness_logs_the_negatives_it_kept():
    lines = []
    E.fitness(log=lines.append)
    assert lines and "tested-and-unchanged" in lines[0]


def test_an_empty_curation_returns_an_empty_frame(monkeypatch):
    monkeypatch.setattr(E, "MEASUREMENTS", ())
    assert E.table().empty
    assert E.fitness(log=lambda *_: None).empty
