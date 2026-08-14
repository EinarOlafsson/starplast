#!/usr/bin/env python3
"""The two claims a map can make, and the guards that stop them being made badly.

Built on a synthetic table where the right answer is known by construction: one cluster of genes
that are all IMC with some unlabelled members hidden among them, one cluster that is one compartment
but two fitness regimes, and a lot of noise. A test that used the real map would be measuring the
map rather than the statistics.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import discovery as D  # noqa: E402


def table(n_noise=300, seed=0):
    """A table with two clusters worth finding and one that should be refused."""
    rng = np.random.default_rng(seed)
    rows, labels = [], []

    # Cluster 0: 40 genes, 30 labelled IMC, 10 unlabelled. The guilt-by-association case.
    for i in range(40):
        known = i < 30
        rows.append({"gene_id": f"IMC{i}", "compartment": "IMC" if known else None,
                     "fitness": float(rng.normal(-1, 0.3)), "n_publications": 0})
        labels.append(0)

    # Cluster 1: 40 rhoptry genes whose fitness falls into two clean groups. The disagreement case.
    for i in range(40):
        rows.append({"gene_id": f"ROP{i}", "compartment": "rhoptry",
                     "fitness": float(rng.normal(-3 if i < 18 else 0.0, 0.25)),
                     "n_publications": 0})
        labels.append(1)

    # Everything else: a mixed bag, at the background rate, and no cluster of its own.
    kinds = ["cytosol", "ER", "nucleus", "mitochondrion", "IMC"]
    for i in range(n_noise):
        rows.append({"gene_id": f"BG{i}", "compartment": kinds[i % len(kinds)],
                     "fitness": float(rng.normal(0, 1)), "n_publications": int(rng.integers(0, 3))})
        labels.append(-1 if i % 3 else 2 + (i // 60))
    return pd.DataFrame(rows), np.array(labels)


# --------------------------------------------------------------------------- guilt, discrete
def test_a_pure_cluster_predicts_its_unlabelled_members():
    nodes, labels = table()
    found = D.guilt(nodes, labels, "compartment")
    imc = found[(found.cluster == 0) & (found.category == "IMC")]
    assert len(imc) == 1, found[["cluster", "category", "purity", "lift", "q"]].to_string()
    row = imc.iloc[0]
    assert row.n_hits == 30 and row.purity == pytest.approx(1.0)
    assert row.n_predicted == 10, "the unlabelled members are the prediction and were miscounted"
    assert row.q < 1e-6 and row.lift > 4
    assert set(row.genes) == {f"IMC{i}" for i in range(30, 40)}


def test_a_cluster_with_nothing_new_in_it_is_not_a_finding():
    """A perfect enrichment over genes that are all already labelled has confirmed something. The
    module exists to produce work, and there is none here."""
    nodes, labels = table()
    nodes.loc[nodes.gene_id.str.startswith("IMC"), "compartment"] = "IMC"
    found = D.guilt(nodes, labels, "compartment")
    assert found[(found.cluster == 0)].empty


def test_a_slice_of_the_proteome_is_refused_however_small_its_p_value():
    """The finding that made these guards necessary: 395 of 1621 genes cytosol against a background
    of 15.9% is a lift of 1.5, a p of 1e-32, and a prediction about nothing. Significance is a
    statement about sample size as much as effect."""
    rng = np.random.default_rng(3)
    n, half = 4000, 2000
    where = np.concatenate([rng.choice(["cytosol", "other", None], half, p=[0.24, 0.56, 0.20]),
                            rng.choice(["cytosol", "other", None], half, p=[0.08, 0.72, 0.20])])
    nodes = pd.DataFrame({"gene_id": [f"G{i}" for i in range(n)], "compartment": where,
                          "n_publications": 0})
    labels = np.array([0] * half + [1] * half)
    found = D.guilt(nodes, labels, "compartment")
    assert found.empty or (found.lift >= D.MIN_LIFT).all()
    assert found.empty or (found.n_cluster <= D.MAX_SHARE * len(labels)).all()


def test_absence_is_not_a_category_to_be_enriched_for():
    """A cluster that is 90% "unassigned" in a proteome that is 60% unassigned has found nothing,
    and would otherwise be the strongest finding in every run."""
    nodes, labels = table()
    nodes.loc[nodes.gene_id.str.startswith("IMC"), "compartment"] = "unassigned"
    found = D.guilt(nodes, labels, "compartment")
    assert "unassigned" not in set(found.category)


# --------------------------------------------------------------------------- guilt, continuous
def test_a_cluster_of_essential_genes_predicts_the_ones_nobody_screened():
    nodes, labels = table()
    nodes.loc[nodes.gene_id.isin([f"IMC{i}" for i in range(35, 40)]), "fitness"] = np.nan
    found = D.guilt(nodes, labels, "fitness")
    row = found[found.cluster == 0].iloc[0]
    assert row.layer_kind == D.CONTINUOUS
    assert row.effect < -D.MIN_EFFECT, "a cluster at -1 against a background of 0 read as no effect"
    assert row.n_predicted == 5, "the unscreened members are the prediction and were miscounted"
    assert row.q < 0.01


def test_a_quantity_that_barely_differs_is_not_a_finding():
    nodes, labels = table()
    nodes["fitness"] = np.random.default_rng(5).normal(0, 1, len(nodes))
    nodes.loc[:10, "fitness"] = np.nan
    found = D.guilt(nodes, labels, "fitness")
    assert found.empty or (found.effect.abs() >= D.MIN_EFFECT).all()


def test_a_layer_that_is_not_in_the_table_finds_nothing_rather_than_raising():
    nodes, labels = table()
    assert D.guilt(nodes, labels, "phlogiston").empty
    assert D.guilt(nodes, labels[:5], "compartment").empty
    assert D.disagreement(nodes, labels, "compartment", "phlogiston").empty


# --------------------------------------------------------------------------- disagreement
def test_one_compartment_and_two_fitness_regimes_is_a_subdivision():
    """The claim nothing about either layer alone can make: these genes agree about where they are
    and disagree about whether the parasite can lose them."""
    nodes, labels = table()
    found = D.disagreement(nodes, labels, "compartment", "fitness")
    row = found[found.cluster == 1].iloc[0]
    assert row.category == "rhoptry"
    assert row.mean_low < -2 and row.mean_high > -1, (row.mean_low, row.mean_high)
    assert row.n_low + row.n_high == 40
    assert row.gap > 2, "two groups three standard deviations apart read as one"


def test_a_cluster_that_agrees_about_both_layers_is_not_a_disagreement():
    """And "agrees" is measured against the spread of the layer across the map, not within the
    cluster. Genes at -3.0 give or take 0.05 split into -2.97 and -3.03 three WITHIN-cluster
    standard deviations apart, and a set that agrees about everything gets reported as a
    subdivision -- which is what this did before the statistic was rescaled."""
    nodes, labels = table()
    nodes.loc[nodes.gene_id.str.startswith("ROP"), "fitness"] = -3.0 + np.random.default_rng(
        1).normal(0, 0.05, 40)
    found = D.disagreement(nodes, labels, "compartment", "fitness")
    assert found.empty or found[found.cluster == 1].empty


def test_a_cluster_that_agrees_about_nothing_is_not_a_disagreement_either():
    """Disagreement needs something to disagree WITH. A cluster of mixed compartments splitting on
    fitness is two facts about a bad cluster, not a subdivision of a category."""
    nodes, labels = table()
    nodes.loc[nodes.gene_id.str.startswith("ROP"), "compartment"] = \
        ["ER", "nucleus", "rhoptry", "cytosol"] * 10
    found = D.disagreement(nodes, labels, "compartment", "fitness")
    assert found.empty or found[found.cluster == 1].empty


def test_two_categorical_layers_disagree_by_splitting():
    nodes, labels = table()
    nodes["phase"] = ["G1"] * len(nodes)
    nodes.loc[nodes.gene_id.isin([f"ROP{i}" for i in range(20)]), "phase"] = "S/M"
    found = D.disagreement(nodes, labels, "compartment", "phase")
    row = found[found.cluster == 1].iloc[0]
    assert row.other_kind == D.DISCRETE
    assert set(row.groups) == {"G1", "S/M"}
    assert sorted(row.counts) == [20, 20]


# --------------------------------------------------------------------------- guards and scoring
def test_a_finding_about_a_layer_that_built_the_map_is_marked_circular():
    """A cluster enriched for localisation in a map built from localisation is arithmetic."""
    nodes, labels = table()
    clean = D.guilt(nodes, labels, "compartment")
    dirty = D.guilt(nodes, labels, "compartment", used_columns={"compartment"})
    assert not clean.circular.any() and dirty.circular.all()
    assert D.yield_score(dirty) == 0.0, "a circular finding was rewarded"
    assert D.yield_score(clean) > 0


def test_a_layer_measured_twice_is_circular_through_either_name():
    """hyperLOPIT's probabilities built the map and `compartment` is what hyperLOPIT calls its
    answer. An exact column match would call that pair independent."""
    from starplast.search import SAME_QUANTITY
    family = next(iter(SAME_QUANTITY.values()))
    assert D._touches(list(family)[0], set(family))


def test_the_score_prefers_many_solid_claims_to_one_overwhelming_one():
    """A configuration that produces one enormous claim about four genes should not beat one that
    produces six solid claims about ninety."""
    one = pd.DataFrame([{"lift": 40.0, "n_predicted": 4, "q": 1e-40, "circular": False}])
    many = pd.DataFrame([{"lift": 4.0, "n_predicted": 15, "q": 1e-4, "circular": False}] * 6)
    assert D.yield_score(many) > D.yield_score(one)


def test_findings_that_fail_correction_count_for_nothing():
    weak = pd.DataFrame([{"lift": 9.0, "n_predicted": 30, "q": 0.4, "circular": False}])
    assert D.yield_score(weak) == 0.0
    assert D.yield_score(pd.DataFrame()) == 0.0
    assert D.yield_score(None) == 0.0


def test_a_quantity_with_no_variation_is_not_an_effect():
    flat = np.zeros(20)
    assert D._effect(flat, flat) == (0.0, 1.0)
    assert D._effect(flat[:2], flat) == (0.0, 1.0)


def test_a_handful_of_values_cannot_be_split_in_two():
    assert D._split(np.array([1.0, 2.0, 3.0]), 1.0)[0] == 0.0
    assert D._split(np.zeros(20), 1.0)[0] == 0.0
    assert D._split(np.array([]), 1.0)[0] == 0.0
    # A scale of zero or nothing at all falls back to 1 rather than dividing by it.
    assert D._split(np.array([0.0, 0, 0, 1, 1, 1, 1]), 0.0)[0] > 0
    assert D._split(np.array([0.0, 0, 0, 1, 1, 1, 1]), np.nan)[0] > 0


# --------------------------------------------------------------------------- the remaining edges
def test_a_category_nobody_has_cannot_be_enriched():
    """Zero successes in the population. scipy would return a number; there is no test to run."""
    assert D._hypergeom(3, 10, 0, 100) == 1.0
    assert D._hypergeom(0, 0, 5, 100) == 1.0
    assert D._hypergeom(3, 10, 5, 0) == 1.0


def test_two_identical_columns_of_numbers_have_no_difference_to_test():
    """Mann-Whitney raises rather than returning 1.0 when every value on both sides is the same,
    which happens whenever a measurement is a constant -- a column of zeros, a flag nobody set."""
    same = np.zeros(10)
    assert D._effect(same, same.copy()) == (0.0, 1.0)


def test_a_split_needs_room_to_cut():
    """Six values is the smallest set this can cut three from each side of, and the guard has to
    hold at exactly that boundary rather than one past it."""
    assert D._split(np.arange(6.0), 1.0)[0] > 0
    assert D._split(np.arange(5.0), 1.0)[0] == 0.0


def test_each_guard_refuses_on_its_own():
    """One test per guard, each with the others satisfied, so a failure names which one broke."""
    nodes, labels = table()
    # Too big a share of the map: the two-clusters-over-everything failure.
    assert D.guilt(nodes, labels, "compartment", max_share=0.01).empty
    # Not enough lift over the background rate.
    assert D.guilt(nodes, labels, "compartment", min_lift=99.0).empty
    # Too few genes carrying the measurement inside the cluster, for a quantity.
    assert D.guilt(nodes, labels, "fitness", min_category=999).empty
    # And a quantity in a cluster that covers too much of the map, which is a separate guard from
    # the categorical one above and was reachable only through this path.
    assert D.guilt(nodes, labels, "fitness", max_share=0.01).empty
    # Too small a cluster to say anything about.
    assert D.guilt(nodes, labels, "compartment", min_cluster=10_000).empty
    # And for a disagreement: not enough measured genes on both sides of a possible split.
    assert D.disagreement(nodes, labels, "compartment", "fitness", min_category=999).empty
    assert D.disagreement(nodes, labels, "compartment", "fitness", min_cluster=10_000).empty


def test_a_cluster_that_is_too_small_or_too_pure_to_predict_is_skipped():
    nodes, labels = table()
    tiny = np.array(labels)
    tiny[:40] = 99                                   # a cluster of 40, all labelled but one
    nodes.loc[nodes.gene_id.str.startswith("IMC"), "compartment"] = "IMC"
    nodes.loc[nodes.gene_id == "IMC0", "compartment"] = None
    found = D.guilt(nodes, tiny, "compartment", min_predictions=3)
    assert found.empty or (found.cluster != 99).all()


def test_a_category_below_the_minimum_count_is_not_a_claim():
    """Three genes agreeing is a coincidence with a p-value attached."""
    nodes, labels = table()
    nodes.loc[nodes.gene_id.isin([f"IMC{i}" for i in range(3)]), "compartment"] = "rare"
    found = D.guilt(nodes, labels, "compartment")
    assert "rare" not in set(found.category)


def test_a_layer_is_circular_through_a_column_that_merely_starts_with_the_family_name():
    """`lopit_prob_er` built the map and `compartment` is what that experiment calls its answer.
    Matching exact column names would call the pair independent."""
    from starplast.search import SAME_QUANTITY
    family = sorted(next(iter(SAME_QUANTITY.values())))
    assert D._touches(family[0], {f"{family[0]}_probability_er"})
    assert not D._touches("something_else_entirely", {"unrelated_column"})
