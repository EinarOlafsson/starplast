#!/usr/bin/env python3
"""Do the scores actually order maps by how good they are?

Every test here builds clusterings whose quality is known BY CONSTRUCTION -- perfect, shuffled,
merged into one, shattered into singletons, corrupted by a controlled fraction -- and asserts that
the numbers put them in the right order. That is the only property of a score that matters for an
optimiser: it is going to climb whatever this returns, so if a worse clustering can score higher,
the optimiser will find that clustering and report it as the answer.

Testing the ordering rather than the values also survives the numbers being tuned, which they will
be. "Perfect beats corrupted beats shuffled" is a claim about the metric; "AUPRC is 0.87" is a claim
about this fixture.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import metrics as M  # noqa: E402


def world(n_per=60, classes=("IMC", "rhoptry", "apicoplast", "nucleus"), spread=0.6, seed=0):
    """Well-separated blobs, one per class, with the truth known exactly."""
    rng = np.random.default_rng(seed)
    centres = np.eye(len(classes)) * 10.0
    coords, truth = [], []
    for i, c in enumerate(classes):
        coords.append(rng.normal(centres[i], spread, size=(n_per, len(classes))))
        truth += [c] * n_per
    return np.vstack(coords), pd.Series(truth)


def corrupt(labels, fraction, rng):
    """Move a share of the genes to a random other cluster."""
    out = np.array(labels)
    move = rng.random(len(out)) < fraction
    out[move] = rng.integers(0, out.max() + 1, int(move.sum()))
    return out


# --------------------------------------------------------------------------- the ordering property
def test_a_perfect_clustering_scores_at_the_top_and_a_shuffled_one_at_chance():
    coords, truth = world()
    perfect = pd.factorize(truth)[0]
    shuffled = np.random.default_rng(1).permutation(perfect)
    good = M.summarise(M.ranking(perfect, truth))
    bad = M.summarise(M.ranking(shuffled, truth))
    assert good["mean_auprc"] > 0.95 and good["mean_auroc"] > 0.95
    # Chance for AUPRC is the class prevalence, which is where lift = 1 -- not zero, and reading it
    # as though it were is how a useless ranking gets reported as "AUPRC 0.25".
    assert bad["mean_lift"] == pytest.approx(1.0, abs=0.35)
    assert good["mean_auprc"] > bad["mean_auprc"]
    # And the AUROC floor is NOT 0.5: leave-one-out drops every positive just below every negative
    # inside an uninformative cluster, so a shuffled clustering lands near 0.36 here. Measured
    # against its own permutation null rather than against a number from a textbook.
    floor = M.chance_level(shuffled, truth, "IMC", n_permutations=8)
    assert bad["mean_auroc"] < 0.55
    assert bad["mean_auroc"] == pytest.approx(floor["auroc"], abs=0.12)
    assert floor["lift"] == pytest.approx(1.0, abs=0.3)


def test_the_score_falls_as_the_clustering_is_corrupted():
    """Monotone in the thing it claims to measure. A score that is high for a perfect clustering and
    high again for a badly corrupted one cannot be climbed toward the good one."""
    coords, truth = world()
    perfect = pd.factorize(truth)[0]
    rng = np.random.default_rng(2)
    scores = [M.summarise(M.ranking(corrupt(perfect, f, rng), truth))["mean_auprc"]
              for f in (0.0, 0.15, 0.35, 0.6, 0.9)]
    assert all(a >= b - 0.02 for a, b in zip(scores, scores[1:])), scores
    assert scores[0] - scores[-1] > 0.4, "corruption barely moved the score"


def test_one_giant_cluster_scores_at_chance_however_pure_the_data_is():
    """The failure this project actually hit: two clusters over eight thousand genes, reported as
    enrichment. A ranking metric has to be immune to it, because the partition metrics are not."""
    coords, truth = world()
    everything = np.zeros(len(truth), dtype=int)
    summary = M.summarise(M.ranking(everything, truth))
    assert summary["mean_lift"] == pytest.approx(1.0, abs=0.05)
    assert M.partition(everything, truth)["largest_share"] == 1.0


def test_singletons_do_not_score_themselves():
    """Leave-one-out, and the whole reason the numbers are honest. Every gene in its own cluster is
    a perfect partition by every purity measure ever written, and predicts nothing at all."""
    coords, truth = world()
    alone = np.arange(len(truth))
    summary = M.summarise(M.ranking(alone, truth))
    assert summary["mean_auprc"] < 0.3, "a clustering of singletons scored as though it predicted"


def test_splitting_a_class_in_two_keeps_the_ranking_and_costs_the_partition():
    """The diagnostic difference between the two families of metric. Over-splitting is nearly free
    for annotation -- the neighbours still carry the right label -- and expensive for any measure of
    agreement between partitions, which is why tuning on ARI alone drives toward merged clusters."""
    coords, truth = world()
    perfect = pd.factorize(truth)[0]
    split = perfect * 2 + (np.arange(len(perfect)) % 2)
    assert M.summarise(M.ranking(split, truth))["mean_auprc"] > 0.9
    assert M.partition(split, truth)["ari"] < M.partition(perfect, truth)["ari"]


def test_noise_is_counted_against_the_map_rather_than_excused():
    """A gene called noise has not been ranked, and dropping it from the denominator hides exactly
    the failure the reader needs to see."""
    coords, truth = world()
    perfect = pd.factorize(truth)[0]
    half_noise = np.where(np.arange(len(perfect)) % 2, perfect, -1)
    assert M.summarise(M.ranking(half_noise, truth))["mean_auprc"] < \
        M.summarise(M.ranking(perfect, truth))["mean_auprc"]
    assert M.partition(half_noise, truth)["noise"] == pytest.approx(0.5, abs=0.01)


# --------------------------------------------------------------------------- edges and absences
def test_a_class_with_no_negatives_has_no_curve_rather_than_a_perfect_one():
    hit = np.ones(20, dtype=bool)
    out = M.curve_scores(np.random.default_rng(0).random(20), hit)
    assert np.isnan(out["auroc"]) and np.isnan(out["auprc"])
    assert np.isnan(M.curve_scores(np.zeros(2), np.array([True, False]))["auroc"])


def test_absent_labels_are_not_a_class_to_score():
    coords, truth = world()
    truth = truth.copy()
    truth[:100] = "unassigned"
    assert "unassigned" not in set(M.ranking(pd.factorize(truth)[0], truth).category)


def test_an_empty_table_summarises_to_zero_rather_than_raising():
    assert M.summarise(pd.DataFrame())["mean_auprc"] == 0.0
    assert M.summarise(None)["n_categories"] == 0


def test_a_map_with_nothing_labelled_scores_nothing():
    coords, _ = world()
    blank = pd.Series([None] * len(coords))
    assert M.knn_purity(coords, blank) == 0.0
    assert M.expected_purity(blank) == 0.0


# --------------------------------------------------------------------------- the map, without a clustering
def test_neighbourhood_purity_beats_chance_on_a_map_that_has_structure():
    coords, truth = world()
    assert M.knn_purity(coords, truth) > 0.9
    assert M.expected_purity(truth) == pytest.approx(0.25, abs=0.01)
    scrambled = np.random.default_rng(3).permutation(coords)
    assert M.knn_purity(scrambled, truth) < M.knn_purity(coords, truth)


def test_the_neighbourhood_ranking_is_reported_beside_the_cluster_one():
    """The gap between them is the diagnostic: structure in the map that the clustering lost."""
    coords, truth = world()
    merged = np.zeros(len(truth), dtype=int)
    table = M.ranking(merged, truth, coords)
    assert (table.nn_auprc > table.auprc).all(), "a lost structure was not visible in the map"


def test_trustworthiness_prefers_a_projection_that_kept_the_neighbourhoods():
    coords, truth = world()
    kept = coords[:, :2] * 1.0
    scrambled = np.random.default_rng(4).normal(size=(len(coords), 2))
    assert M.trustworthiness(coords, kept, k=10) > M.trustworthiness(coords, scrambled, k=10)
    assert np.isnan(M.trustworthiness(coords, coords[:5], k=10)), "mismatched shapes scored anyway"
    assert 0.0 <= M.trustworthiness(coords, kept, k=10, sample=50) <= 1.0


def test_continuous_truth_does_not_fabricate_ranking_classes():
    truth = pd.Series(np.repeat(np.linspace(-2.0, 2.0, 20), 10))
    labels = np.repeat(np.arange(10), 20)
    table = M.ranking(labels, truth)
    report = M.report(labels, truth)
    assert table.empty and table.attrs["continuous_truth"]
    assert np.isnan(report["mean_auprc"])
    assert np.isnan(report["mean_auroc"])


def test_too_few_examples_of_every_category_returns_a_shaped_empty_ranking():
    truth = pd.Series([f"class{i}" for i in range(12)])
    out = M.ranking(np.arange(12) % 2, truth)
    assert out.empty
    assert list(out) == ["category", "n", "auroc", "auprc", "lift", "prevalence"]


def test_the_whole_report_is_one_flat_dict_a_results_table_can_hold():
    coords, truth = world()
    perfect = pd.factorize(truth)[0]
    out = M.report(perfect, truth, coords, X=coords)
    for key in ("mean_auprc", "mean_auroc", "mean_lift", "ari", "ami", "silhouette",
                "knn_purity", "knn_lift", "trustworthiness", "n_clusters", "noise"):
        assert key in out, key
        assert np.isscalar(out[key])
    assert out["knn_lift"] > 3.0
    assert M.report(perfect, truth)["mean_auprc"] > 0.9      # no coordinates: still scores
