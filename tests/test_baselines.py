"""Matched denominators and unknown negatives are required for fair baseline tests."""
from collections import Counter
from itertools import combinations

import numpy as np
import pandas as pd
import pytest

from starplast import organisms as O
from starplast.baselines import (RECIPES, assert_planted_recovery, label_baselines, matched_network_control,
                                positive_only_retrieval, random_ranking, universe_digest, value_baselines)
from starplast.ground_truth import TASKS
from starplast.splits import make_split


def _split():
    ids = [f"gene_{i}" for i in range(40)]
    return make_split(ids, ids, organism=O.TOXOPLASMA, benchmark_id="fixture:v1", group_kind="entity", seed=4)


def test_all_task_recipes_keep_task_specific_denominators_and_limits():
    assert set(RECIPES) == set(TASKS)
    assert all(recipe.baselines and recipe.controls and recipe.cohort_rule and recipe.limitation for recipe in RECIPES.values())


def test_label_prevalence_and_majority_are_train_only_on_the_exact_test_cohort():
    split = _split()
    ids = split.entities("train")
    train = pd.Series(["common"] * (len(ids) - 2) + ["rare"] * 2, index=ids)
    predictions = label_baselines(train, split, seed=5)
    assert tuple(predictions.index) == split.entities("test")
    assert predictions.majority.eq("common").all()
    assert predictions.attrs["prevalence"]["rare"] == 2 / len(ids)
    assert predictions.attrs["evaluation_n"] == len(split.entities("test"))
    assert predictions.attrs["universe_sha256"] == universe_digest(split.entities("test"))
    assert predictions.equals(label_baselines(train, split, seed=5))


def test_central_tendency_never_uses_test_values_and_refuses_cohort_substitution():
    split = _split()
    train = pd.Series(np.arange(len(split.entities("train")), dtype=float), index=split.entities("train"))
    result = value_baselines(train, split)
    assert result["mean"].eq(train.mean()).all() and result["median"].eq(train.median()).all()
    assert tuple(result.index) == split.entities("test")
    poisoned = pd.concat([train, pd.Series(1000000., index=split.entities("test"))])
    with pytest.raises(ValueError, match="Forbidden"):
        value_baselines(poisoned, split)
    with pytest.raises(ValueError, match="complete eligible"):
        value_baselines(train.iloc[:-1], split)
    with pytest.raises(ValueError, match="Missing"):
        value_baselines(train.mask(train.index == train.index[0]), split)


def test_unknown_labels_cannot_enter_majority_or_prevalence():
    split = _split()
    values = pd.Series("known", index=split.entities("train"))
    values.iloc[0] = "unknown"
    with pytest.raises(ValueError, match="Missing/unknown"):
        label_baselines(values, split)


def test_random_rank_and_fixed_depth_preserve_candidate_denominator_and_seed_exclusion():
    candidates = [f"gene_{i}" for i in range(50)]
    result = random_ranking(candidates, seed=2, query_seeds=["seed_outside"])
    assert set(result.entity) == set(candidates) and len(result) == len(candidates)
    assert result.attrs["universe_sha256"] == universe_digest(candidates)
    assert result.equals(random_ranking(candidates, seed=2))
    with pytest.raises(ValueError, match="seeds"):
        random_ranking(candidates, query_seeds=[candidates[0]])
    with pytest.raises(ValueError, match="unique"):
        random_ranking(["a", "a"])


def test_positive_only_recovery_never_reports_unknown_entities_as_negatives():
    result = positive_only_retrieval(["positive", "unknown_1", "unknown_2"], ["positive"], depth=2)
    assert result["known_positive_recall"] == 1
    assert result["observed_positive_fraction"] == .5
    assert result["precision"] is result["auroc"] is result["auprc"] is None
    with pytest.raises(ValueError, match="universe"):
        positive_only_retrieval(["a"], ["outside"], depth=1)


def test_planted_positive_harness_detects_reversed_scores_and_null_fixture():
    ordered = [f"gene_{i}" for i in range(100)]
    positives = ordered[:10]
    assert_planted_recovery(ordered, positives, depth=10, minimum_recall=1)
    with pytest.raises(AssertionError, match="harness failed"):
        assert_planted_recovery(ordered[::-1], positives, depth=10, minimum_recall=.9)
    random = random_ranking(ordered, seed=7)
    with pytest.raises(AssertionError, match="harness failed"):
        assert_planted_recovery(random.entity.tolist(), positives, depth=10, minimum_recall=.9)


def test_network_control_preserves_degrees_and_detection_contacts_without_claiming_mixing():
    edges = [(f"node_{i}", f"node_{(i+1)%20}") for i in range(20)]
    strata = {f"node_{i}": "detected_high" if i % 2 else "detected_low" for i in range(20)}
    result = matched_network_control(edges, strata, swaps=20, seed=8)
    assert result["successful_swaps"] == 20 and result["changed_edges"] > 0
    assert Counter(node for edge in edges for node in edge) == Counter(node for edge in result["edges"] for node in edge)
    contacts = lambda ee: Counter(tuple(sorted((strata[a], strata[b]))) for a, b in ee)
    assert contacts(edges) == contacts(result["edges"])
    neighbours = lambda ee: Counter((node, strata[other]) for a, b in ee for node, other in ((a, b), (b, a)))
    assert neighbours(edges) == neighbours(result["edges"])
    assert result["mixing"] == "not_established" and result["biological_negatives"] == "not_generated"
    assert result == matched_network_control(edges, strata, swaps=20, seed=8)


def test_networks_without_valid_swaps_are_explicitly_untestable_controls():
    nodes = list("abcde")
    edges = list(combinations(nodes, 2))
    result = matched_network_control(edges, dict.fromkeys(nodes, "detected"), swaps=5, max_attempts=30)
    assert result["successful_swaps"] == result["changed_edges"] == 0
    assert result["status"] == "unavailable_no_valid_swaps"
    with pytest.raises(ValueError, match="stratum"):
        matched_network_control(edges, {})
    with pytest.raises(ValueError, match="unweighted"):
        matched_network_control([("a", "b", .9)], {"a": "d", "b": "d"})
