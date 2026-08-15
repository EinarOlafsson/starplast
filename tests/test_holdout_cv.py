import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast.holdout_cv import (ABSTAIN, choose_final_structure, fit_cluster_calls,
                                  inference_table, nested_structure_cv,
                                  predict_cluster_calls, score_predictions)


def test_cluster_calls_are_learned_only_from_the_supplied_training_rows():
    labels = np.array([0, 0, 0, 1, 1, 1])
    truth = np.array(["a", "a", "b", "b", "b", "a"], object)
    calls = fit_cluster_calls(labels, truth, [0, 1, 3, 4], min_support=2, min_probability=.5)
    assert calls[0].category == "a"
    assert calls[1].category == "b"


def test_unsupported_clusters_abstain():
    labels = np.array([0, 1, -1])
    calls = fit_cluster_calls(labels, np.array(["a", "b", "a"]), [0, 1], min_support=2)
    assert predict_cluster_calls(labels, calls, [0, 1, 2]).tolist() == [ABSTAIN] * 3


def test_abstention_reduces_recall_instead_of_disappearing_from_the_score():
    summary, per = score_predictions(["a", "b"], ["a", ABSTAIN], ["a", "b"])
    assert summary["coverage"] == .5
    assert per.set_index("category").loc["b", "recall"] == 0


def test_nested_selection_uses_inner_validation_and_returns_locked_outer_scores():
    truth = pd.Series((["a"] * 10) + (["b"] * 10), name="target")
    good = np.array(([0] * 10) + ([1] * 10))
    bad = np.arange(20) % 2
    outer, per, validation, splits = nested_structure_cv(
        {"good": good, "bad": bad}, truth, folds=2, min_support=1, min_probability=.5)
    assert set(outer.selected_structure) == {"good"}
    assert outer.macro_f1.min() == 1
    assert len(per) == 4 and len(splits) == 2
    assert choose_final_structure(validation) == "good"


def test_inference_rows_carry_inferred_status_and_cross_validated_confidence():
    truth = pd.Series(["a", "a", None], name="target")
    cv = pd.DataFrame({"category": ["a"], "precision": [.8], "recall": [.5], "f1": [.61]})
    out = inference_table(np.array([0, 0, 0]), truth, ["g1", "g2", "g3"], cv,
                          min_support=2, min_probability=.5)
    assert out.loc[0, "gene_id"] == "g3"
    assert out.loc[0, "evidence_status"] == "inferred_from_held_out_structure"
    assert out.loc[0, "nested_cv_precision"] == .8


def test_explicit_unknown_is_inferred_not_scored_as_a_category():
    truth = pd.Series(["a", "a", "unknown"], name="target")
    cv = pd.DataFrame({"category": ["a"], "precision": [1.], "recall": [1.], "f1": [1.]})
    out = inference_table(np.array([0, 0, 0]), truth, ["g1", "g2", "g3"], cv,
                          min_support=2, min_probability=.5)
    assert out.gene_id.tolist() == ["g3"]
