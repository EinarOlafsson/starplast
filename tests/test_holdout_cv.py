import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import holdout_cv as H
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


# --------------------------------------------------------------------------- refusals and skips
def test_noise_and_unlabelled_clusters_never_receive_a_call():
    """A cluster of noise, and a cluster whose training members are all unlabelled, are two ways of
    having nothing to learn from. Neither may produce a call, because a call is what later becomes
    a proposed annotation for a gene nobody has measured."""
    labels = np.array([-1, -1, 0, 0, 0, 1, 1, 1])
    truth = np.array([None, None, "IMC", "IMC", "IMC", None, None, None], dtype=object)
    calls = H.fit_cluster_calls(labels, truth, np.arange(len(labels)), min_support=2)
    assert set(calls) == {0}, "noise or an unlabelled cluster was given a call"
    assert calls[0].category == "IMC"


def test_nested_validation_refuses_a_target_it_cannot_split():
    """Both refusals are about the fold structure rather than the biology, and both are better
    raised than worked around: a stratified fold over a category with two members is a fold whose
    test set is sometimes empty."""
    structures = {"a": np.array([0, 0, 1, 1, 1, 1])}
    one_class = pd.Series(["IMC"] * 6)
    with pytest.raises(ValueError, match="at least two target categories"):
        H.nested_structure_cv(structures, one_class, folds=2)
    too_few = pd.Series(["IMC", "IMC", "ER", "ER", "ER", "ER"])
    with pytest.raises(ValueError, match="at least 5 labels"):
        H.nested_structure_cv(structures, too_few, folds=5)


def test_choosing_from_nothing_is_refused_rather_than_guessed():
    with pytest.raises(ValueError, match="no validation scores"):
        H.choose_final_structure(pd.DataFrame(columns=["structure_id", "macro_f1", "coverage",
                                                       "covered_accuracy"]))


def test_pooled_category_scores_average_the_outer_folds():
    """The locked outer folds are the only honest numbers in the run, and what a reader wants from
    them is the mean and the spread -- a category whose F1 is 0.8 in one fold and 0.2 in another is
    not the same finding as one that scores 0.5 twice."""
    per = pd.DataFrame([
        {"category": "IMC", "precision": 0.8, "recall": 0.6, "f1": 0.7, "support": 10},
        {"category": "IMC", "precision": 0.6, "recall": 0.4, "f1": 0.5, "support": 8},
        {"category": "ER", "precision": 0.9, "recall": 0.9, "f1": 0.9, "support": 12},
    ])
    out = H.pooled_category_scores(per)
    assert list(out.category) == ["ER", "IMC"], "not sorted by f1"
    imc = out[out.category == "IMC"].iloc[0]
    assert imc.f1 == pytest.approx(0.6) and imc.support == 18
    assert imc.f1_sd == pytest.approx(0.1414, abs=1e-3)
    assert pd.isna(out[out.category == "ER"].iloc[0].f1_sd), "one fold has no spread to report"


def test_a_gene_in_a_cluster_with_no_call_gets_no_proposal():
    """The quiet half of the inference table. A gene whose cluster never earned a call is not a
    weak prediction, it is an absence of one, and it must not appear as a row."""
    # Cluster 0 earns a call and holds one unlabelled gene; clusters 1 and 2 hold three more with
    # nothing to learn from. Only the first may produce a row.
    labels = np.array([0, 0, 0, 0, 1, 1, 2])
    truth = pd.Series(["IMC", "IMC", "IMC", None, None, None, None], dtype=object)
    genes = np.array([f"TG{i}" for i in range(7)])
    out = H.inference_table(labels, truth, genes, pd.DataFrame(), min_support=2)
    assert list(out.gene_id) == ["TG3"], f"proposed for genes with no call behind them: {out}"
    assert out.iloc[0].proposed == "IMC"
    assert (out.evidence_status == "inferred_from_held_out_structure").all()
