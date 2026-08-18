#!/usr/bin/env python3
"""Answering a recipe with a classifier instead of a clustering.

The danger here is the opposite of the clustering's. A clustering never sees the label, so it cannot
memorise it; a classifier is fitted to the label and will report memorisation as recovery unless
every prediction it is scored on came from a model that never saw that gene. These tests attack that
directly, and they include the same shuffled-label control the clustering path carries.
"""
import numpy as np
import pandas as pd
import pytest

from starplast import methods as M


def _separable(n=300, noise=0.0, seed=0):
    """Two classes that a linear model can separate, with a knob for how much signal is left."""
    rng = np.random.default_rng(seed)
    y = np.array(["a"] * (n // 2) + ["b"] * (n - n // 2))
    X = rng.normal(size=(n, 6))
    X[:, 0] += np.where(y == "a", 2.5, -2.5) + rng.normal(scale=noise, size=n)
    return X, pd.Series(y, dtype="object")


def test_a_class_too_small_to_learn_is_dropped_rather_than_predicted_never():
    """Kept, a rare class lets a model score well by never predicting it, and the per-label table
    then shows a row of zeros that reads as a failure of the biology."""
    y = pd.Series(["a"] * 100 + ["b"] * 100 + ["rare"] * 5, dtype="object")
    mask, _values = M._usable(y)
    assert mask.sum() == 200 and not mask[-1]


def test_unlabelled_genes_are_not_treated_as_a_class():
    y = pd.Series(["a"] * 60 + ["b"] * 60 + [None] * 40, dtype="object")
    mask, _ = M._usable(y)
    assert mask.sum() == 120


def test_a_classifier_is_scored_only_on_genes_it_never_saw():
    """The mistake this prevents: scoring training predictions reports memorisation as recovery, and
    it is far easier to make here than in the clustering, which never sees the label at all.

    A model with enough capacity fits the training data perfectly. If the reported partition were the
    training prediction, a label with NO signal would still come back perfect.
    """
    rng = np.random.default_rng(0)
    X = rng.normal(size=(240, 40))                     # 40 columns of pure noise
    y = pd.Series(rng.permutation(["a"] * 120 + ["b"] * 120), dtype="object")
    fit = M.logistic(X, y, seed=1)
    agree = np.mean([fit["classes"][p] == v for p, v in zip(fit["partition"], y) if p >= 0])
    assert agree < 0.70, f"noise was 'recovered' at {agree:.2f}; the partition is not out of fold"


def test_a_shuffled_label_does_not_recover_under_the_classifier_either():
    """The negative control the clustering path already carries, applied to this method. If a
    shuffled label scores like the real one, the score is measuring the fit rather than the biology."""
    from starplast.search import score_recovery
    X, y = _separable(300)
    real = M.logistic(X, y, seed=0)
    real_score, _ = score_recovery(np.asarray(real["partition"]), y)
    shuffled = pd.Series(np.random.default_rng(3).permutation(y.to_numpy()), dtype="object")
    fake = M.logistic(X, shuffled, seed=0)
    fake_score, _ = score_recovery(np.asarray(fake["partition"]), shuffled)
    assert fake_score.get("mean_f1", 0) < real_score.get("mean_f1", 1), (
        f"a shuffled label scored {fake_score.get('mean_f1')} against the real "
        f"{real_score.get('mean_f1')}")


def test_real_signal_is_recovered_so_the_control_above_means_something():
    """The other half of the negative control: a test that nothing recovers is passed by a broken
    model."""
    from starplast.search import score_recovery
    X, y = _separable(300)
    fit = M.logistic(X, y, seed=0)
    summary, _ = score_recovery(np.asarray(fit["partition"]), y)
    assert summary["mean_f1"] > 0.8


def test_too_few_labelled_genes_gives_no_partition_rather_than_an_exception():
    X = np.random.default_rng(0).normal(size=(40, 5))
    y = pd.Series([None] * 36 + ["a", "a", "b", "b"], dtype="object")
    fit = M.logistic(X, y)
    assert set(np.unique(fit["partition"])) == {-1} and fit["classes"] == []


def test_unlabelled_genes_are_predicted_by_a_model_trained_on_every_labelled_one():
    """They are the answer, they were in no fold, and holding a fold out of the model that names them
    would only make that answer worse."""
    X, y = _separable(240)
    y = y.copy()
    y.iloc[:40] = None
    fit = M.logistic(np.asarray(X), y, seed=0)
    assert (np.asarray(fit["partition"])[:40] >= 0).all()


def test_the_model_says_which_measurements_it_used():
    """The reason a linear baseline earns its place: a cluster says "these genes go together"; this
    says "this label is predicted by these columns, in this direction"."""
    X, y = _separable(300)
    fit = M.logistic(X, y, seed=0)
    coef = M.coefficients(fit["model"], fit["classes"], [f"col{i}" for i in range(X.shape[1])])
    assert len(coef) and set(coef.columns) == {"label", "feature", "coefficient"}
    # Column 0 is the one carrying the signal, so it must be the strongest.
    assert coef.reindex(coef.coefficient.abs().sort_values(ascending=False).index).iloc[0].feature \
        == "col0"


def test_coefficients_of_an_unfitted_model_are_empty_rather_than_an_error():
    assert M.coefficients(None, [], []).empty


def test_a_model_with_one_coefficient_vector_still_reports(monkeypatch):
    """Not every classifier is one-vs-rest; a plain `coef_` must still be readable."""
    class Plain:
        coef_ = np.array([[0.5, -0.2]])
    out = M.coefficients(Plain(), ["a", "b"], ["x", "y"])
    assert list(out.feature) == ["x", "y"]
