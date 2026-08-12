"""Validating an annotation: hide labels we already have, and see whether they come back.

The point of this module is to attach an error rate to a candidate list. These tests are mostly
about the ways it could produce a flattering number: scoring against genes that were visible when
the cluster was chosen, scoring a label the embedding was built from, or reporting one average
across categories that behave completely differently.
"""
import numpy as np
import pandas as pd
import pytest

from starplast import validate as V


def _labels_and_truth(n_per=40):
    """Three clean clusters, each one category, plus unlabelled genes mixed in."""
    labels = np.repeat([0, 1, 2], n_per)
    truth = pd.Series(["a"] * n_per + ["b"] * n_per + ["c"] * n_per)
    return labels, truth


def test_a_perfect_clustering_recovers_what_was_hidden():
    labels, truth = _labels_and_truth()
    r = V.masked_recovery(labels, truth, "a", folds=3, seed=0)
    assert r.recall == pytest.approx(1.0), "hidden genes were in the right cluster and not found"
    assert r.precision == pytest.approx(1.0)


def test_a_clustering_that_ignores_the_category_scores_badly():
    """The check that the metric can fail: labels assigned at random must not look like recovery."""
    rng = np.random.default_rng(0)
    truth = pd.Series(["a"] * 60 + ["b"] * 60)
    labels = rng.integers(0, 3, size=120)
    r = V.masked_recovery(labels, truth, "a", folds=5, seed=1)
    assert r.precision < 0.6, f"random clusters scored precision {r.precision}"


def test_hidden_genes_do_not_influence_which_cluster_is_chosen():
    """Otherwise the test marks its own homework: the cluster would be picked using the answer."""
    labels = np.array([0] * 10 + [1] * 10)
    truth = pd.Series(["a"] * 10 + ["a"] * 10)
    r = V.masked_recovery(labels, truth, "a", folds=4, seed=0, hold_frac=0.5)
    # Both clusters are all "a", so whichever is chosen from the visible half, the hidden half is
    # split across both -- recall must be well under 1, not 1.
    assert r.recall < 1.0


def test_scoring_a_label_the_embedding_used_is_refused():
    """A cluster matching a feature the map was built from is circular and measures nothing."""
    labels, truth = _labels_and_truth()
    with pytest.raises(ValueError, match="circular"):
        V.masked_recovery(labels, truth, "a", used_columns=["compartment", "a"])


def test_a_category_with_too_few_genes_says_so_rather_than_scoring():
    labels = np.array([0, 0, 1])
    truth = pd.Series(["rare", "rare", "other"])
    r = V.masked_recovery(labels, truth, "rare", folds=3)
    assert r.folds == []
    assert "too few" in r.note
    assert np.isnan(r.f1)


def test_noise_is_never_chosen_as_the_annotating_cluster():
    """"It is in the noise" is not an annotation."""
    labels = np.array([-1] * 30 + [0] * 10)
    truth = pd.Series(["a"] * 40)
    r = V.masked_recovery(labels, truth, "a", folds=3, seed=0)
    assert all(f.cluster >= 0 for f in r.folds)


def test_absence_labels_are_never_treated_as_a_category():
    labels = np.repeat([0, 1], 40)
    truth = pd.Series(["unassigned"] * 40 + ["a"] * 40)
    t = V.validate_all(labels, truth, folds=3, min_size=5)
    assert "unassigned" not in set(t.category)


def test_every_category_is_reported_separately():
    """One global number would hide that a method works for one compartment and not another."""
    labels, truth = _labels_and_truth()
    t = V.validate_all(labels, truth, folds=3, min_size=5)
    assert set(t.category) == {"a", "b", "c"}
    assert list(t.f1) == sorted(t.f1, reverse=True), "not ranked"


def test_a_table_with_nothing_scorable_has_the_right_columns():
    """An empty frame with no columns breaks every caller that reads one."""
    labels = np.array([0, 1])
    truth = pd.Series(["a", "b"])
    t = V.validate_all(labels, truth, min_size=100)
    assert list(t.columns) == ["category", "n_labelled", "n_folds",
                              "precision", "recall", "f1", "note"]


def test_candidates_carry_the_numbers_that_say_how_much_to_believe_them():
    """A candidate must never travel without its cluster's composition."""
    labels = np.array([0] * 10)
    truth = pd.Series(["a"] * 6 + ["b"] * 2 + ["unassigned"] * 2)
    genes = [f"g{i}" for i in range(10)]
    c = V.candidates(labels, truth, "a", 0, genes)
    assert list(c.gene_id) == ["g8", "g9"], "only the unlabelled members are candidates"
    assert c.cluster_frac_category.iloc[0] == pytest.approx(0.6)
    assert c.cluster_frac_contradicting.iloc[0] == pytest.approx(0.2)


def test_the_summary_says_whether_it_re_embedded():
    """A number that did not re-fit is a weaker claim than one that did, and must say so."""
    labels, truth = _labels_and_truth()
    r = V.masked_recovery(labels, truth, "a", folds=2)
    assert "hidden only from scoring" in r.summary()
    r.refit = True
    assert "re-embedded per fold" in r.summary()


def test_missing_values_are_absence_not_a_category():
    labels = np.repeat([0, 1], 30)
    truth = pd.Series([None] * 30 + ["a"] * 30, dtype="object")
    t = V.validate_all(labels, truth, folds=2, min_size=5)
    assert set(t.category) <= {"a"}
