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


# --------------------------------------------------------------------------- propagation (47.2)
def _line_graph(n=120):
    """Two disconnected cliques, so a label seeded in one must not reach the other."""
    a, b = [], []
    for group in (range(0, n // 2), range(n // 2, n)):
        members = list(group)
        for i in members:
            for j in members:
                if i < j:
                    a.append(i); b.append(j)
    return {"L__a": np.array(a), "L__b": np.array(b)}


def test_a_layer_becomes_a_symmetric_normalised_operator():
    """Symmetric rather than row-stochastic: a random walk on the row-stochastic matrix concentrates
    on hubs, and the hubs of the co-mention layers are the genes people write about most."""
    matrix = M.layer_matrix("L", 120, graph=_line_graph())
    assert matrix.shape == (120, 120)
    dense = matrix.toarray()
    assert np.allclose(dense, dense.T), "the operator is not symmetric"
    assert dense.max() <= 1.0 + 1e-9


def test_an_isolated_node_does_not_divide_by_zero():
    graph = {"L__a": np.array([0]), "L__b": np.array([1])}
    dense = M.layer_matrix("L", 5, graph=graph).toarray()
    assert np.isfinite(dense).all() and dense[4].sum() == 0


def test_a_label_does_not_diffuse_across_a_disconnected_component():
    """The claim propagation makes is "these genes are near the labelled ones". If a label reaches a
    component holding none of its seeds, the score is measuring the walk rather than the graph."""
    matrix = M.layer_matrix("L", 120, graph=_line_graph())
    p = M.Propagator(matrix).fit(np.arange(60).reshape(-1, 1), np.zeros(60, dtype=int))
    assert p.scores_[0][:60].sum() > 0
    assert p.scores_[0][60:].sum() == pytest.approx(0.0, abs=1e-9)


def test_propagation_is_scored_only_on_genes_it_was_not_seeded_from():
    """A gene's own seed makes its own score enormous, so scoring a seeded gene reports the seeding.
    Two cliques with opposite labels: the walk must classify held-out members from their neighbours,
    and it must not do so when the label is unrelated to the graph."""
    truth = pd.Series(["a"] * 60 + ["b"] * 60, dtype="object")
    real = M.propagation("L", np.arange(120), truth, 120, graph=_line_graph())
    agree = np.mean([real["classes"][p] == v for p, v in zip(real["partition"], truth) if p >= 0])
    assert agree > 0.95, "the cliques carry the label and it was not recovered"

    shuffled = pd.Series(np.random.default_rng(0).permutation(truth.to_numpy()), dtype="object")
    fake = M.propagation("L", np.arange(120), shuffled, 120, graph=_line_graph())
    agree_fake = np.mean([fake["classes"][p] == v
                          for p, v in zip(fake["partition"], shuffled) if p >= 0])
    assert agree_fake < 0.75, f"a label unrelated to the graph was recovered at {agree_fake:.2f}"


def test_a_node_the_walk_never_reached_carries_no_evidence():
    """It scores zero for every class and takes the first by argmax; the enrichment gate downstream
    is what discards the group, and this pins that it does not crash on the way."""
    graph = {"L__a": np.array([0, 1]), "L__b": np.array([1, 2])}
    matrix = M.layer_matrix("L", 40, graph=graph)
    p = M.Propagator(matrix).fit(np.array([[0], [1]]), np.array([0, 1]))
    assert set(p.predict(np.arange(40).reshape(-1, 1))) <= {0, 1}


def test_a_propagator_reports_no_coefficients_rather_than_raising():
    """Not every method has per-column weights. An empty frame says so; an exception would claim the
    method was broken."""
    matrix = M.layer_matrix("L", 120, graph=_line_graph())
    fitted = M.Propagator(matrix).fit(np.arange(60).reshape(-1, 1), np.zeros(60, dtype=int))
    assert M.coefficients(fitted, ["a"], ["x"]).empty


def test_the_settings_name_the_layer_that_was_walked():
    """A result that does not say which of the thirteen layers produced it cannot be compared with
    another."""
    out = M.propagation("L", np.arange(120), pd.Series(["a"] * 60 + ["b"] * 60, dtype="object"),
                        120, graph=_line_graph())
    assert out["settings"]["method"] == "propagation:L"


def test_the_shipped_graph_loads_and_its_layers_are_walkable():
    """The default path, which every other test here bypasses by handing in a graph. A loader that
    only works on hand-built dictionaries would pass this whole file and fail on the first real run."""
    import os
    from starplast import paths
    if not os.path.exists(paths.cache_file("graph.npz")):
        pytest.skip("built graph not present")
    matrix = M.layer_matrix("xlms", 8140)
    assert matrix.shape == (8140, 8140) and matrix.nnz > 0
    dense_sum = matrix.sum()
    assert np.isfinite(dense_sum)


def test_the_graph_size_is_read_from_the_graph_rather_than_assumed():
    """What lets a caller be told its table is the wrong one, instead of being shown a confident
    answer about the wrong genes."""
    import os
    from starplast import paths
    if not os.path.exists(paths.cache_file("graph.npz")):
        pytest.skip("built graph not present")
    assert M.graph_size() == 8140


def test_no_graph_means_no_size_rather_than_an_exception(monkeypatch, tmp_path):
    from starplast import paths
    monkeypatch.setattr(paths, "cache_file", lambda name: str(tmp_path / name))
    assert M.graph_size() is None
