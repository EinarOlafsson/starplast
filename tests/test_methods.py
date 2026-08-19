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


# --------------------------------------------------------------------------- boosting (47.4)
def test_the_raw_matrix_leaves_missing_values_missing():
    """The reason this method exists. Every other path resolves absence before the model sees it --
    impute, drop the column, drop the gene -- and one of those changes WHICH GENES the labels
    describe. A tree learns a split for the missing branch instead, so "not measured" is evidence."""
    frame = pd.DataFrame({"a": [1.0, np.nan, 3.0], "b": [4.0, 5.0, np.nan],
                          "text": ["x", "y", "z"]})
    X, names = M.raw_matrix(frame, ["a", "b", "text", "absent"])
    assert names == ["a", "b"], "a non-numeric or absent column reached the matrix"
    assert np.isnan(X).sum() == 2, "the gaps were filled by something"


def test_no_numeric_column_gives_an_empty_matrix_rather_than_an_error():
    X, names = M.raw_matrix(pd.DataFrame({"text": ["a", "b"]}), ["text"])
    assert X.shape == (2, 0) and names == []


def test_boosting_learns_from_the_pattern_of_missingness_itself():
    """The sharpest statement of the point: a column whose VALUES carry nothing but whose ABSENCE
    tracks the label. Every imputing policy destroys this signal by construction."""
    n = 300
    y = np.array(["a"] * (n // 2) + ["b"] * (n - n // 2))
    rng = np.random.default_rng(0)
    X = rng.normal(size=(n, 4))
    X[y == "a", 0] = np.nan                       # missing exactly for class a
    fit = M.boosted(X, pd.Series(y, dtype="object"), seed=0)
    agree = np.mean([fit["classes"][p] == v for p, v in zip(fit["partition"], y) if p >= 0])
    assert agree > 0.9, f"missingness carried the label and was not learned ({agree:.2f})"


def test_boosting_is_scored_out_of_fold_like_every_other_supervised_method():
    rng = np.random.default_rng(0)
    X = rng.normal(size=(240, 20))
    y = pd.Series(rng.permutation(["a"] * 120 + ["b"] * 120), dtype="object")
    fit = M.boosted(X, y, seed=1)
    agree = np.mean([fit["classes"][p] == v for p, v in zip(fit["partition"], y) if p >= 0])
    assert agree < 0.70, f"noise was 'recovered' at {agree:.2f}; the partition is not out of fold"


def test_importance_is_measured_by_permutation_not_by_split_count():
    """Split counts are the cheap answer and a misleading one: a column with many distinct values
    offers more places to split and accumulates a high count without predicting anything."""
    X, y = _separable(200)
    fit = M.boosted(X, y, seed=0)
    table = M.importances(fit["model"], [f"col{i}" for i in range(X.shape[1])], X, y, seed=0)
    assert set(table.columns) == {"feature", "importance", "sd"}
    assert table.iloc[0].feature == "col0", "the column carrying the signal is not ranked first"
    assert table.importance.is_monotonic_decreasing


def test_importance_without_a_model_or_data_is_empty_rather_than_an_error():
    assert M.importances(None, ["a"], None, None).empty
    X, y = _separable(60)
    fit = M.boosted(X, y, seed=0)
    assert M.importances(fit["model"], [], X, y).empty
    tiny = pd.Series(["a"] * 8 + [None] * 52, dtype="object")
    assert M.importances(fit["model"], ["c"], X, tiny).empty


# --------------------------------------------------------------------------- multiplex (47.3)
def _two_layer_graph(n=60):
    """Two layers that agree on one split and disagree on another, which is the whole point of a
    consensus: the agreed structure must survive and the contested one must not."""
    half = n // 2
    a1, b1 = [], []
    for i in range(half):
        for j in range(i + 1, half):
            a1.append(i); b1.append(j)
    for i in range(half, n):
        for j in range(i + 1, n):
            a1.append(i); b1.append(j)
    a2, b2 = list(a1), list(b1)                     # layer 2 agrees...
    a2 += [0, 1]; b2 += [n - 1, n - 2]              # ...plus two edges crossing the split
    return {"A__a": np.array(a1), "A__b": np.array(b1),
            "B__a": np.array(a2), "B__b": np.array(b2)}


def test_one_layers_communities_are_found_and_isolates_join_none():
    labels = M.layer_communities("A", 60, graph=_two_layer_graph())
    assert len(set(labels[labels >= 0])) == 2, "the two cliques were not separated"
    lonely = M.layer_communities("A", 70, graph=_two_layer_graph())
    assert (lonely[60:] == -1).all(), "a node with no edge was given a community"


def test_the_consensus_is_over_partitions_rather_than_a_merged_graph():
    """Design decision 2: the thirteen edge types are not one graph and are never merged silently.
    Summing adjacency would let the two attention-biased layers pull every community toward the
    well-published genes."""
    out = M.multiplex_communities(["A", "B"], 60, graph=_two_layer_graph())
    labels = out["partition"]
    assert out["settings"]["method"] == "multiplex"
    assert out["settings"]["layers"] == ["A", "B"]
    groups = {g for g in labels if g >= 0}
    assert len(groups) == 2, f"the agreed split did not survive the consensus ({len(groups)})"
    assert labels[0] == labels[1] and labels[-1] == labels[-2]
    assert labels[0] != labels[-1], "the two cliques were merged"


def test_a_layer_that_cannot_see_a_pair_abstains_rather_than_voting_against_it():
    """Layer sizes here span four orders of magnitude -- `ip_ms` has 64 edges against
    `compartment`'s 118,712 -- so a small layer must not veto what a large one found."""
    graph = _two_layer_graph()
    graph["S__a"], graph["S__b"] = np.array([0]), np.array([1])      # sees two nodes only
    out = M.multiplex_communities(["A", "S"], 60, graph=graph)
    labels = out["partition"]
    assert len({g for g in labels if g >= 0}) == 2, "a tiny layer vetoed the structure"


def test_no_layers_gives_no_communities_rather_than_an_error():
    out = M.multiplex_communities([], 40)
    assert (out["partition"] == -1).all() and out["settings"]["layers"] == []
    assert out["unsupervised"]


def test_a_community_holding_half_the_graph_is_not_a_community():
    """The same guard the clustering path applies: a bisection is not structure."""
    n = 40
    a, b = [], []
    for i in range(n):
        for j in range(i + 1, n):
            a.append(i); b.append(j)                 # one clique of everything
    out = M.multiplex_communities(["W"], n, graph={"W__a": np.array(a), "W__b": np.array(b)})
    assert (out["partition"] == -1).all(), "a single all-inclusive community was kept"


# --------------------------------------------------------------------------- constant columns
def test_a_column_with_one_value_is_dropped_and_nan_is_not_a_value():
    """It cannot produce a split, and since scikit-learn 1.9 it cannot be binned either."""
    X = np.array([[1.0, 5.0, np.nan],
                  [2.0, 5.0, 7.0],
                  [3.0, 5.0, np.nan]])
    step = M._DropConstant().fit(X)
    assert step.keep_ == [0] and step.dropped_ == [1, 2]
    assert step.transform(X).shape == (3, 1)


def test_boosting_survives_a_column_that_is_constant_in_one_fold_only(monkeypatch):
    """Which is why the guard sits inside the estimator rather than over the whole matrix.

    The failure it prevents is not a wrong answer but a crash with no author: scikit-learn's binner
    raises `window shape cannot be larger than input array shape` from inside a joblib worker, three
    frames below anything this project wrote.
    """
    rs = np.random.RandomState(0)
    X = rs.normal(size=(80, 3))
    truth = pd.Series(["a", "b"] * 40)
    # Constant among the first fold's training rows, varying elsewhere.
    X[:40, 2] = 1.0
    fit = M.boosted(X, truth, folds=4)
    assert set(fit["partition"]) <= {0, 1}
    assert fit["settings"]["constant_features_dropped"] >= 0
