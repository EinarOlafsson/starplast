#!/usr/bin/env python3
"""Does the search go uphill, and is uphill the direction of the biology?

Two separate claims, tested separately, because they fail separately.

The first is mechanical: given a scoring function, does the climb find its maximum, avoid wasting
evaluations on coordinates that do nothing, stop when told, and report what it saw. Those are tested
against arithmetic scoring functions with known optima, where any failure is the optimiser's.

The second is the one that matters: **is the score a proxy for the thing we want?** That is tested by
building clusterings whose biological quality is known by construction -- the true grouping, a
merged one, a shattered one, a shuffled one -- and asserting the objectives rank them in the order a
biologist would. An optimiser that climbs perfectly up a score that prefers a merged clustering will
return a merged clustering, and every downstream claim is then a claim about nothing.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import discovery as D, metrics as M, optimize as O  # noqa: E402


def config(**kw):
    base = {"method": "umap", "scaling": "rank", "n_components": 3, "n_neighbors": 25,
            "min_dist": 0.25, "perplexity": 30.0, "algorithm": "hdbscan", "min_cluster_size": 25,
            "min_samples": None, "cluster_selection_epsilon": 0.0,
            "cluster_selection_method": "eom", "n_clusters": 15, "linkage": "ward", "eps": 0.5,
            "blocks": ("expression_summary",)}
    base.update(kw)
    return base


# --------------------------------------------------------------------------- does it climb
def test_the_climb_finds_the_top_of_a_hill_it_can_see():
    """A score with one maximum, in one ordinal coordinate. Anything less than finding it is a bug
    in the walk rather than a hard problem."""
    seen = []

    def evaluate(c):
        seen.append(c["n_neighbors"])
        return -abs(c["n_neighbors"] - 60), {}
    out = O.climb(evaluate, config(n_neighbors=5), restarts=1, max_evaluations=200,
                  log=lambda *_a: None)
    assert out.iloc[0]["n_neighbors"] == 60, out[["n_neighbors", "score"]].head().to_string()


def test_restarts_escape_a_hill_that_is_not_the_highest_one():
    """Two peaks, and the start is on the lower. One restart is the difference between reporting a
    local optimum as the answer and finding the real one."""
    def evaluate(c):
        n = c["n_neighbors"]
        return (3.0 if n == 5 else (10.0 if n == 100 else 0.0)), {}
    lonely = O.climb(evaluate, config(n_neighbors=5), restarts=1, max_evaluations=40,
                     log=lambda *_a: None)
    assert lonely.iloc[0].score == 3.0, "the fixture does not trap a single climb"
    many = O.climb(evaluate, config(n_neighbors=5), restarts=6, max_evaluations=400, seed=3,
                   log=lambda *_a: None)
    assert many.iloc[0].score == 10.0


def test_a_configuration_is_never_evaluated_twice():
    """The neighbourhoods of successive steps overlap heavily. Without the cache, roughly half of a
    real search is spent re-running embeddings it has already seen."""
    calls = []

    def evaluate(c):
        calls.append(O.key(c))
        return float(c["min_cluster_size"]), {}
    O.climb(evaluate, config(), restarts=2, max_evaluations=80, log=lambda *_a: None)
    assert len(calls) == len(set(calls))


def test_coordinates_that_do_nothing_are_not_stepped_along():
    """`min_dist` does nothing under t-SNE and `perplexity` nothing under UMAP. A climber that
    varies them spends most of its evaluations re-testing the configuration it is standing on."""
    for method, dead in (("tsne", "min_dist"), ("umap", "perplexity"), ("pca", "n_neighbors")):
        moves = O.neighbours(config(method=method))
        assert all(m[dead] == config()[dead] for m in moves), f"{method} stepped along {dead}"
    hdb = O.neighbours(config(algorithm="kmeans"))
    assert all(m["min_cluster_size"] == 25 for m in hdb), "kmeans varied an HDBSCAN parameter"
    assert any(m["n_clusters"] != 15 for m in hdb), "kmeans never varied its own cluster count"


def test_the_walk_moves_one_coordinate_at_a_time():
    """What makes the path readable afterwards: each step is a decision with a score behind it."""
    start = config()
    for option in O.neighbours(start, block_pool=["expression_summary", "fitness_screens"]):
        differs = [k for k in start if str(start[k]) != str(option.get(k))]
        assert len(differs) == 1, differs


def test_blocks_are_added_and_removed_but_never_emptied():
    pool = ["expression_summary", "fitness_screens", "protein_features"]
    moves = O.neighbours(config(blocks=("expression_summary",)), block_pool=pool)
    blocks = [m["blocks"] for m in moves if m["blocks"] != ("expression_summary",)]
    assert all(len(b) >= 1 for b in blocks)
    assert any(len(b) == 2 for b in blocks), "no block was ever added"
    two = O.neighbours(config(blocks=("expression_summary", "fitness_screens")), block_pool=pool)
    assert any(len(m["blocks"]) == 1 for m in two), "no block was ever removed"


def test_a_stopped_search_returns_what_it_has_rather_than_nothing():
    """Same contract as the rest of this project: a stop is a decision that enough has been seen."""
    calls = []

    def evaluate(c):
        calls.append(1)
        return float(len(calls)), {}
    out = O.climb(evaluate, config(), restarts=3, max_evaluations=500,
                  should_stop=lambda: len(calls) >= 5, log=lambda *_a: None)
    assert 0 < len(out) <= 6
    assert out.score.is_monotonic_decreasing


def test_the_budget_is_a_budget():
    calls = []

    def evaluate(c):
        calls.append(1)
        return 0.0, {}
    out = O.climb(evaluate, config(), restarts=5, max_evaluations=7, log=lambda *_a: None)
    assert len(calls) == 7 and len(out) == 7


def test_a_random_start_is_a_real_configuration():
    rng = np.random.default_rng(0)
    for _ in range(5):
        c = O.random_start(rng, block_pool=["expression_summary", "fitness_screens"])
        assert set(O.SPACE) <= set(c) and c["blocks"]
        assert O.spec_of(c).method in ("umap", "tsne", "pca")
    assert O.random_start(rng)["blocks"] == ("expression_summary",)


def test_a_search_that_evaluated_nothing_returns_an_empty_table_rather_than_raising():
    out = O.climb(lambda c: (0.0, {}), config(), restarts=1, max_evaluations=0,
                  log=lambda *_a: None)
    assert out.empty


def test_the_recipe_and_the_clustering_arguments_come_apart_cleanly():
    spec = O.spec_of(config(method="tsne", perplexity=12.0, blocks=("a", "b")))
    assert spec.method == "tsne" and spec.perplexity == 12.0 and spec.blocks == ("a", "b")
    assert O.cluster_kw(config(algorithm="kmeans")) == {"n_clusters": 15}
    assert set(O.cluster_kw(config(algorithm="hdbscan"))) == {
        "min_cluster_size", "min_samples", "cluster_selection_epsilon",
        "cluster_selection_method"}
    assert O.cluster_kw(config(algorithm="dbscan")) == {"eps": 0.5, "min_samples": None}
    assert "conjunction" in O.MODES


def test_evaluator_computes_trustworthiness_from_the_matrix(monkeypatch):
    coords, nodes, _labels = biology(n_per=15)
    nodes["expr_fixture"] = coords[:, 0]
    monkeypatch.setattr("starplast.embedding.embed",
                        lambda *_a, **_k: (coords, ["expr_fixture"],
                                          np.ones(len(nodes), dtype=bool), coords.copy()))
    monkeypatch.setattr("starplast.optimize.cluster", lambda X, **_k: np.arange(len(X)) % 4)
    _score, extras = O.evaluator(nodes, mode="auprc", layers=("compartment",),
                                 log=lambda *_a: None)(config())
    assert np.isfinite(extras["trustworthiness"])


# ------------------------------------------------------- is uphill the direction of the biology
def biology(n_per=50, seed=0):
    """A world where the answer is known: four compartments, four blobs, some genes unlabelled.

    Cluster 3 also carries a fitness split, so the same fixture exercises both objectives -- half
    its genes are essential and half are not, which is a subdivision no single layer shows.
    """
    rng = np.random.default_rng(seed)
    classes = ["IMC", "rhoptry", "apicoplast", "nucleus"]
    coords, rows, truth_labels = [], [], []
    for i, c in enumerate(classes):
        centre = np.eye(4)[i] * 12.0
        coords.append(rng.normal(centre, 0.7, size=(n_per, 4)))
        for j in range(n_per):
            rows.append({
                "gene_id": f"{c}{j}",
                # A fifth of each group is unlabelled: those are what guilt by association is for.
                "compartment": None if j % 5 == 0 else c,
                "fitness": float(rng.normal(-3 if (i == 3 and j % 2) else 0.0, 0.3)),
                "n_publications": 0})
            truth_labels.append(i)
    return np.vstack(coords), pd.DataFrame(rows), np.array(truth_labels)


def variants(true_labels, rng):
    """The same world clustered four ways, in the order a biologist would rank them."""
    return {
        "true": np.array(true_labels),
        "split": np.array(true_labels) * 2 + (np.arange(len(true_labels)) % 2),
        "merged": np.zeros(len(true_labels), dtype=int),
        "shuffled": rng.permutation(np.array(true_labels)),
    }


def test_the_yield_score_ranks_the_true_clustering_above_the_wrong_ones():
    """The claim the whole optimiser rests on. If a merged or shuffled clustering can out-score the
    real grouping, the search will find that clustering and every finding after it is noise."""
    _, nodes, true_labels = biology()
    v = variants(true_labels, np.random.default_rng(1))
    score = {k: D.yield_score(D.guilt(nodes, lab, "compartment")) for k, lab in v.items()}
    assert score["true"] > score["shuffled"], score
    assert score["true"] > score["merged"], score
    assert score["merged"] == 0.0, "a single cluster over everything produced findings"
    assert score["shuffled"] == 0.0, "a shuffled clustering produced findings"


def test_the_ranking_metrics_rank_them_the_same_way():
    _, nodes, true_labels = biology()
    v = variants(true_labels, np.random.default_rng(2))
    lift = {k: M.summarise(M.ranking(lab, nodes.compartment))["mean_lift"] for k, lab in v.items()}
    assert lift["true"] > lift["shuffled"] and lift["true"] > lift["merged"]
    # Over-splitting is nearly free for annotation and expensive for agreement. Both are true, and
    # a search told to optimise annotation should not be punished for the split.
    assert lift["split"] > lift["merged"]
    ari = {k: M.partition(lab, nodes.compartment)["ari"] for k, lab in v.items()}
    assert ari["true"] > ari["split"]


def test_disagreement_finds_the_subdivision_and_only_where_it_was_built():
    _, nodes, true_labels = biology()
    found = D.disagreement(nodes, np.array(true_labels), "compartment", "fitness")
    assert len(found) == 1, found[["cluster", "category", "gap"]].to_string()
    assert found.iloc[0].category == "nucleus"
    assert found.iloc[0].gap > 1.0


def test_the_objective_reads_the_map_the_climb_would_be_handed():
    """End to end on the real code path, with PCA and kmeans so it runs in a second: the evaluator
    builds the embedding, clusters it, scores it, and hands back the labels and findings."""
    _, nodes, _ = biology()
    rng = np.random.default_rng(4)
    for i in range(6):
        nodes[f"expr_{i}"] = rng.normal(size=len(nodes))
    evaluate = O.evaluator(nodes, mode="guilt", layers=("compartment",), seed=0)
    score, extras = evaluate(config(method="pca", algorithm="kmeans", n_clusters=4,
                                    blocks=("expression_summary",)))
    assert np.isfinite(score)
    assert len(extras["_labels"]) == len(nodes)
    assert "mean_auprc" in extras and "n_clusters" in extras


def test_every_mode_returns_a_finite_score_on_a_map_that_finds_nothing():
    """A configuration that finds nothing scores zero rather than raising, which is what lets the
    climb walk away from it. This was a real crash: pandas will not concatenate an empty list."""
    _, nodes, _ = biology()
    for i in range(6):
        nodes[f"expr_{i}"] = 0.0
    for mode in O.MODES:
        evaluate = O.evaluator(nodes, mode=mode, layers=("compartment",), against=("fitness",))
        score, extras = evaluate(config(method="pca", algorithm="kmeans", n_clusters=2))
        assert np.isfinite(score), mode


def test_every_configuration_is_streamed_as_it_finishes(qapp=None):
    """Same contract as the rest of this project: a search that can be watched, not one that
    reports at the end. The callback carries the row, the configuration and the artefacts, so a
    caller can draw the map that produced a score the moment the score exists."""
    seen = []

    def evaluate(c):
        return float(c["min_cluster_size"]), {"_labels": np.zeros(3)}
    O.climb(evaluate, config(), restarts=1, max_evaluations=6,
            on_step=lambda row, cfg, extras: seen.append((row, cfg, extras)),
            log=lambda *_a: None)
    assert len(seen) == 6
    row, cfg, extras = seen[0]
    assert row["score"] == cfg["min_cluster_size"] and "_labels" in extras
    assert row["restart"] == 0 and row["step"] == 0


def test_the_block_pool_is_what_this_table_can_actually_offer():
    _, nodes, _ = biology()
    # `n_publications` is in the fixture, and it is what the literature block is made of.
    assert O.block_pool(nodes) == ["literature"]
    for i in range(3):
        nodes[f"expr_{i}"] = 0.0
    nodes["fit_hff"] = 0.0
    assert O.block_pool(nodes) == ["expression_summary", "fitness_screens", "literature"]


def test_the_other_embedding_methods_produce_a_map_of_the_right_shape():
    """t-SNE and PCA are alternatives to UMAP, not decorations: a search may step onto either, and
    a method that raised mid-climb would end the run."""
    from starplast.embedding import EmbeddingSpec, embed
    _, nodes, _ = biology(n_per=20)
    rng = np.random.default_rng(7)
    for i in range(5):
        nodes[f"expr_{i}"] = rng.normal(size=len(nodes))
    for method in ("umap", "tsne", "pca"):
        spec = EmbeddingSpec(blocks=("expression_summary",), method=method, n_components=2,
                             perplexity=8.0, n_neighbors=8)
        coords, names, rows = embed(nodes, spec, log=lambda *_a, **_k: None)
        assert coords.shape == (len(nodes), 2), method
        assert np.isfinite(coords).all(), method


def test_pca_asks_for_no_more_components_than_the_matrix_has():
    """A block set with two columns and a request for three components is an error sklearn raises,
    and raising it mid-search ends the whole run over one configuration."""
    from starplast.embedding import EmbeddingSpec, embed
    _, nodes, _ = biology(n_per=10)
    nodes["expr_only"] = np.arange(len(nodes), dtype=float)
    coords, _names, _rows = embed(nodes, EmbeddingSpec(blocks=("expression_summary",),
                                                       method="pca", n_components=3),
                                  log=lambda *_a, **_k: None)
    assert coords.shape[1] == 1


def test_every_clustering_algorithm_returns_one_label_per_point():
    from starplast.clustering import ALGORITHMS, cluster
    # Blobs rather than a cloud of noise: HDBSCAN is right to call structureless points noise, and
    # a test that fed it noise would be testing that it lies.
    rng = np.random.default_rng(8)
    X = np.vstack([rng.normal(c, 0.5, size=(40, 3)) for c in ([0, 0, 0], [8, 0, 0], [0, 8, 0])])
    for algorithm in ALGORITHMS:
        labels = cluster(X, algorithm=algorithm, n_clusters=4, min_cluster_size=10, eps=2.0,
                         linkage="average")
        assert len(labels) == len(X), algorithm
        assert len(set(labels)) > 1, f"{algorithm} put everything in one group"


def test_hdbscan_takes_the_parameters_that_decide_what_it_returns():
    """min_cluster_size alone cannot express "split this further", which is most of why the
    clustering has been hard to get right."""
    from starplast.clustering import cluster
    rng = np.random.default_rng(9)
    X = np.vstack([rng.normal(c, 0.4, size=(60, 2)) for c in ([0, 0], [1.2, 0], [6, 6])])
    eom = cluster(X, algorithm="hdbscan", min_cluster_size=10, cluster_selection_method="eom")
    leaf = cluster(X, algorithm="hdbscan", min_cluster_size=10, cluster_selection_method="leaf")
    assert len(set(leaf)) >= len(set(eom)), "leaf selection did not split more than eom"
    merged = cluster(X, algorithm="hdbscan", min_cluster_size=10, cluster_selection_epsilon=5.0)
    assert len(set(merged)) <= len(set(eom)), "a large epsilon did not merge anything"
