#!/usr/bin/env python3
"""Searching the space of maps, rather than enumerating it.

`search.search` walks a grid. A grid is honest and it does not scale: adding t-SNE, four clustering
algorithms and HDBSCAN's real parameters takes the existing sweep from a few hundred configurations
to a few hundred thousand, at ten seconds each. The space has to be searched rather than listed.

## Why hill climbing and not gradient descent

There is no gradient. The score is "how much biology does this configuration find", which runs an
embedding and a clustering and a pile of statistics to produce one number -- it is not differentiable
in `min_cluster_size`, and half the coordinates are categorical (which blocks, which algorithm),
where a derivative is not merely unavailable but undefined.

What IS available is that the space is mostly SMOOTH IN THE ORDINAL COORDINATES: 25 neighbours and
30 neighbours give similar maps, and if 30 beat 25 then 40 is worth trying before 5 is. That is the
property hill climbing exploits, and it is the honest analogue of following a gradient here -- take
a step along one coordinate, keep it if the score improved, and give up on that coordinate when it
does not. Categorical coordinates are enumerated instead, which is what "no gradient" means in
practice.

Two things are added because plain hill climbing has two known failures.

**Restarts.** A hill climber finds the top of whatever hill it started on. Several starts from
different corners of the space, keeping the best, is the cheapest defence and the standard one.

**Caching.** The neighbourhoods of successive steps overlap heavily, and the same configuration is
proposed over and over. Every evaluation is stored under its own key, so the second visit is free --
which in practice roughly halves a run. The embedding is cached separately from the clustering,
because most steps change only the clustering and re-running UMAP for those would be the whole cost
of the search for none of the information.

## What it optimises

Whatever it is told to. `evaluator` builds the three that matter: recovery of a held-out label (the
old objective, which measures trust), guilt by association, and layer disagreement (the two from
`discovery`, which measure yield). They pull in different directions and are meant to -- a map tuned
until it recovers compartment perfectly is often a map that has found nothing new, and seeing that
trade-off is the point of being able to optimise for either.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from .clustering import PARAMS, cluster
from .embedding import BLOCKS, METHODS, EmbeddingSpec, columns_for

#: The coordinates of the space, and the values each may take. Ordinal ones are ordered, and that
#: ordering is load-bearing: a step means "the next value along", which is what makes this a climb
#: rather than an enumeration.
SPACE = {
    "method": list(METHODS),
    "scaling": ["rank", "robust", "zscore"],
    "n_components": [2, 3],
    "n_neighbors": [5, 10, 15, 25, 40, 60, 100],
    "min_dist": [0.0, 0.05, 0.1, 0.25, 0.5, 0.8],
    "perplexity": [5.0, 15.0, 30.0, 50.0, 80.0],
    "algorithm": ["hdbscan", "kmeans", "agglomerative", "dbscan"],
    "min_cluster_size": [5, 10, 15, 25, 40, 60, 100, 150],
    "min_samples": [None, 3, 5, 10, 20],
    "cluster_selection_epsilon": [0.0, 0.1, 0.25, 0.5, 1.0],
    "cluster_selection_method": ["eom", "leaf"],
    "n_clusters": [6, 10, 15, 25, 40, 60],
    "linkage": ["ward", "average", "complete"],
    "eps": [0.25, 0.5, 0.75, 1.0, 1.5],
}

#: Which coordinates are ordered. Everything else is categorical and gets enumerated.
ORDINAL = {"n_components", "n_neighbors", "min_dist", "perplexity", "min_cluster_size",
           "min_samples", "cluster_selection_epsilon", "n_clusters", "eps"}

#: Which coordinates belong to the map, and which to the clustering. The split is what lets an
#: embedding be reused across every clustering tried on it.
EMBEDDING_KEYS = ("blocks", "method", "scaling", "n_components", "n_neighbors", "min_dist",
                  "perplexity")

#: Coordinates that do nothing under the current method or algorithm. Varying one of these is a step
#: that cannot change the score, and a climber that takes them wastes most of its evaluations on
#: configurations identical to the one it is standing on.
RELEVANT = {
    "n_neighbors": lambda c: c.get("method") == "umap",
    "min_dist": lambda c: c.get("method") == "umap",
    "perplexity": lambda c: c.get("method") == "tsne",
    "min_cluster_size": lambda c: c.get("algorithm") == "hdbscan",
    "cluster_selection_epsilon": lambda c: c.get("algorithm") == "hdbscan",
    "cluster_selection_method": lambda c: c.get("algorithm") == "hdbscan",
    "min_samples": lambda c: c.get("algorithm") in ("hdbscan", "dbscan"),
    "eps": lambda c: c.get("algorithm") == "dbscan",
    "n_clusters": lambda c: c.get("algorithm") in ("kmeans", "agglomerative"),
    "linkage": lambda c: c.get("algorithm") == "agglomerative",
}


def key(config: dict) -> tuple:
    """A configuration reduced to what actually affects the answer, for the cache.

    Irrelevant coordinates are dropped rather than included: with `min_dist` in the key, a t-SNE at
    min_dist 0.1 and the same t-SNE at 0.5 are two cache entries for one computation, and the cache
    stops working exactly when the search starts moving between methods.
    """
    live = {k: v for k, v in config.items() if RELEVANT.get(k, lambda _c: True)(config)}
    return tuple(sorted((k, tuple(v) if isinstance(v, (list, tuple)) else v)
                        for k, v in live.items()))


def spec_of(config: dict, seed: int = 42) -> EmbeddingSpec:
    """The embedding half of a configuration, as a recipe."""
    return EmbeddingSpec(
        name="search", blocks=tuple(config.get("blocks", ("expression_summary",))),
        scaling=config.get("scaling", "rank"), method=config.get("method", "umap"),
        n_components=int(config.get("n_components", 3)),
        n_neighbors=int(config.get("n_neighbors", 25)),
        min_dist=float(config.get("min_dist", 0.25)),
        perplexity=float(config.get("perplexity", 30.0)), random_state=seed)


def cluster_kw(config: dict) -> dict:
    """The clustering half, as keyword arguments `clustering.cluster` understands."""
    algorithm = config.get("algorithm", "hdbscan")
    return {k: config[k] for k in PARAMS.get(algorithm, ()) if k in config}


def neighbours(config: dict, space: dict = None, block_pool=()) -> list:
    """Every configuration one step from this one.

    One step means one coordinate moved by one place along its ordering, or set to any other value
    where it has no ordering, or one block added to or removed from the map. Moving one coordinate
    at a time is what makes the walk interpretable afterwards -- the path from the start to the
    winner reads as a series of decisions, each with the score that justified it.
    """
    space = space or SPACE
    out = []
    for coordinate, values in space.items():
        if not RELEVANT.get(coordinate, lambda _c: True)(config):
            continue
        current = config.get(coordinate, values[0])
        if coordinate in ORDINAL:
            i = values.index(current) if current in values else 0
            steps = [j for j in (i - 1, i + 1) if 0 <= j < len(values)]
        else:
            steps = [j for j, v in enumerate(values) if v != current]
        for j in steps:
            out.append({**config, coordinate: values[j]})
    blocks = tuple(config.get("blocks", ()))
    for b in block_pool:
        if b in blocks and len(blocks) > 1:
            out.append({**config, "blocks": tuple(x for x in blocks if x != b)})
        elif b not in blocks:
            out.append({**config, "blocks": tuple(blocks) + (b,)})
    return out


def random_start(rng, space: dict = None, block_pool=()) -> dict:
    """Somewhere in the space, chosen without reference to anywhere already looked at."""
    space = space or SPACE
    config = {k: v[int(rng.integers(len(v)))] for k, v in space.items()}
    pool = list(block_pool) or ["expression_summary"]
    n = int(rng.integers(1, min(len(pool), 3) + 1))
    config["blocks"] = tuple(rng.permutation(pool)[:n])
    return config


def climb(evaluate, start: dict, space: dict = None, block_pool=(), restarts: int = 2,
          max_evaluations: int = 120, seed: int = 42, on_step=None, should_stop=None,
          log=print) -> pd.DataFrame:
    """Walk uphill from `start`, restart when stuck, and return everywhere it looked.

    `evaluate` takes a configuration and returns (score, extras) -- extras is whatever the caller
    wants carried into the results table, typically the labels and the findings behind the score.

    Returns the whole visited set rather than the winner, in the order it was visited, with the
    restart and step number. A search that returns only its best configuration cannot be audited:
    the second-best configuration and how far behind it came is the difference between a real
    optimum and a lucky one.
    """
    rng = np.random.default_rng(seed)
    seen, rows, budget = {}, [], int(max_evaluations)

    def look(config, restart, step):
        nonlocal budget
        k = key(config)
        if k in seen:
            return seen[k]
        if budget <= 0 or (should_stop and should_stop()):
            return None
        budget -= 1
        score, extras = evaluate(config)
        seen[k] = score
        row = {"restart": restart, "step": step, "score": float(score),
               **{c: (",".join(map(str, v)) if isinstance(v, tuple) else v)
                  for c, v in config.items() if RELEVANT.get(c, lambda _c: True)(config)},
               **(extras or {})}
        rows.append(row)
        if on_step:
            on_step(row, config, extras)
        return score

    for restart in range(max(int(restarts), 1)):
        config = dict(start) if restart == 0 else random_start(rng, space, block_pool)
        best = look(config, restart, 0)
        if best is None:
            break
        step = 0
        while True:
            step += 1
            options = neighbours(config, space, block_pool)
            rng.shuffle(options)
            moved = False
            for option in options:
                score = look(option, restart, step)
                if score is None:
                    return _finish(rows, log)
                if score > best + 1e-12:
                    best, config, moved = score, option, True
                    break                      # first improvement, not best: cheaper, and on a
                                               # smooth coordinate it converges to the same place
            if not moved:
                log(f"  restart {restart}: settled at {best:.4f} after {step} steps")
                break
    return _finish(rows, log)


def _finish(rows, log) -> pd.DataFrame:
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out = out.sort_values("score", ascending=False).reset_index(drop=True)
    log(f"  {len(out)} configurations evaluated, best {out.score.iloc[0]:.4f}")
    return out


# --------------------------------------------------------------------------- what to optimise
#: What an optimiser may be told to climb. The first three measure yield -- how much the map says
#: about genes nobody has measured -- and the rest measure quality, which is a different question
#: with a different answer, and the gap between them is worth looking at rather than averaging away.
MODES = ("guilt", "disagreement", "both", "recovery", "auprc", "auroc", "knn")


def evaluator(nodes: pd.DataFrame, mode: str = "guilt", layers=(), against=(), seed: int = 42,
              log=print):
    """Build the scoring function for one of the three questions, with its embeddings cached.

    mode "recovery" -- how well a held-out label comes back. Measures whether the map can be
    trusted, and answers nothing about whether it is useful.
    mode "guilt" -- how much the map predicts about genes nobody has measured (`discovery.guilt`).
    mode "disagreement" -- how many categories the map splits along a second measurement
    (`discovery.disagreement`). `against` names the second layers.
    mode "auprc" / "auroc" / "knn" -- the ranking metrics from `metrics`. These ask "how good a
    shortlist would this map give me for each category", which is the question guilt by association
    actually asks, and unlike the partition scores they cannot be won by merging everything into one
    cluster. "knn" scores the EMBEDDING alone, with no clustering involved, which is what to
    optimise first when the clustering is the thing that keeps going wrong.
    mode "both" -- the sum, for when the question is "find me anything".

    Every mode records every other metric alongside its own score. The cost is a few milliseconds
    against an embedding that took thirty seconds, and it means a run optimised for one thing can be
    re-read against another afterwards instead of being run again.
    """
    from . import discovery, metrics
    from .embedding import embed
    from .search import score_recovery
    cache = {}

    def coords_for(config):
        spec = spec_of(config, seed)
        k = tuple(sorted((c, tuple(v) if isinstance(v, (list, tuple)) else v)
                         for c, v in config.items() if c in EMBEDDING_KEYS))
        if k not in cache:
            cache[k] = embed(nodes, spec, log=lambda *_a, **_k: None)
        return cache[k], spec

    def evaluate(config):
        (coords, names, rows), spec = coords_for(config)
        labels = cluster(coords, algorithm=config.get("algorithm", "hdbscan"),
                         **cluster_kw(config))
        used = set(names)
        sub = nodes.loc[rows] if rows is not None and len(rows) == len(coords) else nodes
        found = []
        if mode in ("guilt", "both"):
            for layer in layers:
                found.append(discovery.guilt(sub, labels, layer, used_columns=used))
        if mode in ("disagreement", "both"):
            for layer in layers:
                for other in against:
                    if other != layer:
                        found.append(discovery.disagreement(sub, labels, layer, other,
                                                            used_columns=used))
        truth = sub[layers[0]] if layers and layers[0] in sub.columns else None
        # Every metric, every time. A run optimised for yield that turns out to have poor AUPRC is
        # a run worth knowing about, and computing it later means embedding everything again.
        scored = (metrics.report(labels, truth, coords) if truth is not None
                  else {})
        if mode == "recovery":
            # (summary, per-label table), and the summary is empty when too few genes were both
            # labelled and clustered -- which is exactly what a bad configuration produces.
            summary, _per = score_recovery(labels, truth) if truth is not None else ({}, None)
            score = float(summary.get("mean_f1", 0.0) or 0.0)
            table = pd.DataFrame()
        elif mode in ("auprc", "auroc", "knn"):
            score = float(scored.get({"auprc": "mean_lift", "auroc": "mean_auroc",
                                      "knn": "knn_lift"}[mode], 0.0) or 0.0)
            score = 0.0 if not np.isfinite(score) else score
            table = pd.DataFrame()
        else:
            # A configuration that found nothing is the commonest outcome, not an error: two giant
            # clusters over eight thousand genes produce no claim that survives the guards, and the
            # optimiser needs that back as a score of zero to climb away from it.
            live = [f for f in found if f is not None and not f.empty]
            table = pd.concat(live, ignore_index=True) if live else pd.DataFrame()
            score = discovery.yield_score(table)
        return score, {**{k: v for k, v in scored.items() if k not in ("_labels", "_findings")},
                       "n_findings": int(len(table)), "_labels": labels, "_findings": table}

    return evaluate


def block_pool(nodes: pd.DataFrame) -> list:
    """Which feature blocks this table can actually offer, in a stable order."""
    return [b for b in BLOCKS if columns_for(nodes, EmbeddingSpec(blocks=(b,))).get(b)]
