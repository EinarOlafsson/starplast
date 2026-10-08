"""Shared train-only baselines, matched universes and bounded null-control recipes."""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import hashlib
import json

import numpy as np
import pandas as pd

from . import scorecard as C
from .ground_truth import eligibility_mask


@dataclass(frozen=True)
class BaselineRecipe:
    """The meaning and denominator of a task's baseline/control, not its accuracy."""

    task: str
    baselines: tuple[str, ...]
    controls: tuple[str, ...]
    cohort_rule: str
    limitation: str


RECIPES = {
    C.T_LABEL: BaselineRecipe(C.T_LABEL, ("training majority", "training prevalence"), ("protected-group label null",),
                             "same eligible hidden entities and classes", "No test prevalence used for fitting; unseen classes remain visible"),
    C.T_VALUES: BaselineRecipe(C.T_VALUES, ("training mean", "training median"), ("protected-group value null",),
                              "same eligible finite hidden values with declared units", "No imputation of missing truth or central tendency fitted on test"),
    C.T_RANK: BaselineRecipe(C.T_RANK, ("uniform random ranking",), ("matched candidate-universe null", "degree/detection-stratum network control"),
                            "identical frozen candidate universe, seed exclusions and fixed depth", "Unlabelled entities are not verified negatives; PU precision/AUPRC not identifiable"),
    C.T_SET: BaselineRecipe(C.T_SET, ("uniform fixed-size retrieval",), ("matched candidate-universe null",),
                           "identical frozen candidate universe and set size", "Positive-only recall is limited to known positives"),
    C.T_CLUSTER: BaselineRecipe(C.T_CLUSTER, ("single-cluster recovery",), ("protected-group label null on fixed placements",),
                               "same hidden labelled entities, placements and noise mask", "Geometry/stability do not establish biological module accuracy"),
    C.T_REPL: BaselineRecipe(C.T_REPL, ("discovery-selected null findings",), ("disjoint-validation effect/direction null",),
                            "same discovery selection and independent validation population", "Split-half replication does not establish new-dataset replication"),
}


def _ids(values):
    ids = tuple(values)
    if any(not isinstance(value, str) or not value for value in ids) or len(set(ids)) != len(ids):
        raise ValueError("Candidate entities must be unique explicitly named strings")
    return ids


def universe_digest(entities):
    """Pin the identical ordered candidate universe used by strategy and baseline."""
    return hashlib.sha256(json.dumps(_ids(entities), separators=(",", ":")).encode()).hexdigest()


def _training(values, split, kind):
    values = pd.Series(values)
    _ids(values.index)
    split.guard_fit("model_fit", values.index)
    if set(values.index) != set(split.entities("train")):
        raise ValueError("Baseline must use the same complete eligible training cohort")
    if not eligibility_mask(values, kind).all():
        raise ValueError("Missing/unknown training truth cannot become a baseline label or value")
    return values


def label_baselines(training_values, split, *, role="test", seed=0):
    """Predict training majority and training prevalence on the exact named cohort.

    Prevalence columns are class probabilities from training only; randomized
    calls use a separate frozen seed and do not read evaluation truth.
    """
    train = _training(training_values, split, "categorical").astype(str)
    counts = train.value_counts().sort_index()
    majority = min(counts.index[counts == counts.max()])
    probabilities = counts / counts.sum()
    ids = split.entities(role)
    frame = pd.DataFrame({"majority": majority, "prevalence_call": np.random.default_rng(seed).choice(
        probabilities.index.to_numpy(), size=len(ids), p=probabilities.to_numpy())}, index=ids)
    frame.index.name = "entity"
    frame.attrs.update(universe_sha256=universe_digest(ids), training_n=len(train), evaluation_n=len(ids),
                       prevalence=probabilities.to_dict(), split_identity=split.identity, role=role, seed=seed)
    return frame


def value_baselines(training_values, split, *, role="test"):
    """Predict training mean/median on exactly the strategy's named evaluation cohort."""
    train = pd.to_numeric(_training(training_values, split, "numeric"))
    ids = split.entities(role)
    frame = pd.DataFrame({"mean": float(train.mean()), "median": float(train.median())}, index=ids)
    frame.index.name = "entity"
    frame.attrs.update(universe_sha256=universe_digest(ids), training_n=len(train), evaluation_n=len(ids),
                       split_identity=split.identity, role=role)
    return frame


def random_ranking(candidates, *, seed=0, query_seeds=()):
    """Uniformly rank the exact frozen candidate universe after explicit seed exclusion."""
    candidates = _ids(candidates)
    query_seeds = _ids(query_seeds)
    if set(query_seeds) & set(candidates):
        raise ValueError("Query seeds must be excluded before freezing the candidate universe")
    order = np.random.default_rng(seed).permutation(len(candidates))
    frame = pd.DataFrame({"entity": [candidates[i] for i in order], "rank": np.arange(1, len(candidates) + 1)})
    frame.attrs.update(universe_sha256=universe_digest(candidates), candidate_n=len(candidates), seed=seed, query_seeds=query_seeds)
    return frame


def positive_only_retrieval(ranking, known_positives, *, depth):
    """Report identifiable known-positive recovery without inventing precision/negatives.

    Unlabelled items may be positives. Recall refers only to the supplied known
    positives within the matched candidate universe, not all biological positives.
    """
    candidates = _ids(ranking)
    positives = set(_ids(known_positives))
    if not positives <= set(candidates) or type(depth) is not int or depth < 0 or depth > len(candidates):
        raise ValueError("Positive set and fixed depth must match the candidate universe")
    found = len(positives & set(candidates[:depth]))
    return {"candidate_n": len(candidates), "known_positive_n": len(positives), "depth": depth,
            "recovered_known_positives": found, "known_positive_recall": found / len(positives) if positives else None,
            "observed_positive_fraction": found / depth if depth else None,
            "precision": None, "auroc": None, "auprc": None, "negative_semantics": "unlabelled_is_unknown"}


def assert_planted_recovery(ranking, positives, *, depth, minimum_recall):
    """Fail a synthetic harness when known planted positives are not recovered."""
    if not 0 <= minimum_recall <= 1:
        raise ValueError("Declare a finite minimum synthetic recovery fraction")
    result = positive_only_retrieval(ranking, positives, depth=depth)
    if result["known_positive_recall"] is None or result["known_positive_recall"] < minimum_recall:
        raise AssertionError("Planted-positive harness failed; do not report biological success")
    return result


def matched_network_control(edges, detection_strata, *, seed=0, swaps=100, max_attempts=None):
    """Swap undirected edges while preserving degrees and detection-stratum contacts.

    This bounded control never declares newly rewired pairs biological negatives.
    Successful swaps do not prove a uniformly mixed null; insufficient swaps and
    unestablished mixing stay explicit. Weighted/directed controls need a separate
    recipe rather than silently discarding weights or direction.
    """
    if type(swaps) is not int or swaps <= 0:
        raise ValueError("Request a positive integer number of swaps")
    current = []
    for edge in edges:
        if len(edge) != 2:
            raise ValueError("This recipe supports unweighted undirected endpoint pairs only")
        left, right = edge
        _ids((left, right))
        if left not in detection_strata or right not in detection_strata:
            raise ValueError("Every endpoint requires an explicit detection stratum")
        current.append(tuple(sorted((left, right))))
    if len(set(current)) != len(current):
        raise ValueError("Duplicate undirected edges")
    strata = {node: detection_strata[node] for edge in current for node in edge}
    if any(not isinstance(value, str) or not value for value in strata.values()):
        raise ValueError("Detection strata must be explicit nonempty labels")
    original = list(current)
    degree = Counter(node for edge in current for node in edge)
    contacts = Counter(tuple(sorted((strata[a], strata[b]))) for a, b in current)
    maximum = swaps * 50 if max_attempts is None else max_attempts
    if type(maximum) is not int or maximum < 0:
        raise ValueError("Declare a nonnegative bounded attempt budget")
    rng, edge_set, attempts, successful = np.random.default_rng(seed), set(current), 0, 0
    while successful < swaps and attempts < maximum and len(current) >= 2:
        attempts += 1
        i, j = rng.choice(len(current), 2, replace=False)
        a, b = current[i]
        c, d = current[j]
        if rng.random() < .5:
            a, b = b, a
        if rng.random() < .5:
            c, d = d, c
        if len({a, b, c, d}) != 4 or strata[a] != strata[c] or strata[b] != strata[d]:
            continue
        one, two = tuple(sorted((a, d))), tuple(sorted((c, b)))
        if one in edge_set or two in edge_set:
            continue
        edge_set.remove(current[i])
        edge_set.remove(current[j])
        current[i], current[j] = one, two
        edge_set.update((one, two))
        successful += 1
    assert Counter(node for edge in current for node in edge) == degree
    assert Counter(tuple(sorted((strata[a], strata[b]))) for a, b in current) == contacts
    return {"edges": tuple(sorted(current)), "requested_swaps": swaps, "successful_swaps": successful,
            "attempts": attempts, "degree_preserved": True, "detection_contacts_preserved": True,
            "changed_edges": len(set(current) - set(original)), "mixing": "not_established",
            "status": "generated_requires_mixing_review" if successful else "unavailable_no_valid_swaps",
            "biological_negatives": "not_generated", "seed": seed}
