"""Deliberate group, nested-fit, source and endpoint leakage must be refused."""
from dataclasses import replace
from itertools import combinations
import json

import pandas as pd
import pytest

from starplast import organisms as O
from starplast.splits import (Assignment, PairPartition, ROLES, cold_node_pairs, known_node_pairs,
                             make_exclusions, make_split, read_split, write_split)


def _split(**kwargs):
    ids = [f"gene_{i}" for i in range(32)]
    groups = [f"family_{i // 2}" for i in range(32)]
    return make_split(ids, groups, organism=O.TOXOPLASMA, benchmark_id="fixture:v1", **kwargs)


def test_nested_split_preserves_groups_and_every_entity_exactly_once():
    split = _split(seed=8)
    assert split == _split(seed=8)
    assert split.identity != _split(seed=9).identity
    by_group = {}
    for row in split.assignments:
        by_group.setdefault(row.group, set()).add(row.role)
    assert all(len(roles) == 1 for roles in by_group.values())
    assert {row.entity for row in split.assignments} == {f"gene_{i}" for i in range(32)}
    assert all(split.entities(role) for role in ROLES)
    changed = tuple(replace(row, role="test") if row.entity == split.entities("train")[0] else row for row in split.assignments)
    with pytest.raises(ValueError, match="group crosses"):
        replace(split, assignments=changed)


@pytest.mark.parametrize("stage", ["imputation", "scaling", "representation", "feature_selection", "model_fit"])
def test_preprocessing_and_fitting_cannot_see_any_holdout(stage):
    split = _split()
    split.guard_fit(stage, split.entities("train"))
    for role in ("tune", "calibration", "test"):
        with pytest.raises(ValueError, match="Forbidden"):
            split.guard_fit(stage, split.entities(role))


def test_tuning_refitting_and_calibration_respect_disjoint_roles():
    split = _split()
    split.guard_fit("setting_selection", split.entities("train") + split.entities("tune"))
    split.guard_fit("refit", split.entities("train") + split.entities("tune"))
    split.guard_fit("calibration", split.entities("calibration"))
    for stage, role in (("setting_selection", "test"), ("refit", "calibration"), ("calibration", "train"), ("calibration", "test")):
        with pytest.raises(ValueError, match="Forbidden"):
            split.guard_fit(stage, split.entities(role))
    with pytest.raises(ValueError, match="Unknown"):
        split.guard_fit("unregistered", [])
    with pytest.raises(ValueError, match="Forbidden"):
        split.guard_fit("model_fit", ["outside"])


def test_transductive_graph_access_does_not_permit_test_label_fitting():
    split = _split()
    with pytest.raises(ValueError, match="feature access"):
        split.guard_joint_features(split.entities("test"))
    transductive = replace(split, feature_access="transductive")
    transductive.guard_joint_features(transductive.entities("test"))
    with pytest.raises(ValueError, match="Forbidden"):
        transductive.guard_fit("model_fit", transductive.entities("test"))


def test_target_source_and_derived_layer_contamination_fail_closed():
    split = make_split(list("abcdefgh"), list("abcdefgh"), organism=O.TOXOPLASMA, benchmark_id="fixture:v1",
                       fractions=(.99,.003,.003,.004))
    nodes = pd.DataFrame({"gene_id": split.entities("train")[:4], "compartment": ["n", "c", "n", "c"],
                          "lopit_prob_map": [.9, .8, .95, .85], "unrelated": [2, 5, 0, 3]})
    excluded = make_exclusions(nodes, "compartment", benchmark_id="fixture:v1", split=split)
    excluded.guard_inputs(columns=["unrelated"])
    for column in ("compartment", "lopit_prob_map"):
        with pytest.raises(ValueError, match="contamination"):
            excluded.guard_inputs(columns=[column])
    fitness = nodes.assign(fit_invitro_hff=[1., 2., 3., 4.])
    closure = make_exclusions(fitness, "fit_invitro_hff", benchmark_id="fixture:v1", split=split)
    assert "cofitness" in closure.layers and "structural_hole" in closure.layers
    for layer in ("cofitness", "structural_hole"):
        with pytest.raises(ValueError, match="contamination"):
            closure.guard_inputs(layers=[layer])
    with pytest.raises(ValueError, match="different benchmark"):
        excluded.guard_inputs(benchmark_id="another")
    with pytest.raises(ValueError, match="not in"):
        make_exclusions(nodes, "not_present", benchmark_id="fixture:v1", split=split)
    with pytest.raises(ValueError, match="Forbidden"):
        make_exclusions(nodes.assign(gene_id=list(split.entities("test")) * 4), "compartment", benchmark_id="fixture:v1", split=split)


def test_known_node_and_cold_node_regimes_cannot_be_conflated():
    nodes = [f"node_{i}" for i in range(16)]
    pairs = list(combinations(nodes, 2))
    known = known_node_pairs(pairs, organism=O.TOXOPLASMA, benchmark_id="pairs:v1", seed=4)
    training_nodes = {node for a, b, role in known.pairs if role == "train" for node in (a, b)}
    assert training_nodes == set(nodes)
    assert {role for _a, _b, role in known.pairs} == set(ROLES)
    node_split = make_split(nodes, nodes, organism=O.TOXOPLASMA, benchmark_id="cold:v1", group_kind="entity")
    cold = cold_node_pairs(pairs, node_split)
    train = {node for a, b, role in cold.pairs if role == "train" for node in (a, b)}
    test = {node for a, b, role in cold.pairs if role == "test" for node in (a, b)}
    assert train and test and not train & test
    assert len(cold.pairs) + len(cold.excluded_cross_partition) == len(pairs)
    left = node_split.entities("train")[0]
    right = node_split.entities("test")[0]
    with pytest.raises(ValueError, match="endpoint"):
        PairPartition("cold_node", node_split, ((left, right, "test"),))
    with pytest.raises(ValueError, match="frozen edge"):
        PairPartition("known_node", known.split, ((known.pairs[0][0], known.pairs[0][1],
                                                  "test" if known.pairs[0][2] != "test" else "train"),))


def test_sparse_pair_graphs_and_duplicate_pairs_do_not_invent_holdouts():
    with pytest.raises(ValueError, match="four"):
        known_node_pairs([("a", "b"), ("b", "c")], organism=O.TOXOPLASMA, benchmark_id="sparse")
    with pytest.raises(ValueError, match="Duplicate"):
        known_node_pairs([("a", "b"), ("b", "a")], organism=O.TOXOPLASMA, benchmark_id="duplicate")
    with pytest.raises(ValueError):
        cold_node_pairs([("a", "missing")], _split())


@pytest.mark.parametrize("kind", ["dataset", "context"])
def test_dataset_and_context_groups_are_protected_without_claiming_homology(kind):
    split = make_split([f"observation_{i}" for i in range(16)], [f"source_{i // 2}" for i in range(16)],
                       organism=O.HUMAN, benchmark_id="observations", group_kind=kind)
    assert split.group_kind == kind
    assert len({row.group for row in split.assignments}) == 8


def test_chronological_protocol_keeps_later_groups_out_of_earlier_training():
    dates = [f"2026-01-{i:02d}" for i in range(1, 17)]
    split = make_split(dates, dates, organism=O.HUMAN, benchmark_id="time", group_kind="time", ordered=True)
    assert max(split.entities("train")) < min(split.entities("tune"))
    assert max(split.entities("tune")) < min(split.entities("calibration"))
    assert max(split.entities("calibration")) < min(split.entities("test"))
    with pytest.raises(ValueError, match="chronological"):
        make_split(dates, dates, organism=O.HUMAN, benchmark_id="time", group_kind="time")


def test_split_roundtrip_and_mutation_detection(tmp_path):
    split = _split()
    path = tmp_path / "split.json"
    write_split(path, split)
    assert read_split(path) == (split, None)
    with pytest.raises(FileExistsError):
        write_split(path, split)
    payload = json.loads(path.read_text())
    payload["split"]["assignments"][0]["role"] = "test"
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="identity changed"):
        read_split(path)


def test_exclusion_manifest_is_bound_to_training_population_and_split(tmp_path):
    split = _split()
    train = list(split.entities("train"))
    nodes = pd.DataFrame({"gene_id": train, "compartment": ["n" if i % 2 else "c" for i in range(len(train))]})
    exclusions = make_exclusions(nodes, "compartment", benchmark_id=split.benchmark_id, split=split)
    path = tmp_path / "paired.json"
    write_split(path, split, exclusions)
    assert read_split(path) == (split, exclusions)
    with pytest.raises(ValueError, match="different split"):
        write_split(tmp_path / "invalid.json", replace(split, seed=99), exclusions)
    with pytest.raises(ValueError, match="Forbidden"):
        write_split(tmp_path / "test_selected.json", split, replace(exclusions, fit_entities=split.entities("test")))


def test_direct_time_manifest_cannot_reverse_chronological_roles():
    dates = [f"2026-01-{i:02d}" for i in range(1, 17)]
    split = make_split(dates, dates, organism=O.HUMAN, benchmark_id="time", group_kind="time", ordered=True)
    rows = tuple(replace(row, role="test" if row.role == "train" else "train" if row.role == "test" else row.role)
                 for row in split.assignments)
    with pytest.raises(ValueError, match="Chronological"):
        replace(split, assignments=rows)
