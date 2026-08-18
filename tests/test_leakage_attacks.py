#!/usr/bin/env python3
"""Instruction 42: attack every guard rather than confirm it.

A guard nobody attacks is a guard nobody has checked. Leakage has got into this project three
separate ways, each blind to the others, and each was found by accident rather than by a test that
was looking for it. Every test here tries to get something past a guard and asserts that it fails.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import search as SE
from starplast.embedding import EmbeddingSpec, columns_for


@pytest.fixture(scope="module")
def nodes():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    path = os.path.join(root, "starplast", "data", "nodes.parquet")
    if not os.path.exists(path):
        pytest.skip("built node table not present")
    return pd.read_parquet(path).sample(n=400, random_state=0).reset_index(drop=True)


# --------------------------------------------------------------------------- 1. into the embedding
def test_a_one_hot_categorical_input_cannot_also_be_the_target(nodes):
    """The categorical route into the matrix, which is not a block and so is not covered by the
    block-dropping guard. Fed as a one-hot feature, `compartment` separates the clusters it is then
    scored against by construction."""
    banned = SE.excluded_for(nodes, "compartment", scope="biology")
    spec = EmbeddingSpec(blocks=("Tg_fitness_hff_in_vitro",), categorical=("compartment",))
    without = SE._spec_without(spec, nodes, banned)
    assert "compartment" not in without.categorical, "the target survived as a one-hot feature"


def test_a_categorical_that_is_not_the_target_survives(nodes):
    """The other half: a guard that removed every categorical would pass the test above and be
    useless."""
    banned = SE.excluded_for(nodes, "cellcycle_phase", scope="direct")
    spec = EmbeddingSpec(blocks=("Tg_fitness_hff_in_vitro",), categorical=("compartment",))
    assert "compartment" in SE._spec_without(spec, nodes, banned).categorical


# --------------------------------------------------------------------------- 2. into the score
def test_the_sweep_never_hands_the_battery_a_column_that_is_both_held_out_and_used(nodes):
    """Asserted PER CATEGORY rather than once. `battery` marks anything the map saw as used or
    derived so only held_out counts as evidence -- but if the sweep passes a column in both lists,
    the row is scored as evidence for a feature the map was built from."""
    from starplast.embedding import SLOT_BLOCKS
    usable = [b for b in SLOT_BLOCKS if columns_for(nodes, EmbeddingSpec(blocks=(b,))).get(b)]
    chosen = []
    for _path, blocks in SE.categories_at("evidence", 1, "Tg"):
        chosen.extend([b for b in blocks if b in usable][:1])
    if len(chosen) < 3:
        pytest.skip("this fixture cannot fill blocks in enough categories")
    spec = EmbeddingSpec(blocks=tuple(chosen))
    held_by_category = {" > ".join(path): set(blocks) & set(spec.blocks)
                        for path, blocks in SE.categories_at("evidence", 1, "Tg")}
    seen = []

    def check(row, labels, kept, held_columns):
        # The blocks the sweep BUILT from are everything ticked minus this category's own, which is
        # what `sweep_categories` does. Deriving it from the category string instead -- the first
        # version of this test -- compared block names against a hierarchy path and "found" a leak
        # that was the test's own arithmetic.
        remaining = tuple(b for b in spec.blocks if b not in held_by_category.get(row["category"], set()))
        used = [c for cols in columns_for(nodes, EmbeddingSpec(blocks=remaining)).values()
                for c in cols]
        seen.append((row["category"], set(held_columns) & set(used)))

    SE.sweep_categories(nodes, spec, level=1, on_clustering=check, log=lambda *a: None)
    assert seen, "no category was swept, so this asserts nothing"
    for category, overlap in seen:
        assert not overlap, f"{category}: {sorted(overlap)[:4]} are both held out and used"


# --------------------------------------------------------------------------- 3. through the rows
@pytest.mark.parametrize("policy", ["median", "drop_rows", "zero"])
def test_labels_and_truth_stay_aligned_after_every_missing_value_policy(nodes, policy):
    """A policy that drops genes changes WHICH GENES the labels describe. Scoring a subset against
    the full table silently aligns cluster 3 with the wrong genes -- and every number downstream
    still looks reasonable."""
    from starplast.embedding import NA_POLICIES, embed
    if policy not in NA_POLICIES:
        pytest.skip(f"{policy} is not a policy this build offers")
    spec = EmbeddingSpec(blocks=("Tg_fitness_hff_in_vitro", "Tg_protein_abundance_tachyzoite"),
                         na_policy=policy)
    try:
        coords, _names, kept = embed(nodes, spec, log=lambda *a: None)
    except ValueError:
        pytest.skip(f"{policy} leaves no usable matrix on this fixture")
    sub = nodes.loc[kept] if kept is not None else nodes
    assert len(coords) == len(sub), "the coordinates and the genes they describe disagree"
    truth = sub["compartment"]
    assert len(truth) == len(coords)
    # And the genes really are the ones the policy kept, not the first N of the table.
    if kept is not None and len(sub) < len(nodes):
        assert list(sub.index) != list(range(len(sub))), \
            "the kept rows look like a positional slice rather than a mask"


def test_a_subset_never_scores_against_the_full_table(nodes):
    """The specific shape of the bug: labels of length N scored against a truth column of length M.
    `score_recovery` must not silently broadcast."""
    labels = np.zeros(len(nodes) // 2, dtype=int)
    labels[::3] = 1
    with pytest.raises((ValueError, IndexError)):
        SE.score_recovery(labels, nodes["compartment"])


# --------------------------------------------------------------------------- 6. the facet source
def test_no_slot_carries_another_organisms_life_cycle_stage():
    """Pins a real defect found while implementing 42.6. Deriving facets from the dataset registry
    put `ring` -- a PLASMODIUM stage -- on the Toxoplasma `shared orthogroup`, `shared domain` and
    `fold confidence` slots, and `tachyzoite` on a Plasmodium one, because shared slots declare
    patterns resolving to registry entries that span both organisms. The registry source was refused;
    this is the invariant that would have caught it."""
    import re

    from starplast import slots
    PF_ONLY = {"ring", "trophozoite", "schizont", "ookinete", "gametocyte"}
    TG_ONLY = {"tachyzoite", "bradyzoite"}
    for organism, forbidden in (("Tg", PF_ONLY), ("Pf", TG_ONLY)):
        for slot in slots.all_slots(organism):
            # The FACETS, not the whole address: the last element is the slot key, and this test's
            # first version failed on `tg_resistance_conferring_mutation` -- for the very reason it
            # was written, since "conferring" contains "ring". Whole words, for the same reason.
            facets = [str(part).lower() for part in slots.hierarchy_path(slot, "context")[:-1]]
            for stage in forbidden:
                for facet in facets:
                    assert not re.search(rf"\b{stage}\b", facet), (
                        f"{organism} slot {slot.key} carries the facet {facet!r}, and {stage!r} "
                        f"belongs to the other organism's life cycle")
