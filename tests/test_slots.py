#!/usr/bin/env python3
"""Slots combine answers to biological questions without losing their source."""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import slots


def slot(policy, patterns=("a_", "b_")):
    return slots.Slot("Toxo", "fixture", "transcription", "condition", "gene",
                      tuple(patterns), policy)


def table():
    return pd.DataFrame({"a_one": [1.0, 2.0, np.nan, np.nan],
                         "a_two": [2.0, 3.0, np.nan, np.nan],
                         "b_one": [10.0, np.nan, 30.0, np.nan],
                         "unrelated": [9.0] * 4})


def test_the_catalog_has_both_organisms_and_stable_prefixed_keys():
    catalog = slots.all_slots()
    # A count rather than a property is a test that fails every time the catalog grows and says
    # nothing about whether the growth was right. What matters is that both arms exist, that keys
    # stay prefixed and resolvable, and that the Plasmodium arm covers the same QUESTIONS as the
    # Toxoplasma one -- which is the thing a mirror can silently fail to do.
    assert len(catalog) > 130
    assert {item.organism for item in catalog} == {"Toxo", "Pf"}
    assert all(item.key.startswith(("Toxo_", "Pf_")) for item in catalog)
    assert slots.by_key(slots.all_slots("Toxo")[0].key) is not None

    toxo_axes = {item.axis for item in slots.all_slots("Toxo")}
    pf_axes = {item.axis for item in slots.all_slots("Pf")}
    assert toxo_axes <= pf_axes, f"axes asked of Toxoplasma and not of Plasmodium: {toxo_axes - pf_axes}"


def test_one_uses_only_the_chosen_candidate():
    out = slots.resolve(table(), slot("one"), chosen="b_")
    assert list(out.values) == ["b_one"]
    assert list(out.source) == ["b_", "", "b_", ""]


def test_average_rank_normalizes_before_combining_candidates():
    out = slots.resolve(table(), slot("average"))
    assert out.values.shape == (4, 1)
    # Each candidate is rank-normalised independently before the per-gene mean.
    assert out.values.iloc[0, 0] == pytest.approx(0.0)
    assert out.source.iloc[0] == "a_+b_" and out.source.iloc[1] == "a_"


def test_fill_uses_best_covered_then_records_gap_sources():
    out = slots.resolve(table(), slot("fill"))
    assert out.values.iloc[:3, 0].notna().all()
    assert list(out.source.iloc[:3]) == ["a_", "a_", "b_"]


def test_fill_accepts_nullable_numeric_extension_arrays():
    frame = pd.DataFrame({"a_one": pd.Series([1, pd.NA, 3], dtype="Int64"),
                          "b_one": pd.Series([pd.NA, 2.5, pd.NA], dtype="Float64")})
    out = slots.resolve(frame, slot("fill"))
    assert out.values.dtypes.iloc[0] == np.dtype("float64")
    assert out.values.iloc[:, 0].tolist() == [1.0, 2.5, 3.0]
    assert out.source.tolist() == ["a_", "b_", "a_"]


def test_separate_keeps_columns_but_shares_one_slot():
    out = slots.resolve(table(), slot("separate"))
    assert list(out.values) == ["a_one", "a_two", "b_one"]
    assert out.source.iloc[0] == "a_+b_"
    assert out.source_columns == ("a_one", "a_two", "b_one")


def test_empty_and_invalid_slots_are_explicit():
    assert slots.resolve(table(), slot("one", ("missing_",))).values.empty
    with pytest.raises(ValueError, match="unknown slot policy"):
        slots.resolve(table(), slot("mystery"))
    with pytest.raises(KeyError, match="unknown slot"):
        slots.resolve(table(), "not_a_slot")


def test_a_missing_or_broken_catalog_is_an_explicit_empty_catalog(monkeypatch, tmp_path):
    monkeypatch.setattr(slots, "CATALOG", str(tmp_path / "missing.json"))
    assert slots.all_slots() == ()
    broken = tmp_path / "broken.json"
    broken.write_text("not json")
    monkeypatch.setattr(slots, "CATALOG", str(broken))
    assert slots.all_slots() == ()


def test_slots_are_leaves_in_three_independent_hierarchies():
    measured = next(s for s in slots.all_slots("Toxo") if s.name == "localization · measured")
    topology = next(s for s in slots.all_slots("Toxo") if s.name == "membrane topology")
    assert measured.biology_path == topology.biology_path
    assert measured.evidence_path != topology.evidence_path
    assert measured.key in slots.relationship_tree("Toxo", "biology") \
        ["cell organization"]["localization and topology"]
    assert topology in slots.slots_in_group(("intrinsic and reference", "sequence-derived"),
                                            hierarchy="evidence")


def test_same_target_family_crosses_evidence_branches():
    family = slots.family_slots("subcellular localization")
    assert {s.name for s in family} >= {"localization · measured", "localization · transferred"}
    assert len({s.evidence_path for s in family}) > 1


def test_an_unknown_hierarchy_is_refused_with_the_choices():
    """Three trees exist and a fourth name is a typo, not a request. Naming the choices matters:
    the hierarchies are 'evidence', 'biology' and 'context', which nobody remembers exactly."""
    import pytest
    with pytest.raises(ValueError, match="unknown hierarchy"):
        slots.hierarchy_path(slots.all_slots("Toxo")[0], "taxonomy")
    for good in slots.HIERARCHIES:
        assert isinstance(slots.hierarchy_path(slots.all_slots("Toxo")[0], good), tuple)
