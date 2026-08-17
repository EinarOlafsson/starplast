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
    return slots.Slot("Tg", "fixture", "transcription", "condition", "gene",
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
    assert {item.organism for item in catalog} == {"Tg", "Pf"}
    assert all(item.key.startswith(("Tg_", "Pf_")) for item in catalog)
    assert slots.by_key(slots.all_slots("Tg")[0].key) is not None

    toxo_axes = {item.axis for item in slots.all_slots("Tg")}
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
    measured = next(s for s in slots.all_slots("Tg") if s.name == "localization · measured")
    topology = next(s for s in slots.all_slots("Tg") if s.name == "membrane topology")
    assert measured.biology_path == topology.biology_path
    assert measured.evidence_path != topology.evidence_path
    assert measured.key in slots.relationship_tree("Tg", "biology") \
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
        slots.hierarchy_path(slots.all_slots("Tg")[0], "taxonomy")
    for good in slots.HIERARCHIES:
        assert isinstance(slots.hierarchy_path(slots.all_slots("Tg")[0], good), tuple)


def test_a_slot_that_is_not_measured_per_gene_cannot_be_resolved_against_the_node_table():
    """`unit` decides which table a slot belongs to, and getting that wrong fails silently in the
    worst possible way: a host proteome resolved against parasite genes matches nothing, and
    "matches nothing" looks exactly like "nobody has downloaded this yet". Refused instead."""
    import pytest
    from dataclasses import replace
    gene_slot = slots.all_slots("Tg")[0]
    for unit in ("host_gene", "pair", "ortholog_group"):
        with pytest.raises(ValueError, match="cannot be resolved"):
            slots.resolve(table(), replace(gene_slot, unit=unit))
    assert slots.RESOLVABLE_UNIT == "gene"
    # `metabolite` joined the list when the three metabolism slots were corrected: a metabolite is
    # not a gene, nothing measures a compound's concentration "for" a gene, and getting from one to
    # the other needs a metabolic model. They belong to a metabolite table that does not exist yet,
    # exactly as the host slots belong to host tables.
    assert set(slots.UNITS) == {"gene", "host_gene", "pair", "ortholog_group", "metabolite"}


def test_every_slot_in_the_catalog_declares_a_unit_the_code_knows():
    for item in slots.all_slots():
        assert item.unit in slots.UNITS, f"{item.key} declares unit {item.unit!r}"


# --------------------------------------------------------------------------- is a slot filled
def test_a_pair_slot_is_filled_by_the_graph_and_not_by_a_column():
    """The miscount this function exists to prevent. A slot measured per PAIR -- co-expression,
    co-fitness, shared compartment -- lives in graph.npz and has no node column at all. Asking only
    about columns reported eight Toxoplasma slots as empty while their data was already there, and
    the coverage figure went out twice before anyone checked."""
    import types
    pair = next(s for s in slots.all_slots("Tg") if slots.edge_types(s))
    wanted = slots.edge_types(pair)
    assert wanted, "the fixture is not an edge-backed slot"
    graph = types.SimpleNamespace(files=[f"{wanted[0]}__a", f"{wanted[0]}__b"])
    assert slots.is_filled(pair, table(), graph)
    # Columns alone can never fill it, and no graph means not filled rather than an error.
    assert not slots.is_filled(pair, table(), None)
    assert not slots.is_filled(pair, table(), types.SimpleNamespace(files=["something_else__a"]))


def test_a_column_slot_is_filled_by_the_table_and_ignores_the_graph():
    import types
    column_slot = next(s for s in slots.all_slots("Tg")
                       if s.patterns and not slots.edge_types(s))
    empty = pd.DataFrame({"gene_id": ["a", "b"]})
    assert not slots.is_filled(column_slot, empty, types.SimpleNamespace(files=["anything__a"]))
    assert not slots.is_filled(column_slot, None, None)


def test_coverage_is_answered_in_one_place_so_two_reports_cannot_disagree():
    import numpy as np
    import starplast.paths as P
    nodes = pd.read_parquet(os.path.join(P.data_dir(), "nodes.parquet"))
    graph = np.load(os.path.join(P.data_dir(), "graph.npz"))
    out = slots.coverage("Tg", nodes, graph)
    assert out["n_slots"] == len([s for s in slots.all_slots("Tg") if s.role == "feature"])
    assert 0 < out["filled"] < out["n_slots"], out
    assert len(out["empty"]) == out["n_slots"] - out["filled"]
    # Plasmodium has no table of its own yet, so every one of its slots is empty and says so.
    assert slots.coverage("Pf", nodes, graph)["filled"] == 0


def test_no_plasmodium_slot_claims_a_column_of_the_toxoplasma_table():
    """The structural version of 'Plasmodium is empty', which will outlive that being true.

    Every pattern written in the catalog refers to the Toxoplasma node table, because that is the
    only table there is. A Plasmodium slot that inherits one reads as filled by data about the other
    organism -- and `filled` is what the whole atlas is measured by. `codon_` did this the day three
    Toxoplasma sequence columns were added, because the slot lives in `NEW_SHARED` and both arms are
    built from that list.

    When a Plasmodium node table exists this test still holds: its patterns will name ITS columns,
    and claiming a Toxoplasma one would be the same error it is now.
    """
    import os
    import pandas as pd
    import starplast.paths as P
    nodes = pd.read_parquet(os.path.join(P.data_dir(), "nodes.parquet"))
    leaked = {s.name: sorted(slots.declared_columns(nodes, s))
              for s in slots.all_slots("Pf")
              if s.patterns and len(slots.declared_columns(nodes, s))}
    assert not leaked, f"Plasmodium slots claiming Toxoplasma columns: {leaked}"


def test_a_metabolite_slot_is_refused_against_the_gene_table():
    """The refusal that keeps the two tables apart. Resolving a metabolite slot against genes
    matches nothing, and matching nothing looks exactly like nobody having downloaded it."""
    import pandas as pd
    metab = [s for s in slots.all_slots("Tg") if s.unit == "metabolite"]
    assert metab, "the metabolism axis lost its metabolite slots"
    with pytest.raises(ValueError, match="cannot be resolved against a table of 'gene' rows"):
        slots.resolve(pd.DataFrame(index=[0]), metab[0])


def test_a_metabolite_slot_resolves_against_the_metabolite_table():
    import os
    import pandas as pd
    import starplast.paths as P
    path = os.path.join(P.data_dir(), "metabolites.parquet")
    if not os.path.exists(path):
        pytest.skip("metabolite table not built")
    table = pd.read_parquet(path)
    filled = [s for s in slots.all_slots("Tg")
              if s.unit == "metabolite" and slots.declared_columns(table, s)]
    assert filled, "no metabolite slot resolves against the metabolite table"
    got = slots.resolve(table, filled[0], unit="metabolite")
    assert got.values.notna().to_numpy().sum() > 100


def test_every_unit_with_a_table_names_it():
    """`UNIT_TABLES` is a promise that rows exist for that unit; pair and host_gene are absent
    because edges are not rows and the host tables are not built."""
    assert set(slots.UNIT_TABLES) == {"gene", "metabolite"}
    assert "pair" not in slots.UNIT_TABLES and "host_gene" not in slots.UNIT_TABLES
