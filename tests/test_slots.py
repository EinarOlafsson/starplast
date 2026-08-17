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

    The Plasmodium node table now exists, and the rule survives it in a stronger form. Its columns
    are named the SAME as Toxoplasma's wherever the quantity is the same -- `length` is a protein
    length in both -- because that is what lets the two arms be read side by side. So disjoint
    pattern strings can no longer be what keeps them apart; the table's own accessions are.
    """
    import os
    import pandas as pd
    import starplast.paths as P
    nodes = pd.read_parquet(os.path.join(P.data_dir(), "nodes.parquet"))
    leaked = {s.name: sorted(slots.declared_columns(nodes, s))
              for s in slots.all_slots("Pf")
              if s.patterns and len(slots.declared_columns(nodes, s))}
    assert not leaked, f"Plasmodium slots claiming Toxoplasma columns: {leaked}"
    # And the mirroring is real rather than incidental: at least one Pf pattern IS a Toxoplasma
    # column name, so the assertion above is being held by the species guard rather than by an
    # accident of naming that some later rename would quietly remove.
    shared = {p for s in slots.all_slots("Pf") for p in s.patterns if p in nodes.columns}
    assert shared, "no Pf pattern collides with a Tg column: this test no longer tests anything"


def test_a_plasmodium_slot_is_refused_against_the_toxoplasma_table():
    """The refusal, not merely an empty result -- resolving is where a wrong number would ship."""
    import os
    import pandas as pd
    import starplast.paths as P
    nodes = pd.read_parquet(os.path.join(P.data_dir(), "nodes.parquet"))
    slot = next(s for s in slots.all_slots("Pf") if s.patterns)
    with pytest.raises(ValueError, match="measured in"):
        slots.resolve(nodes, slot)


def test_the_species_of_a_table_is_read_from_its_accessions():
    assert slots.table_organism(pd.DataFrame({"gene_id": ["TGME49_200010"]})) == "Tg"
    assert slots.table_organism(pd.DataFrame({"gene_id": ["PF3D7_0100100"]})) == "Pf"
    assert slots.table_organism(pd.DataFrame({"gene_id": ["something else"]})) is None
    # No id column at all: fall back to the index, which is how several callers hold it.
    assert slots.table_organism(pd.DataFrame(index=["PF3D7_0100100"])) == "Pf"


def test_a_table_of_unknown_species_is_allowed_through():
    """Synthetic frames in tests carry no accessions; refusing those would protect nothing."""
    slot = next(s for s in slots.all_slots("Pf") if s.patterns)
    assert slots.same_species(pd.DataFrame({"length": [1]}), slot)


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
    # `host_gene` joined when the host table gained columns of its own rather than only keys.
    # `pair` stays out: edges and bridges are not rows.
    assert set(slots.UNIT_TABLES) == {"gene", "metabolite", "host_gene"}
    assert "pair" not in slots.UNIT_TABLES


def test_every_empty_toxoplasma_slot_says_why_it_is_empty():
    """A dash says a slot has no data. It must also say whether anyone has looked.

    Without this the slot table cannot tell a question nobody has searched from one that six
    acquisition passes exhausted, and those want opposite next actions. The practical cost of losing
    the distinction is that the next pass re-runs searches that already came back empty -- which
    happened twice during this campaign before the verdicts were written down.

    A new empty slot fails this until someone records what they searched. That is the point: the
    cheapest moment to write down where you looked is immediately after looking.
    """
    import csv
    import os
    path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "instructions", "done", "31_slots_toxoplasma.csv")
    if not os.path.exists(path):
        pytest.skip("slot table not generated")
    with open(path, encoding="utf8") as fh:
        rows = list(csv.DictReader(fh))
    undocumented = [r["slot"] for r in rows
                    if r["grade"] == "-" and not (r.get("blocked_by") or "").strip()]
    assert not undocumented, (
        f"empty slots with no verdict: {undocumented}. Add them to BLOCKED in "
        f"scripts/generate_slot_table.py with what you searched and what would fill them")


def test_an_empty_slot_verdict_distinguishes_missing_from_unreachable():
    """The two blocked states want different things -- an experiment, or a login."""
    import importlib.util
    import os
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    spec = importlib.util.spec_from_file_location(
        "gst", os.path.join(root, "scripts", "generate_slot_table.py"))
    gst = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gst)
    assert gst.BLOCKED, "no empty slot has been triaged"
    for slot, (verdict, searched, wanted) in gst.BLOCKED.items():
        assert verdict in ("missing", "unreachable"), f"{slot} has verdict {verdict!r}"
        assert len(searched) > 60, f"{slot} does not say where anyone looked"
        assert len(wanted) > 30, f"{slot} does not say what would fill it"


def test_a_filled_slot_carries_no_blocked_verdict():
    """Otherwise a stale verdict outlives the gap it described."""
    import csv
    import os
    path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                        "instructions", "done", "31_slots_toxoplasma.csv")
    if not os.path.exists(path):
        pytest.skip("slot table not generated")
    with open(path, encoding="utf8") as fh:
        stale = [r["slot"] for r in csv.DictReader(fh)
                 if r["grade"] != "-" and (r.get("blocked_by") or "").strip()]
    assert not stale, f"filled slots still carrying a blocked verdict: {stale}"


def test_a_plasmodium_pair_slot_is_not_filled_by_the_toxoplasma_graph():
    """Both arms name their layers the same, so a graph alone cannot say whose edges it holds.

    `orthogroup`, `domain` and `coexpression` are the same constructions in both species, which is
    deliberate -- it means a difference between the arms is biology rather than method. The cost is
    that a Plasmodium pair slot handed the Toxoplasma graph finds every layer it asked for. Three
    slots did exactly that the day the Plasmodium graph was built, and the count only moved because
    something else was being checked. The graph is identified by the table it arrives with.
    """
    import os
    import numpy as np
    import pandas as pd
    import starplast.paths as P
    nodes = pd.read_parquet(os.path.join(P.data_dir(), "nodes.parquet"))
    graph = np.load(os.path.join(P.data_dir(), "graph.npz"))
    pf_pairs = [s for s in slots.all_slots("Pf") if slots.edge_types(s)]
    assert pf_pairs, "no Plasmodium slot declares an edge type"
    filled = [s.name for s in pf_pairs if slots.is_filled(s, nodes, graph)]
    assert not filled, f"Plasmodium slots filled by Toxoplasma edges: {filled}"
    # And against its own table and graph they DO fill, or the guard is just refusing everything.
    pf_nodes_path = os.path.join(P.data_dir(), "pf_nodes.parquet")
    pf_graph_path = os.path.join(P.data_dir(), "pf_graph.npz")
    if os.path.exists(pf_nodes_path) and os.path.exists(pf_graph_path):
        pf_nodes = pd.read_parquet(pf_nodes_path)
        pf_graph = np.load(pf_graph_path)
        assert any(slots.is_filled(s, pf_nodes, pf_graph) for s in pf_pairs)


def test_an_edge_slot_with_no_table_to_check_against_is_still_answerable():
    """Callers that hold only a graph are not broken by the species guard."""
    import types
    pair = next(s for s in slots.all_slots("Tg") if slots.edge_types(s))
    wanted = slots.edge_types(pair)[0]
    graph = types.SimpleNamespace(files=[f"{wanted}__a", f"{wanted}__b"])
    assert slots.is_filled(pair, None, graph)
