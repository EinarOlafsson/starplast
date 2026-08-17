#!/usr/bin/env python3
"""The slot tree window: the catalog audit made permanent and visible.

The property that matters here is that the window cannot DISAGREE with the pipeline. It reads
`slots.relationship_tree`, `slots.declared_columns`, `slots.is_filled` and `slots.resolve` rather than
reimplementing any of them, because a second copy of the catalog logic in GUI code is a second thing
to keep in step and it would not be kept in step -- the window would then look authoritative while
being wrong, which is worse than having no window.

Everything here is headless. No test may need a display.
"""
from __future__ import annotations

import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import slots as S  # noqa: E402
from starplast import slot_tree as T  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NODES = os.path.join(ROOT, "starplast", "data", "nodes.parquet")


@pytest.fixture
def nodes():
    if not os.path.exists(NODES):
        pytest.skip("built node table not present")
    return pd.read_parquet(NODES)


# --------------------------------------------------------------------------- it opens with no data
def test_it_opens_with_no_table_at_all(qapp):
    """`nodes=None` is a supported state, not a degenerate one. This is the tool for inspecting an arm
    BEFORE its table exists, and a viewer that needs data to show a slot cannot show the arm that most
    needs looking at."""
    view = T.SlotTree(sources={"Tg": {}, "Pf": {}})
    assert view.tree.topLevelItemCount() > 0, "the catalog did not render without a table"
    assert "slots" in view.alarms.text()


def test_the_window_wraps_the_view(qapp):
    window = T.SlotTreeWindow(sources={"Tg": {}, "Pf": {}})
    assert isinstance(window.view, T.SlotTree)
    assert "slot tree" in window.windowTitle()


# --------------------------------------------------------------------------- it matches the catalog
def test_every_slot_appears_exactly_once_per_hierarchy(qapp):
    """The tree is the catalog re-addressed, not a subset of it. A slot missing from one hierarchy is
    a slot nobody browsing that hierarchy can find."""
    for organism in ("Tg", "Pf"):
        expected = {s.name for s in S.all_slots(organism)}
        for hierarchy in S.HIERARCHIES:
            view = T.SlotTree(sources={organism: {}}, organism=organism)
            view.hierarchy.setCurrentText(hierarchy)
            leaves = view._leaves(view.tree.invisibleRootItem())
            from PyQt6 import QtCore
            names = [i.data(0, QtCore.Qt.ItemDataRole.UserRole) for i in leaves]
            assert len(names) == len(set(names)), f"{organism}/{hierarchy} repeats a slot"
            assert set(names) == expected, f"{organism}/{hierarchy} is missing slots"


def test_a_group_row_counts_the_leaves_beneath_it(qapp):
    view = T.SlotTree(sources={"Tg": {}}, organism="Tg")
    root = view.tree.invisibleRootItem()
    for i in range(root.childCount()):
        group = root.child(i)
        leaves = view._leaves(group)
        assert group.text(1).endswith(f"/{len(leaves)}"), group.text(0)


# --------------------------------------------------------------------------- it cannot disagree
def test_the_empty_count_matches_the_generated_table(qapp):
    """The number this window exists to show, checked against the table the campaign is scored on.

    It caught its own bug on the way in: `slots.is_filled` keys its `tables` argument by UNIT, and
    keying it by file name instead made every metabolite and bridge slot read as empty -- the window
    reported seven empty Toxoplasma slots against the generated table's three.
    """
    csv = os.path.join(ROOT, "instructions", "done", "31_slots.csv")
    if not os.path.exists(csv):
        pytest.skip("the generated slot table is not present")
    table = pd.read_csv(csv)
    table["empty"] = table.grade.astype(str).eq("-")
    sources = T.shipped_sources()
    for organism in ("Tg", "Pf"):
        source = sources.get(organism) or {}
        if source.get("nodes") is None:
            continue
        got = T.audit(organism, source["nodes"], source["tables"], source["graph"])
        expected = int(table[table.organism.eq(organism)]["empty"].sum())
        assert got["empty"] == expected, f"{organism}: window {got['empty']}, table {expected}"


def test_coverage_shown_is_what_resolve_returns(nodes, qapp):
    """The window and the pipeline must not be able to disagree about how many genes a slot covers."""
    slot = next(s for s in S.all_slots("Tg")
                if s.unit == S.RESOLVABLE_UNIT and S.declared_columns(nodes, s))
    row = T._slot_row(slot, nodes, {}, None)
    got = S.resolve(nodes, slot)
    covered = got.values.notna()
    expected = int(covered.any(axis=1).sum()) if covered.ndim > 1 else int(covered.sum())
    assert row["genes"] == expected


# --------------------------------------------------------------------------- the four alarms
def test_a_deliberate_double_claim_turns_the_count_non_zero(qapp, monkeypatch):
    """Zero today, and the count is worth nothing unless it can be made non-zero. A column answering
    two slots means holding one out does not hold the other out."""
    frame = pd.DataFrame({"gene_id": ["TGME49_1"], "shared_column": [1.0]})
    one = S.Slot(organism="Tg", name="first", axis="a", context="c", unit="gene",
                 patterns=("shared_column",), policy="one")
    two = S.Slot(organism="Tg", name="second", axis="a", context="c", unit="gene",
                 patterns=("shared_column",), policy="one")
    monkeypatch.setattr(S, "all_slots", lambda organism=None: (one, two))
    got = T.audit("Tg", frame)
    assert got["double"] == 1
    assert got["orphan"] == 0


def test_a_column_no_slot_claims_is_counted(qapp, monkeypatch):
    frame = pd.DataFrame({"gene_id": ["TGME49_1"], "described": [1.0], "nobody_claims_this": [2.0]})
    slot = S.Slot(organism="Tg", name="only", axis="a", context="c", unit="gene",
                  patterns=("described",), policy="one")
    monkeypatch.setattr(S, "all_slots", lambda organism=None: (slot,))
    assert T.audit("Tg", frame)["orphan"] == 1


def test_a_candidate_with_no_accession_is_counted(qapp, monkeypatch):
    """A citation nobody can download is not a filled slot, and this count exists so it cannot read
    as one."""
    slot = S.Slot(organism="Pf", name="wanted", axis="a", context="c", unit="gene",
                  patterns=(), policy="one",
                  candidates=(("12345", "Nature", ""), ("999", "Cell", "PXD000001")))
    monkeypatch.setattr(S, "all_slots", lambda organism=None: (slot,))
    assert T.audit("Pf", None)["citation"] == 1


def test_the_alarm_line_turns_red_when_a_claim_is_doubled(qapp, monkeypatch):
    frame = pd.DataFrame({"gene_id": ["TGME49_1"], "shared_column": [1.0]})
    one = S.Slot(organism="Tg", name="first", axis="a", context="c", unit="gene",
                 patterns=("shared_column",), policy="one")
    two = S.Slot(organism="Tg", name="second", axis="a", context="c", unit="gene",
                 patterns=("shared_column",), policy="one")
    monkeypatch.setattr(S, "all_slots", lambda organism=None: (one, two))
    view = T.SlotTree(nodes=frame, sources={"Tg": {}}, organism="Tg")
    from starplast import theme as TH
    assert TH.palette_for("dark")["error"] in view.alarms.text()


# --------------------------------------------------------------------------- filtering and details
def test_the_filter_matches_names_columns_and_candidates(nodes, qapp):
    """A group row survives if anything under it matches, so filtering never hides a slot's address."""
    from PyQt6 import QtCore
    view = T.SlotTree(organism="Tg")
    before = len(view._leaves(view.tree.invisibleRootItem()))
    view.filter.setText("bradyzoite")
    after = view._leaves(view.tree.invisibleRootItem())
    assert 0 < len(after) < before
    for item in after:
        name = item.data(0, QtCore.Qt.ItemDataRole.UserRole)
        slot = next(s for s in S.all_slots("Tg") if s.name == name)
        assert view._matches(slot, "bradyzoite"), f"{name} survived a filter it does not match"
    # And a filter nothing answers empties the tree rather than silently showing everything.
    view.filter.setText("no slot mentions this string")
    assert view._leaves(view.tree.invisibleRootItem()) == []


def test_selecting_a_slot_shows_its_address_in_all_three_hierarchies(nodes, qapp):
    from PyQt6 import QtCore
    view = T.SlotTree(organism="Tg")
    leaf = view._leaves(view.tree.invisibleRootItem())[0]
    view.tree.setCurrentItem(leaf)
    html = view.details.toHtml()
    for hierarchy in S.HIERARCHIES:
        assert hierarchy in html, f"{hierarchy} address missing from the details pane"
    assert "axis" in html and "policy" in html


def test_selecting_a_group_says_so_rather_than_showing_nothing(qapp):
    view = T.SlotTree(sources={"Tg": {}}, organism="Tg")
    group = view.tree.topLevelItem(0)
    view.tree.setCurrentItem(group)
    assert "group" in view.details.toHtml().lower()
    assert view._show_details(None) == ""


def test_an_empty_slot_is_coloured_rather_than_left_to_be_noticed(qapp):
    """The empty rows are the point of the window."""
    from PyQt6 import QtCore, QtGui
    from starplast import theme as TH
    view = T.SlotTree(sources={"Pf": {}}, organism="Pf")     # no table: everything reads empty
    leaf = view._leaves(view.tree.invisibleRootItem())[0]
    assert leaf.text(1) == "—"
    assert leaf.foreground(0).color() == QtGui.QColor(TH.palette_for("dark")["warning"])


# --------------------------------------------------------------------------- candidate shapes
def test_a_candidate_is_read_whether_it_is_a_tuple_or_a_dict():
    """Candidates are tuples in the catalog and dicts in the generated table. A viewer that insisted
    on one would show an empty column and look like the data was missing."""
    assert T._accession_of(("1", "Nature", "PXD1")) == "PXD1"
    assert T._accession_of({"accession": "PXD2"}) == "PXD2"
    assert T._accession_of(("1", "Nature")) == ""
    assert T._accession_of("something else") == ""
    assert "Nature" in T._candidate_text(("1", "Nature", "PXD1"))
    assert "Nature" in T._candidate_text({"pmid": "1", "source": "Nature", "title": "t"})
    assert T._candidate_text("plain") == "plain"


def test_shipped_sources_reports_what_is_actually_on_disk():
    got = T.shipped_sources()
    assert set(got) == {"Tg", "Pf"}
    for organism, source in got.items():
        assert set(source) == {"nodes", "graph", "tables"}
        if source["nodes"] is not None:
            assert len(source["nodes"]) > 100, organism
