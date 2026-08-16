#!/usr/bin/env python3
"""Biological-question slots and the policies that combine datasets which answer them.

The generated catalog is package data: the application can use it offline, while
``scripts/generate_slot_table.py`` remains the editable source. A slot is not a filename pattern;
it is a question plus a condition, its candidate measurements, and an explicit combination policy.
"""
from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass

import numpy as np
import pandas as pd

POLICIES = ("one", "average", "fill", "separate")

#: What a slot's rows ARE, which decides which table it can be resolved against. Not decoration: a
#: parasite node table has one row per parasite gene, so resolving a slot of any other unit against
#: it is a category error that would otherwise fail quietly -- a host proteome resolved against
#: Toxoplasma genes matches nothing and returns an empty frame, which reads exactly like a dataset
#: nobody has downloaded yet.
#:
#: gene           -- one parasite gene. The only unit `resolve` accepts, because the node table is
#:                   the only table this application currently holds.
#: host_gene      -- one HOST gene. Lives in a host table keyed by host identifiers; reaches the
#:                   parasite map only through a pair slot, never as a feature column.
#: pair           -- one (parasite gene, host gene) or (parasite gene, parasite gene) observation.
#:                   The bridge between two tables, and the only thing that may cross.
#: ortholog_group -- one orthology group across parasite species. The bridge that carries a
#:                   measurement from one species to another, and the one that must never be
#:                   mistaken for a measurement in the receiving species.
UNITS = ("gene", "host_gene", "pair", "ortholog_group", "metabolite")

#: The unit whose rows are the node table's rows.
RESOLVABLE_UNIT = "gene"
CATALOG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "slots.json")


def _key(text: str) -> str:
    """A stable recipe-safe identifier for a displayed slot name."""
    return re.sub(r"[^a-z0-9]+", "_", str(text).lower()).strip("_")


@dataclass(frozen=True)
class Slot:
    """One biological question and the measurements that can answer it."""
    organism: str
    name: str
    axis: str
    context: str
    unit: str
    patterns: tuple
    policy: str
    candidates: tuple = ()
    # Three independent classifications.  They are deliberately not collapsed into one tree:
    # a membrane-topology slot is sequence evidence about localisation, while an RNA-seq slot can
    # simultaneously concern a life stage and a perturbation.  Keeping the facets separate is what
    # makes both feature selection and leakage exclusion explicit.
    evidence_path: tuple = ()
    biology_path: tuple = ()
    context_path: tuple = ()
    target_family: str = ""
    target_columns: tuple = ()
    role: str = "feature"       # feature | target | metadata | never

    @property
    def key(self) -> str:
        """The identifier stored in embedding recipes."""
        return f"{self.organism}_{_key(self.name)}"


@dataclass
class SlotResult:
    """Resolved values plus per-gene provenance for one slot."""
    values: pd.DataFrame
    source: pd.Series
    source_columns: tuple


def all_slots(organism: str | None = None) -> tuple:
    """The generated catalog, optionally restricted to ``Toxo`` or ``Pf``."""
    try:
        with open(CATALOG, encoding="utf8") as fh:
            records = json.load(fh)
    except (OSError, ValueError):
        return ()
    out = tuple(Slot(organism=str(row["organism"]), name=str(row["name"]),
                     axis=str(row["axis"]), context=str(row["context"]), unit=str(row["unit"]),
                     patterns=tuple(row.get("patterns", ())), policy=str(row["policy"]),
                     candidates=tuple(row.get("candidates", ())),
                     evidence_path=tuple(row.get("evidence_path", ())),
                     biology_path=tuple(row.get("biology_path", ())),
                     context_path=tuple(row.get("context_path", ())),
                     target_family=str(row.get("target_family", "")),
                     target_columns=tuple(row.get("target_columns", ())),
                     role=str(row.get("role", "feature"))) for row in records)
    return tuple(slot for slot in out if organism is None or slot.organism == organism)


def by_key(key: str) -> Slot | None:
    """Resolve a recipe key to its slot."""
    return next((slot for slot in all_slots() if slot.key == key), None)


def source_columns(nodes: pd.DataFrame, slot: Slot) -> tuple:
    """Raw numeric node columns matching a slot, in candidate order without duplicates."""
    out = []
    for pattern in slot.patterns:
        matches = (pattern,) if pattern in nodes.columns else tuple(
            column for column in nodes.columns if str(column).startswith(pattern))
        for column in matches:
            if column not in out and pd.api.types.is_numeric_dtype(nodes[column]):
                out.append(column)
    return tuple(out)


#: A pattern of this form names an EDGE type in the graph rather than a column in the node table.
EDGE_PREFIX = "edge:"


def edge_types(slot: Slot) -> tuple:
    """The graph edge types a slot is filled by, if any."""
    return tuple(str(p)[len(EDGE_PREFIX):] for p in (slot.patterns or ())
                 if str(p).startswith(EDGE_PREFIX))


def is_filled(slot: Slot, nodes: pd.DataFrame = None, graph=None) -> bool:
    """Does this slot have data behind it -- in the node table, OR in the graph?

    Both, and that is the whole reason this is a function rather than a line of whatever script
    happens to be asking. A slot measured per PAIR -- co-expression, co-fitness, shared compartment,
    crosslink MS -- is filled by an edge type in `graph.npz` and has no node column at all. Asking
    `declared_columns` alone reports every one of them as empty, which is exactly what happened:
    eight Toxoplasma slots were counted as unfilled and chased with downloads while their data was
    already in the graph, and the coverage figure went out twice before anyone noticed.

    `graph` is anything with a `files` list, which is what `numpy.load` gives back for an npz.
    """
    wanted = edge_types(slot)
    if wanted:
        if graph is None:
            return False
        present = {str(name).split("__")[0] for name in getattr(graph, "files", ())}
        return all(edge in present for edge in wanted)
    return bool(nodes is not None and declared_columns(nodes, slot))


def coverage(organism: str = "Tg", nodes: pd.DataFrame = None, graph=None) -> dict:
    """How many of an arm's feature slots have data, and which do not.

    The single place this question is answered, so a viewer, a report and a census cannot disagree
    about it -- which they did, before this existed.
    """
    feature = [s for s in all_slots(organism) if s.role == "feature"]
    full = [s for s in feature if is_filled(s, nodes, graph)]
    return {"organism": organism, "n_slots": len(feature), "filled": len(full),
            "empty": tuple(s.name for s in feature if not is_filled(s, nodes, graph))}


def declared_columns(nodes: pd.DataFrame, slot: Slot, numeric_only: bool = False) -> tuple:
    """Every column declared by a slot, including categorical targets and metadata.

    ``source_columns`` is intentionally numeric because it builds matrices.  Leakage closure needs
    the wider declaration: labels, confidence fields and numeric features all travel together.
    """
    out = []
    for pattern in (*slot.patterns, *slot.target_columns):
        matches = (pattern,) if pattern in nodes.columns else tuple(
            column for column in nodes.columns if str(column).startswith(pattern))
        for column in matches:
            if column not in out and (not numeric_only
                                      or pd.api.types.is_numeric_dtype(nodes[column])):
                out.append(column)
    return tuple(out)


HIERARCHIES = ("evidence", "biology", "context")


def hierarchy_path(slot: Slot, hierarchy: str) -> tuple:
    """Return one facet's path from broad group to the slot leaf."""
    if hierarchy not in HIERARCHIES:
        raise ValueError(f"unknown hierarchy {hierarchy!r}; choose from {HIERARCHIES}")
    return tuple(getattr(slot, f"{hierarchy}_path")) + (slot.key,)


def relationship_tree(organism: str = "Tg", hierarchy: str = "evidence") -> dict:
    """A nested, JSON-safe tree whose leaves are slot keys.

    There are three trees rather than a single misleading tree.  A slot is a leaf in every facet;
    callers can therefore omit an assay class, a biological subject, or an experimental context
    without pretending those are the same grouping.
    """
    tree = {}
    for slot in all_slots(organism):
        branch = tree
        for part in hierarchy_path(slot, hierarchy):
            branch = branch.setdefault(part, {})
    return tree


def slots_in_group(path: tuple | list | str, organism: str = "Tg",
                   hierarchy: str = "evidence") -> tuple:
    """Slots below a hierarchy path, useful for selecting or omitting a whole class."""
    wanted = (path,) if isinstance(path, str) else tuple(path)
    return tuple(slot for slot in all_slots(organism)
                 if hierarchy_path(slot, hierarchy)[:len(wanted)] == wanted)


def target_slots(target: str, organism: str = "Tg") -> tuple:
    """Slots that explicitly declare ``target`` as one of their outcomes."""
    return tuple(slot for slot in all_slots(organism) if target in slot.target_columns)


def family_slots(family: str, organism: str = "Tg") -> tuple:
    """All slots estimating the same target quantity, across assays or transfers."""
    return tuple(slot for slot in all_slots(organism) if family and slot.target_family == family)


def _groups(nodes: pd.DataFrame, slot: Slot) -> list:
    """Candidate pattern -> its matching numeric columns."""
    groups = []
    used = set()
    for pattern in slot.patterns:
        matches = (pattern,) if pattern in nodes.columns else tuple(
            column for column in nodes.columns if str(column).startswith(pattern))
        columns = [column for column in matches if column not in used
                   and pd.api.types.is_numeric_dtype(nodes[column])]
        if columns:
            groups.append((pattern, columns))
            used.update(columns)
    return groups


def _rank(frame: pd.DataFrame) -> pd.DataFrame:
    """Column-wise percentile ranks centered on zero, retaining missing values."""
    return frame.rank(method="average", pct=True, na_option="keep") - 0.5


def resolve(nodes: pd.DataFrame, slot: Slot | str, chosen: str | None = None) -> SlotResult:
    """Apply a slot's ``one``, ``average``, ``fill`` or ``separate`` policy.

    ``source`` is one value per gene. For ``fill`` it is the exact candidate supplying that gene;
    for policies that deliberately combine measurements it is a ``+``-joined record of every
    contributing candidate. This is provenance, not a feature, and is never placed in the matrix.
    """
    slot = by_key(slot) if isinstance(slot, str) else slot
    if slot is None:
        raise KeyError("unknown slot")
    if slot.policy not in POLICIES:
        raise ValueError(f"unknown slot policy {slot.policy!r}")
    if slot.unit != RESOLVABLE_UNIT:
        # Refused rather than returned empty. A host-gene or pair slot resolved against a table of
        # parasite genes matches nothing, and "matches nothing" is indistinguishable from "nobody
        # has downloaded this yet" -- so the mistake would be invisible in exactly the place the
        # slot table exists to make visible.
        raise ValueError(
            f"slot {slot.key!r} is measured per {slot.unit!r} and cannot be resolved against a "
            f"table of {RESOLVABLE_UNIT!r} rows; it reaches the map through a bridge slot")
    groups = _groups(nodes, slot)
    empty_source = pd.Series("", index=nodes.index, dtype=object, name=f"{slot.key}_source")
    if not groups:
        return SlotResult(pd.DataFrame(index=nodes.index), empty_source, ())
    raw_columns = tuple(column for _pattern, columns in groups for column in columns)

    if slot.policy == "one":
        selected = next(((pattern, columns) for pattern, columns in groups
                         if chosen in (pattern, *columns)), groups[0])
        pattern, columns = selected
        values = nodes[columns].copy()
        present = values.notna().any(axis=1)
        source = empty_source.mask(present, pattern)
        return SlotResult(values, source, tuple(columns))

    collapsed = [(pattern, nodes[columns].mean(axis=1, skipna=True), columns)
                 for pattern, columns in groups]
    if slot.policy == "average":
        ranked = pd.concat([_rank(nodes[columns]).mean(axis=1, skipna=True)
                            for _pattern, _series, columns in collapsed], axis=1)
        values = ranked.mean(axis=1, skipna=True).to_frame(f"{slot.key}__average")
        sources = pd.DataFrame({pattern: nodes[columns].notna().any(axis=1)
                                for pattern, _series, columns in collapsed})
        source = sources.apply(lambda row: "+".join(row.index[row.to_numpy(bool)]), axis=1)
        source.name = empty_source.name
        return SlotResult(values, source, raw_columns)

    if slot.policy == "fill":
        ordered = sorted(collapsed, key=lambda item: int(item[1].notna().sum()), reverse=True)
        value = pd.Series(np.nan, index=nodes.index, dtype=float)
        source = empty_source.copy()
        for pattern, candidate, _columns in ordered:
            take = value.isna() & candidate.notna()
            # pandas 3 is intentionally strict about assigning nullable extension arrays into a
            # NumPy float Series.  Slots are numeric by construction, so cross that boundary
            # explicitly instead of relying on pandas' version-dependent coercion.
            value.loc[take] = candidate.loc[take].to_numpy(dtype=float, na_value=np.nan)
            source.loc[take] = pattern
        return SlotResult(value.to_frame(f"{slot.key}__fill"), source, raw_columns)

    values = nodes[list(raw_columns)].copy()
    sources = pd.DataFrame({pattern: nodes[columns].notna().any(axis=1)
                            for pattern, _series, columns in collapsed})
    source = sources.apply(lambda row: "+".join(row.index[row.to_numpy(bool)]), axis=1)
    source.name = empty_source.name
    return SlotResult(values, source, raw_columns)
