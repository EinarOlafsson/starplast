#!/usr/bin/env python3
"""A window showing every slot, its place in the tree, and what fills it.

## Why this exists

The slot catalog is 222 entries across three hierarchies and two organisms. Every check run on it so
far was a Python one-liner written for the occasion and thrown away -- which is how the first audit
reported a leak that did not exist, because the one-liner reimplemented pattern matching instead of
calling `slots.declared_columns`.

**A sanity check nobody can see is a sanity check nobody performs.** This is the audit made permanent
and visible: open it and the shape of the whole catalog is in front of you, including the parts that
are empty -- which are the most informative parts, and the reason the empty rows are coloured rather
than hidden.

## The rule this module follows

Nothing here reimplements the catalog. The tree comes from `slots.relationship_tree`, the columns
from `slots.declared_columns`, the coverage from `slots.resolve`, and whether a slot is filled from
`slots.is_filled`. A second copy of that logic in GUI code is a second thing to keep in step, and it
would not be kept in step -- the window would then disagree with the pipeline while looking
authoritative, which is worse than having no window.

## It opens with no data

`nodes=None` is a supported state, not a degenerate one. This is the tool for inspecting an arm
*before* its table exists, and a viewer that needs data to show a slot cannot show the arm that most
needs looking at. With no table, coverage is blank and the catalog is still fully browsable.
"""
from __future__ import annotations

import os

import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets

from . import slots as S
from . import theme as TH

#: What the four alarms mean, shown as tooltips because this window is where somebody meets these
#: concepts for the first time. Each is an audit that used to be run by hand, if it was run at all.
ALARMS = {
    "empty": ("slots with no data",
              "A slot nothing fills. This is the map of what has not been measured yet, and it is "
              "the most useful number in this window -- an empty slot is a question with a known "
              "shape, which is what makes it findable."),
    "double": ("columns claimed twice",
               "A column matched by two slots' patterns. Zero today, and it must stay zero: a column "
               "answering two questions means holding one out does not hold the other out, which is "
               "how a held-out feature leaks back into the map that was built on it."),
    "orphan": ("columns claimed by nobody",
               "A column in the table that no slot's patterns match. Zero today. A column nothing "
               "claims is a measurement the catalog cannot describe, so it can never be held out, "
               "audited or reported."),
    "citation": ("candidates with no accession",
                 "A candidate naming a paper but no downloadable accession. A citation nobody can "
                 "fetch is not a filled slot, and this count exists so that it cannot read as one."),
}

#: The tree's columns. The heading tooltips carry the explanations, per the project rule.
HEADINGS = [
    ("slot", "The question, at its address in this hierarchy. Group rows aggregate what is beneath "
             "them; leaves are slots."),
    ("grade", "Evidence grade of the filling data: A measured in this organism, B measured nearby "
              "and transferred, C predicted or derived. A dash means nothing fills it."),
    ("genes", "How many genes carry a value, computed live from the loaded table by `slots.resolve` "
              "-- not read from a generated CSV, which could be stale."),
    ("coverage", "That count as a share of the table."),
    ("what fills it", "The columns this slot's patterns match in the loaded table, or the candidate "
                      "sources when nothing does."),
]


def _slot_row(slot: S.Slot, nodes, tables, graph=None) -> dict:
    """Everything the tree shows for one slot, computed through the catalog's own functions."""
    filled = S.is_filled(slot, nodes=nodes, graph=graph, tables=tables)
    columns = ()
    genes = None
    if nodes is not None and slot.unit == S.RESOLVABLE_UNIT:
        columns = S.declared_columns(nodes, slot)
        if columns:
            got = S.resolve(nodes, slot)
            values = getattr(got, "values", None)
            if values is not None:
                # A `separate`-policy slot resolves to a FRAME, one column per source kept apart; the
                # others resolve to a Series. Counting genes with any value covers both, and calling
                # int() on the frame's row-wise sum is what raised "cannot convert the series".
                covered = values.notna()
                genes = int(covered.any(axis=1).sum()) if covered.ndim > 1 else int(covered.sum())
    return {"slot": slot, "filled": filled, "columns": columns, "genes": genes,
            "candidates": tuple(slot.candidates or ())}


def shipped_sources() -> dict:
    """Every table and graph the catalog is graded against, keyed by organism.

    Loaded here rather than taken from the window, because this viewer must be able to show an arm
    the open window is not displaying -- that is most of its purpose. Missing files are simply absent
    from the result: an arm that has not been built shows its catalog with no coverage, which is the
    honest picture and the one instruction 40 asks for.
    """
    import numpy as np
    from . import paths
    out = {}
    for organism, (table, graph) in {"Tg": ("nodes.parquet", "graph.npz"),
                                     "Pf": ("pf_nodes.parquet", "pf_graph.npz")}.items():
        nodes = graph_file = None
        path = paths.cache_file(table)
        if os.path.exists(path):
            nodes = pd.read_parquet(path)
        path = paths.cache_file(graph)
        if os.path.exists(path):
            graph_file = np.load(path, allow_pickle=True)
        out[organism] = {"nodes": nodes, "graph": graph_file,
                         "tables": _side_tables(organism)}
    return out


def _side_tables(organism: str) -> dict:
    """The tables a slot of a non-gene unit resolves against, keyed the way `slots.is_filled` reads
    them: by UNIT, plus `"bridge"` for the one bridge table belonging to this arm.

    Keying these by file name instead -- which is the obvious mistake, and the one made here first --
    makes every metabolite and bridge slot read as EMPTY while its data sits on disk. It showed up as
    the window reporting seven empty Toxoplasma slots against the generated table's three, which is
    precisely the disagreement this window exists to catch. It caught its own.
    """
    from . import paths
    out = {}
    metabolites = paths.cache_file(S.UNIT_TABLES.get("metabolite", "metabolites.parquet"))
    if os.path.exists(metabolites):
        out["metabolite"] = pd.read_parquet(metabolites)
    host = paths.cache_file(S.UNIT_TABLES.get("host_gene", "host_proteins.parquet"))
    if os.path.exists(host):
        out["host_gene"] = pd.read_parquet(host)
    # One bridge table per arm, and they must not be swapped: both key their rows `host`, so only the
    # parasite accessions say whose contacts these are. `is_filled` checks that with `same_species`.
    bridge = S.SPECIES_BRIDGE_TABLES.get(organism, {}).get("host")
    if bridge and os.path.exists(paths.cache_file(bridge)):
        out["bridge"] = pd.read_parquet(paths.cache_file(bridge))
    return out


def audit(organism: str, nodes=None, tables=None, graph=None, rows: dict = None) -> dict:
    """The four counts this window exists to make impossible to miss.

    Computed with the same functions the pipeline uses, so the window cannot report a clean catalog
    while the build sees a broken one.
    """
    catalog = S.all_slots(organism)
    # `rows` is the map the window already built. Given it, this walks it instead of asking the
    # catalog the same 119 questions a second time.
    empty = (sum(1 for r in rows.values() if not r["filled"]) if rows else
             sum(1 for s in catalog
                 if not S.is_filled(s, nodes=nodes, graph=graph, tables=tables)))
    citation = sum(1 for s in catalog for c in (s.candidates or ())
                   if not _accession_of(c))
    double = orphan = 0
    if nodes is not None:
        claimed = {}
        for slot in catalog:
            if slot.unit != S.RESOLVABLE_UNIT:
                continue
            columns = rows[slot.name]["columns"] if rows and slot.name in rows else \
                S.declared_columns(nodes, slot)
            for column in columns:
                claimed.setdefault(column, []).append(slot.name)
        double = sum(1 for names in claimed.values() if len(names) > 1)
        # A column claimed by a `role="never"` slot is described -- bookkeeping the catalog knows
        # about -- so it is not an orphan. `gene_id` is the key, not a measurement.
        described = set(claimed) | {"gene_id"}
        orphan = sum(1 for c in nodes.columns if c not in described)
    # The host table gets the same question, because for a long time nothing asked it: the alarm
    # walked `nodes` and only `nodes`, and `pv_enrichment_log2` shipped unclaimed by any slot the
    # whole time. Instruction 48.
    #
    # Both arms' host slots are consulted rather than this organism's. A host column belongs to a
    # TISSUE, not to a parasite -- `rbc_*` is claimed only by a Plasmodium slot and `bmdm_*` only by
    # a Toxoplasma one -- so walking one arm's catalogue would report the other arm's columns as
    # orphans in every window.
    host = (tables or {}).get("host_gene")
    if host is not None and len(host):
        host_claimed = set()
        for slot in S.all_slots():
            if slot.unit == "host_gene":
                host_claimed.update(S.declared_columns(host, slot))
        # `host_id` is the key and `host_name` its label, so neither is a measurement to describe.
        host_described = host_claimed | {"host_id", "host_name"}
        orphan += sum(1 for c in host.columns if c not in host_described)
    return {"empty": empty, "double": double, "orphan": orphan, "citation": citation,
            "n_slots": len(catalog)}


def _accession_of(candidate) -> str:
    """The accession a candidate carries, whatever shape the candidate is.

    Candidates are tuples in the catalog and dicts in the generated table, and this window reads
    both rather than insisting on one -- the alternative is a viewer that shows an empty column and
    looks like the data is missing.
    """
    if isinstance(candidate, dict):
        return str(candidate.get("accession") or "")
    if isinstance(candidate, (tuple, list)):
        return str(candidate[2]) if len(candidate) > 2 else ""
    return ""


def _candidate_text(candidate) -> str:
    """A candidate as one readable line: accession, source, title."""
    if isinstance(candidate, dict):
        parts = [candidate.get("pmid"), candidate.get("source"), candidate.get("title")]
    elif isinstance(candidate, (tuple, list)):
        parts = list(candidate)
    else:
        parts = [candidate]
    return " · ".join(str(p) for p in parts if p)


class SlotTree(QtWidgets.QWidget):
    """The tree, its filters, and the details pane, as one widget so a test can build it headlessly.

    Split from the window for the same reason `build_preferences` is split from `open_preferences`:
    a window that has to be shown to be inspected cannot be tested without a display.
    """

    def __init__(self, nodes: pd.DataFrame = None, tables: dict = None, parent=None,
                 organism: str = None, graph=None, sources: dict = None):
        super().__init__(parent)
        # Given data, use it. Given none, read what is on disk -- this window must be able to show an
        # arm the open map is not displaying, which is most of what it is for.
        self._sources = sources if sources is not None else shipped_sources()
        self.nodes = nodes
        self.graph = graph
        self.tables = tables or {}
        self._explicit = nodes is not None
        self._row_cache = {}
        layout = QtWidgets.QVBoxLayout(self)

        controls = QtWidgets.QHBoxLayout()
        self.hierarchy = QtWidgets.QComboBox()
        self.hierarchy.addItems(list(S.HIERARCHIES))
        self.hierarchy.setToolTip(
            "The same slots, at three different addresses. `evidence` groups by how a thing was "
            "measured, `biology` by what it is about, `context` by where and when. A slot appears "
            "in all three, which is the point of having three: the question 'what do we know about "
            "the bradyzoite' and 'what has been measured by mass spectrometry' are both real.")
        self.hierarchy.currentTextChanged.connect(self.rebuild)
        self.organism = QtWidgets.QComboBox()
        self.organism.addItems(sorted({s.organism for s in S.all_slots()}))
        if organism:
            self.organism.setCurrentText(organism)
        self.organism.setToolTip("Which arm's catalog to show. The two are never merged.")
        self.organism.currentTextChanged.connect(self.rebuild)
        self.filter = QtWidgets.QLineEdit(placeholderText="filter slots, columns, candidates…")
        self.filter.setToolTip(
            "Matches the slot name, the columns it claims and its candidate titles. A group row "
            "survives if anything under it matches, so filtering never hides a slot's address.")
        self.filter.textChanged.connect(self.rebuild)
        controls.addWidget(QtWidgets.QLabel("hierarchy"))
        controls.addWidget(self.hierarchy)
        controls.addWidget(QtWidgets.QLabel("organism"))
        controls.addWidget(self.organism)
        controls.addWidget(self.filter, 1)
        layout.addLayout(controls)

        self.alarms = QtWidgets.QLabel("")
        self.alarms.setWordWrap(True)
        layout.addWidget(self.alarms)

        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        self.tree = QtWidgets.QTreeWidget()
        self.tree.setColumnCount(len(HEADINGS))
        self.tree.setHeaderLabels([h for h, _ in HEADINGS])
        for i, (_, text) in enumerate(HEADINGS):
            self.tree.headerItem().setToolTip(i, TH.tip(text))
        self.tree.currentItemChanged.connect(self._show_details)
        split.addWidget(self.tree)
        self.details = QtWidgets.QTextBrowser()
        self.details.setOpenExternalLinks(True)
        split.addWidget(self.details)
        split.setSizes([620, 380])
        layout.addWidget(split, 1)
        self.rebuild()

    # ------------------------------------------------------------------ building
    def rebuild(self, *_):
        """Redraw the tree for the current hierarchy, organism and filter."""
        organism = self.organism.currentText()
        if not self._explicit:
            source = self._sources.get(organism) or {}
            self.nodes, self.graph = source.get("nodes"), source.get("graph")
            self.tables = source.get("tables") or {}
        hierarchy = self.hierarchy.currentText()
        catalog = S.all_slots(organism)
        # Mapped once per organism, not once per keystroke. A row is `is_filled` plus
        # `declared_columns` plus `resolve` for one slot -- pattern matching over 393 columns and a
        # notna() count -- and none of it depends on which hierarchy is showing or what is typed in
        # the filter. Recomputing all 119 on every rebuild made a rebuild take a second, and the
        # filter rebuilds on every keystroke, so typing "bradyzoite" cost about ten.
        cached = self._row_cache.get(organism)
        if cached is None:
            cached = {s.name: _slot_row(s, self.nodes, self.tables, self.graph) for s in catalog}
            self._row_cache[organism] = cached
        self._rows = cached
        # The key a tree node carries is the LAST element of the slot's own hierarchy path. Rebuilding
        # it from the name instead -- replacing " · " and spaces with underscores -- looked equivalent
        # and was not: 53 of 119 slots failed to match and rendered as empty GROUP rows, which is a
        # window that quietly shows two thirds of the catalog as scaffolding. Ask the catalog.
        self._by_key = {S.hierarchy_path(s, hierarchy)[-1]: s for s in catalog}
        tree = S.relationship_tree(organism, hierarchy)
        self.tree.clear()
        needle = self.filter.text().strip().lower()
        for name, child in tree.items():
            item = self._add(self.tree.invisibleRootItem(), name, child, needle)
            if item is None:
                continue
        self.tree.expandToDepth(0)
        self._show_alarms(organism)
        return self.tree.topLevelItemCount()

    def _add(self, parent, name, child, needle):
        """Add one node, returning it, or None when the filter excludes everything beneath it."""
        slot = self._slot_for(name)
        if slot is not None:
            if needle and not self._matches(slot, needle):
                return None
            return self._add_leaf(parent, slot)
        item = QtWidgets.QTreeWidgetItem([name, "", "", "", ""])
        kept = [self._add(item, key, value, needle) for key, value in (child or {}).items()]
        kept = [k for k in kept if k is not None]
        if not kept and needle:
            return None
        # Accumulated from the children as they are made, rather than re-walking the whole subtree
        # once per group. `_leaves` from every group is quadratic in the catalog, and the catalog is
        # the thing that keeps growing.
        leaves = filled = best = 0
        for child in kept:
            n, f, g = child.data(0, QtCore.Qt.ItemDataRole.UserRole + 3) or (1, 0, 0)
            leaves += n
            filled += f
            best = max(best, g)
        item.setData(0, QtCore.Qt.ItemDataRole.UserRole + 3, (leaves, filled, best))
        item.setText(1, f"{filled}/{leaves}")
        item.setText(4, f"{leaves - filled} empty" if leaves else "")
        item.setText(2, f"max {best:,}" if best else "")
        parent.addChild(item)
        return item

    def _add_leaf(self, parent, slot):
        """One slot as a row, coloured by whether anything fills it."""
        row = self._rows[slot.name]
        genes = row["genes"]
        share = ""
        if genes is not None and self.nodes is not None and len(self.nodes):
            share = f"{100 * genes / len(self.nodes):.1f}%"
        fills = (", ".join(row["columns"][:4]) if row["columns"]
                 else (f"{len(row['candidates'])} candidate(s)" if row["candidates"]
                       else "no candidate"))
        item = QtWidgets.QTreeWidgetItem(
            [f"{slot.organism}_{slot.name}", "A" if row["filled"] else "—",
             f"{genes:,}" if genes else "", share, fills])
        item.setData(0, QtCore.Qt.ItemDataRole.UserRole, slot.name)
        item.setData(0, QtCore.Qt.ItemDataRole.UserRole + 1, row["filled"])
        item.setData(0, QtCore.Qt.ItemDataRole.UserRole + 2, genes or 0)
        # (leaves, filled, best genes) so a parent can add its children up without walking them.
        item.setData(0, QtCore.Qt.ItemDataRole.UserRole + 3,
                     (1, 1 if row["filled"] else 0, genes or 0))
        if not row["filled"]:
            # The empty rows are the point of the window, so they are marked rather than left to be
            # noticed. Colour from the palette: a hex here would stay dark on the light theme.
            colour = QtGui.QColor(TH.palette_for("dark")["warning"])
            for column in range(len(HEADINGS)):
                item.setForeground(column, colour)
        parent.addChild(item)
        return item

    def _slot_for(self, key):
        """The slot a tree key names, or None when the key is a group row."""
        return self._by_key.get(key)

    @staticmethod
    def _leaves(item) -> list:
        """Every leaf beneath a node, so a group row can aggregate what is under it."""
        out = []
        stack = [item]
        while stack:
            node = stack.pop()
            for i in range(node.childCount()):
                child = node.child(i)
                (stack if child.childCount() else out).append(child)
        return out

    def _matches(self, slot, needle: str) -> bool:
        """Whether a slot answers the free-text filter, over its name, columns and candidates."""
        row = self._rows[slot.name]
        hay = " ".join([slot.name, slot.axis or "", slot.context or "",
                        " ".join(row["columns"]),
                        " ".join(_candidate_text(c) for c in row["candidates"])]).lower()
        return needle in hay

    # ------------------------------------------------------------------ reporting
    def _show_alarms(self, organism: str) -> str:
        """The four counts, in words, above the tree."""
        got = audit(organism, self.nodes, self.tables, self.graph, rows=self._rows)
        parts = []
        for key, (label, _tip) in ALARMS.items():
            parts.append(f"<b>{got[key]}</b> {label}")
        text = f"{got['n_slots']} slots · " + " · ".join(parts)
        bad = got["double"] or got["orphan"]
        colour = TH.palette_for("dark")["error" if bad else "fg_muted"]
        self.alarms.setText(f"<span style='color:{colour}'>{text}</span>")
        self.alarms.setToolTip("\n\n".join(f"{label}: {tip}" for label, tip in ALARMS.values()))
        return text

    def _show_details(self, item, _previous=None) -> str:
        """Everything known about the selected slot, including its address in all three hierarchies."""
        if item is None:
            return ""
        name = item.data(0, QtCore.Qt.ItemDataRole.UserRole)
        if not name:
            self.details.setHtml("<p>A group. Select a slot to see what fills it.</p>")
            return ""
        slot = self._by_key.get(f"{self.organism.currentText()}_{name}") or \
            next((s for s in S.all_slots(self.organism.currentText()) if s.name == name), None)
        row = self._rows[name]
        lines = [f"<h3>{slot.organism}_{slot.name}</h3>",
                 f"<p><b>axis</b> {slot.axis} &middot; <b>context</b> {slot.context}<br>"
                 f"<b>unit</b> {slot.unit} &middot; <b>policy</b> {slot.policy} &middot; "
                 f"<b>role</b> {getattr(slot, 'role', '')}<br>"
                 f"<b>target family</b> {getattr(slot, 'target_family', '') or '—'}</p>"]
        lines.append("<p><b>Address</b><br>" + "<br>".join(
            f"{h}: {' &rarr; '.join(S.hierarchy_path(slot, h))}" for h in S.HIERARCHIES) + "</p>")
        if row["columns"]:
            lines.append(f"<p><b>Columns in the loaded table</b> ({len(row['columns'])})<br>"
                         + "<br>".join(row["columns"]) + "</p>")
        else:
            lines.append("<p><b>Columns in the loaded table</b><br>none</p>")
        if row["candidates"]:
            lines.append("<p><b>Candidates</b><br>" + "<br>".join(
                f"{_candidate_text(c)}"
                f"{'' if _accession_of(c) else '  <i>(no accession — cannot be downloaded)</i>'}"
                for c in row["candidates"]) + "</p>")
        html = "".join(lines)
        self.details.setHtml(html)
        return html


class SlotTreeWindow(QtWidgets.QMainWindow):
    """The tree in a window of its own -- shown, never exec'd.

    Same rule as Preferences: a modal here would mean checking the catalog against the map from
    memory, which is exactly the comparison this is for.
    """

    def __init__(self, nodes: pd.DataFrame = None, tables: dict = None, parent=None,
                 organism: str = None, graph=None, sources: dict = None):
        super().__init__(parent)
        self.setWindowTitle("starplast — slot tree")
        self.view = SlotTree(nodes, tables, self, organism=organism, graph=graph, sources=sources)
        self.setCentralWidget(self.view)
        self.resize(1100, 700)
