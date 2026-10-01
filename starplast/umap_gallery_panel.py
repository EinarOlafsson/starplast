#!/usr/bin/env python3
"""The "maps" dock: pregenerated maps to pick from, and how well each label maps onto each.

The gallery ships dozens of maps per organism (`umap_gallery.recipes`), so the list is a TREE grouped
by what the map was built from -- all measurements, all but one kind of evidence, one evidence family,
one individual experiment, and the pairs and triples of families -- with a thumbnail per map colored
by its own clusters. It can be searched, sorted (by the label-free structure score, by how well the
chosen label maps, or by how many genes a map places) and flattened into one best-first list. Each row
states the map's recipe and its structure score, so the list can be read without opening anything.

Clicking a map puts it in the central 3D view through the window's own `use_embedding`, so it is a map
like any other -- genes click, every color mode applies -- and its HDBSCAN clusters become the
window's cluster coloring. Beside it, a sortable score table read in either direction:

* **One label, every map** -- pick a label, see which map it separates on best;
* **One map, every label** -- pick a map, see which labels it organises.

The map on screen is marked in both: with `SHOWN_MARK` on its row in the label view, and by being the
only map in the map view. Picking a map therefore always changes the table, which is the bug this
panel was rewritten around -- see `refresh`.

"Color by label" is one click. "Score the map on screen" runs the same scoring on whatever map the
window shows (a map built in Analysis, say), and adds it to both tables as "on screen".

The panel holds no window reference: it emits what the user asked for, and `install` wires those
requests to a window. That keeps it testable without a full application.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets

from . import umap_gallery as G

#: Edge of a map thumbnail in the list, in pixels.
THUMB = 64
#: Points drawn in a thumbnail. A 64-pixel image cannot show more, and the gallery has dozens of maps
#: to draw: `gallery.thumbnail` paints point by point, so a whole proteome per row is seconds of work.
THUMB_POINTS = 2500
#: The map id the live score of the window's current map is filed under.
ON_SCREEN = "on screen"
LABEL_VIEW, MAP_VIEW = "One label, every map", "One map, every label"
#: Marks the map currently on screen, on its row of the label view's table and in the list.
SHOWN_MARK = "▶ "
#: How the map list can be ordered.
SORT_GROUP = "Group, then gallery order"
SORT_STRUCTURE = "Structure score, best first"
SORT_LABEL = "How well the chosen label maps"
SORT_GENES = "Genes placed, most first"
SORTS = (SORT_GROUP, SORT_STRUCTURE, SORT_LABEL, SORT_GENES)
#: Table columns in order, after the first (Map or Label) column.
TABLE_COLUMNS = ("category_f1", "category_skill", "category_precision", "category_recall",
                 "best_category", "best_n", "best_f1_lower", "best_skill", "coverage", "circular")


def _cell(v) -> QtWidgets.QTableWidgetItem:
    """A table cell storing numbers as numbers, so a column sorts 0.9 above 0.10."""
    item = QtWidgets.QTableWidgetItem()
    if isinstance(v, (bool, np.bool_)):
        item.setText("yes" if v else "")
    elif isinstance(v, (int, float, np.integer, np.floating)) and np.isfinite(v):
        item.setData(QtCore.Qt.ItemDataRole.DisplayRole, float(v))
        item.setText(f"{v:.3f}" if isinstance(v, (float, np.floating)) else str(int(v)))
    elif v is None or (isinstance(v, float) and not np.isfinite(v)):
        item.setText("")
    else:
        item.setText("" if str(v) in ("nan", "None") else str(v))
    item.setFlags(item.flags() & ~QtCore.Qt.ItemFlag.ItemIsEditable)
    return item


def cluster_colors(clusters, theme: str = "dark") -> np.ndarray:
    """RGBA per point: a categorical color per cluster, the theme's unknown gray for noise."""
    from . import theme as TH
    lab = np.asarray(clusters)
    out = np.tile(np.asarray(TH.unknown_color(theme), dtype=float), (len(lab), 1))
    ids = sorted(set(lab[lab >= 0].tolist()))
    for col, k in zip(TH.categorical_colors(max(len(ids), 1), theme), ids):
        out[lab == k, :3] = col[:3]
    return out


def recipe_line(record: dict) -> str:
    """One line naming what a map was built from and how it was settled: the row's second line."""
    blocks, cols = len(record.get("blocks") or ()), len(record.get("columns") or ())
    bits = [f"{cols} columns in {blocks} block(s)",
            f"n_neighbors {record.get('n_neighbors', '?')}",
            f"min_dist {record.get('min_dist', 0.25)}"]
    cl = record.get("clustering") or {}
    if cl:
        bits.append(f"HDBSCAN {cl.get('min_cluster_size')}/{cl.get('min_samples')} "
                    f"{cl.get('cluster_selection_method')}")
    if record.get("min_coverage"):
        bits.append(f"genes measured in ≥{float(record['min_coverage']):.0%} of them")
    bits.append("settings searched" if record.get("tuned") else "settings declared, not searched")
    return " · ".join(str(b) for b in bits)


def structure_line(record: dict) -> str:
    """The label-free structure score of a map and the two things it is made of."""
    if record.get("structure") is None:
        return "structure not recorded"
    return (f"structure {float(record['structure']):.3f} "
            f"(clustered {float(record.get('clustered') or 0):.0%}, "
            f"evenness {float(record.get('evenness') or 0):.2f}, "
            f"separation {float(record.get('separation') or 0):.2f})")


class MapGalleryPanel(QtWidgets.QWidget):
    """Pregenerated maps for one organism, and the label x map score table."""

    #: A gallery map id the user clicked: the window shows it.
    map_chosen = QtCore.pyqtSignal(str)
    #: A label the user wants the map colored by.
    color_label = QtCore.pyqtSignal(str)
    #: Color the map by the shown map's clusters.
    color_clusters = QtCore.pyqtSignal()
    #: Score whatever map the window is showing.
    score_screen = QtCore.pyqtSignal()

    def __init__(self, organism: str, gallery: G.Gallery | None = None, theme: str = "dark",
                 background=(0.06, 0.07, 0.09), parent=None):
        """Build the panel over one organism's part of `gallery` (the shipped one by default)."""
        super().__init__(parent)
        self.organism = organism
        self.gallery = gallery or G.shipped()
        self.theme = theme
        self.background = background
        ok = self.gallery.available() and organism in self.gallery.organisms()
        self.scores = self.gallery.scores(organism) if ok else pd.DataFrame()
        self.records = self.gallery.maps(organism) if ok else []
        self.live = pd.DataFrame()
        self.current_map = None
        self._items: dict = {}          # map id -> its tree item
        self._thumbs: dict = {}         # map id -> QIcon, built once and only when shown

        L = QtWidgets.QVBoxLayout(self)
        L.setContentsMargins(8, 8, 8, 8)
        head = QtWidgets.QLabel(
            f"<b>{len(self.records)} pregenerated maps</b> — click one to show it; the table says "
            f"how well each label maps onto each.")
        head.setWordWrap(True)
        head.setToolTip("Every map is built from measurements only, so no label is scored against a "
                        "map that was built from it. The number is how many ship for this organism.")
        L.addWidget(head)

        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        L.addWidget(split, 1)

        upper = QtWidgets.QWidget()
        uv = QtWidgets.QVBoxLayout(upper)
        uv.setContentsMargins(0, 0, 0, 0)
        find = QtWidgets.QHBoxLayout()
        self.search = QtWidgets.QLineEdit()
        self.search.setPlaceholderText("filter maps…")
        self.search.setClearButtonEnabled(True)
        self.search.setToolTip(
            "Show only the maps whose name, group, evidence family or description contains this "
            "text — try 'fitness', 'transcription', 'all but' or 'pairs'. Space-separated words all "
            "have to match. Clearing the box shows every map again.")
        self.search.textChanged.connect(self._fill_list)
        find.addWidget(self.search, 1)
        self.sort_box = QtWidgets.QComboBox()
        self.sort_box.addItems(SORTS)
        self.sort_box.setToolTip(
            "Order the maps. Group, then gallery order: as the gallery was built. Structure score: "
            "the label-free measure the maps were tuned on — how evenly and separably the map "
            "clusters, with no label involved. How well the chosen label maps: the skill of the "
            "label picked below, so the best map for your question comes first. Genes placed: the "
            "widest maps first, since a map built from a sparse experiment places few genes.")
        self.sort_box.currentTextChanged.connect(self._fill_list)
        find.addWidget(self.sort_box, 1)
        self.group_check = QtWidgets.QCheckBox("Group")
        self.group_check.setChecked(True)
        self.group_check.setToolTip(
            "Keep the maps under headings for what they were built from — all measurements, all but "
            "one kind, one evidence family, one experiment, pairs and triples of families. Unticked, "
            "they become one flat list in the chosen order, which is how to find the single best map.")
        self.group_check.toggled.connect(self._fill_list)
        find.addWidget(self.group_check)
        for w in (self.sort_box,):
            w.setSizeAdjustPolicy(
                QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
            w.setMinimumContentsLength(6)
        uv.addLayout(find)

        self.tree = QtWidgets.QTreeWidget()
        self.tree.setHeaderHidden(True)
        self.tree.setIconSize(QtCore.QSize(THUMB, THUMB))
        self.tree.setRootIsDecorated(True)
        self.tree.setUniformRowHeights(False)
        self.tree.setToolTip(
            "Maps built ahead of time from measurements only, grouped by what each was built from. "
            "Every thumbnail is colored by that map's own clusters (gray is unclustered), and every "
            "row states the map's recipe and its structure score. Click a map to show it in the main "
            "view, where it behaves like any other map; hover a row for the full recipe and the "
            "settings that were searched.")
        self.tree.itemClicked.connect(self._on_tree_item)
        self.tree.itemExpanded.connect(self._draw_thumbs)
        uv.addWidget(self.tree, 1)
        self.shown = QtWidgets.QLabel("No map chosen yet.")
        self.shown.setWordWrap(True)
        self.shown.setToolTip("The map now in the main view, with its recipe and structure score.")
        uv.addWidget(self.shown)
        split.addWidget(upper)

        lower = QtWidgets.QWidget()
        lv = QtWidgets.QVBoxLayout(lower)
        lv.setContentsMargins(0, 0, 0, 0)
        row = QtWidgets.QHBoxLayout()
        self.view_box = QtWidgets.QComboBox()
        self.view_box.addItems([LABEL_VIEW, MAP_VIEW])
        self.view_box.setToolTip(
            "Read the table one of two ways. One label, every map: which map does this label "
            "separate on best? The map on screen is marked. One map, every label: which labels does "
            "the map on screen organise? Sort by any column by clicking its header.")
        self.view_box.currentTextChanged.connect(self.refresh)
        row.addWidget(self.view_box)
        self.label_box = QtWidgets.QComboBox()
        self.label_box.setToolTip(
            "The label to score, to color by, and to rank maps on when the list is sorted by it: a "
            "categorical column such as compartment. In 'one map, every label' view, clicking a row "
            "picks its label here.")
        self.label_box.addItems(self.labels())
        self.label_box.currentTextChanged.connect(self._label_changed)
        row.addWidget(self.label_box, 1)
        self.map_box = QtWidgets.QComboBox()
        self.map_box.setToolTip(
            "The map whose scores the table shows in 'one map, every label' view. Follows the map "
            "you last clicked; 'on screen' appears once the map on screen has been scored.")
        self._fill_map_box()
        self.map_box.currentIndexChanged.connect(self._map_box_changed)
        row.addWidget(self.map_box, 1)
        # Combos sized by a short minimum rather than their longest entry, and buttons allowed to
        # shrink: a dock this wide would otherwise set a minimum window width larger than a small
        # monitor (the window-size clamp test is the check).
        for box in (self.view_box, self.label_box, self.map_box):
            box.setSizeAdjustPolicy(
                QtWidgets.QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
            box.setMinimumContentsLength(6)
        lv.addLayout(row)

        self.table = QtWidgets.QTableWidget()
        self.table.setSortingEnabled(True)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.verticalHeader().setVisible(False)
        self.table.setToolTip(
            "Scores of labels on maps. Categories → clusters: how well the label's categories as a "
            "whole fall into clusters. Best category: the one category that falls into one cluster "
            "best, judged on lower bounds so a tiny category cannot score perfectly by luck. Skill "
            "columns subtract what shuffled labels score on the same map: 0 is chance, 1 perfect. "
            f"The row marked {SHOWN_MARK.strip()} is the map on screen. Double-click a row to show "
            "that map, or to color by that label.")
        self.table.itemSelectionChanged.connect(self._row_selected)
        self.table.cellDoubleClicked.connect(self._row_activated)
        lv.addWidget(self.table, 1)

        buttons = QtWidgets.QHBoxLayout()
        self.color_btn = QtWidgets.QPushButton("Color by label")
        self.color_btn.setToolTip(
            "Color the main map by the label chosen above -- one click. Gray is unlabelled, as "
            "everywhere in this application.")
        self.color_btn.clicked.connect(lambda: self.color_label.emit(self.label_box.currentText()))
        self.clusters_btn = QtWidgets.QPushButton("Color by clusters")
        self.clusters_btn.setToolTip(
            "Color the main map by the HDBSCAN clusters of the map on screen -- the clusters the "
            "scores were computed against. Gray is unclustered.")
        self.clusters_btn.clicked.connect(self.color_clusters.emit)
        self.score_btn = QtWidgets.QPushButton("Score map on screen")
        self.score_btn.setToolTip(
            "Score every label on the map the main view is showing now -- for example one you built "
            "in Analysis. Its clusters are used if it has them; otherwise it is clustered with "
            "HDBSCAN exactly as the gallery maps were. The result is added to the table as "
            "'on screen'.")
        self.score_btn.clicked.connect(self.score_screen.emit)
        for b in (self.color_btn, self.clusters_btn, self.score_btn):
            b.setMinimumWidth(40)
            b.setSizePolicy(QtWidgets.QSizePolicy.Policy.Ignored,
                            QtWidgets.QSizePolicy.Policy.Fixed)
            buttons.addWidget(b, 1)
        lv.addLayout(buttons)
        self.note = QtWidgets.QLabel("")
        self.note.setWordWrap(True)
        self.note.setToolTip("What the selected row says, in words.")
        lv.addWidget(self.note)
        split.addWidget(lower)
        split.setStretchFactor(0, 2)
        split.setStretchFactor(1, 3)

        if not self.records:
            self.note.setText("No pregenerated maps ship for this organism.")
        if "compartment" in self.labels():
            self.label_box.setCurrentText("compartment")
        self._fill_list()
        self.refresh()

    # ---------------------------------------------------------------- data
    def labels(self) -> list:
        """Labels scored for this organism, in table order."""
        frames = [f for f in (self.scores, self.live) if not f.empty]
        if not frames:
            return []
        return list(dict.fromkeys(pd.concat(frames).label.astype(str)))

    def all_scores(self) -> pd.DataFrame:
        """Shipped scores plus the live score of the map on screen, if any."""
        frames = [f for f in (self.scores, self.live) if not f.empty]
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    def record(self, map_id: str) -> dict:
        """One map's manifest record, or an empty dict for the map on screen."""
        for r in self.records:
            if r["id"] == map_id:
                return r
        return {}

    def title(self, map_id: str) -> str:
        """A map's readable name."""
        if map_id == ON_SCREEN:
            return "Map on screen (scored live)"
        return self.record(map_id).get("title", map_id)

    def label_skill(self) -> dict:
        """map id -> the chosen label's `category_skill` on it, for sorting and for the rows."""
        s = self.all_scores()
        if s.empty:
            return {}
        s = s[s.label == self.label_box.currentText()]
        return {str(r["map"]): float(r["category_skill"])
                for _, r in s.iterrows() if np.isfinite(r["category_skill"])}

    # ---------------------------------------------------------------- the list
    def _matches(self, record: dict) -> bool:
        """Whether a map passes the filter box: every word of it, anywhere in the map's text."""
        words = self.search.text().lower().split()
        if not words:
            return True
        hay = " ".join(str(record.get(k, "")) for k in
                       ("id", "title", "group", "family", "description")).lower()
        return all(w in hay for w in words)

    def _thumb(self, map_id: str) -> QtGui.QIcon | None:
        """The map's thumbnail, drawn once. A thumbnail is a convenience, not the map."""
        if map_id in self._thumbs:
            return self._thumbs[map_id]
        icon = None
        try:
            m = self.gallery.load(self.organism, map_id)
            from .gallery import thumbnail
            xyz, lab = m["xyz"], m["clusters"]
            if len(xyz) > THUMB_POINTS:              # a 64-pixel image cannot show more
                take = np.linspace(0, len(xyz) - 1, THUMB_POINTS).astype(int)
                xyz, lab = xyz[take], lab[take]
            img = thumbnail(xyz, cluster_colors(lab, self.theme), size=THUMB,
                            background=self.background, point=1.2)
            icon = QtGui.QIcon(QtGui.QPixmap.fromImage(img))
        except Exception:
            icon = None
        self._thumbs[map_id] = icon
        return icon

    def _draw_thumbs(self, parent=None):
        """Give every map row under `parent` (or every visible row) its thumbnail."""
        items = ([parent.child(i) for i in range(parent.childCount())] if parent is not None
                 else list(self._items.values()))
        for it in items:
            map_id = it.data(0, QtCore.Qt.ItemDataRole.UserRole)
            if not map_id or not it.icon(0).isNull():
                continue
            icon = self._thumb(map_id)
            if icon is not None:
                it.setIcon(0, icon)

    def _row_text(self, r: dict, skill: dict) -> str:
        """The two lines of a map's row: its name, then its size, structure and label skill."""
        bits = [f"{int(r.get('n_genes') or 0):,} genes", f"{r.get('n_clusters', '?')} clusters"]
        if r.get("structure") is not None:
            bits.append(f"structure {float(r['structure']):.2f}")
        sk = skill.get(r["id"])
        if sk is not None:
            bits.append(f"{self.label_box.currentText()} skill {sk:+.3f}")
        mark = SHOWN_MARK if r["id"] == self.current_map else ""
        return f"{mark}{r.get('title', r['id'])}\n" + " · ".join(bits)

    def _fill_list(self, *_):
        """Rebuild the map tree for the current filter, sort order and grouping."""
        self.tree.clear()
        self._items = {}
        skill = self.label_skill()
        recs = [r for r in self.records if self._matches(r)]
        key = self.sort_box.currentText()
        if key == SORT_STRUCTURE:
            recs.sort(key=lambda r: -float(r.get("structure") or 0.0))
        elif key == SORT_LABEL:
            recs.sort(key=lambda r: -skill.get(r["id"], -np.inf))
        elif key == SORT_GENES:
            recs.sort(key=lambda r: -int(r.get("n_genes") or 0))
        grouped = self.group_check.isChecked()
        parents: dict = {}
        for r in recs:
            item = QtWidgets.QTreeWidgetItem([self._row_text(r, skill)])
            item.setData(0, QtCore.Qt.ItemDataRole.UserRole, r["id"])
            searched = r.get("settings_searched") or []
            item.setToolTip(0, "\n\n".join(x for x in (
                r.get("description", ""), recipe_line(r), structure_line(r),
                (f"{len(searched)} settings searched; this one won on structure alone, with no "
                 f"label involved." if searched else ""),
                ("Sorted here by " + key.lower() + ".")) if x))
            if grouped:
                g = r.get("group") or "Maps"
                p = parents.get(g)
                if p is None:
                    p = QtWidgets.QTreeWidgetItem([g])
                    p.setFlags(QtCore.Qt.ItemFlag.ItemIsEnabled)
                    p.setToolTip(0, f"{g}: a heading, not a map. Click the arrow to open it.")
                    self.tree.addTopLevelItem(p)
                    parents[g] = p
                p.addChild(item)
            else:
                self.tree.addTopLevelItem(item)
            self._items[r["id"]] = item
        for g, p in parents.items():
            p.setText(0, f"{g} — {p.childCount()} map(s)")
        # The first heading opens, and the one holding the map on screen: dozens of thumbnails are
        # drawn point by point, so the rest wait until their heading is expanded.
        if grouped:
            tops = [self.tree.topLevelItem(i) for i in range(self.tree.topLevelItemCount())]
            for p in tops[:1]:
                p.setExpanded(True)
            cur = self._items.get(self.current_map)
            if cur is not None and cur.parent() is not None:
                cur.parent().setExpanded(True)
        else:
            self._draw_thumbs(None)
        cur = self._items.get(self.current_map)
        if cur is not None:
            self.tree.blockSignals(True)
            self.tree.setCurrentItem(cur)
            self.tree.blockSignals(False)
            self.tree.scrollToItem(cur)

    def items(self) -> list:
        """The map rows in display order: (map id, tree item). Headings are not rows."""
        out = []

        def walk(item):
            map_id = item.data(0, QtCore.Qt.ItemDataRole.UserRole)
            if map_id:
                out.append((map_id, item))
            for i in range(item.childCount()):
                walk(item.child(i))
        for i in range(self.tree.topLevelItemCount()):
            walk(self.tree.topLevelItem(i))
        return out

    def _on_tree_item(self, item, _col=0):
        map_id = item.data(0, QtCore.Qt.ItemDataRole.UserRole)
        if not map_id:                              # a heading, not a map
            item.setExpanded(not item.isExpanded())
            return
        self._on_item(map_id)

    def _on_item(self, item_or_id):
        """A map was chosen in the list: show it and follow it everywhere."""
        map_id = (item_or_id if isinstance(item_or_id, str)
                  else item_or_id.data(0, QtCore.Qt.ItemDataRole.UserRole))
        self.select_map(map_id)
        self.map_chosen.emit(map_id)

    def select_map(self, map_id: str):
        """Mark a map as the current one (list, map box, table) without asking the window."""
        self.current_map = map_id
        cur = self._items.get(map_id)
        if cur is not None:
            self.tree.blockSignals(True)
            self.tree.setCurrentItem(cur)
            self.tree.blockSignals(False)
        i = self.map_box.findData(map_id)
        if i >= 0:
            # Blocked, because `refresh` runs once below either way: letting the combo's own signal
            # through refreshed the table only when the index actually moved, which is how choosing
            # a second map used to leave the first map's scores on screen.
            self.map_box.blockSignals(True)
            self.map_box.setCurrentIndex(i)
            self.map_box.blockSignals(False)
        self._relabel_rows()
        self.refresh()

    def _relabel_rows(self):
        """Move `SHOWN_MARK` to the current map's row without rebuilding the whole tree."""
        skill = self.label_skill()
        for map_id, item in self._items.items():
            r = self.record(map_id)
            if r:
                item.setText(0, self._row_text(r, skill))

    def _fill_map_box(self):
        current = self.map_box.currentData()
        self.map_box.blockSignals(True)
        self.map_box.clear()
        ids = [r["id"] for r in self.records] + ([ON_SCREEN] if not self.live.empty else [])
        for map_id in ids:
            self.map_box.addItem(self.title(map_id), map_id)
        i = self.map_box.findData(current)
        self.map_box.setCurrentIndex(max(i, 0))
        self.map_box.blockSignals(False)

    def _map_box_changed(self, *_):
        """The map combo was used directly: that IS choosing a map, so follow it as a click would."""
        map_id = self.map_box.currentData()
        if map_id and map_id != self.current_map:
            self.current_map = map_id
            self._relabel_rows()
        self.refresh()

    # ---------------------------------------------------------------- the table
    def refresh(self, *_):
        """Refill the table for the chosen view, label and map.

        Both views depend on the map that is current, which is what the rewrite fixed. The label view
        lists every map for one label, so its ROWS do not change when another map is picked -- but
        the marked row does, the selected row does, and the sentence under the table does, and those
        are what told the user nothing had happened. The map view is filtered to the current map, so
        picking a second map replaces its contents outright.
        """
        by_label = self.view_box.currentText() == LABEL_VIEW
        self.map_box.setVisible(not by_label)
        rec = self.record(self.current_map) if self.current_map else {}
        if self.current_map is None:
            self.shown.setText("No map chosen yet — click one above.")
        else:
            self.shown.setText(f"<b>{self.title(self.current_map)}</b><br>{recipe_line(rec)}"
                               f"<br>{structure_line(rec)}" if rec else
                               f"<b>{self.title(self.current_map)}</b>")
        s = self.all_scores()
        if s.empty:
            self.table.setRowCount(0)
            return
        if by_label:
            s = s[s.label == self.label_box.currentText()].copy()
            first = "map"
        else:
            s = s[s["map"] == self.map_box.currentData()].copy()
            first = "label"
        s = s.sort_values("category_skill", ascending=False).reset_index(drop=True)
        # The row's identity (a map id or a label) rides on its first cell, so sorting by a
        # column cannot attribute a row to the wrong map; the map's readable title is displayed.
        keys = s[first].astype(str).tolist()
        if by_label:
            s["map"] = [(SHOWN_MARK if m == self.current_map else "") + self.title(m)
                        for m in s["map"]]
        if "circular" in s:
            # Stored as True/False (or empty when unknown, for the map on screen); shown as a word.
            s["circular"] = ["yes" if str(v) == "True" else "no" if str(v) == "False" else ""
                             for v in s["circular"]]
        self.note.setText("")
        cols = [first] + [c for c in TABLE_COLUMNS if c in s.columns]
        self.table.setSortingEnabled(False)
        self.table.clear()
        self.table.setColumnCount(len(cols))
        self.table.setRowCount(len(s))
        for j, c in enumerate(cols):
            name, tip = G.SCORE_COLUMNS.get(c, (c, ""))
            h = QtWidgets.QTableWidgetItem(name)
            h.setToolTip(tip)
            self.table.setHorizontalHeaderItem(j, h)
        for i, r in enumerate(s[cols].itertuples(index=False)):
            for j, v in enumerate(r):
                item = _cell(v)
                if j == 0:
                    item.setData(QtCore.Qt.ItemDataRole.UserRole, keys[i])
                    item.setToolTip(f"{self.title(keys[i])} — the map on screen."
                                    if keys[i] == self.current_map else str(item.text()))
                self.table.setItem(i, j, item)
        self.table.resizeColumnsToContents()
        self.table.setSortingEnabled(True)
        if by_label and self.current_map is not None:
            self._select_key(self.current_map)

    def _select_key(self, key: str):
        """Select and scroll to the row whose identity is `key`, after any sorting."""
        for i in range(self.table.rowCount()):
            if self._row_key(i) == key:
                self.table.selectRow(i)
                self.table.scrollToItem(self.table.item(i, 0))
                return

    def _row_key(self, row: int):
        item = self.table.item(row, 0)
        return None if item is None else item.data(QtCore.Qt.ItemDataRole.UserRole)

    def _row_selected(self):
        rows = {i.row() for i in self.table.selectedItems()}
        if not rows:
            return
        row = min(rows)
        key = self._row_key(row)
        if key is None:
            return
        s = self.all_scores()
        if self.view_box.currentText() == LABEL_VIEW:
            hit = s[(s.label == self.label_box.currentText()) & (s["map"] == key)]
        else:
            hit = s[(s.label == key) & (s["map"] == self.map_box.currentData())]
            self.label_box.blockSignals(True)
            self.label_box.setCurrentText(str(key))
            self.label_box.blockSignals(False)
        if not hit.empty:
            self.note.setText(sentence(hit.iloc[0], self.title(hit.iloc[0]["map"])))

    def _row_activated(self, row: int, _col: int):
        key = self._row_key(row)
        if key is None:
            return
        if self.view_box.currentText() == LABEL_VIEW:
            if key != ON_SCREEN:
                self._on_item(key)
        else:
            self.label_box.setCurrentText(str(key))
            self.color_label.emit(str(key))

    def _label_changed(self, *_):
        """A new label: the rows' skill figures and, if it is the sort key, the whole order change."""
        self._fill_list()
        if self.view_box.currentText() == LABEL_VIEW:
            self.refresh()

    # ---------------------------------------------------------------- live
    def set_live(self, scores: pd.DataFrame):
        """File the scores of the map on screen under 'on screen' and show them."""
        self.live = scores.assign(map=ON_SCREEN) if not scores.empty else pd.DataFrame()
        known = set(self.labels())
        for lab in self.live.get("label", []):
            if self.label_box.findText(lab) < 0 and lab in known:
                self.label_box.addItem(lab)
        self._fill_map_box()
        i = self.map_box.findData(ON_SCREEN)
        if i >= 0:
            self.map_box.blockSignals(True)
            self.map_box.setCurrentIndex(i)
            self.map_box.blockSignals(False)
            self.current_map = ON_SCREEN
        self.view_box.setCurrentText(MAP_VIEW)
        self._relabel_rows()
        self.refresh()


def sentence(r, title: str) -> str:
    """One row of the table, in words."""
    def f(v):
        try:
            return f"{float(v):.2f}"
        except (TypeError, ValueError):
            return "–"
    circ = " The map was built from a column that restates this label, so this is expected." \
        if str(r.get("circular")) == "True" else ""
    return (f"On {title}, {r['label']}'s categories map to clusters at F1 {f(r['category_f1'])} "
            f"(skill {f(r['category_skill'])}; precision {f(r['category_precision'])}, recall "
            f"{f(r['category_recall'])}). Best single category: {r['best_category']} "
            f"(n = {r['best_n']}), size-aware score {f(r['best_f1_lower'])}, skill "
            f"{f(r['best_skill'])}.{circ}")


# --------------------------------------------------------------------------- the window
def organism_of(window) -> str | None:
    """The window's space code."""
    from . import organisms
    try:
        return organisms.by_species(window.species).code
    except (KeyError, AttributeError):
        return None


def show_map(window, panel: MapGalleryPanel, map_id: str):
    """Put a gallery map in the central view, with its clusters as the window's clusters."""
    m = panel.gallery.aligned(panel.organism, map_id, window.nodes["gene_id"].astype(str))
    window.use_embedding(m["xyz"], m["rows"])
    full = np.full(window.n, -1, dtype=int)
    full[m["rows"]] = m["clusters"]
    window.cluster_labels = full
    window.redraw()
    info = m["info"]
    window.statusBar().showMessage(
        f"map: {info['title']} -- {len(m['rows']):,} genes, {info['n_clusters']} clusters; "
        f"'Color by label' in the maps panel colors it by the label chosen there")


def color_by_label(window, label: str):
    """Color the main map by a label: the color-by box and the category color mode, together."""
    if window.category_box.findText(label) < 0:
        window.statusBar().showMessage(f"{label} is not a color source in this table")
        return
    window.category_box.setCurrentText(label)
    window.set_color_mode("compartment")


def score_on_screen(window, panel: MapGalleryPanel) -> pd.DataFrame:
    """Score every label on the map the window shows; clusters it first when it has none."""
    from . import strategies as S
    from .clustering import cluster
    placed = (np.ones(window.n, bool) if window.placed is None
              else np.asarray(window.placed, bool))
    rows = np.flatnonzero(placed)
    lab = window.cluster_labels
    if lab is None or len(lab) != window.n:
        lab = np.full(window.n, -1, dtype=int)
        lab[rows] = cluster(np.asarray(window.xyz)[rows], "hdbscan", **G.CLUSTERING)
        window.cluster_labels = lab
    ctx = S.Context(window.nodes, graph={}, organism=panel.organism)
    labels = [c for c in panel.labels() if c in window.nodes.columns]
    scores = G.score_map(ctx, rows, np.asarray(lab)[rows], labels=labels or None,
                         map_id=ON_SCREEN, columns=None)
    panel.set_live(scores)
    window.statusBar().showMessage(f"scored {len(scores)} labels on the map on screen "
                                   f"({len(rows):,} genes)")
    return scores


def install(window) -> MapGalleryPanel | None:
    """Add the maps dock to a window, tabbed beside Analysis. None when no gallery ships for it."""
    from . import theme as TH
    code = organism_of(window)
    gallery = G.shipped()
    if code is None or not gallery.available() or code not in gallery.organisms():
        return None
    try:
        from . import help_index
        help_index.DOCK_HELP.setdefault(
            "maps", "Pregenerated maps to pick from, and how well each label maps onto each.")
    except Exception:                                   # the help search is optional here
        pass
    bg = TH.rgbf(TH.palette_for(window.theme)["bg"])[:3]
    panel = MapGalleryPanel(code, gallery=gallery, theme=window.theme, background=bg)
    panel.map_chosen.connect(lambda map_id: show_map(window, panel, map_id))
    panel.color_label.connect(lambda label: color_by_label(window, label))
    panel.color_clusters.connect(lambda: window.set_color_mode("clusters"))
    panel.score_screen.connect(lambda: score_on_screen(window, panel))
    d = QtWidgets.QDockWidget("maps", window)
    d.setFeatures(QtWidgets.QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
    d.setWidget(panel)
    window.addDockWidget(QtCore.Qt.DockWidgetArea.RightDockWidgetArea, d)
    anchor = getattr(window, "analysis_dock", None) or getattr(window, "right_dock", None)
    if anchor is not None:
        window.tabifyDockWidget(anchor, d)
    window.maps_dock = d
    window.maps_panel = panel
    return panel
