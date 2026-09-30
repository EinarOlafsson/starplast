#!/usr/bin/env python3
"""The "maps" dock: pregenerated maps to pick from, and how well each label maps onto each.

A row per shipped map (`umap_gallery.Gallery`) with a thumbnail colored by its own clusters; clicking
one puts it in the central 3D view through the window's own `use_embedding`, so it is a map like any
other -- genes click, every color mode applies -- and its HDBSCAN clusters become the window's
cluster coloring. Beside it, a sortable score table read in either direction:

* **One label, every map** -- pick a label, see which map it separates on best;
* **One map, every label** -- pick a map, see which labels it organises.

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
#: The map id the live score of the window's current map is filed under.
ON_SCREEN = "on screen"
LABEL_VIEW, MAP_VIEW = "One label, every map", "One map, every label"
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

        L = QtWidgets.QVBoxLayout(self)
        L.setContentsMargins(8, 8, 8, 8)
        head = QtWidgets.QLabel(
            "<b>Pregenerated maps</b> — click one to show it; the table says how well each label "
            "maps onto each.")
        head.setWordWrap(True)
        L.addWidget(head)

        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        L.addWidget(split, 1)

        self.list = QtWidgets.QListWidget()
        self.list.setIconSize(QtCore.QSize(THUMB, THUMB))
        self.list.setSpacing(2)
        self.list.setWordWrap(True)
        self.list.setToolTip(
            "Maps built ahead of time from measurements only -- all of them at three neighbourhood "
            "sizes, all but localization, and one per kind of evidence. Each thumbnail is colored "
            "by that map's own clusters (gray is unclustered). Click a map to show it in the main "
            "view, where it behaves like any other map.")
        self.list.itemClicked.connect(self._on_item)
        split.addWidget(self.list)

        lower = QtWidgets.QWidget()
        lv = QtWidgets.QVBoxLayout(lower)
        lv.setContentsMargins(0, 0, 0, 0)
        row = QtWidgets.QHBoxLayout()
        self.view_box = QtWidgets.QComboBox()
        self.view_box.addItems([LABEL_VIEW, MAP_VIEW])
        self.view_box.setToolTip(
            "Read the table one of two ways. One label, every map: which map does this label "
            "separate on best? One map, every label: which labels does this map organise? Sort by "
            "any column by clicking its header.")
        self.view_box.currentTextChanged.connect(self.refresh)
        row.addWidget(self.view_box)
        self.label_box = QtWidgets.QComboBox()
        self.label_box.setToolTip(
            "The label to score and to color by: a categorical column such as compartment. In "
            "'one map, every label' view, clicking a row picks its label here.")
        self.label_box.addItems(self.labels())
        self.label_box.currentTextChanged.connect(self._label_changed)
        row.addWidget(self.label_box, 1)
        self.map_box = QtWidgets.QComboBox()
        self.map_box.setToolTip(
            "The map whose scores the table shows in 'one map, every label' view. Follows the "
            "map you last clicked; 'on screen' appears once the map on screen has been scored.")
        self._fill_map_box()
        self.map_box.currentIndexChanged.connect(self.refresh)
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
            "Double-click a row to show that map, or to color by that label.")
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
            "in Analysis. Its clusters are used if it has them; otherwise it is clustered with HDBSCAN "
            "exactly as the gallery maps were (25 / 5, leaf). The result is "
            "added to the table as 'on screen'.")
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
        self._fill_list()
        if "compartment" in self.labels():
            self.label_box.setCurrentText("compartment")
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

    def title(self, map_id: str) -> str:
        """A map's readable name."""
        if map_id == ON_SCREEN:
            return "Map on screen (scored live)"
        for r in self.records:
            if r["id"] == map_id:
                return r["title"]
        return map_id

    # ---------------------------------------------------------------- the list
    def _fill_list(self):
        self.list.clear()
        for r in self.records:
            text = (f"{r['title']}\n{r['n_genes']:,} genes · {r['n_clusters']} clusters · "
                    f"{r['noise_fraction']:.0%} unclustered")
            item = QtWidgets.QListWidgetItem(text)
            item.setData(QtCore.Qt.ItemDataRole.UserRole, r["id"])
            item.setToolTip(f"{r['description']}\n\n{len(r['columns'])} measurement columns in "
                            f"{len(r['blocks'])} block(s); n_neighbors {r['n_neighbors']}; "
                            f"HDBSCAN min cluster size {G.CLUSTERING['min_cluster_size']}, "
                            f"min samples {G.CLUSTERING['min_samples']}, "
                            f"{G.CLUSTERING['cluster_selection_method']} selection.")
            try:
                m = self.gallery.load(self.organism, r["id"])
                from .gallery import thumbnail
                img = thumbnail(m["xyz"], cluster_colors(m["clusters"], self.theme), size=THUMB,
                                background=self.background, point=1.2)
                item.setIcon(QtGui.QIcon(QtGui.QPixmap.fromImage(img)))
            except Exception:                       # a thumbnail is a convenience, not the map
                pass
            self.list.addItem(item)

    def _on_item(self, item):
        map_id = item.data(QtCore.Qt.ItemDataRole.UserRole)
        self.select_map(map_id)
        self.map_chosen.emit(map_id)

    def select_map(self, map_id: str):
        """Mark a map as the current one (list, map box, table) without asking the window."""
        self.current_map = map_id
        for i in range(self.list.count()):
            it = self.list.item(i)
            if it.data(QtCore.Qt.ItemDataRole.UserRole) == map_id:
                self.list.blockSignals(True)
                self.list.setCurrentItem(it)
                self.list.blockSignals(False)
        i = self.map_box.findData(map_id)
        if i >= 0:
            self.map_box.setCurrentIndex(i)
        self.refresh()

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

    # ---------------------------------------------------------------- the table
    def refresh(self, *_):
        """Refill the table for the chosen view, label and map."""
        by_label = self.view_box.currentText() == LABEL_VIEW
        self.map_box.setVisible(not by_label)
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
            s["map"] = [self.title(m) for m in s["map"]]
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
                self.table.setItem(i, j, item)
        self.table.resizeColumnsToContents()
        self.table.setSortingEnabled(True)

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
                self.select_map(key)
                self.map_chosen.emit(key)
        else:
            self.label_box.setCurrentText(str(key))
            self.color_label.emit(str(key))

    def _label_changed(self, *_):
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
            self.map_box.setCurrentIndex(i)
        self.view_box.setCurrentText(MAP_VIEW)
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
