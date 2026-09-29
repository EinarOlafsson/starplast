"""The star map: one gene in the middle, the genes linked to it around it, and who says so.

The 3D map answers "what is near what" in one embedding. This view answers a different question --
"what links THIS gene to others, and through which evidence or inference" -- by drawing the chosen
gene (the "star") at the centre and its linked genes on rings around it. Every edge is a link from
:mod:`starplast.star_edges` and carries its provenance:

* colour  = its source: each measured layer, each strategy, and one colour for your own runs;
* line    = solid for measured, dashed for inferred;
* width   = its strength (the score's percentile within its run);
* hover   = the gene's id and product, or the link's layer / strategy / run, setting and score.

Click a neighbour to re-centre on it (Back returns), double-click it to select it on the 3D map.
Depth is one or two hops. Each source keeps only its strongest links around each gene
(``per source``), so a dense layer cannot hide a sparse one. The list on the right is the legend,
the per-source toggle and the per-source link count at once.

Runs made in the Strategies tab arrive through :meth:`StarMapPanel.add_result` and appear at once, in
the "your runs" colour; they are also kept on disk (`star_edges.UserEdgeStore`) and come back next
session. :func:`install` puts the panel in the main window as a dock beside Strategies.
"""
from __future__ import annotations

import math
import os

import numpy as np
import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets

from . import star_edges as E
from . import theme as TH

#: How many links each source keeps around each expanded gene, by default.
DEFAULT_PER_SOURCE = 6
#: The most genes one star can hold.
MAX_NODES = 160
#: Ring radii, in scene units, before zoom.
RING1, RING2 = 190.0, 360.0

TIPS = {
    "gene": "Type a gene id (or pick one from the list) and press Enter to put it at the centre. "
            "Its linked genes are drawn around it.",
    "centre": "Put the typed gene at the centre of the star map.",
    "back": "Return to the gene that was at the centre before this one.",
    "depth": "How far out to draw: 1 hop shows the genes linked to the centre; 2 hops also shows "
             "the genes linked to those (each keeping a third as many links, so the map stays "
             "readable).",
    "per_source": "How many links each source (a measured layer, a strategy, your runs) keeps "
                  "around each gene -- the strongest ones. A dense layer such as shared "
                  "compartment joins hundreds of genes; this keeps it from hiding a sparse one "
                  "such as crosslinks.",
    "follow": "When on, selecting a gene anywhere else in the application (the 3D map, search, a "
              "results table) also centres it here.",
    "map": "Select the centre gene on the 3D map and fly to it, with its evidence in the side "
           "panel. Double-clicking a gene here does the same for that gene.",
    "all": "Show links from every source.",
    "none": "Hide links from every source; then tick the ones to compare.",
    "measured": "Show only the measured layers (the data), hiding every inferred source.",
    "sources": "The legend and the switches: one row per source with its colour, and how many "
               "links the centre gene has in it. Untick a row to hide that source. Measured "
               "layers are drawn solid, strategy inferences dashed, and your own runs in gold.",
    "view": "The star map. The centre gene is the star; linked genes sit on the rings. Hover a "
            "gene or a link to read where it comes from; click a gene to re-centre on it; "
            "double-click to select it on the 3D map; scroll to zoom; drag to pan.",
    "info": "What the pointer is over: a gene's id and product, or a link's source, run, "
            "setting and score.",
}


def _qcolor(rgb, alpha: float = 1.0) -> QtGui.QColor:
    c = QtGui.QColor.fromRgbF(float(rgb[0]), float(rgb[1]), float(rgb[2]), float(alpha))
    return c


class _Node(QtWidgets.QGraphicsEllipseItem):
    """A gene on the star map: hover to read it, click to re-centre, double-click to select."""

    def __init__(self, panel, gene: int, r: float, fill: QtGui.QColor, ring: QtGui.QColor):
        super().__init__(-r, -r, 2 * r, 2 * r)
        self.panel, self.gene = panel, int(gene)
        self.setBrush(QtGui.QBrush(fill))
        pen = QtGui.QPen(ring)
        pen.setWidthF(1.6)
        pen.setCosmetic(True)
        self.setPen(pen)
        self.setAcceptHoverEvents(True)
        self.setCursor(QtCore.Qt.CursorShape.PointingHandCursor)
        self.setZValue(10)

    def hoverEnterEvent(self, ev):
        self.panel._hover_gene(self.gene)
        super().hoverEnterEvent(ev)

    def mousePressEvent(self, ev):
        if ev.button() == QtCore.Qt.MouseButton.LeftButton:
            # Deferred: re-centring clears the scene, and this item is still inside its own event.
            QtCore.QTimer.singleShot(0, lambda g=self.gene: self.panel.centre_on(g))
            ev.accept()
            return
        super().mousePressEvent(ev)

    def mouseDoubleClickEvent(self, ev):
        QtCore.QTimer.singleShot(0, lambda g=self.gene: self.panel.select_gene_on_map(g))
        ev.accept()


class _Edge(QtWidgets.QGraphicsPathItem):
    """One link: its pen says its source and strength; hovering it says its provenance."""

    def __init__(self, panel, row, path: QtGui.QPainterPath, pen: QtGui.QPen):
        super().__init__(path)
        self.panel, self.row = panel, row
        self.setPen(pen)
        self.setAcceptHoverEvents(True)
        self.setZValue(1)
        self._pen = QtGui.QPen(pen)

    def shape(self):
        s = QtGui.QPainterPathStroker()
        s.setWidth(8)                     # a hairline is hard to hover; the target is wider
        return s.createStroke(self.path())

    def hoverEnterEvent(self, ev):
        pen = QtGui.QPen(self._pen)
        pen.setWidthF(self._pen.widthF() + 2.0)
        self.setPen(pen)
        self.setZValue(5)
        self.panel._hover_edge(self.row)
        super().hoverEnterEvent(ev)

    def hoverLeaveEvent(self, ev):
        self.setPen(self._pen)
        self.setZValue(1)
        super().hoverLeaveEvent(ev)


class _View(QtWidgets.QGraphicsView):
    """The canvas: a faint starfield and the two rings behind the graph; wheel zooms, drag pans."""

    def __init__(self, scene, panel):
        super().__init__(scene)
        self.panel = panel
        self.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        self.setDragMode(QtWidgets.QGraphicsView.DragMode.ScrollHandDrag)
        self.setTransformationAnchor(QtWidgets.QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        self.setMinimumSize(260, 260)
        self.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        rng = np.random.default_rng(7)
        self._stars = rng.uniform(-1, 1, (160, 2)) * 900, rng.uniform(0.3, 1.2, 160)

    def resizeEvent(self, ev):
        super().resizeEvent(ev)
        self.panel.fit()

    def wheelEvent(self, ev):
        f = 1.15 if ev.angleDelta().y() > 0 else 1 / 1.15
        self.scale(f, f)

    def drawBackground(self, painter, rect):
        pal = self.panel.palette_colors()
        painter.fillRect(rect, QtGui.QColor(pal["bg"]))
        dim = QtGui.QColor(pal["fg_dim"])
        dim.setAlphaF(0.35)
        painter.setPen(QtCore.Qt.PenStyle.NoPen)
        painter.setBrush(dim)
        for (x, y), s in zip(*self._stars):
            painter.drawEllipse(QtCore.QPointF(x, y), s, s)
        ring = QtGui.QColor(pal["border"])
        pen = QtGui.QPen(ring)
        pen.setCosmetic(True)
        pen.setStyle(QtCore.Qt.PenStyle.DotLine)
        painter.setPen(pen)
        painter.setBrush(QtCore.Qt.BrushStyle.NoBrush)
        for r in self.panel.rings:
            painter.drawEllipse(QtCore.QPointF(0, 0), r, r)


class StarMapPanel(QtWidgets.QWidget):
    """A navigable network centred on one gene, drawing measured and inferred links by source."""

    #: The gene now at the centre.
    centre_changed = QtCore.pyqtSignal(str)
    #: A gene to select on the 3D map.
    select_on_map = QtCore.pyqtSignal(str)
    status = QtCore.pyqtSignal(str)

    def __init__(self, nodes: pd.DataFrame, organism: str, ctx=None, graph=None,
                 store: E.UserEdgeStore | None = None, theme="dark",
                 shipped_path: str | None = None, parent=None):
        """Build the panel over one organism's table. Links are loaded on first use.

        :param ctx: a `strategies.Context` over the same table; needed to turn a user run into
            links (and its graph is reused for the measured layers).
        :param graph: measured layers (`layer__a/b/w` arrays); default: `ctx.graph`, else the
            space's own graph file.
        :param store: where the user's runs are kept; `None` keeps them for this session only.
        :param theme: a theme name, or a callable returning the current one.
        """
        super().__init__(parent)
        self.nodes = nodes.reset_index(drop=True)
        self.organism = organism
        self.ctx, self._graph, self.store = ctx, graph, store or E.UserEdgeStore(None)
        self._theme, self._shipped = theme, shipped_path
        self.gene_ids = self.nodes["gene_id"].astype(str).to_numpy()
        from .embedding import as_text
        self.products = (as_text(self.nodes["product"]).to_numpy() if "product" in self.nodes
                         else np.array([""] * len(self.nodes), dtype=object))
        self._pos = {g.upper(): i for i, g in enumerate(self.gene_ids)}
        self._index: E.EdgeIndex | None = None
        self.colors: dict = {}
        self.labels: dict = {}
        self.kinds: dict = {}
        self.enabled: set | None = None      # None: every group, including ones added later
        self.centre: int | None = None
        self.history: list = []
        self.rings: tuple = ()
        self._pending: str | None = None
        self.last_nodes = pd.DataFrame(columns=["gene", "hop", "parent"])
        self.last_edges = E.empty()

        self.gene_edit = QtWidgets.QLineEdit()
        self.gene_edit.setPlaceholderText("gene id…")
        comp = QtWidgets.QCompleter(list(self.gene_ids), self)
        comp.setCaseSensitivity(QtCore.Qt.CaseSensitivity.CaseInsensitive)
        comp.setFilterMode(QtCore.Qt.MatchFlag.MatchContains)
        self.gene_edit.setCompleter(comp)
        self.gene_edit.returnPressed.connect(self._centre_typed)
        self.centre_btn = QtWidgets.QPushButton("Centre")
        self.centre_btn.clicked.connect(self._centre_typed)
        self.back_btn = QtWidgets.QPushButton("◂ Back")
        self.back_btn.clicked.connect(self.back)
        self.back_btn.setEnabled(False)
        self.depth = QtWidgets.QComboBox()
        self.depth.addItems(["1 hop", "2 hops"])
        self.depth.currentIndexChanged.connect(lambda _i: self.redraw())
        self.per_source = QtWidgets.QSpinBox()
        self.per_source.setRange(1, 40)
        self.per_source.setValue(DEFAULT_PER_SOURCE)
        self.per_source.setPrefix("per source ")
        self.per_source.valueChanged.connect(lambda _v: self.redraw())
        self.follow_box = QtWidgets.QCheckBox("Follow selection")
        self.follow_box.setChecked(True)
        self.map_btn = QtWidgets.QPushButton("Select on 3D map")
        self.map_btn.clicked.connect(lambda: self.select_gene_on_map(self.centre))
        self.all_btn = QtWidgets.QPushButton("All")
        self.all_btn.clicked.connect(lambda: self.set_enabled(None))
        self.none_btn = QtWidgets.QPushButton("None")
        self.none_btn.clicked.connect(lambda: self.set_enabled(set()))
        self.measured_btn = QtWidgets.QPushButton("Measured")
        self.measured_btn.clicked.connect(self.only_measured)
        self.sources = QtWidgets.QListWidget()
        self.sources.itemChanged.connect(self._source_toggled)
        self.sources.setMinimumWidth(190)
        self.sources.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.sources.setTextElideMode(QtCore.Qt.TextElideMode.ElideMiddle)
        self.legend = QtWidgets.QLabel()
        self.legend.setWordWrap(True)
        self.headline = QtWidgets.QLabel("No gene at the centre yet.")
        self.headline.setWordWrap(True)
        self.info = QtWidgets.QLabel("Hover a gene or a link.")
        self.info.setWordWrap(True)
        self.info.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.info.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        self.info.setMinimumHeight(54)
        self.scene = QtWidgets.QGraphicsScene(self)
        self.scene.setSceneRect(-900, -900, 1800, 1800)
        self.view = _View(self.scene, self)
        for w, key in ((self.gene_edit, "gene"), (self.centre_btn, "centre"),
                       (self.back_btn, "back"), (self.depth, "depth"),
                       (self.per_source, "per_source"), (self.follow_box, "follow"),
                       (self.map_btn, "map"), (self.all_btn, "all"), (self.none_btn, "none"),
                       (self.measured_btn, "measured"), (self.sources, "sources"),
                       (self.view, "view"), (self.info, "info"), (self.legend, "sources"),
                       (self.headline, "view")):
            w.setToolTip(TH.tip(TIPS[key]))

        top = QtWidgets.QHBoxLayout()
        top.addWidget(self.back_btn)
        top.addWidget(self.gene_edit, 1)
        top.addWidget(self.centre_btn)
        top.addWidget(self.depth)
        top.addWidget(self.per_source)
        bottom = QtWidgets.QHBoxLayout()
        bottom.addWidget(self.follow_box)
        bottom.addStretch(1)
        bottom.addWidget(self.map_btn)
        side = QtWidgets.QVBoxLayout()
        toggles = QtWidgets.QHBoxLayout()
        for b in (self.all_btn, self.none_btn, self.measured_btn):
            toggles.addWidget(b)
        side.addLayout(toggles)
        side.addWidget(self.sources, 1)
        side.addWidget(self.legend)
        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        left = QtWidgets.QWidget()
        ll = QtWidgets.QVBoxLayout(left)
        ll.setContentsMargins(0, 0, 0, 0)
        ll.addWidget(self.headline)
        ll.addWidget(self.view, 1)
        ll.addWidget(self.info)
        right = QtWidgets.QWidget()
        right.setLayout(side)
        split.addWidget(left)
        split.addWidget(right)
        split.setStretchFactor(0, 3)
        split.setStretchFactor(1, 1)
        lay = QtWidgets.QVBoxLayout(self)
        lay.addLayout(top)
        lay.addWidget(split, 1)
        lay.addLayout(bottom)
        QtGui.QShortcut(QtGui.QKeySequence(QtCore.Qt.Key.Key_Backspace), self.view,
                        activated=self.back)

    # ------------------------------------------------------------------ data
    def palette_colors(self) -> dict:
        """The current theme's palette."""
        name = self._theme() if callable(self._theme) else self._theme
        return TH.palette_for(name or "dark")

    def theme_name(self) -> str:
        """The current theme's name."""
        return (self._theme() if callable(self._theme) else self._theme) or "dark"

    @property
    def index(self) -> E.EdgeIndex:
        """Every link of this table (measured, shipped and the user's), built on first use."""
        if self._index is None:
            idx = E.EdgeIndex(len(self.gene_ids))
            graph = self._graph
            if graph is None and self.ctx is not None:
                graph = self.ctx.graph
            if graph is None:
                graph = E.load_graph(self.organism, self.gene_ids)
            idx.add(E.measured_edges(graph))
            idx.add(E.shipped_edges(self.organism, self.gene_ids, self._shipped))
            idx.add(self.store.load(self.organism, self.gene_ids))
            self._index = idx
            self._assign_colors()
            self._fill_sources()
        return self._index

    def _group_order(self) -> list:
        g = self.index.groups() if self._index is not None else pd.DataFrame()
        if g.empty:
            return []
        measured = sorted(g.loc[g["kind"] == E.MEASURED, "group"].unique())
        inferred = [x for x in g.loc[g["kind"] == E.INFERRED, "group"].unique()
                    if x != E.YOUR_RUNS]
        order = {k: i for i, k in enumerate(E.SHIPPED_STRATEGIES)}
        inferred = sorted(inferred, key=lambda k: (order.get(k, 999), k))
        mine = [E.YOUR_RUNS] if (g["group"] == E.YOUR_RUNS).any() else []
        return measured + inferred + mine

    def _assign_colors(self):
        groups = self._group_order()
        known = [g for g in groups if g != E.YOUR_RUNS]
        cols = TH.categorical_colors(max(len(known), 1), self.theme_name())
        self.colors = {g: _qcolor(c) for g, c in zip(known, cols)}
        self.colors[E.YOUR_RUNS] = QtGui.QColor(self.palette_colors()["warning"])
        g = self._index.groups()
        self.kinds = dict(zip(g["group"], g["kind"])) if len(g) else {}
        kinds = self.kinds
        self.labels = {}
        for g in groups:
            if kinds.get(g) == E.MEASURED:
                self.labels[g] = f"{g} · {E.layer_text(g)}"
            elif g == E.YOUR_RUNS:
                self.labels[g] = "your runs"
            else:
                self.labels[g] = E.strategy_title(g)

    def _fill_sources(self):
        """The legend list: a swatch, the source, and the centre's links in it."""
        _ = self.index                  # built (and filled once) before this fill starts
        self.sources.blockSignals(True)
        self.sources.clear()
        totals = self.index.groups().groupby("group")["links"].sum().to_dict()
        here = self.index.counts(self.centre) if self.centre is not None else {}
        for g in self._group_order():
            kind = E.MEASURED if self._is_measured(g) else E.INFERRED
            n = int(here.get(g, 0))
            it = QtWidgets.QListWidgetItem(f"{self._short(g)}   {n:,}")
            it.setData(QtCore.Qt.ItemDataRole.UserRole, g)
            it.setIcon(self._swatch(g, kind))
            it.setFlags(it.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            on = self.enabled is None or g in self.enabled
            it.setCheckState(QtCore.Qt.CheckState.Checked if on else QtCore.Qt.CheckState.Unchecked)
            it.setToolTip(TH.tip(f"{self.labels.get(g, g)}. {kind.capitalize()}. "
                                 f"{n:,} links of the centre gene; {int(totals.get(g, 0)):,} in all."))
            self.sources.addItem(it)
        self.sources.blockSignals(False)
        self.legend.setText(self._legend_html())

    def _is_measured(self, group: str) -> bool:
        _ = self.index
        return self.kinds.get(group) == E.MEASURED

    def _short(self, g: str) -> str:
        if g == E.YOUR_RUNS:
            return "★ your runs"
        text = self.labels.get(g, g)
        return text.split(":")[0] if ":" in text else g

    def _swatch(self, g: str, kind: str) -> QtGui.QIcon:
        pm = QtGui.QPixmap(28, 12)
        pm.fill(QtCore.Qt.GlobalColor.transparent)
        p = QtGui.QPainter(pm)
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        pen = QtGui.QPen(self.colors.get(g, QtGui.QColor("#888")))
        pen.setWidthF(3.0)
        if kind == E.INFERRED:
            pen.setStyle(QtCore.Qt.PenStyle.DashLine)
        p.setPen(pen)
        p.drawLine(2, 6, 26, 6)
        p.end()
        return QtGui.QIcon(pm)

    def _legend_html(self) -> str:
        pal = self.palette_colors()
        gold = self.colors.get(E.YOUR_RUNS, QtGui.QColor(pal["warning"])).name()
        return (f"<span style='color:{pal['fg_muted']}'>── measured layer &nbsp; ╌╌ inferred "
                f"(strategy) &nbsp; <span style='color:{gold}'>╌╌ your runs</span><br>"
                f"width = strength within its run &nbsp; ● size = hop</span>")

    # ------------------------------------------------------------------ user runs
    def add_result(self, result) -> int:
        """Add the links of a strategy run made in the application; returns how many.

        They are kept by the store (on disk when it has a root) and drawn at once under "your runs".
        A run that relates no genes to each other adds nothing and says so.
        """
        if self.ctx is None:
            self.status.emit("star map: this run cannot be linked without the strategy context")
            return 0
        _ = self.index
        frame = self.store.add(result, self.ctx, self.organism)
        if not len(frame):
            self.status.emit(f"star map: {result.strategy} states no gene-gene links")
            return 0
        self.index.add(frame)
        self._assign_colors()
        if self.enabled is not None:
            self.enabled.add(E.YOUR_RUNS)
        self._fill_sources()
        self.redraw()
        self.status.emit(f"star map: {len(frame):,} links from your {result.strategy} run "
                         f"({frame['run'].iloc[0]})")
        return len(frame)

    # ------------------------------------------------------------------ navigation
    def position(self, gene) -> int | None:
        """A gene id (or position) as a position in this table, or None."""
        if gene is None:
            return None
        if isinstance(gene, (int, np.integer)):
            return int(gene) if 0 <= int(gene) < len(self.gene_ids) else None
        return self._pos.get(str(gene).strip().upper())

    def centre_on(self, gene, remember: bool = True) -> bool:
        """Put a gene at the centre and redraw. False if the gene is not in this table."""
        pos = self.position(gene)
        if pos is None:
            self.status.emit(f"star map: no gene {gene!r} in this table")
            return False
        if remember and self.centre is not None and self.centre != pos:
            self.history.append(self.centre)
        self.centre = pos
        self.back_btn.setEnabled(bool(self.history))
        self.gene_edit.setText(self.gene_ids[pos])
        self._fill_sources()
        self.redraw()
        self.centre_changed.emit(self.gene_ids[pos])
        return True

    def back(self) -> bool:
        """Return to the previous centre."""
        if not self.history:
            return False
        prev = self.history.pop()
        return self.centre_on(prev, remember=False)

    def follow(self, gene) -> None:
        """Centre a gene selected elsewhere, if following is on (deferred while hidden)."""
        if not self.follow_box.isChecked():
            return
        pos = self.position(gene)
        if pos is None or pos == self.centre:
            return
        if self.isVisible():
            self.centre_on(pos)
        else:
            self._pending = self.gene_ids[pos]
            self.gene_edit.setText(self._pending)

    def showEvent(self, ev):
        """Draw a gene selected while the panel was hidden."""
        super().showEvent(ev)
        if self._pending:
            pending, self._pending = self._pending, None
            self.centre_on(pending)

    def select_gene_on_map(self, gene) -> None:
        """Ask the main window to select this gene on the 3D map."""
        pos = self.position(gene)
        if pos is not None:
            self.select_on_map.emit(self.gene_ids[pos])

    def _centre_typed(self):
        text = self.gene_edit.text().strip()
        if text:
            self.centre_on(text)

    # ------------------------------------------------------------------ sources
    def groups_shown(self) -> list:
        """The groups whose links are drawn."""
        order = self._group_order()
        return order if self.enabled is None else [g for g in order if g in self.enabled]

    def set_enabled(self, groups) -> None:
        """Show exactly these groups (None: all of them, including ones added later)."""
        self.enabled = None if groups is None else set(groups)
        self._fill_sources()
        self.redraw()

    def only_measured(self) -> None:
        """Show the measured layers only."""
        self.set_enabled({g for g in self._group_order() if self._is_measured(g)})

    def _source_toggled(self, item):
        g = item.data(QtCore.Qt.ItemDataRole.UserRole)
        on = item.checkState() == QtCore.Qt.CheckState.Checked
        if self.enabled is None:
            self.enabled = set(self._group_order())
        (self.enabled.add if on else self.enabled.discard)(g)
        self.redraw()

    # ------------------------------------------------------------------ drawing
    def layout_positions(self, nodes: pd.DataFrame, edges: pd.DataFrame) -> dict:
        """gene -> (x, y): the centre at the origin, hop 1 on the inner ring grouped by source,
        hop 2 on the outer ring beside the gene that reached it."""
        pos = {}
        if nodes.empty:
            return pos
        centre = int(nodes.loc[nodes["hop"] == 0, "gene"].iloc[0])
        pos[centre] = (0.0, 0.0)
        order = {g: i for i, g in enumerate(self._group_order())}
        first = nodes[nodes["hop"] == 1]["gene"].tolist()
        best = {}
        for row in edges.itertuples(index=False):
            for g, o in ((row.a, row.b), (row.b, row.a)):
                if o == centre and g in first:
                    cur = best.get(g)
                    key = (order.get(row.group, 999), -float(row.strength))
                    if cur is None or key < cur:
                        best[g] = key
        first.sort(key=lambda g: (best.get(g, (999, 0.0)), g))
        n1 = len(first)
        r1 = max(RING1, min(300.0, n1 * 7.0))
        angle = {}
        for i, g in enumerate(first):
            a = -math.pi / 2 + 2 * math.pi * i / max(n1, 1)
            angle[g] = a
            pos[g] = (r1 * math.cos(a), r1 * math.sin(a))
        second = nodes[nodes["hop"] == 2]
        if len(second):
            r2 = r1 + max(RING2 - RING1, min(260.0, len(second) * 1.6))
            # Each hop-1 gene owns the arc facing it; the genes it reached share that arc.
            sector = 2 * math.pi / max(n1, 1)
            kids = {}
            for r in second.itertuples(index=False):
                kids.setdefault(int(r.parent), []).append(int(r.gene))
            for parent, genes in kids.items():
                base = angle.get(parent, -math.pi / 2)
                m = len(genes)
                for j, g in enumerate(sorted(genes)):
                    a = base + sector * 0.9 * ((j + 0.5) / m - 0.5)
                    rr = r2 + (j % 3) * 22.0           # staggered, so a crowded arc stays legible
                    pos[g] = (rr * math.cos(a), rr * math.sin(a))
            self.rings = (r1, r2)
        else:
            self.rings = (r1,)
        return pos

    def redraw(self) -> None:
        """Recompute the star around the centre and draw it."""
        self.scene.clear()
        if self.centre is None:
            self.rings = ()
            return
        nodes, edges = E.star(self.index, self.centre, groups=self.groups_shown(),
                              depth=self.depth.currentIndex() + 1,
                              per_group=self.per_source.value(), max_nodes=MAX_NODES)
        self.last_nodes, self.last_edges = nodes, edges
        pos = self.layout_positions(nodes, edges)
        pal = self.palette_colors()
        hop = dict(zip(nodes["gene"], nodes["hop"]))
        # Parallel links between the same two genes are fanned out, so each can be seen and hovered.
        fan = {}
        for row in edges.itertuples(index=False):
            fan.setdefault((min(row.a, row.b), max(row.a, row.b)), []).append(row)
        for (a, b), rows in fan.items():
            if a not in pos or b not in pos:
                continue
            (x1, y1), (x2, y2) = pos[a], pos[b]
            dx, dy = x2 - x1, y2 - y1
            length = math.hypot(dx, dy) or 1.0
            nx, ny = -dy / length, dx / length
            for i, row in enumerate(rows):
                off = (i - (len(rows) - 1) / 2) * 9.0
                path = QtGui.QPainterPath(QtCore.QPointF(x1, y1))
                path.quadTo(QtCore.QPointF((x1 + x2) / 2 + nx * off * 2, (y1 + y2) / 2 + ny * off * 2),
                            QtCore.QPointF(x2, y2))
                col = QtGui.QColor(self.colors.get(row.group, QtGui.QColor(pal["fg_muted"])))
                # Links reaching the outer ring are drawn faint: the centre's own links are what the
                # view is about, and the second hop is context around them.
                outer = max(hop.get(a, 2), hop.get(b, 2)) >= 2
                s = float(row.strength)
                col.setAlphaF((0.15 + 0.3 * s) if outer else (0.55 + 0.4 * s))
                pen = QtGui.QPen(col)
                pen.setCosmetic(True)
                pen.setWidthF((0.6 + 1.2 * s) if outer else (0.8 + 3.2 * s))
                if row.kind == E.INFERRED:
                    pen.setStyle(QtCore.Qt.PenStyle.DashLine)
                pen.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
                self.scene.addItem(_Edge(self, row, path, pen))
        font = QtGui.QFont(self.font())
        small = QtGui.QFont(font)
        small.setPointSizeF(max(6.5, font.pointSizeF() * 0.8))
        n1 = int((nodes["hop"] == 1).sum())
        # Each gene's ring colour is the source of its strongest link on the map.
        ring = {}
        for row in edges.sort_values("strength", ascending=True).itertuples(index=False):
            ring[row.a] = ring[row.b] = row.group
        for row in nodes.itertuples(index=False):
            g, hop = int(row.gene), int(row.hop)
            if g not in pos:
                continue
            r = (13.0, 7.0, 4.5)[hop]
            if hop == 0:
                fill = QtGui.QColor(pal["accent"])
                edge = QtGui.QColor(pal["accent_hi"])
                halo = QtGui.QRadialGradient(0, 0, 46)
                glow = QtGui.QColor(pal["accent"])
                glow.setAlphaF(0.35)
                halo.setColorAt(0.0, glow)
                glow.setAlphaF(0.0)
                halo.setColorAt(1.0, glow)
                h = self.scene.addEllipse(-46, -46, 92, 92, QtGui.QPen(QtCore.Qt.PenStyle.NoPen),
                                          QtGui.QBrush(halo))
                h.setZValue(0)
            else:
                fill = QtGui.QColor(pal["surface_hi"])
                edge = QtGui.QColor(self.colors.get(ring.get(g), QtGui.QColor(pal["fg_muted"])))
            item = _Node(self, g, r, fill, edge)
            item.setPos(*pos[g])
            item.setToolTip(self._gene_html(g))
            self.scene.addItem(item)
            if hop == 0 or (hop == 1 and n1 <= 48):
                label = self.scene.addSimpleText(self.gene_ids[g], font if hop == 0 else small)
                label.setBrush(QtGui.QColor(pal["fg"] if hop == 0 else pal["fg_muted"]))
                br = label.boundingRect()
                x, y = pos[g]
                if hop == 0:
                    label.setPos(-br.width() / 2, 18)
                else:
                    d = math.hypot(x, y) or 1.0
                    ux, uy = x / d, y / d
                    lx = x + ux * 12 - (br.width() if ux < -0.2 else (br.width() / 2 if abs(ux) <= 0.2 else 0))
                    ly = y + uy * 12 - br.height() / 2
                    label.setPos(lx, ly)
                label.setZValue(11)
        shown = edges.groupby("group").size().to_dict() if len(edges) else {}
        total = sum(self.index.counts(self.centre).values())
        self.headline.setText(
            f"<b>{self.gene_ids[self.centre]}</b> · {self.products[self.centre] or 'no product'}"
            f" &nbsp; <span style='color:{pal['fg_muted']}'>{len(edges):,} links drawn of "
            f"{total:,} (enabled sources, strongest {self.per_source.value()} per source), "
            f"{len(nodes) - 1:,} genes, {len(shown)} sources</span>")
        self.fit()

    def fit(self) -> None:
        """Frame the whole star in the view."""
        if self.scene.items():
            rect = self.scene.itemsBoundingRect().adjusted(-40, -40, 40, 40)
            self.view.fitInView(rect, QtCore.Qt.AspectRatioMode.KeepAspectRatio)

    # ------------------------------------------------------------------ hover text
    def _gene_html(self, g: int) -> str:
        counts = self.index.counts(g)
        top = ", ".join(f"{self._short(k)} {v}" for k, v in
                        sorted(counts.items(), key=lambda kv: -kv[1])[:6])
        return (f"<b>{self.gene_ids[g]}</b><br>{self.products[g] or 'no product'}"
                f"<br><span style='color:#9aa6b4'>{sum(counts.values()):,} links: {top}</span>")

    def edge_html(self, row) -> str:
        """A link's provenance as the hover text shows it."""
        a, b = self.gene_ids[int(row.a)], self.gene_ids[int(row.b)]
        if row.kind == E.MEASURED:
            src = f"measured layer <b>{row.source}</b>: {E.layer_text(row.source)}"
            run = "the space's graph (data)"
        else:
            src = f"inferred by <b>{E.strategy_title(row.source)}</b>"
            run = ("your run " if row.origin == E.ORIGIN_USER else "shipped run ") + str(row.run)
            if row.created:
                run += f" · {row.created}"
        parts = [f"<b>{a}</b> — <b>{b}</b>", src, f"{E.HOW.get(row.how, row.how)}; {run}"]
        if row.setting:
            parts.append(f"setting: {row.setting}")
        parts.append(f"score {float(row.score):.4g} · strength {float(row.strength):.2f}")
        if row.note:
            parts.append(str(row.note))
        return "<br>".join(parts)

    def _hover_gene(self, g: int):
        self.info.setText(self._gene_html(g))

    def _hover_edge(self, row):
        self.info.setText(self.edge_html(row))


def install(window) -> StarMapPanel | None:
    """Add the star map to the main window as a dock beside Strategies, and wire it up.

    Strategy runs feed it (`result_ready`); "Select on 3D map" goes through the window's own search;
    selecting a gene in the window centres it here when following is on (the window calls
    :meth:`StarMapPanel.follow` from its pick handler).
    """
    from . import organisms, paths
    strategies = getattr(window, "strategy_panel", None)
    ctx = getattr(strategies, "ctx", None)
    code = ctx.organism if ctx is not None else organisms.detect(window.nodes["gene_id"])
    if code is None:
        return None
    store = E.UserEdgeStore(os.path.join(paths.user_cache_dir(), "star_edges"))
    panel = StarMapPanel(window.nodes, code, ctx=ctx, store=store,
                         theme=lambda: getattr(window, "theme", "dark"))
    panel.status.connect(lambda m: window.statusBar().showMessage(m))
    panel.select_on_map.connect(window._workflow_gene)
    if strategies is not None:
        strategies.result_ready.connect(panel.add_result)
    dock = QtWidgets.QDockWidget("star map", window)
    dock.setFeatures(QtWidgets.QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
    dock.setWidget(panel)
    window.addDockWidget(QtCore.Qt.DockWidgetArea.RightDockWidgetArea, dock)
    anchor = getattr(window, "strategies_dock", None) or getattr(window, "analysis_dock", None)
    if anchor is not None:
        window.tabifyDockWidget(anchor, dock)
    if getattr(window, "right_dock", None) is not None:
        window.right_dock.raise_()
    window.star_map, window.star_map_dock = panel, dock
    return panel
