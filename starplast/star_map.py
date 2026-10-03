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

**The larger network.** The mode box chooses how much is drawn: one gene (the star above), the
*neighbourhood* of that gene (breadth-first, strongest links first, up to `E.MAX_HOPS` hops and
`E.NEIGHBOURHOOD_MAX` genes), or the *whole network* (up to `E.OVERVIEW_MAX` genes and
`E.LARGE_EDGE_MAX` links, with the largest clusters coloured). Both larger views degrade by keeping
the strongest links and the best-linked genes, and print the sentence that says what they left out.
They are drawn by :class:`_NetworkCanvas`, which builds one path per source and strength bucket once
and paints it into a cached pixmap, so panning, zooming and hovering rebuild nothing.

**What counts as a link** is the user's to decide: which sources may contribute, a minimum strength,
a cap on how many links one gene may bring, and "at least N different sources must agree"
(:class:`star_edges.Definition`). The resulting network is counted live, and a definition can be
saved and reloaded by name (:class:`star_edges.DefinitionStore`, under `paths.user_cache_dir()`).

**Nothing in this module moves the graph when the pointer does.** A hover changes a colour, a
highlight and the text in a fixed-height box, and nothing else: the layout is cached per centre,
mode and definition, the view is framed by the rectangle that layout fixed, and hover text is
rebuilt on a timer. The history of why each of those is necessary is in
`instructions/open/58_star_map_steady_larger_defined.md`.

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
#: The fixed margin an edge item adds around its path, so its bounding rectangle never depends on
#: whether it is hovered.
HOVER_MARGIN = 3.0
#: Milliseconds a pointer movement waits before the hover text is rebuilt. Sweeping across a dense
#: graph fires hundreds of hover events a second; without this the panel spends its time on text.
HOVER_DEBOUNCE_MS = 40
#: Milliseconds a change to the connection rules waits before the network is rebuilt, so dragging
#: the strength slider recomputes once and not a hundred times.
RULE_DEBOUNCE_MS = 180

#: The three views. The star is one gene; the other two are the larger network.
MODE_STAR, MODE_NEIGHBOURHOOD, MODE_NETWORK = "star", "neighbourhood", "network"
MODE_LABELS = ((MODE_STAR, "one gene (star)"), (MODE_NEIGHBOURHOOD, "neighbourhood"),
               (MODE_NETWORK, "whole network"))
#: Default hops and gene cap for the neighbourhood view.
DEFAULT_HOPS = 2
DEFAULT_SIZE_CAP = 600
#: The minimap's side, in pixels.
MINIMAP = 132
#: How close (in pixels) the pointer must come to a gene in the large views to pick it up.
PICK_RADIUS = 9.0
#: Strength buckets the large views batch their line drawing into.
WIDTH_BUCKETS = 3
#: How many clusters of the larger views get a colour of their own, and how many genes a cluster
#: needs before it is worth one.
CLUSTER_COLOURS = 12
MIN_CLUSTER = 4
#: Fixed heights, in pixels, for the two labels that change their text as the pointer moves. They
#: are fixed on purpose: a label that grows with its text resizes the view, and a view that resizes
#: re-frames the graph.
INFO_HEIGHT = 58
HEADLINE_HEIGHT = 22

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
            "setting and score. This box keeps its height whatever it says, so reading it never "
            "moves the graph.",
    "mode": "How much of the network to draw. 'one gene (star)' puts the chosen gene in the middle "
            "with its links on rings. 'neighbourhood' grows outwards from it for as many hops as "
            "you ask, up to the gene cap. 'whole network' lays out every gene the connection rules "
            "below leave, with its clusters visible.",
    "hops": "How far the neighbourhood reaches from the centre gene: 1 is its direct links, 2 also "
            "their links, and so on. The walk takes the strongest links first, so the gene cap "
            "keeps the best-supported part of the neighbourhood.",
    "size_cap": "The most genes the larger views will draw. Raising it shows more and costs more "
            "time to lay out; lowering it keeps the picture readable. Whatever is left out is "
            "stated on screen.",
    "min_strength": "A link is drawn only when its strength -- its score's percentile within its "
            "own run -- is at least this. 0 keeps every link; raise it to keep only each source's "
            "best-supported ones.",
    "max_per_gene": "The most links any one gene may contribute, its strongest first. This is what "
            "stops a hub (a gene in hundreds of shared-compartment links) from filling the picture "
            "by itself. 'no cap' lets every link through.",
    "min_sources": "Draw a pair of genes only when at least this many DIFFERENT sources say they "
            "are linked -- two independent sources agreeing is much better evidence than one. 1 "
            "asks for no agreement at all.",
    "counts": "The network these rules define, right now: how many links survive them, between how "
            "many genes, from how many sources -- and what each rule removed.",
    "definition": "A saved set of connection rules, by name. Pick one to load it, or type a new "
            "name and press Save to keep the rules now set (the sources ticked on the right, the "
            "minimum strength, the per-gene cap and the agreement rule). Saved definitions are "
            "kept with your other state and come back next session.",
    "save_def": "Save the connection rules now set under the name in the box beside this, so you "
            "can come back to exactly this definition of the network.",
    "delete_def": "Forget the saved definition named in the box beside this.",
    "canvas": "The larger network. Drag to pan, scroll to zoom, hover a gene to read it (the graph "
            "does not move), click a gene to make it the centre, double-click to select it on the "
            "3D map. The small picture in the corner is the whole layout with your view marked; "
            "click it to jump.",
    "fit_btn": "Frame the whole network in the view again, at the zoom that shows all of it.",
    "note": "What the caps left out of this picture, if anything.",
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
        self.panel.hover_gene(self.gene)
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
    """One link: its pen says its source and strength; hovering it says its provenance.

    Hovering changes the pen's COLOUR only. It once widened the pen by two units, which grew the
    item's bounding rectangle, which grew the scene's, which made the next fit re-scale the whole
    graph: the map twitched under the pointer. Nothing here may change geometry.
    """

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

    def boundingRect(self):
        """The path's rectangle, widened by a FIXED margin -- never by the hover state."""
        return super().boundingRect().adjusted(-HOVER_MARGIN, -HOVER_MARGIN,
                                               HOVER_MARGIN, HOVER_MARGIN)

    def hoverEnterEvent(self, ev):
        pen = QtGui.QPen(self._pen)                    # same width: geometry must not move
        col = QtGui.QColor(self._pen.color())
        col.setAlphaF(1.0)
        pen.setColor(col)
        self.setPen(pen)
        self.panel.hover_edge(self.row)
        super().hoverEnterEvent(ev)

    def hoverLeaveEvent(self, ev):
        self.setPen(self._pen)
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
        """Re-frame the graph when the VIEW really changed size, and never otherwise.

        This handler is why the map used to twitch. The hover text goes into a label in the same
        layout; a longer line made the label taller, the view shorter, and this handler re-fitted
        the graph -- so moving the pointer rescaled everything. Qt also sends a resize event with an
        unchanged size, so the size is compared before anything is re-framed.
        """
        super().resizeEvent(ev)
        if ev.oldSize() != ev.size() and self.panel.mode() == MODE_STAR:
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


class _NetworkCanvas(QtWidgets.QWidget):
    """A few thousand genes, drawn once into paths and then panned and zoomed over.

    The star view is a `QGraphicsScene` with one item per gene and per link, which is right for a
    hundred of them and wrong for a few thousand. This widget keeps the whole picture as one
    `QPainterPath` per (source, strength bucket) -- built once, in scene coordinates, when the layout
    is set -- and paints those paths through the current transform into a cached pixmap. Panning,
    zooming and hovering do not rebuild anything:

    * the layout is fixed by the caller and never recomputed here, least of all in `paintEvent`;
    * hovering repaints the cached pixmap plus a small overlay, so no geometry moves;
    * the pointer's gene is looked up by distance over the frozen positions, and reported through a
      debounce timer, so sweeping across the graph does not flood the panel with text.
    """

    #: A gene was clicked (make it the centre).
    picked = QtCore.pyqtSignal(int)
    #: A gene was double-clicked (select it on the 3D map).
    opened = QtCore.pyqtSignal(int)
    #: The gene under the pointer, or -1 for none. Debounced.
    hovered = QtCore.pyqtSignal(int)

    def __init__(self, panel, parent=None):
        super().__init__(parent)
        self.panel = panel
        self.setMouseTracking(True)
        self.setMinimumSize(260, 260)
        self.setFocusPolicy(QtCore.Qt.FocusPolicy.StrongFocus)
        self._genes = np.zeros(0, dtype=np.int64)
        self._pos = np.zeros((0, 2), dtype=np.float64)
        self._sizes = np.zeros(0, dtype=np.float64)
        self._fills: list = []
        self._paths: list = []               # (QPainterPath, QColor, width, dashed)
        self._labels: dict = {}
        self._rect = QtCore.QRectF(-500, -500, 1000, 1000)
        self._scale = 1.0
        self._eye = QtCore.QPointF(0, 0)     # the scene point at the middle of the widget
        self._cache: QtGui.QPixmap | None = None
        self._cache_for: tuple | None = None
        self._mini: QtGui.QPixmap | None = None
        self._hover = -1
        self._drag: QtCore.QPoint | None = None
        self._mini_at: QtCore.QRect | None = None
        self._moved = False          # the user has panned or zoomed since the last fit
        self._pending = -1
        self._timer = QtCore.QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(HOVER_DEBOUNCE_MS)
        self._timer.timeout.connect(self._settle)

    # ------------------------------------------------------------------ the frozen picture
    def set_graph(self, genes, pos, sizes, fills, paths, labels) -> None:
        """Take a finished layout: positions, node colours and sizes, edge paths, labels.

        Everything is in scene coordinates and is treated as frozen from here on.
        """
        self._genes = np.asarray(genes, dtype=np.int64)
        self._pos = np.asarray(pos, dtype=np.float64).reshape(-1, 2)
        self._sizes = np.asarray(sizes, dtype=np.float64).reshape(-1)
        self._fills = list(fills)
        self._paths = list(paths)
        self._labels = dict(labels or {})
        self._hover = self._pending = -1
        if len(self._pos):
            lo, hi = self._pos.min(0), self._pos.max(0)
            pad = max(40.0, 0.06 * float(max(hi[0] - lo[0], hi[1] - lo[1], 1.0)))
            self._rect = QtCore.QRectF(lo[0] - pad, lo[1] - pad,
                                       (hi[0] - lo[0]) + 2 * pad, (hi[1] - lo[1]) + 2 * pad)
        else:
            self._rect = QtCore.QRectF(-500, -500, 1000, 1000)
        self._mini = self._mini_at = None
        self.fit()

    def clear(self) -> None:
        """Forget the picture."""
        self.set_graph(np.zeros(0), np.zeros((0, 2)), np.zeros(0), [], [], {})

    def fit(self) -> None:
        """Frame the whole layout, at the zoom that shows all of it."""
        w, h = max(self.width(), 1), max(self.height(), 1)
        self._scale = min(w / max(self._rect.width(), 1e-6), h / max(self._rect.height(), 1e-6))
        self._eye = self._rect.center()
        self._moved = False
        self._invalidate()

    def _invalidate(self):
        self._cache = None
        self.update()

    def resizeEvent(self, ev):
        """A real resize re-frames the layout, unless the user has panned or zoomed since.

        The layout itself never changes here -- only how much of it the widget shows.
        """
        super().resizeEvent(ev)
        if ev.oldSize() == ev.size():
            return
        if self._moved:
            self._invalidate()
        else:
            self.fit()

    # ------------------------------------------------------------------ coordinates
    def to_device(self, xy) -> np.ndarray:
        """Scene points (n, 2) as widget pixels."""
        xy = np.asarray(xy, dtype=np.float64).reshape(-1, 2)
        return np.column_stack([(xy[:, 0] - self._eye.x()) * self._scale + self.width() / 2.0,
                                (xy[:, 1] - self._eye.y()) * self._scale + self.height() / 2.0])

    def to_scene(self, point) -> QtCore.QPointF:
        """A widget point as a scene point."""
        return QtCore.QPointF((point.x() - self.width() / 2.0) / self._scale + self._eye.x(),
                              (point.y() - self.height() / 2.0) / self._scale + self._eye.y())

    def _transform(self) -> QtGui.QTransform:
        t = QtGui.QTransform()
        t.translate(self.width() / 2.0, self.height() / 2.0)
        t.scale(self._scale, self._scale)
        t.translate(-self._eye.x(), -self._eye.y())
        return t

    # ------------------------------------------------------------------ painting
    def _render(self) -> QtGui.QPixmap:
        pal = self.panel.palette_colors()
        ratio = self.devicePixelRatioF()
        pm = QtGui.QPixmap(int(self.width() * ratio), int(self.height() * ratio))
        pm.setDevicePixelRatio(ratio)
        pm.fill(QtGui.QColor(pal["bg"]))
        p = QtGui.QPainter(pm)
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        p.setTransform(self._transform())
        for path, colour, width, dashed in self._paths:
            pen = QtGui.QPen(colour)
            pen.setWidthF(width)
            pen.setCosmetic(True)                 # the transform must not thicken the lines
            if dashed:
                pen.setStyle(QtCore.Qt.PenStyle.DashLine)
            p.setPen(pen)
            p.setBrush(QtCore.Qt.BrushStyle.NoBrush)
            p.drawPath(path)
        p.resetTransform()
        if len(self._pos):
            dev = self.to_device(self._pos)
            by: dict = {}
            for i in range(len(dev)):
                by.setdefault((self._fills[i].rgba(), round(float(self._sizes[i]), 1)),
                              []).append(QtCore.QPointF(dev[i, 0], dev[i, 1]))
            for (rgba, size), points in by.items():
                pen = QtGui.QPen(QtGui.QColor.fromRgba(rgba))
                pen.setWidthF(max(1.5, size))
                pen.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
                p.setPen(pen)
                p.drawPoints(QtGui.QPolygonF(points))
        p.end()
        return pm

    def _draw_minimap(self, p: QtGui.QPainter):
        pal = self.panel.palette_colors()
        if self._mini is None and len(self._pos):
            pm = QtGui.QPixmap(MINIMAP, MINIMAP)
            pm.fill(QtGui.QColor(pal["surface"]))
            mp = QtGui.QPainter(pm)
            mp.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
            s = min(MINIMAP / max(self._rect.width(), 1e-6), MINIMAP / max(self._rect.height(), 1e-6))
            mt = QtGui.QTransform()
            mt.translate(MINIMAP / 2.0, MINIMAP / 2.0)
            mt.scale(s, s)
            mt.translate(-self._rect.center().x(), -self._rect.center().y())
            mp.setTransform(mt)
            for path, colour, _w, _d in self._paths:
                col = QtGui.QColor(colour)
                col.setAlphaF(min(1.0, col.alphaF() * 0.8))
                pen = QtGui.QPen(col)
                pen.setWidthF(0.7)
                pen.setCosmetic(True)
                mp.setPen(pen)
                mp.setBrush(QtCore.Qt.BrushStyle.NoBrush)
                mp.drawPath(path)
            mp.end()
            self._mini = pm
        if self._mini is None:
            return
        x, y = self.width() - MINIMAP - 8, self.height() - MINIMAP - 8
        self._mini_at = QtCore.QRect(x, y, MINIMAP, MINIMAP)
        p.setOpacity(0.92)
        p.drawPixmap(x, y, self._mini)
        p.setOpacity(1.0)
        pen = QtGui.QPen(QtGui.QColor(pal["border"]))
        pen.setWidthF(1.0)
        p.setPen(pen)
        p.setBrush(QtCore.Qt.BrushStyle.NoBrush)
        p.drawRect(self._mini_at)
        # where the view is, inside the whole layout
        s = min(MINIMAP / max(self._rect.width(), 1e-6), MINIMAP / max(self._rect.height(), 1e-6))
        vw = self.width() / self._scale * s
        vh = self.height() / self._scale * s
        cx = x + MINIMAP / 2.0 + (self._eye.x() - self._rect.center().x()) * s
        cy = y + MINIMAP / 2.0 + (self._eye.y() - self._rect.center().y()) * s
        pen.setColor(QtGui.QColor(pal["accent_hi"]))
        pen.setWidthF(1.4)
        p.setPen(pen)
        p.drawRect(QtCore.QRectF(cx - vw / 2, cy - vh / 2, vw, vh))

    def paintEvent(self, ev):
        if self._cache is None or self._cache_for != (self.size(), self.devicePixelRatioF()):
            self._cache = self._render()
            self._cache_for = (self.size(), self.devicePixelRatioF())
        p = QtGui.QPainter(self)
        p.drawPixmap(0, 0, self._cache)
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        pal = self.panel.palette_colors()
        if 0 <= self._hover < len(self._pos):
            # The only thing hover draws: a ring and a name, over the cached picture. The picture
            # itself, and every position in it, is untouched.
            d = self.to_device(self._pos[self._hover:self._hover + 1])[0]
            pen = QtGui.QPen(QtGui.QColor(pal["accent_hi"]))
            pen.setWidthF(2.0)
            p.setPen(pen)
            p.setBrush(QtCore.Qt.BrushStyle.NoBrush)
            r = max(7.0, float(self._sizes[self._hover]))
            p.drawEllipse(QtCore.QPointF(d[0], d[1]), r, r)
            name = self.panel.gene_ids[int(self._genes[self._hover])]
            p.setPen(QtGui.QColor(pal["fg"]))
            p.drawText(QtCore.QPointF(d[0] + r + 4, d[1] - r - 2), name)
        for gene, text in self._labels.items():
            i = int(gene)
            if not 0 <= i < len(self._pos):
                continue
            d = self.to_device(self._pos[i:i + 1])[0]
            p.setPen(QtGui.QColor(pal["fg"]))
            p.drawText(QtCore.QPointF(d[0] + 8, d[1] - 6), str(text))
        self._draw_minimap(p)
        p.end()

    # ------------------------------------------------------------------ pointer
    def _at(self, point) -> int:
        if not len(self._pos):
            return -1
        dev = self.to_device(self._pos)
        d = np.hypot(dev[:, 0] - point.x(), dev[:, 1] - point.y())
        i = int(d.argmin())
        return i if d[i] <= PICK_RADIUS + max(2.0, float(self._sizes[i])) else -1

    def gene_at(self, point) -> int:
        """The gene id at a widget point, or -1."""
        i = self._at(point)
        return int(self._genes[i]) if i >= 0 else -1

    def _settle(self):
        if self._pending != self._hover:
            self._hover = self._pending
            self.update()
        self.hovered.emit(int(self._genes[self._hover]) if self._hover >= 0 else -1)

    def mouseMoveEvent(self, ev):
        if self._drag is not None:
            delta = ev.position().toPoint() - self._drag
            self._drag = ev.position().toPoint()
            self._eye -= QtCore.QPointF(delta.x() / self._scale, delta.y() / self._scale)
            self._moved = True
            self._invalidate()
            return
        i = self._at(ev.position())
        if i != self._pending:
            self._pending = i
            self._timer.start()            # coalesced: a sweep costs one update, not hundreds

    def leaveEvent(self, ev):
        self._pending = -1
        self._timer.start()
        super().leaveEvent(ev)

    def mousePressEvent(self, ev):
        if ev.button() != QtCore.Qt.MouseButton.LeftButton:
            return super().mousePressEvent(ev)
        if getattr(self, "_mini_at", None) is not None and \
                self._mini_at.contains(ev.position().toPoint()):
            s = min(MINIMAP / max(self._rect.width(), 1e-6),
                    MINIMAP / max(self._rect.height(), 1e-6))
            self._eye = QtCore.QPointF(
                self._rect.center().x() + (ev.position().x() - self._mini_at.x() - MINIMAP / 2) / s,
                self._rect.center().y() + (ev.position().y() - self._mini_at.y() - MINIMAP / 2) / s)
            self._moved = True
            self._invalidate()
            return
        i = self._at(ev.position())
        if i >= 0:
            self.picked.emit(int(self._genes[i]))
            return
        self._drag = ev.position().toPoint()
        self.setCursor(QtCore.Qt.CursorShape.ClosedHandCursor)

    def mouseReleaseEvent(self, ev):
        self._drag = None
        self.unsetCursor()
        super().mouseReleaseEvent(ev)

    def mouseDoubleClickEvent(self, ev):
        i = self._at(ev.position())
        if i >= 0:
            self.opened.emit(int(self._genes[i]))
            return
        super().mouseDoubleClickEvent(ev)

    def wheelEvent(self, ev):
        f = 1.15 if ev.angleDelta().y() > 0 else 1 / 1.15
        before = self.to_scene(ev.position())
        self._scale = float(min(max(self._scale * f, 0.02), 60.0))
        after = self.to_scene(ev.position())
        self._eye += before - after                 # zoom around the pointer
        self._moved = True
        self._invalidate()


class StarMapPanel(QtWidgets.QWidget):
    """A navigable network centred on one gene, drawing measured and inferred links by source."""

    #: The gene now at the centre.
    centre_changed = QtCore.pyqtSignal(str)
    #: A gene to select on the 3D map.
    select_on_map = QtCore.pyqtSignal(str)
    status = QtCore.pyqtSignal(str)

    def __init__(self, nodes: pd.DataFrame, organism: str, ctx=None, graph=None,
                 store: E.UserEdgeStore | None = None, theme="dark",
                 shipped_path: str | None = None, definitions: "E.DefinitionStore | None" = None,
                 parent=None):
        """Build the panel over one organism's table. Links are loaded on first use.

        :param ctx: a `strategies.Context` over the same table; needed to turn a user run into
            links (and its graph is reused for the measured layers).
        :param graph: measured layers (`layer__a/b/w` arrays); default: `ctx.graph`, else the
            space's own graph file.
        :param store: where the user's runs are kept; `None` keeps them for this session only.
        :param theme: a theme name, or a callable returning the current one.
        :param definitions: where the user's named connection definitions are kept; `None` keeps
            them for this session only.
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
        self.last_counts: dict = {}
        self.definitions = definitions or E.DefinitionStore(None)
        # Everything below is a cache, and every one of them is keyed by what it depends on, so a
        # picture is computed once and then held still. `_frame_rect` is the rectangle `fit` frames:
        # it comes from the layout and NEVER from the scene's current item bounds, which change when
        # an item is hovered.
        self._defined_cache: tuple | None = None
        self._layout_cache: dict = {}
        self._frame_rect: QtCore.QRectF | None = None
        self._hover_pending = None
        self._note = ""
        self._ready = False              # nothing is drawn (and no links are loaded) until built

        self.gene_edit = QtWidgets.QLineEdit()
        self.gene_edit.setPlaceholderText("gene id…")
        self.gene_edit.setToolTip(TH.tip(
            "The gene at the centre of the map. Type an identifier -- it completes on any part of "
            "one -- and press Enter to put that gene in the middle with everything linked to it "
            "around it. Clicking a neighbour moves the centre there, so you can walk the network."))
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
        self.mode_box = QtWidgets.QComboBox()
        for key, label in MODE_LABELS:
            self.mode_box.addItem(label, key)
        self.mode_box.currentIndexChanged.connect(self._mode_changed)
        self.depth = QtWidgets.QComboBox()
        self.depth.addItems(["1 hop", "2 hops"])
        self.depth.currentIndexChanged.connect(lambda _i: self.redraw())
        self.per_source = QtWidgets.QSpinBox()
        self.per_source.setRange(1, 40)
        self.per_source.setValue(DEFAULT_PER_SOURCE)
        self.per_source.setPrefix("per source ")
        self.per_source.valueChanged.connect(lambda _v: self.redraw())
        self.hops = QtWidgets.QSpinBox()
        self.hops.setRange(1, E.MAX_HOPS)
        self.hops.setValue(DEFAULT_HOPS)
        self.hops.setPrefix("hops ")
        self.hops.valueChanged.connect(lambda _v: self.redraw())
        self.size_cap = QtWidgets.QSpinBox()
        self.size_cap.setRange(50, max(E.OVERVIEW_MAX, E.NEIGHBOURHOOD_MAX))
        self.size_cap.setSingleStep(100)
        self.size_cap.setValue(DEFAULT_SIZE_CAP)
        self.size_cap.setPrefix("up to ")
        self.size_cap.setSuffix(" genes")
        self.size_cap.valueChanged.connect(lambda _v: self.redraw())
        self.fit_btn = QtWidgets.QPushButton("Fit")
        self.fit_btn.clicked.connect(self.fit)
        self.follow_box = QtWidgets.QCheckBox("Follow selection")
        self.follow_box.setChecked(True)
        self.map_btn = QtWidgets.QPushButton("Select on 3D map")
        self.map_btn.clicked.connect(lambda: self.select_gene_on_map(self.centre))

        # ---- what counts as a connection
        self.min_strength = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.min_strength.setRange(0, 100)
        self.min_strength.setValue(0)
        self.min_strength.setMinimumWidth(110)
        self.min_strength.valueChanged.connect(self._rules_changed)
        self.strength_label = QtWidgets.QLabel("0.00")
        self.strength_label.setMinimumWidth(34)
        self.max_per_gene = QtWidgets.QSpinBox()
        self.max_per_gene.setRange(0, 200)
        self.max_per_gene.setValue(0)
        self.max_per_gene.setPrefix("≤ ")
        self.max_per_gene.setSuffix(" per gene")
        self.max_per_gene.setSpecialValueText("no cap per gene")
        self.max_per_gene.valueChanged.connect(self._rules_changed)
        self.min_sources = QtWidgets.QSpinBox()
        self.min_sources.setRange(1, 8)
        self.min_sources.setValue(1)
        self.min_sources.setPrefix("≥ ")
        self.min_sources.setSuffix(" sources must agree")
        self.min_sources.setSpecialValueText("any one source")
        self.min_sources.valueChanged.connect(self._rules_changed)
        self.definition_box = QtWidgets.QComboBox()
        self.definition_box.setEditable(True)
        self.definition_box.setInsertPolicy(QtWidgets.QComboBox.InsertPolicy.NoInsert)
        self.definition_box.setMinimumWidth(150)
        self.definition_box.lineEdit().setPlaceholderText("definition name…")
        self.definition_box.activated.connect(self._definition_chosen)
        self.save_def_btn = QtWidgets.QPushButton("Save")
        self.save_def_btn.clicked.connect(self.save_definition)
        self.delete_def_btn = QtWidgets.QPushButton("Forget")
        self.delete_def_btn.clicked.connect(self.forget_definition)
        self.counts_label = QtWidgets.QLabel("no links yet")
        self.counts_label.setWordWrap(False)
        self.counts_label.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.note_label = QtWidgets.QLabel("")
        self.note_label.setWordWrap(False)
        self.note_label.setTextFormat(QtCore.Qt.TextFormat.RichText)

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
        self.headline.setWordWrap(False)
        self.headline.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.info = QtWidgets.QLabel("Hover a gene or a link.")
        self.info.setWordWrap(True)
        self.info.setTextFormat(QtCore.Qt.TextFormat.RichText)
        self.info.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        self.info.setAlignment(QtCore.Qt.AlignmentFlag.AlignTop | QtCore.Qt.AlignmentFlag.AlignLeft)
        # A FIXED height, and a size policy that ignores what it says. This is the fix for the
        # twitch the user reported: the hover text used to change this label's height, which resized
        # the view, which re-fitted the graph -- so moving the pointer moved the map.
        self.info.setFixedHeight(INFO_HEIGHT)
        self.info.setSizePolicy(QtWidgets.QSizePolicy.Policy.Ignored,
                                QtWidgets.QSizePolicy.Policy.Fixed)
        for w in (self.headline, self.counts_label, self.note_label):
            w.setSizePolicy(QtWidgets.QSizePolicy.Policy.Ignored,
                            QtWidgets.QSizePolicy.Policy.Fixed)
        self.headline.setFixedHeight(HEADLINE_HEIGHT)

        self.scene = QtWidgets.QGraphicsScene(self)
        self.scene.setSceneRect(-900, -900, 1800, 1800)
        self.view = _View(self.scene, self)
        self.canvas = _NetworkCanvas(self)
        self.canvas.picked.connect(lambda g: self.centre_on(int(g)))
        self.canvas.opened.connect(lambda g: self.select_gene_on_map(int(g)))
        self.canvas.hovered.connect(self._canvas_hover)
        self.stack = QtWidgets.QStackedWidget()
        self.stack.addWidget(self.view)
        self.stack.addWidget(self.canvas)
        # Hover text is rebuilt on a timer, not on every mouse-move event: a sweep across a dense
        # graph fires hundreds a second, and each one used to rebuild the label's rich text.
        self._hover_timer = QtCore.QTimer(self)
        self._hover_timer.setSingleShot(True)
        self._hover_timer.setInterval(HOVER_DEBOUNCE_MS)
        self._hover_timer.timeout.connect(self._show_hover)
        # And the connection rules are applied on a timer, so dragging the slider rebuilds once.
        self._rule_timer = QtCore.QTimer(self)
        self._rule_timer.setSingleShot(True)
        self._rule_timer.setInterval(RULE_DEBOUNCE_MS)
        self._rule_timer.timeout.connect(self._apply_rules)

        for w, key in ((self.gene_edit, "gene"), (self.centre_btn, "centre"),
                       (self.back_btn, "back"), (self.depth, "depth"),
                       (self.per_source, "per_source"), (self.follow_box, "follow"),
                       (self.map_btn, "map"), (self.all_btn, "all"), (self.none_btn, "none"),
                       (self.measured_btn, "measured"), (self.sources, "sources"),
                       (self.view, "view"), (self.info, "info"), (self.legend, "sources"),
                       (self.headline, "view"), (self.mode_box, "mode"), (self.hops, "hops"),
                       (self.size_cap, "size_cap"), (self.min_strength, "min_strength"),
                       (self.strength_label, "min_strength"),
                       (self.max_per_gene, "max_per_gene"), (self.min_sources, "min_sources"),
                       (self.counts_label, "counts"), (self.definition_box, "definition"),
                       (self.save_def_btn, "save_def"), (self.delete_def_btn, "delete_def"),
                       (self.canvas, "canvas"), (self.fit_btn, "fit_btn"),
                       (self.note_label, "note")):
            w.setToolTip(TH.tip(TIPS[key]))

        top = QtWidgets.QHBoxLayout()
        top.addWidget(self.back_btn)
        top.addWidget(self.gene_edit, 1)
        top.addWidget(self.centre_btn)
        top.addWidget(self.mode_box)
        top.addWidget(self.depth)
        top.addWidget(self.per_source)
        top.addWidget(self.hops)
        top.addWidget(self.size_cap)
        top.addWidget(self.fit_btn)
        rules = QtWidgets.QHBoxLayout()
        rules.addWidget(QtWidgets.QLabel("a link counts when"))
        rules.addWidget(QtWidgets.QLabel("strength ≥"))
        rules.addWidget(self.min_strength)
        rules.addWidget(self.strength_label)
        rules.addWidget(self.max_per_gene)
        rules.addWidget(self.min_sources)
        rules.addStretch(1)
        rules.addWidget(self.definition_box)
        rules.addWidget(self.save_def_btn)
        rules.addWidget(self.delete_def_btn)
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
        ll.addWidget(self.stack, 1)
        ll.addWidget(self.note_label)
        ll.addWidget(self.info)
        right = QtWidgets.QWidget()
        right.setLayout(side)
        split.addWidget(left)
        split.addWidget(right)
        split.setStretchFactor(0, 3)
        split.setStretchFactor(1, 1)
        lay = QtWidgets.QVBoxLayout(self)
        # The two control rows hold a dozen fixed-width widgets between them, which asked for a
        # 995-pixel panel: docked and tabified that became the whole WINDOW's minimum width, so the
        # window could not be made to fit a small screen. In a scroll area the rows keep their own
        # width and the panel can be narrow, scrolling the controls instead of the window.
        for row in (top, rules):
            holder = QtWidgets.QWidget()
            holder.setLayout(row)
            strip = QtWidgets.QScrollArea()
            strip.setWidget(holder)
            strip.setWidgetResizable(True)
            strip.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
            strip.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            strip.setFixedHeight(holder.sizeHint().height() + 2)
            lay.addWidget(strip)
        lay.addWidget(self.counts_label)
        lay.addWidget(split, 1)
        lay.addLayout(bottom)
        # The map itself may be any width; without this the splitter's own children would put the
        # 995 back.
        for child in (left, right):
            child.setMinimumWidth(120)
        QtGui.QShortcut(QtGui.QKeySequence(QtCore.Qt.Key.Key_Backspace), self,
                        activated=self.back)
        self._refresh_definitions()
        self._mode_changed()
        self._ready = True

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
        self._defined_cache = None            # new links: the rules must be applied again
        self._layout_cache.clear()
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


    # ------------------------------------------------------------------ what counts as a connection
    def mode(self) -> str:
        """Which view is showing: `MODE_STAR`, `MODE_NEIGHBOURHOOD` or `MODE_NETWORK`."""
        return self.mode_box.currentData() or MODE_STAR

    def set_mode(self, mode: str) -> None:
        """Show one of the three views."""
        for i in range(self.mode_box.count()):
            if self.mode_box.itemData(i) == mode:
                self.mode_box.setCurrentIndex(i)
                return

    def _mode_changed(self, *_a):
        mode = self.mode()
        big = mode != MODE_STAR
        self.stack.setCurrentWidget(self.canvas if big else self.view)
        self.depth.setVisible(mode == MODE_STAR)
        self.per_source.setVisible(mode == MODE_STAR)
        self.hops.setVisible(mode == MODE_NEIGHBOURHOOD)
        self.size_cap.setVisible(big)
        self.fit_btn.setVisible(big)
        cap = E.OVERVIEW_MAX if mode == MODE_NETWORK else E.NEIGHBOURHOOD_MAX
        if self.size_cap.value() > cap:
            self.size_cap.blockSignals(True)
            self.size_cap.setValue(cap)
            self.size_cap.blockSignals(False)
        self.size_cap.setMaximum(cap)
        self.redraw()

    def current_definition(self) -> E.Definition:
        """The connection definition the controls now describe."""
        return E.Definition(name=self.definition_box.currentText().strip(),
                            groups=None if self.enabled is None else tuple(self.enabled),
                            min_strength=self.min_strength.value() / 100.0,
                            max_per_gene=self.max_per_gene.value(),
                            min_sources=self.min_sources.value())

    def set_definition(self, definition: E.Definition) -> None:
        """Put a definition into the controls and redraw once."""
        for w in (self.min_strength, self.max_per_gene, self.min_sources, self.sources):
            w.blockSignals(True)
        self.min_strength.setValue(int(round(definition.min_strength * 100)))
        self.max_per_gene.setValue(definition.max_per_gene)
        self.min_sources.setValue(definition.min_sources)
        for w in (self.min_strength, self.max_per_gene, self.min_sources, self.sources):
            w.blockSignals(False)
        self.enabled = None if definition.groups is None else set(definition.groups)
        self.strength_label.setText(f"{definition.min_strength:.2f}")
        self._fill_sources()
        self.redraw()

    def _rules_changed(self, *_a):
        """A rule moved: show its value at once, rebuild the network on the timer."""
        self.strength_label.setText(f"{self.min_strength.value() / 100.0:.2f}")
        self._rule_timer.start()

    def _apply_rules(self):
        self.redraw()

    def defined(self) -> tuple:
        """The links this definition calls edges: (`EdgeIndex`, frame, counts).

        Cached on the definition's signature, so moving the pointer, changing the view or
        re-centring never re-applies the rules.
        """
        definition = self.current_definition()
        key = definition.signature()
        if self._defined_cache is not None and self._defined_cache[0] == key:
            return self._defined_cache[1:]
        frame, counts = E.apply_definition(self.index.edges, definition, len(self.gene_ids))
        idx = E.EdgeIndex(len(self.gene_ids))
        idx.add(frame)
        _ = idx.edges                            # build the lookup now, off any paint path
        self._defined_cache = (key, idx, frame, counts)
        self._layout_cache.clear()
        return idx, frame, counts

    def _counts_html(self, counts: dict) -> str:
        pal = self.palette_colors()
        dropped = [(n, counts.get(k, 0)) for k, n in
                   (("dropped_source", "source off"), ("dropped_strength", "too weak"),
                    ("dropped_agreement", "too few sources"), ("dropped_cap", "over the per-gene cap"))
                   if counts.get(k, 0)]
        tail = ("; ".join(f"{v:,} {n}" for n, v in dropped)) if dropped else "nothing removed"
        return (f"<b>{E.counts_text(counts)}</b> "
                f"<span style='color:{pal['fg_muted']}'>— {tail} "
                f"(of {counts.get('from_sources', 0)} sources in this table)</span>")

    # ------------------------------------------------------------------ saved definitions
    def _refresh_definitions(self, select: str = "") -> None:
        text = select or self.definition_box.currentText()
        self.definition_box.blockSignals(True)
        self.definition_box.clear()
        self.definition_box.addItems(self.definitions.names())
        self.definition_box.setCurrentText(text)
        self.definition_box.blockSignals(False)

    def _definition_chosen(self, _index: int):
        saved = self.definitions.get(self.definition_box.currentText())
        if saved is not None:
            self.set_definition(saved)
            self.status.emit(f"star map: connection definition “{saved.name}” — {saved.describe()}")

    def save_definition(self) -> str:
        """Keep the rules now set under the name in the box; returns the name, or ""."""
        name = self.definition_box.currentText().strip()
        if not name:
            self.status.emit("star map: type a name for this connection definition first")
            return ""
        self.definitions.save(self.current_definition(), name)
        self._refresh_definitions(name)
        self.status.emit(f"star map: saved connection definition “{name}”")
        return name

    def forget_definition(self) -> bool:
        """Forget the saved definition named in the box."""
        name = self.definition_box.currentText().strip()
        if not self.definitions.remove(name):
            self.status.emit(f"star map: no saved connection definition “{name}”")
            return False
        self._refresh_definitions("")
        self.status.emit(f"star map: forgot connection definition “{name}”")
        return True

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

    def _layout_key(self) -> tuple:
        """Everything a layout depends on. Hover is deliberately not part of it."""
        definition = self.current_definition()
        # The whole-network layout is the same picture whichever gene is at the centre, so the
        # centre is left out of its key: re-centring there recolours, it does not re-lay-out.
        centre = None if self.mode() == MODE_NETWORK else self.centre
        return (self.mode(), centre, definition.signature(), self.depth.currentIndex(),
                self.per_source.value(), self.hops.value(), self.size_cap.value())

    def redraw(self) -> None:
        """Rebuild the view for the centre, the mode and the connection rules -- and only for those.

        Nothing here runs on a hover. The positions it computes are cached under
        :meth:`_layout_key`, so asking for the same picture twice gives the same coordinates, and
        `fit` frames the rectangle they define rather than whatever the scene currently measures.
        """
        if not self._ready:
            return
        idx, frame, counts = self.defined()
        self.last_counts = counts
        self.counts_label.setText(self._counts_html(counts))
        if self.mode() != MODE_STAR:
            self._redraw_large(idx, frame)
            return
        self.canvas.clear()
        self.scene.clear()
        self.note_label.setText("")
        if self.centre is None:
            self.rings = ()
            self._frame_rect = None
            return
        nodes, edges = E.star(idx, self.centre, groups=self.groups_shown(),
                              depth=self.depth.currentIndex() + 1,
                              per_group=self.per_source.value(), max_nodes=MAX_NODES)
        self.last_nodes, self.last_edges = nodes, edges
        key = self._layout_key()
        cached = self._layout_cache.get(key)
        if cached is None:
            pos = self.layout_positions(nodes, edges)
            self._layout_cache[key] = (pos, self.rings)
        else:
            pos, self.rings = cached
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
        # The frame is fixed HERE, from the positions, once per layout. `fit` uses nothing else.
        xs = [p[0] for p in pos.values()] or [0.0]
        ys = [p[1] for p in pos.values()] or [0.0]
        pad = 70.0
        self._frame_rect = QtCore.QRectF(min(xs) - pad, min(ys) - pad,
                                         max(xs) - min(xs) + 2 * pad,
                                         max(ys) - min(ys) + 2 * pad)
        self.fit()

    # ------------------------------------------------------------------ the larger network
    def _redraw_large(self, idx: E.EdgeIndex, frame: pd.DataFrame) -> None:
        """Lay out and draw a neighbourhood or the whole network on the cached canvas.

        The expensive parts -- choosing the genes, finding the clusters, running the spring layout --
        happen here, once per :meth:`_layout_key`, and the result is cached. The canvas then turns
        them into one path per (source, strength bucket) and paints those into a pixmap, so panning,
        zooming and hovering never come back to this method.
        """
        self.scene.clear()
        self.rings = ()
        key = self._layout_key()
        cached = self._layout_cache.get(key)
        if cached is None:
            if self.mode() == MODE_NETWORK:
                nodes, edges, note = E.overview(frame, len(self.gene_ids),
                                                max_nodes=self.size_cap.value(),
                                                max_edges=E.LARGE_EDGE_MAX)
            elif self.centre is None:
                nodes, edges, note = (pd.DataFrame(columns=["gene", "hop", "parent"]),
                                      E.empty(), "")
            else:
                nodes = E.neighbourhood(frame, len(self.gene_ids), [self.centre],
                                        hops=self.hops.value(), max_nodes=self.size_cap.value())
                keep = set(nodes["gene"].tolist())
                edges = frame[frame["a"].isin(keep) & frame["b"].isin(keep)].reset_index(drop=True)
                reached = int(len(nodes))
                note = ("" if reached < self.size_cap.value() else
                        f"stopped at the {self.size_cap.value():,}-gene cap — "
                        f"the neighbourhood reaches further")
            genes = nodes["gene"].to_numpy(dtype=np.int64)
            cluster = E.clusters(genes, edges)
            xy = E.layout(genes, edges, cluster)
            cached = (genes, cluster, xy, edges, note)
            self._layout_cache[key] = cached
        genes, cluster, xy, edges, note = cached
        self.last_nodes = pd.DataFrame({"gene": genes, "hop": np.zeros(len(genes), dtype=np.int64),
                                        "parent": np.full(len(genes), -1, dtype=np.int64)})
        self.last_edges = edges
        self._note = note
        pal = self.palette_colors()
        self.note_label.setText(
            f"<span style='color:{pal['warning']}'>{note}</span>" if note else "")
        at = {int(g): i for i, g in enumerate(genes)}
        paths = self._edge_paths(edges, xy, at)
        sizes, fills = self._node_style(genes, cluster, edges, at)
        labels = {}
        if self.centre is not None and int(self.centre) in at and self.mode() != MODE_NETWORK:
            labels[at[int(self.centre)]] = self.gene_ids[int(self.centre)]
        self.canvas.set_graph(genes, xy, sizes, fills, paths, labels)
        head = (f"<b>{self.gene_ids[self.centre]}</b> · neighbourhood, {self.hops.value()} hop"
                f"{'' if self.hops.value() == 1 else 's'}"
                if self.mode() == MODE_NEIGHBOURHOOD and self.centre is not None
                else "<b>the whole network</b>")
        self.headline.setText(
            f"{head} &nbsp; <span style='color:{pal['fg_muted']}'>{len(genes):,} genes, "
            f"{len(edges):,} links, {len(set(cluster.tolist()))} clusters</span>")

    def _edge_paths(self, edges: pd.DataFrame, xy: np.ndarray, at: dict) -> list:
        """One `QPainterPath` per (source, strength bucket): what the canvas paints, built once."""
        out = []
        if edges is None or edges.empty:
            return out
        pal = self.palette_colors()
        s = edges["strength"].to_numpy(dtype=np.float64)
        bucket = np.minimum((s * WIDTH_BUCKETS).astype(int), WIDTH_BUCKETS - 1)
        ia = np.array([at.get(int(v), -1) for v in edges["a"].to_numpy()], dtype=np.int64)
        ib = np.array([at.get(int(v), -1) for v in edges["b"].to_numpy()], dtype=np.int64)
        ok = (ia >= 0) & (ib >= 0)
        groups = edges["group"].to_numpy()
        kinds = edges["kind"].to_numpy()
        for group in pd.unique(groups):
            for b in range(WIDTH_BUCKETS):
                pick = ok & (groups == group) & (bucket == b)
                if not pick.any():
                    continue
                path = QtGui.QPainterPath()
                for u, v in zip(ia[pick], ib[pick]):
                    path.moveTo(xy[u, 0], xy[u, 1])
                    path.lineTo(xy[v, 0], xy[v, 1])
                col = QtGui.QColor(self.colors.get(group, QtGui.QColor(pal["fg_muted"])))
                col.setAlphaF(0.22 + 0.22 * b)
                dashed = str(kinds[pick][0]) == E.INFERRED
                out.append((path, col, 0.7 + 0.7 * b, dashed))
        return out

    def _node_style(self, genes, cluster, edges: pd.DataFrame, at: dict) -> tuple:
        """Sizes (by how many links a gene has) and fills (by cluster) for the large views."""
        pal = self.palette_colors()
        degree = np.zeros(len(genes), dtype=np.float64)
        if edges is not None and len(edges):
            for column in ("a", "b"):
                idx = np.array([at.get(int(v), -1) for v in edges[column].to_numpy()])
                np.add.at(degree, idx[idx >= 0], 1.0)
        top = max(degree.max(), 1.0) if len(degree) else 1.0
        sizes = 3.0 + 6.0 * np.sqrt(degree / top)
        # Only the biggest clusters get a colour of their own. A network of a few thousand genes
        # falls into hundreds of little ones, and colouring every one of them turns the picture into
        # confetti that says nothing; the small ones are drawn in the muted colour instead.
        cluster = np.asarray(cluster, dtype=np.int64)
        wheel = TH.categorical_colors(CLUSTER_COLOURS, self.theme_name())
        rest = QtGui.QColor(pal["fg_muted"])
        rest.setAlphaF(0.75)
        fills = [rest] * len(cluster)
        if len(cluster):
            counts = np.bincount(cluster)
            big = np.argsort(-counts, kind="stable")[:CLUSTER_COLOURS]
            seat = {int(c): i for i, c in enumerate(big) if counts[c] >= MIN_CLUSTER}
            fills = [(_qcolor(wheel[seat[int(c)]], 0.95) if int(c) in seat else rest)
                     for c in cluster]
        if self.centre is not None and int(self.centre) in at:
            i = at[int(self.centre)]
            fills[i] = QtGui.QColor(pal["accent_hi"])
            sizes[i] = max(sizes[i], 11.0)
        return sizes, fills

    def fit(self) -> None:
        """Frame the whole graph in the view, from the rectangle the LAYOUT fixed.

        The rectangle comes from the frozen positions, not from `scene.itemsBoundingRect()`: that
        one answers differently while an item is hovered, and framing the view by it made the graph
        jump under the pointer.
        """
        if self.mode() != MODE_STAR:
            self.canvas.fit()
            return
        if self._frame_rect is not None and self.scene.items():
            self.view.fitInView(self._frame_rect, QtCore.Qt.AspectRatioMode.KeepAspectRatio)

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

    def hover_gene(self, g: int) -> None:
        """The pointer is over a gene: queue its text (nothing is drawn or laid out)."""
        self._hover_pending = ("gene", int(g))
        self._hover_timer.start()

    def hover_edge(self, row) -> None:
        """The pointer is over a link: queue its provenance."""
        self._hover_pending = ("edge", row)
        self._hover_timer.start()

    def _canvas_hover(self, gene: int):
        if int(gene) < 0:
            self._hover_pending = ("none", -1)
            self._hover_timer.start()
        else:
            self.hover_gene(int(gene))

    def _show_hover(self):
        """Put the queued hover text in the box. This is the ONLY thing a hover changes."""
        pending, self._hover_pending = self._hover_pending, None
        if pending is None:
            return
        what, value = pending
        if what == "gene":
            self.info.setText(self._gene_html(int(value)))
        elif what == "edge":
            self.info.setText(self.edge_html(value))
        else:
            self.info.setText("Hover a gene or a link.")

    # Kept for callers of the old private names.
    _hover_gene = hover_gene
    _hover_edge = hover_edge


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
    definitions = E.DefinitionStore(os.path.join(paths.user_cache_dir(), "star_edges"))
    panel = StarMapPanel(window.nodes, code, ctx=ctx, store=store, definitions=definitions,
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
