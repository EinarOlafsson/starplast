#!/usr/bin/env python3
"""The gallery: every embedding a walk produces, visible while the walk is still running.

A hyperparameter walk builds one map per configuration, scores it, and -- until this existed -- threw
the map away and reported a number. That is the wrong way round for the question being asked. Nobody
sweeps 288 configurations to find out that trustworthiness was 0.91; they sweep them to find a map
where a compartment falls out as its own island, and no scalar says whether that happened. So the
maps are kept and shown, and the scores travel with them as captions rather than as the whole answer.

Two modes, because two different things are being done:

* **grid** -- a wall of thumbnails, for finding the one configuration that looks unlike the others.
  Comparison is the point, so every thumbnail is rendered at one size, from one projection, on one
  background: a difference between two tiles is then a difference between two maps.
* **scroll** -- one large view with a slider stepping through the configurations in order, for
  watching what a single hyperparameter does as it changes. A grid shows which is different; a scroll
  shows how it got that way.

**Both fill in as the walk runs.** `add()` takes one `tuning.WalkStep` and is called the moment that
configuration finishes, which is the difference between a tool and a progress bar -- run 12 can be
opened and clicked into while run 13 is still computing. In scroll mode the view follows the newest
map only while the slider is already at the end, the way a log tail does; once the user steps back to
look at something, the walk stops moving the view out from under them.

**Thumbnails are painted, not grabbed from the GL view.** Rendering 288 offscreen framebuffers needs
a live GL context per thumbnail, is slow, and fails on exactly the headless machines the test suite
runs on. A thumbnail is a small orthographic projection of the same coordinates -- deterministic,
fast enough to keep up with a walk, and testable without a display. It is a picture; the map it opens
into is the application's own 3D view, where a gene can be clicked like any other.

Color comes from the caller (`colour_fn`), so a thumbnail is colored by whatever the main map is
colored by. Gray means unknown there and means unknown here.
"""
from __future__ import annotations

import numpy as np
from PyQt6 import QtCore, QtGui, QtWidgets

#: Edge of a grid thumbnail, in pixels. Large enough that a cluster reads as a cluster, small enough
#: that a 288-configuration walk is a wall rather than a scroll.
THUMB = 150
#: Edge of the single large view in scroll mode.
LARGE = 460
#: Fallback point colour when no `colour_fn` is given: an unsaturated blue-grey that reads on both
#: the dark and the light themes without claiming to encode anything.
DEFAULT_POINT = (0.55, 0.72, 0.90, 0.85)


def project(coords: np.ndarray, size: int, margin: float = 0.06):
    """Coordinates to pixels: an orthographic view down the third axis, framed to the box.

    Framed on the data rather than on a fixed extent, because an embedding's scale depends on its
    feature set -- a shared fixed frame would show one walk's maps as specks and another's clipped.
    Both axes are scaled by the SAME factor, so a map that is genuinely elongated still looks
    elongated; scaling each axis to fill the box would make every configuration look equally
    isotropic, which is precisely the difference the gallery exists to show.
    """
    xy = np.asarray(coords, dtype=float)[:, :2]
    lo, hi = xy.min(axis=0), xy.max(axis=0)
    span = float(np.max(hi - lo))
    box = size * (1.0 - 2.0 * margin)
    scale = box / span if span > 1e-9 else 1.0
    centre = (lo + hi) / 2.0
    px = (xy - centre) * scale + size / 2.0
    # y is flipped: the embedding's axes point up, a raster's point down, and without this every
    # thumbnail is a mirror image of the map it opens into.
    return px[:, 0], size - px[:, 1]


def thumbnail(coords, colours=None, size: int = THUMB, background=(0.06, 0.07, 0.09),
              point: float = 2.0) -> QtGui.QImage:
    """One configuration's map as a small image.

    Points are drawn back to front along the third component so that near points cover far ones, the
    same way the 3D view occludes rather than blending -- overlapping points that sum would turn a
    dense thumbnail white and every configuration would look alike.
    """
    img = QtGui.QImage(size, size, QtGui.QImage.Format.Format_ARGB32_Premultiplied)
    img.fill(QtGui.QColor.fromRgbF(*background[:3]))
    coords = np.asarray(coords, dtype=float)
    if coords.ndim != 2 or len(coords) == 0:
        return img
    x, y = project(coords, size)
    c = np.asarray(colours, dtype=float) if colours is not None else None
    if c is None or c.ndim != 2 or len(c) != len(coords):
        c = np.tile(np.asarray(DEFAULT_POINT, dtype=float), (len(coords), 1))
    if c.shape[1] == 3:
        c = np.hstack([c, np.ones((len(c), 1))])
    order = (np.argsort(coords[:, 2]) if coords.shape[1] > 2 else np.arange(len(coords)))
    p = QtGui.QPainter(img)
    p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
    pen = QtGui.QPen()
    pen.setWidthF(point)
    pen.setCapStyle(QtCore.Qt.PenCapStyle.RoundCap)
    for i in order:
        pen.setColor(QtGui.QColor.fromRgbF(*np.clip(c[i], 0.0, 1.0)))
        p.setPen(pen)
        p.drawPoint(QtCore.QPointF(float(x[i]), float(y[i])))
    p.end()
    return img


def caption(step) -> str:
    """The configuration and what it scored, in one line.

    The hyperparameters lead because they are what distinguishes one tile from the next; the scores
    follow because they are the reason to look closer, not the answer. Missing scores are left out
    rather than printed as nan -- a walk run without the cluster check has no cluster count, and
    "clusters=nan" reads as a failure rather than as a question that was not asked.
    """
    row = getattr(step, "row", {}) or {}
    bits = []
    if row.get("blocks"):
        bits.append(str(row["blocks"]))
    bits += [f"nn={row.get('n_neighbors', '?')}", f"md={row.get('min_dist', '?')}"]
    t = row.get("trustworthiness")
    if t is not None and np.isfinite(t):
        bits.append(f"trust {t:.3f}")
    # A scored configuration leads with what it scored; a plain walk has no score to show.
    for key, name in (("mean_f1", "mean F1"), ("best_f1", "best F1")):
        v = row.get(key)
        if v is not None and np.isfinite(v):
            bits.append(f"{name} {v:.3f}")
    k = row.get("n_clusters_hdbscan", row.get("n_clusters"))
    if k is not None:
        noise = row.get("noise_frac")
        bits.append(f"{int(k)} clusters" + (f", {noise:.0%} noise" if noise is not None else ""))
    return "  ·  ".join(bits)


class GalleryPanel(QtWidgets.QWidget):
    """Thumbnails of every configuration a walk has finished, in grid or scroll form.

    Holds the `WalkStep`s themselves, not only their pictures, so that clicking a tile can hand the
    real coordinates back to the window rather than a bitmap. `chosen` is what the window connects to.
    """

    #: A configuration the user asked to see, as a `tuning.WalkStep`. The window shows it in the
    #: central 3D view, where it behaves like any other map -- clickable genes, colour modes, edges.
    chosen = QtCore.pyqtSignal(object)

    def __init__(self, colour_fn=None, background=(0.06, 0.07, 0.09), parent=None):
        super().__init__(parent)
        #: Called with a step's boolean gene mask; returns one RGBA row per coordinate. Supplied by
        #: the window so a thumbnail is coloured by whatever the main map is coloured by.
        self.colour_fn = colour_fn
        self.background = background
        self.steps = []

        bar = QtWidgets.QHBoxLayout()
        self.mode_box = QtWidgets.QComboBox()
        self.mode_box.addItems(["grid", "scroll"])
        self.mode_box.setToolTip(
            "grid puts every configuration side by side at one size, which is how you spot the one "
            "map that is not like the others. scroll steps through them in order in a single large "
            "view, which is how you see what changing one hyperparameter actually does. The same "
            "maps either way.")
        self.mode_box.currentTextChanged.connect(self.set_mode)
        self.count = QtWidgets.QLabel("no maps yet")
        self.count.setToolTip(
            "How many configurations have finished, out of how many the walk expects to run. Maps "
            "appear here as they are computed, so this is a walk in progress rather than a report "
            "on one that ended.")
        bar.addWidget(QtWidgets.QLabel("show as"))
        bar.addWidget(self.mode_box)
        bar.addStretch(1)
        bar.addWidget(self.count)

        self.list = QtWidgets.QListWidget()
        self.list.setViewMode(QtWidgets.QListView.ViewMode.IconMode)
        self.list.setIconSize(QtCore.QSize(THUMB, THUMB))
        self.list.setResizeMode(QtWidgets.QListView.ResizeMode.Adjust)
        self.list.setMovement(QtWidgets.QListView.Movement.Static)
        self.list.setSpacing(8)
        self.list.setWordWrap(True)
        self.list.setToolTip(
            "Every configuration this walk has finished. Click one to open it in the main view, "
            "where it is a map like any other -- click a gene for its evidence, color by anything, "
            "draw edges. The caption is the configuration and what it scored; a map that looks "
            "structured is a lead to check, not a result.")
        self.list.itemClicked.connect(self._on_item)

        page = QtWidgets.QWidget()
        pv = QtWidgets.QVBoxLayout(page)
        pv.setContentsMargins(0, 0, 0, 0)
        self.big = QtWidgets.QLabel()
        self.big.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.big.setMinimumSize(240, 240)
        self.big.setToolTip(
            "One configuration at a time, in the order the walk ran them. Stepping through in order "
            "is how a hyperparameter's effect becomes visible: n_neighbors rising turns local detail "
            "into global shape, and no score reports that.")
        self.caption = QtWidgets.QLabel("")
        self.caption.setAlignment(QtCore.Qt.AlignmentFlag.AlignCenter)
        self.caption.setWordWrap(True)
        self.slider = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.slider.setMinimum(0)
        self.slider.setMaximum(0)
        self.slider.setToolTip(
            "Step through the configurations in the order the walk ran them. While it sits at the "
            "newest map it follows the walk; move it back and it stays where you put it, so a walk "
            "that is still running cannot pull the view off what you are looking at.")
        self.slider.valueChanged.connect(self.show_index)
        self.open_btn = QtWidgets.QPushButton("open this map in the main view")
        self.open_btn.setToolTip(
            "Show this configuration in the 3D view, replacing what is there. It becomes the map: "
            "genes are clickable, every color mode applies, and the walk keeps running behind it. "
            "Genes outside the walk's subsample have no position in this embedding and are hidden "
            "rather than drawn at the origin.")
        self.open_btn.clicked.connect(self.expand_current)
        pv.addWidget(self.big, 1)
        pv.addWidget(self.caption)
        pv.addWidget(self.slider)
        pv.addWidget(self.open_btn)

        self.stack = QtWidgets.QStackedWidget()
        self.stack.addWidget(self.list)
        self.stack.addWidget(page)

        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.addLayout(bar)
        self.hint = QtWidgets.QLabel(
            "<i>Run <b>walk hyperparameters</b> on the Map tab. Maps appear here one at a time, as "
            "each is computed.</i>")
        self.hint.setWordWrap(True)
        lay.addWidget(self.hint)
        lay.addWidget(self.stack, 1)

    # ------------------------------------------------------------------ contents
    def clear(self):
        """Drop every map. Called when a new walk starts, so two walks cannot be read as one."""
        self.steps = []
        self.list.clear()
        self.big.clear()
        self.caption.setText("")
        self.slider.setMaximum(0)
        self.count.setText("no maps yet")
        self.hint.show()

    #: Colours for a step that carries its own clustering. Built once, here, rather than taken from
    #: the window's palette, because a thumbnail of a scored configuration has to show the
    #: clustering that was scored -- the map's current colour mode is about something else.
    CLUSTER_COLOURS = [
        (0.95, 0.75, 0.20), (0.35, 0.70, 0.95), (0.45, 0.85, 0.45), (0.95, 0.45, 0.55),
        (0.70, 0.55, 0.95), (0.30, 0.85, 0.80), (0.95, 0.60, 0.30), (0.60, 0.80, 0.35),
    ]
    #: Unclustered points. Grey, for the same reason grey means unknown everywhere else: HDBSCAN
    #: calling a gene noise is a finding about that gene, not a gap in the drawing.
    NOISE_COLOUR = (0.45, 0.45, 0.48, 0.55)

    def cluster_colours(self, labels) -> np.ndarray:
        """One color per point from a clustering, noise in gray."""
        labels = np.asarray(labels)
        out = np.tile(np.asarray(self.NOISE_COLOUR, dtype=float), (len(labels), 1))
        ids = sorted({int(v) for v in labels if v >= 0})
        for k, cid in enumerate(ids):
            out[labels == cid, :3] = self.CLUSTER_COLOURS[k % len(self.CLUSTER_COLOURS)]
            out[labels == cid, 3] = 0.9
        return out

    def colours_for(self, step):
        """The point colors for one step: its own clustering if it has one, else the window's."""
        labels = getattr(step, "labels", None)
        if labels is not None and len(labels) == len(step.coords):
            return self.cluster_colours(labels)
        if self.colour_fn is None:
            return None
        try:
            c = self.colour_fn(step.genes)
        except Exception as exc:
            # A gallery is a viewer. A colouring that fails -- a column dropped, a clustering of the
            # wrong length -- must cost the colour, not the picture, and must say so once rather than
            # raising out of a signal handler where Qt will simply print it and continue.
            print(f"starplast: gallery coloring unavailable ({type(exc).__name__}: {exc})")
            return None
        c = np.asarray(c, dtype=float) if c is not None else None
        return c if (c is not None and len(c) == len(step.coords)) else None

    def add(self, step):
        """Add one finished configuration, and show it if the view is at the end.

        Following the newest map only when the slider is already at the last one is the log-tail
        rule: a walk that yanked the view to run 13 while run 12 was being read would make the
        gallery unusable during exactly the period it exists for.
        """
        at_end = self.slider.value() >= self.slider.maximum()
        self.steps.append(step)
        img = thumbnail(step.coords, self.colours_for(step), size=THUMB,
                        background=self.background)
        it = QtWidgets.QListWidgetItem(QtGui.QIcon(QtGui.QPixmap.fromImage(img)), caption(step))
        it.setData(QtCore.Qt.ItemDataRole.UserRole, len(self.steps) - 1)
        it.setToolTip(caption(step) + "\nClick to open this map in the main view.")
        self.list.addItem(it)
        total = getattr(step, "total", 0) or 0
        self.count.setText(f"{len(self.steps)} of {total} computed" if total
                           else f"{len(self.steps)} computed")
        self.hint.setVisible(False)
        self.slider.blockSignals(True)
        self.slider.setMaximum(len(self.steps) - 1)
        self.slider.blockSignals(False)
        if at_end:
            self.slider.setValue(len(self.steps) - 1)
            self.show_index(len(self.steps) - 1)
        return it

    # ------------------------------------------------------------------ modes
    def set_mode(self, mode: str):
        """Switch between the wall of thumbnails and the single stepped view.

        The chooser is set from here as well as reading into it, so that a mode changed in code --
        by a test, by a keyboard shortcut, by anything but the combo box itself -- does not leave the
        control claiming one mode while the panel shows the other.
        """
        self.stack.setCurrentIndex(1 if mode == "scroll" else 0)
        if self.mode_box.currentText() != mode:
            self.mode_box.blockSignals(True)          # or it calls straight back into here
            self.mode_box.setCurrentText(mode)
            self.mode_box.blockSignals(False)
        if mode == "scroll":
            self.show_index(self.slider.value())

    def show_index(self, i: int):
        """Render one configuration into the large view."""
        if not (0 <= i < len(self.steps)):
            return
        step = self.steps[i]
        img = thumbnail(step.coords, self.colours_for(step), size=LARGE,
                        background=self.background, point=2.6)
        self.big.setPixmap(QtGui.QPixmap.fromImage(img))
        self.caption.setText(f"<b>{i + 1} of {len(self.steps)}</b> &nbsp; {caption(step)}")

    def current(self):
        """The step the scroll view is showing, or None if nothing has arrived yet."""
        i = self.slider.value()
        return self.steps[i] if 0 <= i < len(self.steps) else None

    # ------------------------------------------------------------------ opening one
    def expand_current(self):
        """Ask for the scroll view's configuration to be shown in the central 3D view."""
        step = self.current()
        if step is not None:
            self.chosen.emit(step)

    def _on_item(self, item):
        i = item.data(QtCore.Qt.ItemDataRole.UserRole)
        if i is not None and 0 <= i < len(self.steps):
            # Keep the scroll view in step with the grid, so switching modes after clicking a tile
            # shows the tile that was clicked rather than wherever the slider happened to be.
            self.slider.setValue(i)
            self.chosen.emit(self.steps[i])
