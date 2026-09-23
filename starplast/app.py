#!/usr/bin/env python3
"""Desktop browser for gene evidence in Toxoplasma gondii and Plasmodium falciparum.

Each point represents a gene in a precomputed or user-built embedding. Colours show
annotations or measurements; independently selectable edge layers show the type of
relationship between genes. Proximity in the embedding suggests similar input
features, but does not establish a physical interaction or shared function.

The window also hosts table import, feature selection, clustering, held-out recovery,
and export tools. Long analyses run through the job runner so the interface stays
responsive. Use :mod:`starplast.discover` for discovery searches without a display.
"""
from __future__ import annotations

import inspect
import json
import os
import sys

import numpy as np
import pandas as pd

# pyqtgraph binds to whichever Qt it finds in sys.modules first. If anything imported PySide6 earlier,
# its half of the GL widget comes from PySide6 while ours comes from PyQt6, and the app fails to import.
# This package depends on PyQt6, so say so rather than depending on import order.
os.environ.setdefault("PYQTGRAPH_QT_LIB", "PyQt6")

from PyQt6 import QtCore, QtGui, QtWidgets  # noqa: E402

from . import theme as TH  # noqa: E402
import pyqtgraph as pg  # noqa: E402
import pyqtgraph.opengl as gl  # noqa: E402

from . import paths  # noqa: E402
from . import sprite as _sprite  # noqa: E402

# Before any GL widget is built, and at import rather than in `main`, because the window is
# constructed directly by the suite and by anything embedding it. pyqtgraph 0.14 refuses a
# GLViewWidget whose surface format reports less than OpenGL 2.1, and Qt's untouched default reports
# 2.0 on every driver -- so without this the map fails to draw on hardware that exceeds the
# requirement, and says the driver is at fault while doing it.
_sprite.ensure_gl_format()
from .chat import ChatPanel  # noqa: E402
from .console import ConsolePanel  # noqa: E402
from .jobs import FAILED, JobRunner  # noqa: E402
from . import logging_util  # noqa: E402
# The same absence labels the held-out search excludes. Shared rather than restated, so "unassigned"
# cannot come to mean one thing in the scoring and another in the browser.
from .search import ABSENCE_LABELS as ABSENCE  # noqa: E402

# Resolved rather than assumed. The cache used to be located as "the directory above the package",
# which is true in a source checkout and false in an installed wheel.
DATA = paths.data_dir()

EDGE_TYPES = [
    ("comention", "co-mention (33,924 abstracts)"),
    ("comention_ft", "co-mention (open-access full texts)"),
    ("orthogroup", "shared orthogroup"),
    ("coexpression", "co-expression (stage series)"),
    ("compartment", "shared compartment (hyperLOPIT)"),
    ("cofitness", "co-fitness (7 CRISPR screens)"),
    ("domain", "shared InterPro domain"),
    ("xlms", "crosslink MS — measured physical proximity"),
    ("ip_ms", "IP-MS — replicated pulldown of a tagged bait"),
    ("struct", "structural similarity (Foldseek TM ≥ 0.7)"),
    ("structural_hole", "structural hole (biology links them, literature does not)"),
    ("unwritten_interaction", "measured to bind, never written about"),
]
# Both co-mention types are attention-biased and both carry a corrected residual, so the attention toggle
# governs each of them. They are kept separate because they are different populations: every abstract in
# the field, versus only the papers a publisher deposited open access.
COMENTION = ("comention", "comention_ft")
FIT = ["fit_invitro_hff", "fit_invivo_PE", "fit_invivo_lung", "fit_invivo_liver",
       "fit_invivo_spleen", "fit_naive_bmdm", "fit_ifng"]

EDGE_CAP = 20000        # per type, on drawing only. Stated in the tooltip rather than applied silently.

COLOR_MODES = ["compartment", "compartment (incl. transferred)", "clusters", "in vitro fitness",
                "publications", "depth of attention", "structure confidence (pLDDT)",
                "cyst / tachyzoite expression", "annotations"]

# None means "follow the point style". The rest are absolute pixel sizes.
POINT_SIZES = [("Automatic", None), ("Tiny (2 px)", 2.0), ("Small (4 px)", 4.0),
               ("Medium (7 px)", 7.0), ("Large (11 px)", 11.0), ("Huge (16 px)", 16.0)]

#: The point size a fresh install starts at. `Automatic` -- follow whatever the point style says --
#: was the old default and it made the two size controls argue: the style said 5 px, the menu said
#: "Automatic", and neither told a reader which number was in force. A stated size is one number, and
#: 7 px is the one that reads at 8,140 genes without the cloud closing up.
DEFAULT_POINT_SIZE = 7.0

#: Help shared by the display menu and its preference controls.
DISPLAY_HELP = {
    "Theme": "Choose interface colours. Paper gives a light background suitable for exported figures.",
    "Colour map": "Choose how numerical values or categories are coloured. Auto matches the data type.",
    "Point style": "Change point size and opacity. An explicit point size in View overrides the style size.",
    "Overlap": "Control how overlapping genes are drawn. Additive blending emphasizes density but can obscure colours.",
    "Point render mode": "Choose flat discs or shaded spheres. This changes appearance without changing the analysis.",
    "Light render mode": "Soft lighting adds depth; ray traced lighting adds shadows from dense regions and costs more rendering time.",
    "Light target": "Aim the light with the cursor or the selected gene. Lighting does not select genes or alter measurements.",
    "Pointer beam": "Change the width and softness of cursor-controlled lighting. Gene-centred lights are unaffected.",
    "Pointer response": "Choose how quickly the light follows the cursor. Smooth and cinematic add a delay.",
    "Target marker": "Show a visual marker at the light target. The marker represents no biological measurement.",
    "Light mood": "Choose neutral, cool, or warm lighting. Neutral preserves category colours most directly.",
    "Background": "Show or hide the animated panel background. Turn it off for a simpler display.",
    "Window size": "Choose the initial window size. Sizes larger than the monitor are reduced to fit.",
    "Full screen": "Use the whole monitor. Turn this off to restore a movable window with its title bar.",
    "Fade with distance": "Fade distant genes for depth. Turn this off when comparing the exact colours of points.",
    "Reference grid": "Show a plane below the map to make rotation easier to follow. Its spacing has no biological unit.",
}

EDGE_EXPLANATION = (
    "They are kept as separate layers rather than added together into one 'interaction' edge.\n\n"
    "These layers are not measurements of the same thing. A crosslink-MS edge is a "
    "measured physical contact. A co-mention edge is two genes appearing in one abstract, which "
    "happens to popular genes far more than to related ones. A co-expression edge is a correlation "
    "across a stage series.\n\n"
    "Combining them would let the weakest inference borrow the credibility of the strongest "
    "measurement: a merged score of 0.8 cannot tell you whether two proteins were measured touching "
    "or merely mentioned together, and there is no honest weighting that recovers the difference "
    "afterwards. So they are drawn in separate colors and toggled independently, and the answer to "
    "'are these two related' is 'by which evidence'.")

MAP_EXPLANATION = (
    "Position represents similarity under the saved feature recipe, reduced to three dimensions. "
    "The packaged maps use the same balanced embedding builder as interactive maps. Their archives "
    "record the selected features, ordered genes, data hashes and algorithm that actually ran.\n\n"
    "The current default uses numeric biological features with median imputation, robust scaling "
    "and balanced blocks. Categorical compartment labels are not one-hot encoded. Numeric localization "
    "confidence can still be an input, so the display is not independent validation of localization. "
    "No literature column is an input, but uneven measurement coverage can still shape the map.\n\n"
    "Proximity suggests a question; it does not establish function, interaction or a probability. "
    "Use Tools > Predict a trait for separate held-out evaluation, calibration and unsupported-call "
    "abstention. Its models fit their own inputs and exclude the target and derived evidence. "
    "Grey means unknown, not a measured zero.")

# How many distinct values a column may have and still be offered as a filter category. Above this it
# is an identifier rather than a class -- orthogroup has 7,331 values, and a list that long is not a
# filter, it is a scrolling exercise.
MAX_CATEGORY_VALUES = 60
#: How far the interface text can be scaled. Below 0.8 the compartment counts stop being readable
#: at a glance; above 1.6 the six analysis tabs no longer fit a laptop screen.
UI_SCALE_RANGE = (0.8, 1.6)
#: How tall the caption under the cell diagram is, whatever it says. Three lines at the default
#: font: enough for "the X color is on a shape shared with Y -- click it to step through them", and
#: fixed so that a longer note scrolls instead of moving the drawing.
DIAGRAM_NOTE_HEIGHT = 54

#: The startup window size, chosen in Preferences. `MATCH_SCREEN` is the default and means the work
#: area of whichever monitor the window opens on -- the screen minus its taskbar and panels, which is
#: the largest size that actually fits. It is not the raw screen resolution: on a 1080p monitor a
#: 1920x1080 window is taller than the space available to it, and the title bar goes off the top.
#:
#: The fixed sizes are offered because this program is run on two machines with different monitors and
#: a figure made at one size should be reproducible at that size on the other. Every one of them is
#: CLAMPED to the work area before use, so choosing 3840x2160 on a 1080p screen gives a window that
#: fits rather than one that hangs off the edge.
MATCH_SCREEN = "match screen"
WINDOW_SIZES = (MATCH_SCREEN, "1280 × 720", "1366 × 768", "1600 × 900", "1920 × 1080",
                "2560 × 1440", "3200 × 1800", "3840 × 2160")
# Columns that are categorical by dtype but meaningless as a filter: provenance flags and internals.
CATEGORY_DENYLIST = {"gene_id", "product", "symbol", "orthogroup"}


def as_text(s):
    """A plain string Series with missing values as "".

    `.astype(str)` is not enough. Under pandas 3 a string column has dtype `str` and keeps NA through
    the cast, so the result mixes `str` with `nan`: sorting it raises TypeError, and the NaNs compare
    unequal to everything including themselves. Under pandas 2 the same cast produced the literal
    "nan". Neither is what the caller wants, so missing becomes "", which ABSENCE already treats as
    absence -- the same convention the held-out search uses.
    """
    return s.astype("object").where(s.notna(), "").astype(str)


#: Prefix marking a color source that is a kept clustering rather than a column.
RUN_PREFIX = "clustering: "
#: Prefix marking a numeric column shown as bins rather than as a ramp.
BIN_PREFIX = "binned: "


def numeric_columns(nodes) -> list[str]:
    """Every column that is a quantity worth binning, best first.

    A quantity binned behaves like a category, which is what makes it comparable with a clustering
    -- the comparison the color-by panel exists for. Identifiers and flags are left out: binning
    `gene_id` is not a question anybody has.
    """
    from pandas.api import types as pdt
    out = []
    for c in nodes.columns:
        if c in CATEGORY_DENYLIST:
            continue
        s = nodes[c]
        if not pdt.is_numeric_dtype(s) or pdt.is_bool_dtype(s):
            continue
        if s.notna().sum() < 30 or s.nunique(dropna=True) < 5:
            continue
        out.append(c)
    front = [c for c in ("fit_invitro_hff", "n_publications", "mean_plddt", "expr_tachy",
                         "expr_cyst") if c in out]
    return front + sorted(c for c in out if c not in front)


def category_columns(nodes) -> list[str]:
    """Every column that can serve as a filter class, best first.

    The compartment field used to be the only one, hardcoded, which meant the map could be filtered by
    the single worst-recovered property in it and by nothing else. Any column with a small number of
    repeated values is a legitimate way to slice the proteome, so the panel offers all of them and
    lets the user decide which question they are asking.
    """
    from pandas.api import types as pdt
    out = []
    for c in nodes.columns:
        if c in CATEGORY_DENYLIST:
            continue
        s = nodes[c]
        # Asked as "is it NOT a quantity", not as "is the dtype object". pandas 3 gives string columns
        # dtype `str` where pandas 2 gave `object`, so testing for object found nothing but the one
        # boolean column -- every filter category silently disappeared in a newer pandas while every
        # test went on passing against the older one.
        if pdt.is_numeric_dtype(s) and not pdt.is_bool_dtype(s):
            continue
        if pdt.is_datetime64_any_dtype(s) or pdt.is_timedelta64_dtype(s):
            continue
        n = s.nunique(dropna=True)
        if 2 <= n <= MAX_CATEGORY_VALUES:
            out.append(c)
    # compartment first because it is the historical default and what most users arrive looking for,
    # then the other measured classes, then everything else alphabetically.
    front = [c for c in ("compartment", "compartment_best", "cellcycle_phase",
                         "stage_enriched_derived", "lit_tier", "attention_depth") if c in out]
    return front + sorted(c for c in out if c not in front)

# Qualitative palette; "unassigned" is deliberately grey, because a missing hyperLOPIT call means
# unknown (assignment tracks abundance) and must not read as a 27th compartment.
# Kept as a fallback: the categorical color map chosen in the UI supersedes it.
PALETTE = [
    (0.90, 0.24, 0.24), (0.20, 0.55, 0.90), (0.25, 0.75, 0.35), (0.95, 0.65, 0.15),
    (0.65, 0.35, 0.85), (0.15, 0.80, 0.78), (0.95, 0.45, 0.70), (0.55, 0.75, 0.20),
    (0.85, 0.35, 0.10), (0.35, 0.45, 0.85), (0.10, 0.65, 0.50), (0.80, 0.80, 0.20),
    (0.60, 0.20, 0.45), (0.30, 0.70, 0.95), (0.75, 0.55, 0.30), (0.45, 0.35, 0.70),
    (0.95, 0.80, 0.45), (0.20, 0.40, 0.35), (0.85, 0.55, 0.55), (0.40, 0.60, 0.50),
    (0.70, 0.70, 0.90), (0.55, 0.45, 0.20), (0.30, 0.85, 0.60), (0.90, 0.40, 0.45),
    (0.50, 0.50, 0.95), (0.65, 0.85, 0.35),
]
GREY = (0.45, 0.45, 0.48)   # fallback; the live value comes from TH.unknown_color(theme)

# Roughly how many edges can be drawn at full alpha before they stop being separable lines and become a
# filled region. Beyond it, alpha is scaled down rather than edges being dropped -- density stays
# visible as brightness instead of being silently truncated.
EDGE_INK_TARGET = 1500

# Depth of attention is categorical (see literature.DEPTH_OF): named in a title / in an abstract / only in
# a body or caption. Distinct hues rather than a ramp, because the tiers are not a measured quantity.
DEPTH_COLOR = {"focal": (0.98, 0.86, 0.30),          # the paper is about this gene
                "substantive": (0.35, 0.70, 0.95),    # a stated part of the paper's claims
                "incidental": (0.55, 0.35, 0.60)}     # mentioned in passing, or listed in a table


#: The parasite arms this window can open, and the cache each lives in. Instruction 39's rule is that
#: nothing is merged: one table per species, and anything crossing a boundary crosses through a bridge
#: that can be inspected and disbelieved. So this is a CHOICE of table, never a union -- one combined
#: node table would put two id spaces in one index and make "cluster 5 is 71% IMC" a claim about a
#: mixture of organisms.
#:
#: Until this existed, `load` opened `nodes.parquet` by name, and the entire Plasmodium arm -- node
#: table, graph, host bridge, 41 filled slots -- was data the browser could not open.
SPECIES = {
    "Toxoplasma gondii": {"nodes": "nodes.parquet", "graph": "graph.npz", "code": "Tg"},
    "Plasmodium falciparum": {"nodes": "pf_nodes.parquet", "graph": "pf_graph.npz", "code": "Pf"},
}
DEFAULT_SPECIES = "Toxoplasma gondii"


def available_species() -> list:
    """The species whose cache is actually present, in menu order.

    Checked rather than assumed: a checkout that has never built the Plasmodium arm should offer one
    species, not offer two and fail on the second.
    """
    return [name for name, where in SPECIES.items()
            if os.path.exists(os.path.join(DATA, where["nodes"]))
            and os.path.exists(os.path.join(DATA, where["graph"]))]


def load(species: str = DEFAULT_SPECIES):
    """Load one species' built cache: node table, coordinates, edge layers and stored models.

    Raises with the command that builds the cache if it is absent, rather than letting pandas raise
    a file-not-found from inside a constructor, which tells the user nothing about what to do.
    """
    where = SPECIES.get(species) or SPECIES[DEFAULT_SPECIES]
    npz, pq = os.path.join(DATA, where["graph"]), os.path.join(DATA, where["nodes"])
    if not (os.path.exists(npz) and os.path.exists(pq)):
        raise SystemExit("No cached graph. Run:  python -m starplast.build_graph")
    nodes = pd.read_parquet(pq)
    if where["code"] == "Tg":
        from .structure_catalog import attach_features
        nodes = attach_features(nodes, os.path.join(DATA, "af3_features.parquet"))
    z = np.load(npz)
    if "gene_ids" in z.files and not np.array_equal(z["gene_ids"].astype(str), nodes.gene_id.astype(str)):
        raise ValueError("graph gene order does not match the node table; rebuild the graph cache")
    xyz = z["xyz"].astype(np.float32)
    edges = {}
    for k, _ in EDGE_TYPES:
        if f"{k}__a" in z.files:
            # `r` is the attention-corrected residual, and only the co-mention layers have one --
            # every other layer's weight IS its strength. The Toxoplasma builder happens to write `r`
            # for all of them, so indexing it unconditionally worked there and raised KeyError on the
            # first graph built by anything else, which is what the Plasmodium arm is.
            edges[k] = {"a": z[f"{k}__a"], "b": z[f"{k}__b"], "w": z[f"{k}__w"],
                        "r": z[f"{k}__r"] if f"{k}__r" in z.files else z[f"{k}__w"]}
    # Optional: the crosslink model table, which says how a measured interaction is thought to happen.
    mp = os.path.join(DATA, "crosslink_models.parquet")
    models = pd.read_parquet(mp) if os.path.exists(mp) else pd.DataFrame()
    return nodes, xyz, edges, models


#: What the left mouse button does. Two modes rather than a held modifier, because each is done in
#: long stretches -- navigate for a while, then select for a while -- and a key held down for a
#: minute is a worse control than a mode set once.
INTERACTION_MODES = ("navigate", "select")
#: Rotation constraints. A free orbit never returns to the same view twice, which is exactly wrong
#: for comparing two maps; constrained to an axis, the view is one number you can come back to.
NAVIGATE_AXES = ("free", "x", "y", "z")
#: The two gate shapes. They answer different questions -- see `Map3D.finish_gate`.
GATE_SHAPES = ("lasso (2D)", "brush (3D)")


class _GateOverlay(QtWidgets.QWidget):
    """The gate being drawn, painted over the GL view.

    A transparent child widget rather than painting inside `paintGL`: mixing QPainter into a
    QOpenGLWidget's GL painting is driver-dependent, and this project already has one class of bug
    that appears only on someone else's machine. The overlay takes no mouse events, so the view
    underneath behaves exactly as it did.
    """

    def __init__(self, parent):
        """Create a transparent overlay for drawing selection outlines over the map."""
        super().__init__(parent)
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_TransparentForMouseEvents)
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_NoSystemBackground)
        self.setAttribute(QtCore.Qt.WidgetAttribute.WA_TranslucentBackground)
        self.points: list = []          # widget-space points of a lasso
        self.circle = None              # (x, y, radius in pixels) of a brush

    def paintEvent(self, ev):
        """Draw whichever gate is in progress. Nothing at all when none is."""
        if not self.points and self.circle is None:
            return
        p = QtGui.QPainter(self)
        self.draw(p)
        p.end()

    def draw(self, p: QtGui.QPainter):
        """The drawing itself, onto any painter.

        Separated from `paintEvent` so it can be checked without a window: rendering a widget into
        an image brings the widget background with it, and a test that cannot tell "drew a lasso"
        from "filled the rectangle" is not checking the lasso.
        """
        p.setRenderHint(QtGui.QPainter.RenderHint.Antialiasing, True)
        pen = QtGui.QPen(QtGui.QColor(255, 255, 255, 220))
        pen.setWidthF(1.6)
        p.setPen(pen)
        p.setBrush(QtGui.QBrush(QtGui.QColor(255, 255, 255, 28)))
        if self.circle is not None:
            x, y, r = self.circle
            p.drawEllipse(QtCore.QPointF(x, y), r, r)
        elif len(self.points) > 1:
            p.drawPolygon(QtGui.QPolygonF([QtCore.QPointF(*xy) for xy in self.points]))

    def show_lasso(self, points):
        """Show a lasso through these widget-space points."""
        self.points, self.circle = list(points), None
        self.update()

    def show_brush(self, x, y, r):
        """Show a brush of radius `r` pixels centered at (x, y)."""
        self.points, self.circle = [], (float(x), float(y), float(r))
        self.update()

    def clear(self):
        """Remove whatever was being drawn."""
        self.points, self.circle = [], None
        self.update()


def inside_polygon(x, y, poly) -> np.ndarray:
    """Which of the points (x, y) fall inside a polygon. Ray casting, vectorised over the points.

    Written out rather than taken from matplotlib: it runs over 8,140 points as the mouse moves, and
    the dependency is not otherwise needed at run time. Points with NaN coordinates -- genes behind
    the camera -- are outside by construction, which is the answer that matches what is on screen.
    """
    x, y = np.asarray(x, dtype=float), np.asarray(y, dtype=float)
    inside = np.zeros(x.shape, dtype=bool)
    poly = [(float(a), float(b)) for a, b in poly]
    if len(poly) < 3:
        return inside
    ok = np.isfinite(x) & np.isfinite(y)
    xs = np.where(ok, x, 0.0)
    ys = np.where(ok, y, 0.0)
    n = len(poly)
    for i in range(n):
        x1, y1 = poly[i]
        x2, y2 = poly[(i + 1) % n]
        # Points whose rightward ray crosses this edge. The `y2 != y1` guard is the horizontal-edge
        # case, which would otherwise divide by zero and count a whole row of points as crossings.
        straddles = (y1 > ys) != (y2 > ys)
        with np.errstate(divide="ignore", invalid="ignore"):
            xint = (x2 - x1) * (ys - y1) / np.where(y2 != y1, y2 - y1, np.nan) + x1
        inside ^= straddles & np.isfinite(xint) & (xs < xint)
    return inside & ok


class Map3D(gl.GLViewWidget):
    """GLViewWidget plus click-picking and gated selection, neither of which pyqtgraph provides."""

    picked = QtCore.pyqtSignal(int)
    #: The genes inside a gate the user drew, as an index array. An empty gate is a real answer --
    #: it selected nothing -- and is emitted rather than swallowed.
    gated = QtCore.pyqtSignal(object)

    def __init__(self, xyz):
        """Create a navigable OpenGL view for an (n_genes, 3) coordinate array."""
        super().__init__()
        self.xyz = xyz
        self._proj_kind = None      # resolved once, see _projection_matrix
        # Which points may be picked. None means all of them; a mask is set when the displayed
        # embedding covers only some of the genes, so a hidden gene cannot be selected by clicking
        # where it would have been.
        self.pickable = None
        #: The pointer, as a direction in this widget's frame. None until the mouse has been in it.
        self.pointer = None
        self.pointer_px = None
        # Without this, Qt delivers a move event only while a BUTTON IS HELD. So the light that
        # follows the pointer -- the default source -- moved only while the map was being dragged,
        # and the drag is also what orbits the camera: the light appeared to be stuck to the cloud
        # rather than to the pointer, and sat frozen the rest of the time. Reported as "I can't see
        # the mouse light". pyqtgraph's own handler tests ev.buttons() and ignores a hover, so
        # tracking costs nothing beyond the events.
        self.setMouseTracking(True)
        self.mode = "navigate"
        self.axis = "free"
        self.gate_shape = GATE_SHAPES[0]
        self._gate = None               # points of the lasso, or (anchor index, radius) of a brush
        self.overlay = _GateOverlay(self)
        self.overlay.setGeometry(self.rect())
        self.fit_view()

    def resizeEvent(self, ev):
        """Keep the gate overlay the size of the view it is drawn on."""
        super().resizeEvent(ev)
        self.overlay.setGeometry(self.rect())

    def fit_view(self, margin=1.35):
        """Frame the data rather than assuming a fixed distance.

        The embedding is renormalized on every build and its extent changes with the feature set, so a
        hard-coded 170 left the map as a small island in a large empty viewport.
        """
        if len(self.xyz):
            centre = self.xyz.mean(0)
            radius = float(np.linalg.norm(self.xyz - centre, axis=1).max())
            self.setCameraPosition(pos=pg.Vector(*centre), distance=max(radius * margin, 10.0))
        else:
            self.setCameraPosition(distance=170)

    def animate_distance(self, target: float, ms: int = 420):
        """Ease the camera to a new distance instead of cutting to it.

        A cut between levels of detail loses the viewer's place: the whole picture changes at once and
        nothing carries over, so the eye has to re-find the region it was looking at. Easing keeps the
        correspondence between what was on screen and what is now on screen.
        """
        start = float(self.opts.get("distance", target))
        if abs(target - start) < 1e-6:
            return
        if getattr(self, "_cam_timer", None) is not None:
            self._cam_timer.stop()                # a second change mid-flight retargets, never queues
        steps = max(int(ms / 16), 1)
        self._cam_step = 0

        def tick():
            """Advance the camera animation using smooth interpolation and stop at its target."""
            self._cam_step += 1
            f = min(self._cam_step / steps, 1.0)
            f = f * f * (3.0 - 2.0 * f)           # smoothstep: no jerk at either end
            self.setCameraPosition(distance=start + (target - start) * f)
            if f >= 1.0:
                self._cam_timer.stop()

        self._cam_timer = QtCore.QTimer(self)
        self._cam_timer.timeout.connect(tick)
        self._cam_timer.start(16)

    def data_radius(self) -> float:
        """Radius of the point cloud about its center, for framing and grid sizing."""
        if not len(self.xyz):
            return 100.0
        return float(np.linalg.norm(self.xyz - self.xyz.mean(0), axis=1).max())

    def _projection_matrix(self):
        """`projectionMatrix()` across pyqtgraph versions.

        The signature has changed twice: no arguments, then `(region=None)`, and in current releases
        `(region, viewport)` with both required. Calling the old way on a new pyqtgraph raises TypeError
        inside the mouse handler, which is why clicking a gene printed a traceback and selected nothing.

        The convention is decided by INSPECTING the signature, not by catching TypeError. A try/except
        chain cannot tell "you called me wrongly" from "something inside me raised TypeError", so it
        would silently retry a call that failed for an unrelated reason and report the wrong cause.
        """
        fn = self.projectionMatrix
        if self._proj_kind is None:
            try:
                params = [p for name, p in inspect.signature(fn).parameters.items()
                          if p.kind in (p.POSITIONAL_OR_KEYWORD, p.KEYWORD_ONLY)]
                required = [p for p in params if p.default is inspect.Parameter.empty]
                names = {p.name for p in params}
                if "viewport" in names and any(p.name == "viewport" for p in required):
                    self._proj_kind = "region_viewport"
                elif "region" in names:
                    self._proj_kind = "region"
                else:
                    self._proj_kind = "none"
            except (TypeError, ValueError):          # C extension without introspection
                self._proj_kind = "none"

        dpr = self.devicePixelRatioF() if hasattr(self, "devicePixelRatioF") else 1.0
        # Device pixels, not logical ones: on a HiDPI display a logical viewport puts the picking ray
        # in the wrong place, which reads as "clicking is slightly off" rather than as a failure.
        rect = (0, 0, int(self.width() * dpr), int(self.height() * dpr))
        if self._proj_kind == "region_viewport":
            # BOTH arguments are indexed by pyqtgraph 0.14 -- `region[0]`, `viewport` unpacked into
            # four names -- so neither may be None. Passing None for the region raised
            # "'NoneType' object is not subscriptable" from inside the mouse handler, which is exactly
            # the crash this whole compat layer exists to prevent, arriving by a different door.
            # The region for picking is the whole widget, so it IS the viewport.
            return fn(rect, rect)
        if self._proj_kind == "region":
            # Here `region=None` is the documented "use the whole viewport" default.
            return fn(None)
        return fn()

    def _mvp(self):
        """Return the 4 × 4 world-to-clip projection matrix used for point picking."""
        m = self._projection_matrix() * self.viewMatrix()
        return np.array([[r.x(), r.y(), r.z(), r.w()]
                         for r in (m.row(i) for i in range(4))], dtype=float)

    def project(self):
        """Node positions in widget (logical) pixels; NaN where behind the camera."""
        mvp = self._mvp()
        h = np.hstack([self.xyz, np.ones((len(self.xyz), 1), dtype=float)])
        p = h @ mvp.T
        w = p[:, 3].copy()
        w[np.abs(w) < 1e-9] = np.nan
        ndc = p[:, :3] / w[:, None]
        sx = (ndc[:, 0] + 1.0) * 0.5 * self.width()
        sy = (1.0 - ndc[:, 1]) * 0.5 * self.height()
        sx[w <= 0] = np.nan
        return sx, sy

    # ------------------------------------------------------------------ gating
    def begin_gate(self, x: float, y: float):
        """Start a gate at a widget-space point. The shape decides what that means.

        Split from the mouse handler so the whole gating path can be driven -- and tested -- without
        a window manager. Every display claim in this project has to be checkable headlessly, and a
        gate that returns the wrong genes is a display claim.
        """
        if self.gate_shape.startswith("brush"):
            i = self.nearest(x, y)
            self._gate = ("brush", i, x, y, 0.0)
            if i is not None:
                self.overlay.show_brush(x, y, 0.0)
        else:
            self._gate = ("lasso", [(x, y)])
            self.overlay.show_lasso([(x, y)])

    def extend_gate(self, x: float, y: float):
        """Continue the gate to a new point."""
        if not self._gate:
            return
        if self._gate[0] == "brush":
            _, i, cx, cy, _ = self._gate
            r = float(np.hypot(x - cx, y - cy))
            self._gate = ("brush", i, cx, cy, r)
            self.overlay.show_brush(cx, cy, r)
        else:
            pts = self._gate[1]
            # Sub-pixel moves are dropped: a freehand lasso otherwise accumulates thousands of
            # points, and the polygon test is linear in them.
            if not pts or np.hypot(x - pts[-1][0], y - pts[-1][1]) >= 2.0:
                pts.append((x, y))
                self.overlay.show_lasso(pts)

    def finish_gate(self):
        """Close the gate, emit the genes inside it, and return them.

        **2D lasso**: a polygon in screen space, so it takes everything behind it as well. That is
        what "grab that visual cluster" means -- the cluster you can see is a screen-space object.

        **3D brush**: a ball in world space around the gene under the press, with the drag setting
        its radius. Deep in the cloud a lasso also catches the far side, which looks like a
        selection of one structure and is a selection of two; a world-space ball cannot.
        """
        gate, self._gate = self._gate, None
        self.overlay.clear()
        if not gate:
            return np.array([], dtype=int)
        try:
            sx, sy = self.project()
        except Exception as exc:            # a mouse handler is the wrong place to raise
            print(f"starplast: gating unavailable ({type(exc).__name__}: {exc})")
            self.gated.emit(np.array([], dtype=int))
            return np.array([], dtype=int)
        if gate[0] == "brush":
            _, i, cx, cy, r_px = gate
            hit = np.zeros(len(self.xyz), dtype=bool)
            if i is not None and r_px > 2:
                # Pixels to world units, measured on this view rather than assumed: the conversion
                # depends on the camera distance, and a fixed factor would make the brush grow and
                # shrink as you zoom while the circle on screen did not.
                scale = self._world_per_pixel(i, cx, cy, sx, sy)
                d = np.linalg.norm(self.xyz - self.xyz[i], axis=1)
                hit = d <= r_px * scale
        else:
            hit = inside_polygon(sx, sy, gate[1])
        if self.pickable is not None and len(self.pickable) == len(hit):
            hit &= self.pickable
        idx = np.flatnonzero(hit)
        self.gated.emit(idx)
        return idx

    def _world_per_pixel(self, i: int, cx: float, cy: float, sx, sy) -> float:
        """World units per screen pixel near gene `i`, measured from the projection itself.

        Taken from the spread of the genes actually near that point on screen rather than from the
        camera parameters, so it holds whatever projection convention pyqtgraph is using this
        release -- the same reason `_projection_matrix` inspects rather than assumes.
        """
        near = np.flatnonzero(np.isfinite(sx) & (np.hypot(sx - cx, sy - cy) < 120))
        if len(near) < 8:
            return float(self.data_radius() / max(self.width(), 1)) * 2.0
        d_px = np.hypot(sx[near] - cx, sy[near] - cy)
        d_world = np.linalg.norm(self.xyz[near] - self.xyz[i], axis=1)
        ok = d_px > 1e-6
        return float(np.median(d_world[ok] / d_px[ok])) if ok.any() else 1.0

    def nearest(self, x: float, y: float, within: float = 1e9):
        """The gene nearest a widget point, or None. Shared by picking and by the brush anchor."""
        try:
            sx, sy = self.project()
        except Exception as exc:
            print(f"starplast: picking unavailable ({type(exc).__name__}: {exc})")
            return None
        d = np.hypot(sx - x, sy - y)
        if self.pickable is not None and len(self.pickable) == len(d):
            d = np.where(self.pickable, d, np.nan)
        if np.all(np.isnan(d)):
            return None
        i = int(np.nanargmin(d))
        return i if d[i] < within else None

    def under_pointer(self):
        """Where the gene under the pointer IS, in world coordinates, or None.

        Deliberately the nearest gene on screen rather than the nearest along a ray: what a reader
        means by "that cluster" is the one they can see under the cursor, and the two answers differ
        only where a nearer gene projects closer to the pointer than the one being looked at.
        """
        if self.pointer_px is None or not len(self.xyz):
            return None
        i = self.nearest(*self.pointer_px, within=max(self.width(), self.height()) * 0.25)
        return None if i is None else np.asarray(self.xyz[i], dtype=float)

    # ------------------------------------------------------------------ mouse
    def mousePressEvent(self, ev):
        """Start a selection gate or pass navigation input to the OpenGL view."""
        if self.mode == "select" and ev.button() == QtCore.Qt.MouseButton.LeftButton:
            p = ev.position()
            self.begin_gate(p.x(), p.y())
            return
        super().mousePressEvent(ev)

    def camera_basis(self):
        """(eye, right, up, forward) in world coordinates, read off the camera.

        Everything screen-relative needs this. A light asked for "top left" or "where the pointer
        is" is a statement about the SCREEN, and putting it at a fixed world position instead means
        it stops being top-left the moment the map is orbited -- which is exactly what "the corner
        lights do not seem to move" was.
        """
        eye = np.array(self.cameraPosition(), dtype=float)
        m = self.viewMatrix()
        # A view matrix's rows are the camera's axes in world space: row 0 right, row 1 up, row 2
        # points AWAY from what is being looked at, so forward is its negative.
        right = np.array([m.row(0).x(), m.row(0).y(), m.row(0).z()], dtype=float)
        up = np.array([m.row(1).x(), m.row(1).y(), m.row(1).z()], dtype=float)
        forward = -np.array([m.row(2).x(), m.row(2).y(), m.row(2).z()], dtype=float)
        return eye, right, up, forward

    def mouseMoveEvent(self, ev):
        # Where the pointer is, as a direction in the view's own frame: x right, y UP (Qt counts
        # down), z toward the eye. Kept here rather than computed from a world position because the
        # map rotates -- a light fixed in world space swings away from the pointer the moment
        # anything moves, which reads as the light being broken.
        """Update the mouse light and active gate, or rotate the map in navigation mode."""
        pos = ev.position()
        w, h = max(self.width(), 1), max(self.height(), 1)
        # Clamped to the widget, because a drag that leaves it keeps delivering moves: without this
        # the pointer reads several widths out, and the light it places goes with it -- off the side
        # of the map, lighting nothing, until the mouse comes back.
        clamp = lambda v: float(min(max(v, -1.0), 1.0))
        self.pointer = (clamp(2.0 * pos.x() / w - 1.0), clamp(1.0 - 2.0 * pos.y() / h), 1.0)
        # And in pixels, because the gene under the pointer is found by projecting the map onto the
        # widget: the light that follows the pointer needs the DEPTH of what is being pointed at,
        # and normalised coordinates have thrown that away.
        self.pointer_px = (float(pos.x()), float(pos.y()))
        if self.mode == "select" and self._gate:
            p = ev.position()
            self.extend_gate(p.x(), p.y())
            return
        if self.axis == "free" or ev.buttons() != QtCore.Qt.MouseButton.LeftButton:
            super().mouseMoveEvent(ev)
            return
        # Constrained orbit. The camera is spherical, so each axis is one of its two angles held
        # still: about z is azimuth alone, about x or y is elevation with azimuth pinned to that
        # axis. Two maps compared from "x" are compared from the same place.
        pos = ev.position()
        last = getattr(self, "mousePos", pos)
        self.mousePos = pos
        diff = pos - last
        if self.axis == "z":
            self.orbit(-diff.x(), 0)
        else:
            self.setCameraPosition(azimuth=0.0 if self.axis == "x" else 90.0)
            self.orbit(0, diff.y())

    def mouseReleaseEvent(self, ev):
        """Complete a selection gate or select the gene under a navigation click."""
        if self.mode == "select" and ev.button() == QtCore.Qt.MouseButton.LeftButton:
            self.finish_gate()
            return
        super().mouseReleaseEvent(ev)
        if ev.button() != QtCore.Qt.MouseButton.LeftButton:
            return
        p = ev.position()
        try:
            sx, sy = self.project()
        except Exception as exc:            # a mouse handler is the wrong place to raise
            print(f"starplast: picking unavailable ({type(exc).__name__}: {exc})")
            return
        d = np.hypot(sx - p.x(), sy - p.y())
        if self.pickable is not None and len(self.pickable) == len(d):
            # NaN rather than a large number: nanargmin ignores it, and there is no distance at
            # which an unplaced gene should win.
            d = np.where(self.pickable, d, np.nan)
        if np.all(np.isnan(d)):
            return
        i = int(np.nanargmin(d))
        if d[i] < 14:
            self.picked.emit(i)


class Window(QtWidgets.QMainWindow):
    """The application window: a 3D map, the panels around it, and everything they can do.

    Holds the view state that the menus mutate -- level of detail, coloring, point size, which edge
    layers are active, the filter category -- and redraws from it. State lives here as plain
    attributes rather than being read back out of widgets, so a headless caller can drive the whole
    interface without a window manager, which is how the tests exercise it.
    """
    def __init__(self, species: str = None):
        """Load one organism and construct its map, evidence docks, menus, and analysis panel."""
        super().__init__()
        self.theme = 'dark'
        self.point_style = TH.DEFAULT_POINT_STYLE
        self.point_mode = 'occlude'
        self._spin_speed = 0.35
        self.depth_cue = True
        self.show_ground = True
        self.cmap_name = None
        self.apply_log_settings()
        # Which parasite table this window is showing. One species per window, never a union -- see
        # `SPECIES`. Remembered, so the arm someone works in is the one that opens next time.
        stored = str(self.settings().value("data/species", DEFAULT_SPECIES) or DEFAULT_SPECIES)
        offered = available_species() or [DEFAULT_SPECIES]
        self.species = species or (stored if stored in offered else offered[0])
        self.nodes, self.xyz, self.edges, self.models = load(self.species)
        self.n = len(self.nodes)
        self.sel = None
        self.setWindowTitle(f"starplast — {self.species} knowledge map")
        self.setWindowIcon(QtGui.QIcon(os.path.join(os.path.dirname(__file__), "data", "icons", "starplast.svg")))
        # Was a flat resize(1580, 950), which is not the reason the window used to open too large --
        # that was a 2,897 px minimum height, fixed in `analysis_panel._scrolled` -- but it did mean
        # the window was the same size on a 1080p monitor and a 4K one. Applied last, after the size
        # is known to be reachable.
        self.apply_window_size()

        from .runs import RunStore
        self.categories = category_columns(self.nodes)
        self.numerics = numeric_columns(self.nodes)
        #: How many bins a numeric color source is cut into. Quantile bins, so the choice is about
        #: how fine a distinction to draw rather than about the shape of the distribution.
        self.bins = 5
        #: Kept clusterings. Beside the cache rather than in it: these are the user's runs and must
        #: survive a rebuild of the data.
        # Display settings, read once from where they were left. The lighting timer exists whether
        # or not it is running: creating it lazily inside a handler is how a timer ends up owned by
        # whichever thread happened to touch it first.
        s = QtCore.QSettings("starplast", "starplast")
        from . import ambient as _ambient, lighting as _lighting
        self._ambient_mode = str(s.value("display/ambient", "none"))
        self._ambient = {
            "speed": float(s.value("display/ambient_speed", _ambient.DEFAULT_SPEED, type=float)),
            "size": float(s.value("display/ambient_size", _ambient.DEFAULT_SIZE, type=float)),
            "density": float(s.value("display/ambient_density", _ambient.DEFAULT_DENSITY,
                                     type=float)),
        }
        self._ambient_widget = None
        old_mode = str(s.value("display/lighting", _lighting.DEFAULT_MODE))
        mode = "soft" if old_mode == "lit" else old_mode
        old_points = s.value("display/light_point_mode",
                             s.value("display/light_finish", _lighting.DEFAULT_POINT_MODE))
        self._lighting = {
            "mode": mode if mode in _lighting.MODES else "off",
            "source": _lighting.normalize_source(s.value("display/light_source",
                                                          _lighting.DEFAULT_SOURCE)),
            "point_mode": _lighting.normalize_point_mode(old_points),
            "mood": (str(s.value("display/light_mood", _lighting.DEFAULT_MOOD))
                     if str(s.value("display/light_mood", _lighting.DEFAULT_MOOD))
                     in _lighting.LIGHT_MOODS else _lighting.DEFAULT_MOOD),
            "pointer_mode": (str(s.value("display/light_pointer_mode",
                                         _lighting.DEFAULT_POINTER_MODE))
                             if str(s.value("display/light_pointer_mode",
                                            _lighting.DEFAULT_POINTER_MODE))
                             in _lighting.POINTER_MODES else _lighting.DEFAULT_POINTER_MODE),
            "response": (str(s.value("display/light_response", _lighting.DEFAULT_RESPONSE))
                         if str(s.value("display/light_response", _lighting.DEFAULT_RESPONSE))
                         in _lighting.RESPONSES else _lighting.DEFAULT_RESPONSE),
            "target_marker": (str(s.value("display/light_target_marker",
                                          _lighting.DEFAULT_TARGET_MARKER))
                              if str(s.value("display/light_target_marker",
                                             _lighting.DEFAULT_TARGET_MARKER))
                              in _lighting.TARGET_MARKERS else _lighting.DEFAULT_TARGET_MARKER),
        }
        #: How opaque the panels are over the drifting background. 1.0 is the old look.
        self._container_opacity = float(s.value("display/container_opacity", 1.0, type=float))
        self._light_t = 0.0
        #: Text size for the whole interface, as a multiplier on the font this desktop asked for.
        #: Set on the application rather than on each widget: every layout then measures its own
        #: contents at the new size, which is what keeps a longer label from being clipped instead
        #: of merely smaller.
        self._ui_scale = float(s.value("display/ui_scale", 1.0, type=float))
        #: The desktop's own font, kept so scaling is always applied to it rather than to whatever
        #: the last scaling produced -- compounding 1.2 three times is 1.7, and the text creeps.
        app = QtWidgets.QApplication.instance()
        self._base_font = QtGui.QFont(app.font()) if app is not None else QtGui.QFont()
        self._base_colors = self._base_sizes = None
        self._sprite_state = None          # (point render mode, light direction) for the sprite
        self._grid = None                     # where the map is solid, for shadows and bounces
        self._grid_for = None                 # the coordinates that grid was built from
        self._smoothed_pointer = None         # continuous screen ray, eased without gene snapping
        self._light_timer = QtCore.QTimer(self)
        self._light_timer.timeout.connect(self._light_tick)
        self.runs = RunStore(os.path.join(paths.user_cache_dir(), "runs"))
        self.runs.load_all()
        self.category = "compartment" if "compartment" in self.categories else self.categories[0]
        comps = sorted(as_text(self.nodes[self.category]).unique())
        absent = {str(x).lower() for x in ABSENCE}
        self.comps = ([c for c in comps if c.lower() not in absent]
                      + [c for c in comps if c.lower() in absent])
        self.color_of = dict(zip(self.comps,
                                  TH.categorical_colors(len(self.comps), self.theme)))
        for c in self.comps:
            if c.lower() in absent:
                self.color_of[c] = GREY

        # Coloring and edge state live in the menus, so they are plain
        # attributes with menu actions over them rather than widgets read out of a side panel.
        self.color_mode = "compartment"
        self.point_size = DEFAULT_POINT_SIZE   # None would mean "whatever the point style says"
        self.edge_on = {k: (k in ("comention", "cofitness") and k in self.edges)
                        for k, _ in EDGE_TYPES}
        self.attn_on = True                   # a correctness default, not a preference
        self.all_edges_on = False
        # A clustering from the analysis panel, so the map can be colored by it. Looking at
        # structure beside a held-out variable is what this application is for, and until this
        # existed the Clusters tab computed labels and discarded them.
        self.cluster_labels = None
        # Which genes the displayed embedding has coordinates for. None means "the built cache",
        # which covers all of them; a mask arrives with any map built over a subsample.
        self.placed = None
        # A gated set of genes: the natural input to annotation, and the only way to select a group
        # without clicking 8,140 times. None means no gate has been drawn; an EMPTY array means one
        # was and it caught nothing, which is a different thing and is reported as such.
        self.gated = None
        #: Every import this session, with the choices that produced it. An imported column whose
        #: provenance is a memory of which dropdowns were set cannot be defended three weeks later.
        self.imports = []
        # Which genes carry a saved annotation. Its own mask and its own color, because an
        # annotation is a fourth thing beside measurement, inference and absence, and reading as any
        # of the three is the failure this application is built to prevent.
        self.annotated = None
        self._spin_home = None                # where spin started, so it can be put back

        self.view = Map3D(self.xyz)
        self.view.picked.connect(self.on_pick)
        self.view.gated.connect(self.on_gated)
        from . import sprite as SP
        self.scatter = SP.ShadedScatter(pos=self.xyz, size=5.0, pxMode=True)
        # GLScatterPlotItem blends additively by default, which sums the colors of overlapping points.
        # With 8,140 genes in dense UMAP clusters every mode rendered as one white blob and the color
        # encoding -- the thing the map is for -- was invisible. Translucent blending with depth testing
        # makes nearer points occlude farther ones instead of adding to them.
        self.scatter.setGLOptions("translucent")
        self.view.addItem(self.scatter)
        self.halo_item = None
        self.emitter_item = None
        self.grid_item = None
        self._depth_ctx = None
        self.edge_items = []

        self.setCentralWidget(self.view)
        self.jobs = JobRunner(self)
        self.addDockWidget(QtCore.Qt.DockWidgetArea.LeftDockWidgetArea, self._left())
        self.addDockWidget(QtCore.Qt.DockWidgetArea.RightDockWidgetArea, self._right())
        self._analysis()
        self._tool_docks()
        self.status = self.statusBar()
        self._progress()
        self._menus()
        # Right-click belongs on the thing being configured. The spin button used to be a permanent
        # control the size of a paragraph for something toggled once a session.
        self.view.setContextMenuPolicy(QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
        self.view.customContextMenuRequested.connect(self._context_menu)
        self.apply_theme(self.theme)
        self.set_ui_scale(self._ui_scale)
        TH.wrap_tooltips(self)
        self.redraw()
        if self._ambient_mode != "none":
            self._apply_ambient()
        if self._lighting["mode"] != "off":
            self._light_timer.start(60)

    # ------------------------------------------------------------------ importing
    def import_data(self, path: str = ""):
        """Read a user's table, offer the preprocessing, and add the columns to this session."""
        from .importer import READABLE
        if not path:
            path, _ = QtWidgets.QFileDialog.getOpenFileName(
                self, "Import data", "",
                "Tables (" + " ".join(f"*{e}" for e in READABLE) + ");;All files (*)")
        if not path:
            return None
        d = self.build_import_dialog(path)
        if d is None or not d.exec():
            return None
        return self.apply_import(d)

    def build_import_dialog(self, path: str):
        """The import dialog, built but not shown.

        Split from `import_data` for the reason every dialog here is: `exec` blocks until a human
        closes it, so a test that called it would hang rather than fail.
        """
        from .importer import (DUPLICATES, QUANTIFICATIONS, describe_columns, numeric_columns,
                               read_any, sheets_in, suggest_quantification)
        from .embedding import NA_POLICIES, SCALINGS
        from .tuning import suggest_gene_column
        try:
            sheets = sheets_in(path)
            table = read_any(path, sheet=sheets[0] if sheets else 0)
        except Exception as exc:
            QtWidgets.QMessageBox.warning(self, "Import", f"{path} could not be read:\n{exc}")
            return None

        d = QtWidgets.QDialog(self)
        d.setWindowTitle(f"Import {os.path.basename(path)}")
        d.table, d.path = table, path
        lay = QtWidgets.QVBoxLayout(d)
        form = QtWidgets.QFormLayout()

        d.sheet = QtWidgets.QComboBox()
        d.sheet.addItems([str(s) for s in sheets] or ["(single table)"])
        d.sheet.setEnabled(bool(sheets))
        d.sheet.setToolTip("Which sheet holds the data. Published supplements routinely put the "
                           "table on sheet 3, behind a legend and a blank sheet.")

        d.gene_column = QtWidgets.QComboBox()
        d.gene_column.addItems([str(c) for c in table.columns])
        guess = suggest_gene_column(table)[0]
        if guess is not None:
            d.gene_column.setCurrentText(str(guess))
        d.gene_column.setToolTip(
            "The column holding gene identifiers. Malformed ones are repaired rather than refused: "
            "`TgME49.208830` and `gene|TGME49_208830|v2` are both real formats from published "
            "supplements and both match nothing unrepaired. Previous and strain accessions are "
            "resolved forward -- one 2019 screen uses pre-2012 ids for every gene and contributed "
            "nothing at all until that was done.")

        quant, why = suggest_quantification(table, numeric_columns(table, str(guess or "")))
        d.quantification = QtWidgets.QComboBox()
        d.quantification.addItems(QUANTIFICATIONS)
        d.quantification.setCurrentText(quant)
        d.quantification.setToolTip(
            "What the numbers ARE, which decides whether a log is taken and whether each column is "
            "centred. THE RANGE DECIDES, NOT THE FILENAME: the same GEO series ships an FPKM file "
            "reaching 16,520 and derived columns still called FPKM that stop at 9.7 because they "
            "were logged upstream. Log it twice and real variation compresses to nothing; skip it "
            "and one gene dominates every distance. The table below is the evidence.")

        d.scaling = QtWidgets.QComboBox()
        d.scaling.addItems(SCALINGS)
        d.scaling.setCurrentText("rank")
        d.scaling.setToolTip(
            "rank is the safe default: published screens carry inverted sign conventions, ~64x "
            "differences in spread and heavy tails that z-scoring does not tame.")

        d.na_policy = QtWidgets.QComboBox()
        d.na_policy.addItems(NA_POLICIES)
        d.na_policy.setCurrentText("median")
        d.na_policy.setToolTip(
            "What to do about missing values. `indicator` is left to the embedding, which adds "
            "missingness as its own weighted block -- doing it here as well would count absence "
            "twice.")

        d.duplicates = QtWidgets.QComboBox()
        d.duplicates.addItems(DUPLICATES)
        d.duplicates.setToolTip(
            "How to combine several rows for one gene. A hit table often lists one row per guide or "
            "per peptide, and which of them is the claim is your call rather than this program's.")

        d.flip = QtWidgets.QCheckBox("flip the sign")
        d.flip.setToolTip(
            "Some screens are inverted relative to others -- naive-BMDM and IFN-gamma against the "
            "in vitro screen, for instance. Flip if yours disagrees, and rank-normalize before "
            "pooling it with another, because that is where the 64x spread bites.")

        d.prefix = QtWidgets.QLineEdit("imported_")
        d.prefix.setToolTip(
            "Every imported column keeps this prefix, so it cannot be mistaken for a measurement "
            "that shipped with the cache. Nothing is written to the cache either way.")

        form.addRow("sheet", d.sheet)
        form.addRow("identifiers", d.gene_column)
        form.addRow("these numbers are", d.quantification)
        form.addRow("scaling", d.scaling)
        form.addRow("missing values", d.na_policy)
        form.addRow("duplicate genes", d.duplicates)
        form.addRow("direction", d.flip)
        form.addRow("column prefix", d.prefix)
        lay.addLayout(form)

        note = QtWidgets.QLabel(f"<b>{quant}</b> — {why}")
        note.setWordWrap(True)
        lay.addWidget(note)
        d.note = note

        d.preview = QtWidgets.QTableWidget()
        d.preview.setAlternatingRowColors(True)
        d.preview.setToolTip("Each column's range, before anything is done to it. `max` is what "
                             "gives the quantification away: a column reaching 16,520 has not been "
                             "logged, one stopping at 9.7 has.")
        described = describe_columns(table, numeric_columns(table, str(guess or "")))
        d.preview.setRowCount(min(len(described), 50))
        d.preview.setColumnCount(len(described.columns))
        d.preview.setHorizontalHeaderLabels([str(c) for c in described.columns])
        for i, (_, r) in enumerate(described.head(50).iterrows()):
            for j, v in enumerate(r):
                item = QtWidgets.QTableWidgetItem(f"{v:.4g}" if isinstance(v, float) else str(v))
                d.preview.setItem(i, j, item)
        d.preview.resizeColumnsToContents()
        lay.addWidget(d.preview, 1)

        # Re-read the sheet, and re-judge the numbers, when the sheet changes: the legend sheet and
        # the data sheet are different tables and a guess made on one is nonsense about the other.
        d.sheet.currentTextChanged.connect(lambda s: self._reload_import_sheet(d, s))

        buttons = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(d.accept)
        buttons.rejected.connect(d.reject)
        lay.addWidget(buttons)
        d.resize(760, 560)
        return d

    def _reload_import_sheet(self, dialog, sheet: str):
        """Read another sheet of the same workbook into an open import dialog."""
        from .importer import read_any, suggest_quantification
        try:
            dialog.table = read_any(dialog.path, sheet=sheet)
        except Exception as exc:
            self.status.showMessage(f"sheet {sheet!r} could not be read: {exc}")
            return
        dialog.gene_column.clear()
        dialog.gene_column.addItems([str(c) for c in dialog.table.columns])
        quant, why = suggest_quantification(dialog.table)
        dialog.quantification.setCurrentText(quant)
        dialog.note.setText(f"<b>{quant}</b> — {why}")

    def apply_import(self, dialog):
        """Run the import the dialog describes and put its columns into this session."""
        from .importer import merge_into, preprocess
        from .identity import GeneIndex
        try:
            resolve = GeneIndex.load().resolve
        except Exception:
            # No identity tables: the import still works on current accessions, and says it will not
            # reach the older ones rather than pretending it did.
            resolve = None
            self.status.showMessage("identity tables unavailable -- previous and strain accessions "
                                    "will not be resolved")
        try:
            imported, record = preprocess(
                dialog.table, gene_column=dialog.gene_column.currentText(),
                quantification=dialog.quantification.currentText(),
                scaling=dialog.scaling.currentText(),
                na_policy=dialog.na_policy.currentText(),
                duplicates=dialog.duplicates.currentText(),
                flip=dialog.flip.isChecked(),
                prefix=dialog.prefix.text().strip() or "imported_",
                resolve=resolve, log=print)
        except ValueError as exc:
            QtWidgets.QMessageBox.warning(self, "Import", str(exc))
            self.status.showMessage(f"import refused: {exc}")
            return None
        self.nodes = merge_into(self.nodes, imported)
        self.imports.append(record)
        self.numerics = numeric_columns(self.nodes)
        self.categories = category_columns(self.nodes)
        self._refresh_sources()
        if getattr(self, "panel", None) is not None:
            self.panel.nodes = self.nodes
            self.panel.add_imported(list(imported.columns))
        self.status.showMessage(
            f"imported {len(imported.columns)} column(s) for {record['genes']:,} genes "
            f"({record['quantification']}, {record['scaling']}) -- tick 'imported' on the Data tab "
            f"to build a map from them")
        return record

    # ------------------------------------------------------------------ logging
    def settings(self):
        """Where preferences persist. One place, so a new setting cannot invent its own file."""
        return QtCore.QSettings("starplast", "starplast")

    # ------------------------------------------------------------------ window size
    def work_area(self) -> QtCore.QRect:
        """The space a window may occupy on the monitor it is on.

        `availableGeometry`, not `geometry`: the difference is the taskbar, and it is the difference
        between a window that fits and one whose title bar is off the top of the screen. Falls back to
        the primary screen because during `__init__` the window has not been shown and has no screen
        of its own yet.
        """
        screen = self.screen() or QtGui.QGuiApplication.primaryScreen()
        return screen.availableGeometry() if screen else QtCore.QRect(0, 0, 1280, 720)

    def window_settings(self) -> dict:
        """The stored window preferences, with the defaults that apply on a fresh install."""
        s = self.settings()
        size = str(s.value("window/size", MATCH_SCREEN, type=str) or MATCH_SCREEN)
        return {"size": size if size in WINDOW_SIZES else MATCH_SCREEN,
                "fullscreen": bool(s.value("window/fullscreen", False, type=bool))}

    def apply_window_size(self, **changes) -> str:
        """Resize to the chosen size, or go full screen, and say what happened.

        Every branch ends up clamped to `work_area`, including the fixed sizes, so no choice here can
        produce a window larger than the display. That is the whole point of the control: the previous
        behaviour was one hardcoded size for every monitor, and a size too large for a screen is not
        recoverable by dragging when the title bar is off the top.

        Returns the sentence shown in Preferences, because a setting whose effect was silently clamped
        needs to say so -- otherwise choosing 4K on a 1080p screen looks like the setting was ignored.
        """
        cfg = {**self.window_settings(), **changes}
        s = self.settings()
        s.setValue("window/size", cfg["size"])
        s.setValue("window/fullscreen", bool(cfg["fullscreen"]))
        area = self.work_area()
        if cfg["fullscreen"]:
            self.showFullScreen()
            return f"full screen on a {area.width()} × {area.height()} work area"
        if self.isFullScreen():
            self.showNormal()
        if cfg["size"] == MATCH_SCREEN:
            want = QtCore.QSize(area.width(), area.height())
        else:
            w, h = (int(p.strip()) for p in cfg["size"].replace("×", "x").split("x"))
            want = QtCore.QSize(w, h)
        fitted = want.boundedTo(QtCore.QSize(area.width(), area.height()))
        self.resize(fitted)
        # Centre it: a clamped window pinned at its old top-left can still sit half off the screen.
        self.move(area.x() + max(0, (area.width() - fitted.width()) // 2),
                  area.y() + max(0, (area.height() - fitted.height()) // 2))
        note = f"{fitted.width()} × {fitted.height()}"
        if fitted != want:
            note += f" (clamped from {want.width()} × {want.height()} to fit this monitor)"
        return note

    def log_settings(self) -> dict:
        """The stored logging preferences, with the defaults that apply on a fresh install.

        Off, and WARNING to the console. Off because writing files to someone's disk uninvited is
        how a tool loses trust; WARNING rather than silent because a warning nobody enabled a log to
        see is a silent failure with extra steps.
        """
        s = self.settings()
        return {"enabled": s.value("logging/enabled", False, type=bool),
                "file_level": str(s.value("logging/file_level", "DEBUG")),
                "console_level": str(s.value("logging/console_level", "WARNING")),
                "directory": str(s.value("logging/directory", logging_util.log_dir()))}

    def apply_log_settings(self, **changes) -> str:
        """Apply the stored settings, with any changes, and persist what was applied.

        Persisted only after `configure` has accepted it, so a directory that could not be written
        does not come back on the next launch as though it had worked.
        """
        cfg = {**self.log_settings(), **changes}
        path = logging_util.configure(**cfg)
        self.log = logging_util.get_logger(__name__)
        if changes:
            # Written only when something was actually chosen. Persisting on every startup would
            # mean the application rewrites the user's settings file merely for having been opened.
            s = self.settings()
            for k, v in cfg.items():
                s.setValue(f"logging/{k}", v)
            where = f"writing to {path}" if path else "not writing a file"
            self.statusBar().showMessage(
                f"logging: {where}; console shows {cfg['console_level']} and above")
            if hasattr(self, "log_path"):
                self.log_path.setText(path or "no file is being written")
        return path

    # ------------------------------------------------------------------ appearance
    def _cmap_for(self, values):
        """The color map for a column: the user's choice if it suits the data, else the right default.

        A diverging map on a strictly positive quantity invents a midpoint, and a sequential map on a
        residual hides its sign, so the column's kind has the final say over an inappropriate choice.
        """
        kind = TH.kind_for_column(values)
        if self.cmap_name and TH.CMAPS.get(self.cmap_name, (None,))[0] == kind:
            return TH.CMAPS[self.cmap_name][1]
        return TH.DEFAULT_CMAP[kind]

    def apply_theme(self, name: str):
        """Repaint everything from the palette -- widgets, GL background, and the compartment colors."""
        self.theme = name
        app = QtWidgets.QApplication.instance()
        if app is not None:
            app.setStyleSheet(TH.stylesheet(name, self._container_opacity, self._ui_scale))
        self.view.setBackgroundColor(pg.mkColor(TH.palette_for(name)["bg"]))
        # Recolor the classes for this ground, then restore the deliberate grey for "unknown".
        self.color_of = dict(zip(self.comps,
                                  TH.categorical_colors(len(self.comps), name, self.cmap_name)))
        self.color_of["unassigned"] = TH.unknown_color(name)[:3]
        if getattr(self, "diagram", None) is not None:
            # Recolored with the map, from the same dict, so a theme change moves both together.
            self._refresh_diagram()
        if hasattr(self, "gallery"):
            # Thumbnails already painted keep the ground they were painted on; the large view is
            # re-rendered on the spot so at least the one being looked at follows the theme.
            self.gallery.background = TH.rgbf(TH.palette_for(name)["bg"])[:3]
            self.gallery.show_index(self.gallery.slider.value())
        if hasattr(self, "theme_box") and self.theme_box.currentText() != name:
            self.theme_box.blockSignals(True)
            self.theme_box.setCurrentText(name)
            self.theme_box.blockSignals(False)
        self.redraw()

    def apply_point_style(self):
        """Re-apply the point style and redraw."""
        st = TH.POINT_STYLES[self.point_style]
        self.scatter.setGLOptions(TH.gl_options(self.point_mode))
        self.scatter.setData(pos=self.xyz, color=self.colors(self.visible_mask()),
                             size=st["size"])
        self.redraw()

    def _analysis(self):
        """The analysis dock: data selection, tuning, clustering, the battery and the search.

        Tabbed with the evidence panel rather than given its own area, so the map keeps the width. Import
        failures are caught: the panel needs scikit-learn and umap-learn, and the browser must still open
        without them.
        """
        try:
            from .analysis_panel import AnalysisPanel
            from .annotations import AnnotationStore
            from .tuning import EmbeddingStore
        except Exception as e:
            self.statusBar().showMessage(f"analysis panel unavailable: {e}")
            return
        d = QtWidgets.QDockWidget("analysis")
        d.setFeatures(QtWidgets.QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
        # Shares the window's runner, so an analysis appears in the Jobs panel and can be stopped
        # there like anything else. With its own private thread it was invisible and unstoppable,
        # and it also blocked every other tab for the several minutes a search takes.
        # Annotations live beside the cache, not inside it: the cache is built and can be rebuilt,
        # while these are the user's own proposals and must survive `build_graph`.
        self.annotations = AnnotationStore(
            os.path.join(paths.user_cache_dir(), "annotations.csv"))
        # Saved embeddings sit inside the package, beside the cache they were computed from -- unless
        # the state directory is overridden, which is how the suite keeps its own walks out of the
        # user's store. Read here rather than at import time so setting the variable takes effect.
        embeddings = os.path.join(os.environ.get(paths.ENV_STATE) or DATA, "embeddings")
        panel = AnalysisPanel(self.nodes, store=EmbeddingStore(embeddings),
                              runner=self.jobs, annotations=self.annotations)
        panel.status.connect(lambda m: self.statusBar().showMessage(m))
        panel.embedding_ready.connect(self.use_embedding)
        panel.clusters_ready.connect(self.use_clusters)
        panel.annotations_changed.connect(self.refresh_annotations)
        d.setWidget(panel)
        self.addDockWidget(QtCore.Qt.DockWidgetArea.RightDockWidgetArea, d)
        self.tabifyDockWidget(self.right_dock, d)
        self.right_dock.raise_()
        self.analysis_dock = d
        self.panel = panel
        self._gallery()

    def _open_workflows(self, tab=0):
        """Open task-oriented exploration, prediction and measured-screen comparison."""
        from .workflows import WorkflowDialog
        if not hasattr(self, 'workflows_dialog'):
            self.workflows_dialog = WorkflowDialog(self.nodes, runner=self.jobs, parent=self)
            self.workflows_dialog.gene_selected.connect(self._workflow_gene)
        self.workflows_dialog.tabs.setCurrentIndex(tab)
        self.workflows_dialog.show()
        self.workflows_dialog.raise_()

    def _workflow_gene(self, gene):
        """Keep a guided-workflow selection aligned with the map and evidence panel."""
        self.search.setText(gene)
        self.do_search()

    def _gallery(self):
        """The walk gallery, along the bottom where a wall of thumbnails has room to be a wall.

        Fed from the analysis panel's per-configuration signal, so a thumbnail appears as each map is
        computed rather than when the walk ends. Clicking one shows it in the central view, which is
        why this is a dock beside the map rather than a window over it: the gallery picks, the map
        displays, and what is displayed is the application's own 3D view with everything it can do.
        """
        from .gallery import GalleryPanel
        bg = TH.rgbf(TH.palette_for(self.theme)["bg"])[:3]
        self.gallery = GalleryPanel(color_fn=self.colors_for_genes, background=bg)
        self.gallery.chosen.connect(self.show_walk_map)
        self.gallery_dock = QtWidgets.QDockWidget("gallery", self)
        self.gallery_dock.setWidget(self.gallery)
        self.addDockWidget(QtCore.Qt.DockWidgetArea.BottomDockWidgetArea, self.gallery_dock)
        self.gallery_dock.hide()
        self.panel.walk_started.connect(self._walk_started)
        self.panel.walk_step.connect(self.gallery.add)
        # A search's configurations go to the same wall, each colored by the clustering that was
        # scored. The two walks produce the same kind of thing -- a map with a number attached --
        # and looking at them in two different places would be an accident of implementation.
        self.panel.search_step.connect(self.gallery.add)

    def _walk_started(self):
        """Clear the gallery and show it, so the first thumbnail lands somewhere visible."""
        self.gallery.clear()
        self.gallery_dock.show()
        self.gallery_dock.raise_()

    def colors_for_genes(self, mask):
        """The current coloring, restricted to a subset of genes -- what the gallery paints with.

        Taken from `colors` rather than reimplemented so a thumbnail is colored by exactly what the
        map is colored by, including the rule that gray means unknown. The visibility filter is
        deliberately not applied: a thumbnail showing only the filtered classes would look like a
        different embedding rather than the same one seen through a filter.
        """
        c = self.colors(np.ones(self.n, bool))
        m = np.asarray(mask)
        return c[m.astype(bool)] if m.dtype == bool else c[m.astype(int)]

    def use_embedding(self, coords, rows):
        """Swap the displayed map for one the analysis panel just built.

        `rows` says which genes this embedding has coordinates for, and the rest are recorded as
        unplaced rather than drawn. A walk embeds a seeded subsample, so most of the proteome has no
        position in one of its maps -- and every unplaced gene used to be left at the origin, where
        6,000 of them formed a dense lump in the middle of the map that could be clicked, filtered
        and counted like real structure. That is absence rendered as a value, which is the one thing
        this application exists not to do.
        """
        rows = np.asarray(rows)
        placed = (rows.astype(bool) if rows.dtype == bool
                  else np.isin(np.arange(self.n), rows.astype(int)))
        coords = np.asarray(coords, dtype=np.float32)
        # Unplaced genes sit at the centre of the cloud rather than at the origin, so that framing,
        # the horizon grid and the depth cue are computed off the real extent. They are never drawn.
        full = np.repeat(coords.mean(0)[None, :], self.n, axis=0).astype(np.float32)
        full[placed] = coords
        self.xyz = full
        self.view.xyz = full
        self.placed = placed
        # Picking works off the same mask: an invisible point at the centre of the map is still the
        # nearest point to a click there, so hiding without this selects a gene that is not shown.
        self.view.pickable = placed
        if self.sel is not None and not placed[self.sel]:
            self.sel = None
        self.scatter.setData(pos=self.xyz)
        self.view.fit_view()
        self.redraw()
        n = int(placed.sum())
        self.statusBar().showMessage(
            f"showing a rebuilt map over {n:,} genes"
            + (f"; the other {self.n - n:,} are not in this embedding and are hidden"
               if n < self.n else ""))

    def show_walk_map(self, step):
        """Show one configuration from the gallery in the central view.

        The expanded map is not a picture of a map: it goes through the same path as "build this
        map", so genes are clickable, every color mode applies and edges draw as usual.
        """
        self.use_embedding(step.coords, step.genes)
        self.statusBar().showMessage(
            f"walk configuration {step.index} of {step.total}: {step.label}"
            f"  ·  {int(np.sum(step.genes)):,} genes of {self.n:,}; the rest are not in this "
            f"embedding")

    # ------------------------------------------------------------------ menus
    def _menus(self):
        """Every setting, in the menu bar. The panel keeps only what is used continuously."""
        mb = self.menuBar()

        # ---- File
        f = mb.addMenu("&File")
        a = f.addAction("Export image…")
        a.setShortcut("Ctrl+E")
        a.triggered.connect(self.export_image)
        f.addAction("Export visible genes (CSV)…").triggered.connect(self.export_visible)
        f.addAction("Export gated selection (CSV)…").triggered.connect(self.export_gated)
        a = f.addAction("Import data…")
        a.setShortcut("Ctrl+I")
        a.setToolTip("Read your own table -- CSV, TSV, Excel or parquet -- resolve its identifiers "
                     "through the identity layer, and choose how it is normalized. Every choice is "
                     "recorded with the imported columns, and nothing is written into the shipped "
                     "cache.")
        a.triggered.connect(self.import_data)
        f.addSeparator()
        # One window per species. Switching opens a new window rather than swapping the table under
        # this one: nearly every panel here was built from the node table -- the category list, the
        # colour map, the feature blocks, the saved runs -- and rebuilding all of that in place is a
        # much larger change than opening the arm someone asked for.
        species = f.addMenu("Species")
        species.setToolTipsVisible(True)
        self.species_group = QtGui.QActionGroup(self)
        for name in available_species():
            act = species.addAction(name)
            act.setCheckable(True)
            act.setChecked(name == self.species)
            act.setToolTip(
                "One table per species, never a union. Merging them would put two identifier spaces "
                "in one index and make a statement like 'cluster 5 is 71% IMC' a claim about a "
                "mixture of organisms. What crosses between them is a bridge slot, which can be "
                "inspected and disbelieved.")
            self.species_group.addAction(act)
            act.triggered.connect(lambda _c=False, n=name: self.open_species(n))
        f.addSeparator()
        a = f.addAction("Preferences…")
        a.setShortcut("Ctrl+,")
        a.setMenuRole(QtWidgets.QMenu.__mro__ and QtGui.QAction.MenuRole.PreferencesRole)
        a.triggered.connect(self.open_preferences)
        f.addSeparator()
        # Results are the expensive thing this program produces -- a search is minutes to hours --
        # and until now they lived only until the window closed.
        a = f.addAction("Save all analysis results…")
        a.setToolTip("Every tab's table in one file, loadable back into the tabs it came from.")
        a.triggered.connect(lambda: self.panel.save_all_results()
                            if hasattr(self, "panel") else None)
        a = f.addAction("Load analysis results…")
        a.setToolTip("Load a saved bundle back into every tab it names. A loaded row is as "
                     "clickable as a computed one: it carries the recipe that rebuilds its map.")
        a.triggered.connect(lambda: self.panel.load_all_results()
                            if hasattr(self, "panel") else None)
        f.addAction("Export relationships (CSV)…").triggered.connect(self.export_relationships)
        f.addAction("Export graph (GraphML)…").triggered.connect(self.export_graphml)
        f.addSeparator()
        q = f.addAction("Quit")
        q.setShortcut("Ctrl+Q")
        q.triggered.connect(self.close)

        # ---- View
        v = mb.addMenu("&View")
        col = v.addMenu("Color by")
        self.color_group = QtGui.QActionGroup(self)
        for name in COLOR_MODES:
            act = col.addAction(name)
            act.setCheckable(True)
            act.setChecked(name == self.color_mode)
            self.color_group.addAction(act)
            act.triggered.connect(lambda _c, n=name: self.set_color_mode(n))

        ps = v.addMenu("Point size")
        self.size_group = QtGui.QActionGroup(self)
        for label, val in POINT_SIZES:
            act = ps.addAction(label)
            act.setCheckable(True)
            act.setChecked(val == DEFAULT_POINT_SIZE)
            act.setData(val)
            self.size_group.addAction(act)
            act.triggered.connect(lambda _c, s=val: self.set_point_size(s))

        mouse = v.addMenu("Left mouse button")
        mouse.setToolTipsVisible(True)
        self.mode_group = QtGui.QActionGroup(self)
        for name in INTERACTION_MODES:
            act = mouse.addAction(name.capitalize())
            act.setCheckable(True)
            act.setChecked(name == self.view.mode)
            act.setToolTip({"navigate": "Drag to rotate the map; click a gene to select it.",
                            "select": "Drag to draw a gate and take the genes inside it. Clicking "
                                      "one gene at a time is not a way to select a cluster."}[name])
            self.mode_group.addAction(act)
            act.triggered.connect(lambda _c, n=name: self.set_interaction_mode(n))
        mouse.addSeparator()
        self.axis_group = QtGui.QActionGroup(self)
        for name in NAVIGATE_AXES:
            act = mouse.addAction(f"Rotate: {name}")
            act.setCheckable(True)
            act.setChecked(name == self.view.axis)
            act.setToolTip(
                "Free orbit never returns to the same view twice, which is exactly wrong for "
                "comparing two maps. Constrained to an axis, the view is one number you can come "
                "back to." if name == "free" else
                f"Rotate about {name} only, so the view is reproducible.")
            self.axis_group.addAction(act)
            act.triggered.connect(lambda _c, n=name: self.set_navigate_axis(n))
        mouse.addSeparator()
        self.gate_group = QtGui.QActionGroup(self)
        for name in GATE_SHAPES:
            act = mouse.addAction(f"Gate: {name}")
            act.setCheckable(True)
            act.setChecked(name == self.view.gate_shape)
            act.setToolTip(
                "A polygon in screen space, so it takes everything behind it too -- which is what "
                "'grab that visual cluster' means." if name.startswith("lasso") else
                "A ball in world space around the gene you press on, sized by the drag. Deep in "
                "the cloud a lasso also catches the far side: that looks like one structure and is "
                "two.")
            self.gate_group.addAction(act)
            act.triggered.connect(lambda _c, n=name: self.set_gate_shape(n))

        v.addSeparator()
        self.spin_act = v.addAction("Spin")
        self.spin_act.setCheckable(True)
        self.spin_act.setShortcut("Ctrl+R")
        self.spin_act.setToolTip(
            "Rotate the map continuously. Depth in a 3D scatter reads only when it moves: a still "
            "frame of 8,140 points is a flat disc however good the shading is, and parallax is the "
            "only cue that separates a near point from a far one.")
        self.spin_act.toggled.connect(self.toggle_spin)
        self.spin_busy_act = v.addAction("Spin while a job runs")
        self.spin_busy_act.setCheckable(True)
        self.spin_busy_act.setChecked(True)
        self.spin_busy_act.setToolTip(
            "Spin while something is running and stop where it started, so the motion says 'working' "
            "without also moving the camera you had set up.")
        v.addSeparator()
        v.addAction("Reset view / clear filters").triggered.connect(self.reset)

        # ---- Edges
        self.edge_menu = mb.addMenu("&Edges")
        self.edge_menu.setToolTipsVisible(True)
        self.edge_act = {}
        for k, label in EDGE_TYPES:
            act = self.edge_menu.addAction(label)
            act.setCheckable(True)
            act.setEnabled(k in self.edges)
            act.setChecked(self.edge_on.get(k, False) and k in self.edges)
            n = len(self.edges[k]["a"]) if k in self.edges else 0
            act.setToolTip(f"{n:,} edges" if n else "not present in this build")
            act.toggled.connect(lambda on, key=k: self.set_edge(key, on))
            self.edge_act[k] = act
        self.edge_menu.addSeparator()
        self.all_edges_act = self.edge_menu.addAction("Draw all active edges")
        self.all_edges_act.setCheckable(True)
        self.all_edges_act.setChecked(self.all_edges_on)
        self.all_edges_act.setToolTip(
            f"Off, edges are drawn only for the gene you have selected -- the only view you can "
            f"actually trace. On, every active edge is drawn, up to the strongest {EDGE_CAP:,} per "
            f"type, and alpha falls as the count rises so a dense layer reads as brightness rather "
            f"than as a solid sheet over the map.")
        self.all_edges_act.toggled.connect(lambda on: (setattr(self, "all_edges_on", on),
                                                       self.redraw()))
        self.attn_act = self.edge_menu.addAction("Attention-corrected co-mention")
        self.attn_act.setCheckable(True)
        self.attn_act.setChecked(self.attn_on)
        self.attn_act.setToolTip(
            "Raw co-mention counts track how often a gene is studied, not how related two genes are. "
            "Corrected shows log2 observed/expected given each gene's own publication count.")
        self.attn_act.toggled.connect(lambda on: (setattr(self, "attn_on", on), self.redraw()))
        self.edge_menu.addSeparator()
        self.edge_menu.addAction("Why are these never combined?").triggered.connect(self.explain_edges)

        # ---- Tools
        t = mb.addMenu("&Tools")
        t.setToolTipsVisible(True)
        for index, title in enumerate(('Explore a gene…', 'Predict a trait…', 'Compare a screen…')):
            action = t.addAction(title)
            action.setToolTip('Open a guided workflow with source evidence, held-out evaluation and portable exports.')
            action.triggered.connect(lambda _checked=False, tab=index: self._open_workflows(tab))
        t.addSeparator()
        # A toggle beside the panel toggles, not a "…" that only ever opens. It sits with the other
        # things you show and hide, which is where a reader looks for it.
        self.slot_tree_act = t.addAction("Slot tree")
        self.slot_tree_act.setCheckable(True)
        self.slot_tree_act.setShortcut("Ctrl+T")
        self.slot_tree_act.setToolTip(
            "Every slot, its address in each of the three hierarchies, and what fills it -- with the "
            "empty ones coloured, because they are the map of what has not been measured.\n\n"
            "It reads the catalog's own functions rather than a copy of them, so it cannot report a "
            "clean catalog while the build sees a broken one, and it shows BOTH arms whichever one "
            "this window is displaying.")
        self.slot_tree_act.toggled.connect(self.toggle_slot_tree)
        t.addSeparator()
        for dock in (self.console_dock, self.jobs_dock, self.chat_dock,
                     getattr(self, "analysis_dock", None), getattr(self, "gallery_dock", None),
                     self.right_dock):
            if dock is not None:
                t.addAction(dock.toggleViewAction())

        # The same action in both menus, deliberately. It was asked for in Tools and found in View, and
        # one QAction in two places keeps a single tick rather than two that can disagree.
        v.addSeparator()
        v.addAction(self.slot_tree_act)

        h = mb.addMenu("&Help")
        h.setToolTipsVisible(True)
        for label, page in (("User guide", "guide.html"), ("Python API", "API.html")):
            action = h.addAction(label)
            action.setToolTip("Open the Starplast documentation in your web browser. Requires an internet connection.")
            action.triggered.connect(lambda _checked=False, target=page:
                QtGui.QDesktopServices.openUrl(QtCore.QUrl("https://einarolafsson.github.io/starplast/" + target)))
        h.addSeparator()
        h.addAction("What this map does and does not show").triggered.connect(self.explain_map)
        h.addAction("Precision, recall, and how each can be gamed").triggered.connect(
            self.explain_scoring)
        h.addSeparator()
        a = h.addAction("About starplast")
        a.setMenuRole(QtGui.QAction.MenuRole.AboutRole)
        a.triggered.connect(self.about)

    def _context_menu(self, pos):
        """Right-click on the map. Shows the menu; `build_context_menu` makes it."""
        m = self.build_context_menu()
        m.exec(self.view.mapToGlobal(pos))
        return m

    def display_choices(self) -> list:
        """Every display setting that is a CHOICE, as (label, options, current, apply).

        One table, so the right-click menu and Preferences cannot drift apart -- they used to be able
        to: a control added to the dialog was not in the menu, and the menu is where a reader already
        is when they want to change how the map looks, because they are right-clicking the map.

        Numeric settings are deliberately absent. A slider is not a menu item, and text size, spin
        speed, panel opacity and the three blob numbers stay in Preferences where they can be dragged.
        """
        from . import lighting as L
        auto = "auto (match the data)"
        return [
            ("Theme", list(TH.THEMES), self.theme, self.apply_theme),
            ("Colour map", [auto] + list(TH.CMAPS), self.cmap_name or auto, self._on_cmap),
            ("Point style", list(TH.POINT_STYLES), self.point_style, self._on_point_style),
            ("Overlap", list(TH.POINT_MODES), self.point_mode, self._on_point_mode),
            ("Point render mode", list(L.POINT_MODES), self._lighting["point_mode"],
             lambda v: self.set_lighting_option("point_mode", v)),
            ("Light render mode", list(L.MODES), self._lighting["mode"], self.set_lighting),
            ("Light target", list(L.SOURCES), self._lighting["source"],
             lambda v: self.set_lighting_option("source", v)),
            ("Pointer beam", list(L.POINTER_MODES), self._lighting["pointer_mode"],
             lambda v: self.set_lighting_option("pointer_mode", v)),
            ("Pointer response", list(L.RESPONSES), self._lighting["response"],
             lambda v: self.set_lighting_option("response", v)),
            ("Target marker", list(L.TARGET_MARKERS), self._lighting["target_marker"],
             lambda v: self.set_lighting_option("target_marker", v)),
            ("Light mood", list(L.LIGHT_MOODS), self._lighting["mood"],
             lambda v: self.set_lighting_option("mood", v)),
            ("Background", ["none", "blobs"], self._ambient_mode, self.set_ambient),
            ("Window size", list(WINDOW_SIZES), self.window_settings()["size"],
             lambda v: self._on_window_change(size=v)),
        ]

    def display_toggles(self) -> list:
        """The display settings that are booleans, as (label, state, apply)."""
        return [
            ("Full screen", self.isFullScreen(),
             lambda on: self._on_window_change(fullscreen=bool(on))),
            ("Fade with distance", self.depth_cue, self.set_depth_cue),
            ("Reference grid", self.show_ground, self.set_show_ground),
        ]

    def set_depth_cue(self, on: bool) -> bool:
        """Whether distance fades and shrinks a point. Named so the checkbox in Preferences and the
        item in the right-click menu drive one function rather than two equivalent lambdas."""
        self.depth_cue = bool(on)
        self.redraw()
        return self.depth_cue

    def set_show_ground(self, on: bool) -> bool:
        """Whether the horizon grid is drawn."""
        self.show_ground = bool(on)
        self.redraw()
        return self.show_ground

    def _add_display_menu(self, parent) -> QtWidgets.QMenu:
        """Add the whole display section to a menu, built from the two tables above.

        One builder rather than the lists written out again: every option list here is the same object
        Preferences offers, so a finish retired or a mode added appears in both places without anyone
        having to remember there is a second place.
        """
        menu = parent.addMenu("Display")
        menu.setToolTipsVisible(True)
        self._display_groups = []
        for label, options, current, apply in self.display_choices():
            sub = menu.addMenu(label)
            sub.setToolTipsVisible(True)
            sub.menuAction().setToolTip(DISPLAY_HELP[label])
            group = QtGui.QActionGroup(self)
            group.setExclusive(True)
            for name in options:
                act = sub.addAction(str(name))
                act.setCheckable(True)
                act.setChecked(str(name) == str(current))
                act.setData(str(name))
                act.setToolTip(f"{name}: {DISPLAY_HELP[label]}")
                act.setStatusTip(DISPLAY_HELP[label])
                group.addAction(act)
                # `apply` and `name` bound as defaults. A closure over the loop variables would give
                # every action the LAST option in the list, which is the classic form of this bug.
                act.triggered.connect(lambda _checked=False, fn=apply, v=str(name): fn(v))
            self._display_groups.append((label, group))
        menu.addSeparator()
        for label, state, apply in self.display_toggles():
            act = menu.addAction(label)
            act.setToolTip(DISPLAY_HELP[label])
            act.setCheckable(True)
            act.setChecked(bool(state))
            act.toggled.connect(lambda on, fn=apply: fn(on))
        menu.addSeparator()
        menu.addAction("All display settings…").triggered.connect(self.open_preferences)
        return menu

    def build_context_menu(self):
        """Construct the menu without showing it.

        Split from `_context_menu` for the same reason `build_preferences` is split from
        `open_preferences`: exec() enters a modal loop and does not return until a human closes the
        menu, so a test that called it hung forever instead of failing.
        """
        m = QtWidgets.QMenu(self)
        m.addAction(self.spin_act)
        m.addSeparator()
        mouse = m.addMenu("Left mouse button")
        for act in self.mode_group.actions():
            mouse.addAction(act)
        mouse.addSeparator()
        for act in self.gate_group.actions():
            mouse.addAction(act)
        if self.gated is not None and len(self.gated):
            m.addAction(f"Export the {len(self.gated):,} gated genes (CSV)…").triggered.connect(
                self.export_gated)
            m.addAction("Clear the gate").triggered.connect(self.clear_gate)
        m.addSeparator()
        ps = m.addMenu("Point size")
        for act in self.size_group.actions():
            ps.addAction(act)
        self._add_display_menu(m)
        m.addSeparator()
        m.addAction("Export image…").triggered.connect(self.export_image)
        m.addAction("Export visible genes (CSV)…").triggered.connect(self.export_visible)
        m.addAction("Export relationships (CSV)…").triggered.connect(self.export_relationships)
        m.addAction("Export graph (GraphML)…").triggered.connect(self.export_graphml)
        m.addSeparator()
        m.addAction("Reset view / clear filters").triggered.connect(self.reset)
        return m

    # ------------------------------------------------------------------ menu state
    def set_color_mode(self, name: str):
        """Change what color encodes. Each mode is a different claim about the data."""
        self.color_mode = name
        self.redraw()

    def set_point_size(self, size):
        """Set an absolute point size, or None to follow the point style."""
        self.point_size = size
        self.redraw()

    def set_edge(self, key: str, on: bool):
        """Turn one edge layer on or off. Layers are never merged."""
        self.edge_on[key] = bool(on)
        self.redraw()

    def explain_edges(self):
        """Explain, in words, why relationship types are kept separate."""
        QtWidgets.QMessageBox.information(self, "Why edge types are kept separate", EDGE_EXPLANATION)

    def explain_map(self):
        """Explain what the map is and what held-out testing says it does not support."""
        QtWidgets.QMessageBox.information(self, "What this map shows", MAP_EXPLANATION)

    def about(self):
        """What this is, what it is standing on, and who drew the cell."""
        QtWidgets.QMessageBox.about(self, "About starplast", self.about_text())
        return self.about_text()

    def about_text(self) -> str:
        """The About text, built rather than written down, so it cannot go stale.

        Everything in it is read at the moment it is asked for: the version, what is doing the
        computing, how big the shipped cache is and how many genes it holds. A dialog that claims a
        version the program is not running is worse than no dialog.
        """
        from . import __version__, gpu, paths
        b = gpu.backend()
        credit = getattr(self.diagram, "credit", {}) if self.diagram is not None else {}
        return (
            f"<h3>starplast {__version__}</h3>"
            f"<p>A 3D browser for the <i>Toxoplasma gondii</i> knowledge map: {len(self.nodes):,} "
            f"genes, positioned by the saved balanced feature recipe. "
            f"The layout is exploratory; guided prediction uses separate held-out evaluation.</p>"
            f"<p><b>Computing:</b> UMAP {b['umap']}, clustering {b['cluster']}.<br>"
            f"<b>Cache:</b> {paths.data_dir()}</p>"
            f"<p>Proximity here is a hypothesis to check, never evidence on its own. Grey means "
            f"unknown: never zero, never a category.</p>"
            + (f"<p>Cell drawing by {credit.get('creator')} (SwissBioPics, SIB), "
               f"<a href='{credit.get('license')}'>CC BY 4.0</a>.</p>" if credit.get("creator")
               else "")
            + "<p><a href='https://github.com/EinarOlafsson/starplast'>"
              "github.com/EinarOlafsson/starplast</a></p>")

    def explain_scoring(self):
        """The precision/recall explainer, read from objectives.py so it cannot drift from the code.

        Shown in a monospace, selectable box: it contains a table, and a table reflowed into a
        proportional font is unreadable.
        """
        from .objectives import EXPLANATION
        d = QtWidgets.QDialog(self)
        d.setWindowTitle("Precision, recall, and how each can be gamed")
        lay = QtWidgets.QVBoxLayout(d)
        view = QtWidgets.QPlainTextEdit(EXPLANATION)
        view.setReadOnly(True)
        view.setLineWrapMode(QtWidgets.QPlainTextEdit.LineWrapMode.NoWrap)
        view.setFont(QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.SystemFont.FixedFont))
        view.setMinimumSize(760, 520)
        lay.addWidget(view)
        b = QtWidgets.QDialogButtonBox(QtWidgets.QDialogButtonBox.StandardButton.Close)
        b.rejected.connect(d.reject)
        b.accepted.connect(d.accept)
        lay.addWidget(b)
        return d

    # ------------------------------------------------------------------ interaction modes
    def set_interaction_mode(self, name: str):
        """Switch what the left mouse button does, and say so -- a silent mode change is a trap."""
        self.view.mode = name if name in INTERACTION_MODES else "navigate"
        self.status.showMessage(
            "select: drag to gate a set of genes; the gate's composition appears on the right"
            if self.view.mode == "select" else
            "navigate: drag to rotate, click a gene to select it")

    def set_navigate_axis(self, name: str):
        """Constrain the orbit to one axis, or free it."""
        self.view.axis = name if name in NAVIGATE_AXES else "free"
        self.status.showMessage(
            "free orbit -- it will not return to this view exactly" if self.view.axis == "free"
            else f"rotating about {self.view.axis} only, so this view can be returned to")

    def set_gate_shape(self, name: str):
        """Choose the 2D lasso or the 3D brush."""
        self.view.gate_shape = name if name in GATE_SHAPES else GATE_SHAPES[0]

    def on_gated(self, idx):
        """Take a gated set of genes: mark them, report what is in it, and offer it for export."""
        idx = np.asarray(idx, dtype=int)
        self.gated = idx
        self.redraw()
        if not len(idx):
            self.status.showMessage("the gate caught no genes -- it is empty, not broken")
            return
        self.detail.setHtml(self.describe_gate(idx))
        self.status.showMessage(
            f"{len(idx):,} genes gated  ·  File ▸ Export gated selection, or right-click the map")

    def clear_gate(self):
        """Drop the gated set and go back to showing the whole map."""
        self.gated = None
        self.detail.setHtml("<p style='color:#888'>Click a gene.</p>")
        self.redraw()
        self.status.showMessage("gate cleared")

    def describe_gate(self, idx) -> str:
        """What is in a gated set: its composition by the current category, and what carries none.

        Composition rather than a list, because the question a gate is drawn to answer is "what is
        this clump". The unlabeled count leads, because those genes are the candidates a gate
        exists to produce and most of this proteome is among them.
        """
        idx = np.asarray(idx, dtype=int)
        vals = self.category_values().to_numpy()[idx]
        absent = {str(x).lower() for x in ABSENCE}
        counts = pd.Series(vals).value_counts()
        named = [(v, n) for v, n in counts.items() if str(v).lower() not in absent]
        unlabelled = int(sum(n for v, n in counts.items() if str(v).lower() in absent))
        rows = "".join(
            f"<tr><td>{v}</td><td align='right'>{n}</td>"
            f"<td align='right' style='color:#888'>{n / len(idx):.0%}</td></tr>"
            for v, n in named[:15])
        return (f"<h3>{len(idx):,} genes gated</h3>"
                f"<p style='color:#888'>Composition by <b>{self.category}</b>. "
                f"<b>{unlabelled:,}</b> carry no value for it — those are the candidates a gate is "
                f"drawn to find, and the map is a hypothesis about them, not evidence.</p>"
                f"<table width='100%'>{rows}</table>")

    def export_gated(self, path: str = "", columns=None):
        """Write the gated genes to CSV, with their coordinates.

        A separate action from "export visible" rather than a mode of it: the two answer different
        questions -- everything passing the filter, versus the set someone drew a gate around -- and
        one export that silently meant whichever had happened last would be worse than either.
        """
        if self.gated is None or not len(self.gated):
            self.status.showMessage("no gated selection -- set the left button to Select and drag "
                                    "a gate over some genes")
            return None
        path = path or self._ask_path("Export gated genes", "CSV (*.csv)", "starplast_gated.csv")
        if not path:
            return None
        if columns is None:
            d = self.choose_export_columns()
            if not d.exec():
                return None
            columns = self._ticked(d)
        cols = ["gene_id"] + [c for c in columns if c != "gene_id" and c in self.nodes.columns]
        out = self.nodes.iloc[self.gated][cols].copy()
        out["x"], out["y"], out["z"] = (self.xyz[self.gated, 0], self.xyz[self.gated, 1],
                                        self.xyz[self.gated, 2])
        out.to_csv(path, index=False)
        self.status.showMessage(f"wrote {len(out):,} gated genes and {len(cols)} columns to {path}")
        return path

    def refresh_annotations(self):
        """Re-read the annotations file and redraw, switching to the color that shows them.

        Re-read rather than tracked: the file is meant to be shared and hand-edited, and a window
        holding its own idea of what is in it would disagree with the file the moment anyone did.
        """
        try:
            self.annotated = self.annotations.mask(self.nodes.gene_id)
        except Exception as exc:                       # a bad row must not take the window down
            print(f"starplast: annotations unavailable ({type(exc).__name__}: {exc})")
            return
        n = int(self.annotated.sum())
        if n:
            self.set_color_mode("annotations")
        else:
            self.redraw()
        self.status.showMessage(
            f"{n:,} genes carry an annotation -- drawn in their own color, which nothing else uses; "
            f"they are proposals, not measurements")

    def use_clusters(self, labels):
        """Take a clustering from the analysis panel, KEEP it, and color the map by it.

        Kept rather than only drawn: the second run used to replace the first with no way back, so
        the one comparison this application is for -- does this structure survive different settings
        -- could not be made by looking. It arrives in the color-by panel under a timestamp name.

        The coloring switches automatically, because a user who has just pressed "cluster this map"
        wants to see the clusters and leaving it on compartment made the button look inert.
        """
        self.cluster_labels = np.asarray(labels)
        recipe = {}
        panel = getattr(self, "panel", None)
        if panel is not None:
            try:
                recipe = {"spec": panel.spec().to_dict(), "algorithm": panel.algo.currentText(),
                          "min_cluster_size": panel.mcs.value(), "eps": panel.eps.value()}
            except Exception:      # the panel is optional; a run without its recipe still beats none
                recipe = {}
        self.keep_run(self.cluster_labels, recipe=recipe)
        self.set_color_mode("clusters")
        n = len(set(self.cluster_labels[self.cluster_labels >= 0]))
        self.status.showMessage(f"coloring by {n} clusters, kept as a run you can come back to; "
                                f"gray is unclustered, which is a real answer and not a missing one")

    def _window_tab(self):
        """How large the window opens, and whether it opens full screen.

        Its own tab rather than a row in Display, because these two are the only settings that take
        effect before anything is drawn, and they are the ones to reach for when the window is the
        wrong size for the monitor -- which has to be findable without scrolling a tab that is mostly
        lighting controls.
        """
        cfg = self.window_settings()
        w = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(w)

        self.window_size_box = QtWidgets.QComboBox()
        self.window_size_box.addItems(WINDOW_SIZES)
        self.window_size_box.setCurrentText(cfg["size"])
        self.window_size_box.setToolTip(
            "The size the window opens at. 'match screen' is the default and means the work area of "
            "whichever monitor it opens on -- the screen minus its taskbar, which is the largest size "
            "that actually fits.\n\n"
            "The fixed sizes are for reproducing a figure at the same size on a different monitor. "
            "Every one of them is clamped to the monitor before use, so picking a size larger than "
            "the screen gives a window that fits rather than one whose title bar is off the top.")
        self.window_size_box.currentTextChanged.connect(
            lambda text: self._on_window_change(size=text))

        # A switch rather than a checkbox, to match the GPU control and spacr.
        self.fullscreen_switch = TH.Switch("", checked=cfg["fullscreen"])
        self.fullscreen_switch.setToolTip(
            "Fill the whole screen with no title bar. Off by default: the map is usually read beside "
            "something else, and a window with no title bar cannot be moved.\n\n"
            "While this is on the size above has no effect; turning it off restores that size.")
        self.fullscreen_switch.toggled.connect(
            lambda on: self._on_window_change(fullscreen=bool(on)))

        self.window_note = QtWidgets.QLabel("")
        self.window_note.setWordWrap(True)
        self._refresh_window_note()

        form.addRow("window size", self.window_size_box)
        form.addRow("full screen", self.fullscreen_switch)
        form.addRow("", self.window_note)
        return w

    def _on_window_change(self, **changes) -> str:
        """Apply a window setting and report the size actually used."""
        note = self.apply_window_size(**changes)
        self._refresh_window_note(note)
        self.status.showMessage(f"window: {note}")
        return note

    def _refresh_window_note(self, note: str = "") -> str:
        """Say what the window is and what the monitor allows, since the two can differ."""
        area = self.work_area()
        text = note or f"{self.width()} × {self.height()}"
        self.window_note.setText(f"now {text} — this monitor allows up to "
                                 f"{area.width()} × {area.height()}")
        return text

    def _display_tab(self):
        """Scenery and lighting: the two settings that change how the map LOOKS rather than what it
        says. Kept apart from Appearance because everything here costs frames, and a control that
        costs frames should be somewhere a person can find it again to turn it off."""
        from . import ambient, lighting
        w = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(w)

        self.ambient_box = QtWidgets.QComboBox()
        self.ambient_box.addItems(["none", "blobs"])
        self.ambient_box.setCurrentText(self._ambient_mode)
        self.ambient_box.setToolTip(
            "Soft colour blobs drifting behind the panels -- spaCR's own background, with the same "
            "controls. Scenery: it never covers a control and never takes a click.")
        self.ambient_box.currentTextChanged.connect(self.set_ambient)

        def slider(bounds, value, step=0.05):
            """Create a numeric preference control with explicit bounds, step, and initial value."""
            s = QtWidgets.QDoubleSpinBox()
            s.setRange(*bounds)
            s.setSingleStep(step)
            s.setValue(value)
            return s

        self.ambient_speed = slider(ambient.SPEED_RANGE, self._ambient["speed"])
        self.ambient_size = slider(ambient.SIZE_RANGE, self._ambient["size"])
        self.ambient_density = slider(ambient.DENSITY_RANGE, self._ambient["density"])
        for box, key in ((self.ambient_speed, "speed"), (self.ambient_size, "size"),
                         (self.ambient_density, "density")):
            box.setToolTip({"speed": "Speed of the background animation. Higher values move the blobs faster.",
                            "size": "Size of background blobs. This changes the background only.",
                            "density": "Density of background blobs. Higher values add more visual detail."}[key])
            box.valueChanged.connect(lambda v, k=key: self.set_ambient_option(k, v))

        self.light_box = QtWidgets.QComboBox()
        self.light_box.addItems(list(lighting.MODES))
        self.light_box.setCurrentText(self._lighting["mode"])
        self.light_box.setToolTip(
            "How light reaches the points.\n\n"
            "soft uses broad fill. ray traced lets dense clusters shadow genes behind them. In "
            "glossy and metallic modes, an OpenGL vertex shader marches GPU-accelerated rays "
            "through a 3D density texture; flat mode and old drivers use the CPU fallback. These "
            "are volumetric density rays—not Vulkan/path tracing or ray-core triangle tracing.")
        self.light_box.currentTextChanged.connect(self.set_lighting)

        self.light_source = QtWidgets.QComboBox()
        self.light_source.addItems(list(lighting.SOURCES))
        self.light_source.setCurrentText(self._lighting["source"])
        self.light_source.setToolTip(
            "What interaction positions the light.\n\n"
            "mouse flashlight casts continuously from the camera through the cursor—it never snaps "
            "to a gene or changes depth when overlapping points cross. selected gene emits from the "
            "clicked point. selected gene and its edges adds emitters at visible-edge neighbors.")
        self.light_source.currentTextChanged.connect(lambda v: self.set_lighting_option("source", v))

        self.pointer_mode = QtWidgets.QComboBox()
        self.pointer_mode.addItems(list(lighting.POINTER_MODES))
        self.pointer_mode.setCurrentText(self._lighting["pointer_mode"])
        self.pointer_mode.setToolTip(
            "The mouse beam shape. Broad, focused and soft flashlights differ in cone width and "
            "edge softness; parallel wash behaves like a distant studio lamp. This affects mouse "
            "targeting only—clicked-gene lights always emanate from the clicked coordinates.")
        self.pointer_mode.currentTextChanged.connect(
            lambda v: self.set_lighting_option("pointer_mode", v))

        self.light_response = QtWidgets.QComboBox()
        self.light_response.addItems(list(lighting.RESPONSES))
        self.light_response.setCurrentText(self._lighting["response"])
        self.light_response.setToolTip(
            "How quickly the flashlight follows the cursor. Direct is immediate; smooth removes "
            "small hand motion; cinematic trails slowly. All are continuous screen rays and none "
            "selects a data point.")
        self.light_response.currentTextChanged.connect(
            lambda v: self.set_lighting_option("response", v))

        self.target_marker = QtWidgets.QComboBox()
        self.target_marker.addItems(list(lighting.TARGET_MARKERS))
        self.target_marker.setCurrentText(self._lighting["target_marker"])
        self.target_marker.setToolTip(
            "A purely visual marker showing where the light is aimed. None draws nothing; halo, "
            "beacon and pulse are comparison candidates and never encode genes or enter analysis.")
        self.target_marker.currentTextChanged.connect(
            lambda v: self.set_lighting_option("target_marker", v))

        self.light_mood = QtWidgets.QComboBox()
        self.light_mood.addItems(list(lighting.LIGHT_MOODS))
        self.light_mood.setCurrentText(self._lighting["mood"])
        self.light_mood.setToolTip(
            "The soft light's color temperature. Neutral preserves category colors most directly; "
            "cool blue and warm are visibly colored fills whose peak intensity is capped so they "
            "do not turn dense regions white.")
        self.light_mood.currentTextChanged.connect(lambda v: self.set_lighting_option("mood", v))

        self.point_render = QtWidgets.QComboBox()
        self.point_render.addItems(list(lighting.POINT_MODES))
        self.point_render.setCurrentText(self._lighting["point_mode"])
        self.point_render.setToolTip(
            "How each gene is drawn. flat is a simple data-color disc with gentle diffuse light "
            "but no spherical highlight. glossy 3D reconstructs a sphere normal for every fragment, "
            "uses a GGX highlight, and writes the curved surface depth. metallic 3D reflects a soft "
            "studio environment with its data color, a darker body, and a strong Fresnel rim. These "
            "are OpenGL material shaders, not renamed textures; an old driver falls back safely.")
        self.point_render.currentTextChanged.connect(
            lambda v: self.set_lighting_option("point_mode", v))

        self.container_opacity = QtWidgets.QDoubleSpinBox()
        self.container_opacity.setRange(0.35, 1.0)
        self.container_opacity.setSingleStep(0.05)
        self.container_opacity.setValue(self._container_opacity)
        self.container_opacity.setToolTip(
            "How solid the panels are over the drifting background.\n\n"
            "Only the containers take it. A translucent FIELD would put moving colour behind text "
            "somebody is reading, and the point of a background is that it stays behind things. It "
            "stops at 0.35 for the same reason: below that the compartment list is legible only "
            "while the blobs happen to be elsewhere.")
        self.container_opacity.valueChanged.connect(self.set_container_opacity)


        # A slider, because this is a thing people drag until it looks right rather than a number
        # they know in advance -- and the label beside it says what the number is, since a slider
        # with no read-out cannot be set back to where it was.
        zoom_row = QtWidgets.QWidget()
        zoom_lay = QtWidgets.QHBoxLayout(zoom_row)
        zoom_lay.setContentsMargins(0, 0, 0, 0)
        self.zoom_box = QtWidgets.QSlider(QtCore.Qt.Orientation.Horizontal)
        self.zoom_box.setRange(int(UI_SCALE_RANGE[0] * 100), int(UI_SCALE_RANGE[1] * 100))
        self.zoom_box.setSingleStep(5)
        self.zoom_box.setPageStep(10)
        self.zoom_box.setValue(int(round(self._ui_scale * 100)))
        self.zoom_label = QtWidgets.QLabel(f"{self._ui_scale:.2f}×")
        self.zoom_box.setToolTip(
            "Text size for the whole interface, as a multiple of the size this desktop asked for.\n\n"
            "It moves the application's font rather than a stylesheet, so every layout re-measures "
            "its own contents: a longer label makes a wider row instead of being cut off at the "
            "old width. Nothing here is ever clipped -- if a panel cannot fit its text it grows, "
            "and if the window cannot fit the panel the panel scrolls.")
        self.zoom_box.valueChanged.connect(lambda v: self.set_ui_scale(v / 100.0))
        zoom_lay.addWidget(self.zoom_box, 1)
        zoom_lay.addWidget(self.zoom_label)

        form.addRow("text size", zoom_row)
        form.addRow("background", self.ambient_box)
        form.addRow("blob speed", self.ambient_speed)
        form.addRow("blob size", self.ambient_size)
        form.addRow("blob density", self.ambient_density)
        form.addRow("panel opacity", self.container_opacity)
        form.addRow(QtWidgets.QLabel(""))
        form.addRow("light render mode", self.light_box)
        form.addRow("light target", self.light_source)
        form.addRow("pointer beam", self.pointer_mode)
        form.addRow("pointer response", self.light_response)
        form.addRow("target marker", self.target_marker)
        form.addRow("light mood", self.light_mood)
        form.addRow("point render mode", self.point_render)
        return w

    def set_ui_scale(self, scale: float) -> float:
        """Scale every piece of text in the program, and let the layouts follow.

        On the QApplication's font rather than on a stylesheet: a stylesheet font-size does not
        change what a widget reports as its size hint, so the text grew and the boxes did not, and
        labels were cut off at the old width. Changing the application font invalidates every
        layout, which is the whole point -- a wider label makes a wider row.
        """
        scale = float(min(max(scale, UI_SCALE_RANGE[0]), UI_SCALE_RANGE[1]))
        self._ui_scale = scale
        QtCore.QSettings("starplast", "starplast").setValue("display/ui_scale", scale)
        if getattr(self, "zoom_label", None) is not None:
            self.zoom_label.setText(f"{scale:.2f}×")
        app = QtWidgets.QApplication.instance()
        if app is not None:
            font = QtGui.QFont(self._base_font)
            size = self._base_font.pointSizeF()
            if size > 0:
                font.setPointSizeF(size * scale)
            else:
                font.setPixelSize(max(int(self._base_font.pixelSize() * scale), 1))
            # `setFont` on the application and nothing else. Qt propagates it to every widget that
            # has not been given a font of its own and re-lays them out; walking allWidgets() to set
            # it by hand touched widgets that were mid-deletion and segfaulted the interpreter.
            app.setFont(font)
            # And the stylesheet, which carries font sizes of its own -- a stylesheet font-size
            # BEATS the application font, so without this the setting changed the font and the
            # sheet immediately overrode it on every widget. That is why it appeared to do nothing.
            app.setStyleSheet(TH.stylesheet(self.theme, self._container_opacity, scale))
            self.updateGeometry()
        return scale

    def set_ambient(self, mode: str) -> str:
        """Turn the drifting background on or off, and remember the choice."""
        self._ambient_mode = "blobs" if str(mode) == "blobs" else "none"
        QtCore.QSettings("starplast", "starplast").setValue("display/ambient", self._ambient_mode)
        self._apply_ambient()
        return self._ambient_mode

    def set_ambient_option(self, key: str, value) -> None:
        """One of speed, size or density."""
        self._ambient[key] = float(value)
        QtCore.QSettings("starplast", "starplast").setValue(f"display/ambient_{key}", float(value))
        if self._ambient_widget is not None:
            self._ambient_widget.configure(**{key: float(value)})

    def _apply_ambient(self) -> None:
        """Build the background widget on first use, and show or hide it."""
        from .ambient import AmbientWidget
        if self._ambient_mode == "none":
            if self._ambient_widget is not None:
                self._ambient_widget.hide()
            return
        if self._ambient_widget is None:
            pal = TH.palette_for(self.theme)
            # Parented to the WINDOW, not to one panel: it has to be behind every container -- the
            # find panel, the compartment list, the cell, the analysis tabs -- which is what the
            # panels' own opacity then lets through. Behind a single panel it was scenery for one
            # corner of the screen.
            self._ambient_widget = AmbientWidget(
                colors=[pal["accent"], pal["accent_lo"], pal.get("info", pal["accent_hi"])],
                background=pal["bg"], parent=self)
            self._ambient_widget.lower()
            self.installEventFilter(self)
        self._ambient_widget.configure(**self._ambient)
        self._ambient_widget.setGeometry(self.rect())
        self._ambient_widget.show()
        self._ambient_widget.lower()

    def eventFilter(self, obj, ev):
        """Keep the background the size of the panel it sits behind."""
        if (obj is self and self._ambient_widget is not None
                and ev.type() == QtCore.QEvent.Type.Resize):
            self._ambient_widget.setGeometry(self.rect())
        return super().eventFilter(obj, ev)

    def set_lighting(self, mode: str) -> str:
        """Choose off, soft illumination, or one of the density-ray shadow depths."""
        from .lighting import MODES
        mode = "soft" if mode == "lit" else mode       # migrate old recipes/settings
        self._lighting["mode"] = mode if mode in MODES else "off"
        QtCore.QSettings("starplast", "starplast").setValue("display/lighting",
                                                            self._lighting["mode"])
        if self._lighting["mode"] == "off":
            self._light_timer.stop()
            self._refresh_sprite()
            self.redraw()
        else:
            self._light_timer.start(60)
        return self._lighting["mode"]

    def set_lighting_option(self, key: str, value) -> str:
        """Set the interaction target, light mood, or point render mode and apply it now."""
        from .lighting import (LIGHT_MOODS, POINTER_MODES, POINT_MODES, RESPONSES, SOURCES,
                               TARGET_MARKERS, normalize_point_mode, normalize_source)
        key = "point_mode" if key == "finish" else key   # public compatibility, not a UI control
        if key == "source":
            value = (normalize_source(value) if str(value) in (*SOURCES, "mouse")
                     else self._lighting["source"])
        elif key == "point_mode":
            known = str(value) in (*POINT_MODES, "2D", "matt", "satin", "glossy", "metallic")
            value = normalize_point_mode(value) if known else self._lighting["point_mode"]
        elif key == "mood":
            value = value if value in LIGHT_MOODS else self._lighting["mood"]
        elif key == "pointer_mode":
            value = value if value in POINTER_MODES else self._lighting["pointer_mode"]
        elif key == "response":
            value = value if value in RESPONSES else self._lighting["response"]
        elif key == "target_marker":
            value = value if value in TARGET_MARKERS else self._lighting["target_marker"]
        else:
            return ""
        self._lighting[key] = value
        QtCore.QSettings("starplast", "starplast").setValue(f"display/light_{key}",
                                                            self._lighting[key])
        if key == "point_mode":
            self._refresh_sprite()
        if self._lighting["mode"] != "off":
            self._light_tick()
        return str(self._lighting[key])

    def set_container_opacity(self, value: float) -> float:
        """How much of the background shows through the panels."""
        self._container_opacity = float(min(max(float(value), 0.35), 1.0))
        QtCore.QSettings("starplast", "starplast").setValue("display/container_opacity",
                                                            self._container_opacity)
        self.apply_theme(self.theme)
        return self._container_opacity

    def edge_neighbours(self, index: int) -> list:
        """Genes joined to this one by an edge type that is currently DRAWN.

        Currently drawn, not merely present: the question the light is answering is "what is this
        gene connected to in the picture in front of me", and a light on a neighbour whose edge type
        is switched off would point at a relationship the reader cannot see.
        """
        out = set()
        for k, _ in EDGE_TYPES:
            if not self.edge_on.get(k) or k not in self.edges:
                continue
            e = self.edges[k]
            a, b = np.asarray(e["a"]), np.asarray(e["b"])
            out.update(b[a == index].tolist())
            out.update(a[b == index].tolist())
        out.discard(index)
        return sorted(out)

    def frame_lights(self) -> list:
        """The lights for this frame, from whichever source is chosen."""
        from .lighting import RESPONSES, light_at, mood_color
        try:
            basis = self.view.camera_basis()
        except Exception:                       # no GL context yet -- offscreen, or mid-startup
            basis = None
        mood = mood_color(self._lighting["mood"])
        pointer = getattr(self.view, "pointer", None)
        if self._lighting["source"] == "mouse flashlight":
            target = np.asarray(pointer or (0.0, 0.0, 1.0), dtype=float)
            response = float(RESPONSES[self._lighting["response"]])
            if self._smoothed_pointer is None:
                self._smoothed_pointer = target
            else:
                self._smoothed_pointer += response * (target - self._smoothed_pointer)
            pointer = tuple(self._smoothed_pointer)
        else:
            self._smoothed_pointer = None
        lit = light_at(self.xyz, self._lighting["source"], self._light_t,
                        1, 0.0, self._light_radius(),
                        pointer=pointer, basis=basis, selected=self.sel,
                        neighbours=(self.edge_neighbours(self.sel)
                                    if self.sel is not None
                                    and self._lighting["source"].endswith("edges") else None),
                        color=mood, pointer_mode=self._lighting["pointer_mode"],
                        fov=float(self.view.opts.get("fov", 60.0)),
                        aspect=max(float(self.view.width()), 1.0) /
                               max(float(self.view.height()), 1.0))

        ray_mode = "ray traced" in self._lighting["mode"]
        grid = self.occupancy() if ray_mode else None
        absorb = 0.86 if self._lighting["mode"] == "deep ray traced" else 0.44
        gpu_rays = self.scatter.set_scene(self._lighting["point_mode"], lit, mood, grid, absorb)
        return lit if gpu_rays else self.cast_rays(lit)

    def _refresh_sprite(self, lit=None) -> bool:
        """Rebuild the ball each gene is drawn as, if the finish or the light has moved.

        The sprite is where a gene stops being a disc of colour and becomes a lit sphere, so the
        highlight on it has to come from the same direction as the light on the cloud -- otherwise
        every ball is lit from the top left while the map is lit from the right, and the eye reads
        the two as different scenes. Rebuilt only when the direction has actually moved, because
        this runs on the light timer: a 64x64 array is cheap and an upload every frame for a light
        that has turned by a thousandth of a degree is still waste.
        """
        from . import sprite as SP
        finish = self._lighting["point_mode"]
        where = SP.DEFAULT_LIGHT
        try:
            basis = self.view.camera_basis()
        except Exception:
            basis = None
        lit = self.frame_lights() if lit is None and self._lighting["mode"] != "off" else lit
        if lit is None:
            from .lighting import mood_color
            self.scatter.set_scene(finish, [], mood_color(self._lighting["mood"]), None)
        # The GPU shader receives world-space lights every frame; rebuilding and uploading a CPU
        # texture as they move is both wasted work and a source of one-frame highlight jumps. Only
        # the old-context fallback bakes the light direction into a texture.
        gpu_material = self.scatter.gpu_material_enabled(finish)
        if not gpu_material and basis is not None and lit:
            centre = self.xyz.mean(axis=0) if len(self.xyz) else np.zeros(3)
            where = SP.to_screen(np.asarray(lit[0]["pos"], dtype=float) - centre, basis)
        was = self._sprite_state
        moved = was is None or was[0] != finish or float(np.dot(
            np.asarray(where) / max(float(np.linalg.norm(where)), 1e-9),
            np.asarray(was[1]) / max(float(np.linalg.norm(was[1])), 1e-9))) < 0.995
        if not moved:
            return False
        self._sprite_state = (finish, where)
        self.scatter.set_sprite(SP.texture(finish, where))
        return True

    def _eye(self):
        """Where the camera is, in world coordinates, or None before there is one."""
        try:
            return self.view.camera_basis()[0]
        except Exception:
            return None

    def occupancy(self):
        """The map binned into a coarse box, for asking what is in the way of what.

        Cached against the coordinates themselves: switching embedding or filtering genes changes
        the cloud, and everything else -- orbiting, relighting, sixty frames a second -- does not.
        """
        from . import rays
        if self._grid is None or self._grid_for is not self.xyz:
            self._grid = rays.build(self.xyz)
            self._grid_for = self.xyz
        return self._grid

    def cast_rays(self, lit: list) -> list:
        """Trace volumetric shadow rays when the explicit ray-traced mode is selected.

        Each output light carries one transmittance per gene, computed by marching the segment from
        that gene to the light through the cached density volume. Soft mode deliberately skips this
        cost and cannot cast shadows.
        """
        if "ray traced" not in self._lighting["mode"] or not lit or not len(self.xyz):
            return lit
        grid = self.occupancy()
        strength = 0.86 if self._lighting["mode"] == "deep ray traced" else 0.44
        return [dict(light, shadow=np.power(grid.transmittance(self.xyz, light["pos"]),
                                            strength / 0.55)) for light in lit]

    def _draw_emitter(self, lit: list) -> None:
        """Draw the selected experimental target marker without encoding data."""
        marker = self._lighting["target_marker"]
        if marker == "none" or not lit:
            if self.emitter_item is not None:
                self.emitter_item.setVisible(False)
            return
        # One target per light, and `lit` is already known non-empty above, so this list is too --
        # no guard beneath it, for the same reason the GSE129869 loader has none: a check that
        # cannot fire tells a reader the case is handled when nothing handles it.
        pos = np.asarray([np.asarray(light.get("target", light["pos"]), dtype=np.float32)
                          for light in lit], dtype=np.float32)
        phase = 0.5 + 0.5 * np.sin(self._light_t * 2.4)
        size = {"halo": 15.0, "beacon": 10.0, "pulse": 12.0 + 12.0 * phase}[marker]
        alpha = {"halo": 0.34, "beacon": 0.72, "pulse": 0.28 + 0.38 * phase}[marker]
        if marker == "beacon":
            lifted = pos + np.array([0.0, 0.0, self._light_radius() * 0.08], np.float32)
            pos = np.concatenate((pos, lifted), axis=0)
        color = np.tile(np.array([1.0, 0.86, 0.42, alpha], np.float32), (len(pos), 1))
        if self.emitter_item is None:
            self.emitter_item = gl.GLScatterPlotItem(pos=pos, color=color, size=size, pxMode=True)
            self.emitter_item.setGLOptions("translucent")
            self.view.addItem(self.emitter_item)
        else:
            self.emitter_item.setData(pos=pos, color=color, size=size, pxMode=True)
            self.emitter_item.setVisible(True)

    def _light_radius(self) -> float:
        """How far out the lights orbit: outside the cloud, so they light it rather than sit in it."""
        return float(np.abs(self.xyz).max() or 1.0) * 1.6

    def _light_ground(self, lit) -> None:
        """Light the grid with the same lights as the points.

        A lit cloud over an unlit grid reads as two pictures: the horizon is what the eye uses to
        judge where the light is coming from, and leaving it flat throws that away.

        Computed here rather than through `shade`, which infers a point's normal from the direction
        out of the cloud's centre. That is right for a point cloud and wrong for a floor: the probe
        sits BELOW the centre, so the inferred normal pointed down and every light above the map
        gave the grid exactly zero. A plane has a real normal, and it is up.
        """
        if self.grid_item is None or not len(self.xyz):
            return
        p = TH.palette_for(self.theme)
        base = np.array(TH.rgbf(p["border"])[:3], dtype=float)
        centre = self.xyz.mean(0).astype(float)
        r = self.view.data_radius() or 1.0
        probe = np.array([centre[0], centre[1], float(self.xyz[:, 2].min()) - r * 0.08])
        up = np.array([0.0, 0.0, 1.0])

        # Accumulate RGB rather than one brightness scalar. Otherwise a warm/cool light changes
        # the genes but leaves the horizon neutral, which makes the two look like separate scenes.
        total = np.full(3, 0.35, dtype=float)
        for light in lit:
            to_light = np.asarray(light["pos"], dtype=float) - probe
            dist = float(np.linalg.norm(to_light)) or 1.0
            strength = max(float(np.dot(up, to_light / dist)), 0.0) / (1.0 + dist / (2.0 * r))
            total += strength * np.asarray(light.get("color", np.ones(3)), dtype=float)
        rgb = np.clip(base * np.minimum(total, 2.0) * 1.5, 0.0, 1.0)
        # The grid is drawn faint, so colour alone is not visible: the light moves its ALPHA too,
        # between the theme's own value and about three times it. That is what makes a lit horizon
        # read as lit rather than as a slightly different grey.
        flat = 70 if TH.is_light(p) else 40
        brightness = float(np.mean(total))
        self.grid_item.setColor((int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255),
                                 int(np.clip(flat * min(brightness, 2.2), 15, 210))))

    def _light_tick(self) -> None:
        """One frame of moving light: reshade the points from the flat colours redraw computed.

        From `_base_colors` rather than from whatever is on screen, because shading an
        already-shaded array darkens it a little more every frame until the map goes black.
        """
        from .lighting import POINT_MODES, shade
        if self._lighting["mode"] == "off" or self._base_colors is None:
            return
        self._light_t += 0.06
        lit = self.frame_lights()
        if self.scatter.gpu_material_enabled(self._lighting["point_mode"]):
            # The fragment shader needs the unlit category color as its albedo. Feeding it the old
            # CPU highlight made every sphere reflect an already-white point and washed the map out.
            colors = np.array(self._base_colors, copy=True)
        else:
            colors = shade(self.xyz, self._base_colors, lit,
                           point_mode=self._lighting["point_mode"], eye=self._eye())
        size_scale = POINT_MODES[self._lighting["point_mode"]]["size"]
        self.scatter.setData(pos=self.xyz, color=colors, size=self._base_sizes * size_scale)
        self._refresh_sprite(lit)
        self._draw_emitter(lit)
        self._light_ground(lit)

    def open_preferences(self):
        """Appearance settings, in a window of their own.

        Shown rather than exec'd, and kept on `self`: a modal dialog freezes the map behind it, so
        every setting had to be judged from memory of what the map looked like a moment ago. This
        one stays open beside the window, is moved independently of it, and changes take effect
        under it while it sits there.
        """
        if getattr(self, "_prefs", None) is None:
            self._prefs = self.build_preferences()
            self._prefs.setModal(False)
            self._prefs.setWindowFlag(QtCore.Qt.WindowType.Window, True)
        self._prefs.show()
        self._prefs.raise_()
        self._prefs.activateWindow()
        return self._prefs

    def toggle_slot_tree(self, on: bool):
        """Show or hide the slot tree. The menu item is a checkbox, so it has to be able to close it."""
        window = self.open_slot_tree() if on else getattr(self, "_slot_tree", None)
        if window is not None:
            window.setVisible(bool(on))
        return window

    def open_slot_tree(self):
        """The slot tree, in a window of its own. Shown rather than exec'd, and kept on `self`.

        A modal here would mean checking the catalog against the map from memory, which is the exact
        comparison it is for.
        """
        from .slot_tree import SlotTreeWindow
        if getattr(self, "_slot_tree", None) is None:
            # Opens on the arm this window is showing. Defaulting to whichever organism code sorts
            # first meant the tree opened on Plasmodium while the map showed Toxoplasma, and a filter
            # typed against what was on screen returned nothing.
            self._slot_tree = SlotTreeWindow(
                organism=SPECIES[self.species]["code"], parent=self)
        self._slot_tree.show()
        self._slot_tree.raise_()
        return self._slot_tree

    def open_species(self, name: str):
        """Open another species' map in a window of its own, and remember the choice.

        Returns the new window. The old one is closed only once the new one exists: building it reads
        a cache and lays out several thousand genes, and closing first would leave nothing on screen
        if that failed.
        """
        if name == self.species or name not in SPECIES:
            return self
        # Written HERE, not in the constructor. Persisting on every construction meant that merely
        # opening a window changed what opens next time -- and in the suite, one test opening the
        # Plasmodium arm silently moved eighteen later tests onto a table whose genes they do not
        # contain. Remembering is a consequence of CHOOSING, which is this method.
        self.settings().setValue("data/species", name)
        other = Window(species=name)
        other.show()
        # Kept on the new window so Python does not collect the object that owns the running event
        # filters and timers the moment this method returns.
        other._opened_from = self
        self.close()
        return other

    def _gpu_wanted(self) -> bool:
        """Whether the GPU switch is on, from settings so it survives a restart.

        Where nobody has touched the switch, the answer is whether a backend is THERE -- the same
        default `gpu.enabled()` applies. It used to hardcode False here, so on a machine with CUDA
        installed the program ran its searches on the GPU while Preferences showed the switch off:
        the setting disagreed with the behaviour, and the switch was the one lying.
        """
        from . import gpu
        s = QtCore.QSettings("starplast", "starplast")
        if s.contains("compute/gpu"):
            return bool(s.value("compute/gpu", type=bool))
        return any(gpu.available()[k] for k in ("cuml", "cupy", "torch"))

    def _on_gpu(self, on: bool) -> str:
        """Remember the choice and say what it will actually do.

        Reported rather than assumed: a switch that silently does nothing because no backend is
        installed is worse than no switch, and this is the sentence that tells the difference.
        """
        from . import gpu
        QtCore.QSettings("starplast", "starplast").setValue("compute/gpu", bool(on))
        note = gpu.describe()
        if hasattr(self, "gpu_note"):
            self.gpu_note.setText(note)
        self.status.showMessage(note)
        return note

    def compare_backends(self, sample: int = 2000):
        """Run the same embedding both ways and show the two maps beside each other.

        In a job, because it is two UMAPs and the window should stay usable; the dialog is built by
        `build_comparison` so a test can inspect it without entering a modal loop.
        """
        from .benchmark import compare
        self.status.showMessage("building the same map on the CPU and on the GPU…")
        spec = self.panel.spec() if getattr(self, "panel", None) is not None else None
        job = self.run_job(lambda: compare(self.nodes, spec, sample=sample,
                                           log=lambda m: self.status.showMessage(str(m))),
                           "CPU vs GPU")

        def show(jid: int, ok: bool):
            # One-shot, and only for THIS job: the runner's signal carries every job's completion,
            # and a comparison dialog opening because some unrelated search finished would be a
            # window appearing for no reason the user can connect to anything they did.
            """Display the backend comparison when its own job completes successfully."""
            if jid != job.id:
                return
            self.jobs.finished.disconnect(show)
            if ok and isinstance(job.result, dict):
                self.build_comparison(job.result).exec()
            else:
                self.status.showMessage(f"comparison failed: {job.error or 'no result'}")

        self.jobs.finished.connect(show)
        return job

    def build_comparison(self, result: dict):
        """The two maps, the two clocks, and how much they agree -- built, not shown."""
        from .gallery import thumbnail
        d = QtWidgets.QDialog(self)
        d.setWindowTitle("CPU vs GPU")
        lay = QtWidgets.QVBoxLayout(d)
        maps = QtWidgets.QHBoxLayout()
        for key, title in (("cpu", "CPU"), ("gpu", "GPU")):
            run = result.get(key)
            box = QtWidgets.QVBoxLayout()
            label = QtWidgets.QLabel()
            if run is not None:
                img = thumbnail(run["coords"], size=260, background=TH.rgbf(
                    TH.palette_for(self.theme)["bg"])[:3])
                label.setPixmap(QtGui.QPixmap.fromImage(img))
                text = f"<b>{title}</b><br>{run['backend']}<br>{run['seconds']:.1f} s"
            else:
                label.setText("not available")
                text = f"<b>{title}</b><br>—"
            box.addWidget(label)
            box.addWidget(QtWidgets.QLabel(text))
            maps.addLayout(box)
        lay.addLayout(maps)
        note = QtWidgets.QLabel(result.get("note", ""))
        note.setWordWrap(True)
        lay.addWidget(note)
        close = QtWidgets.QPushButton("close")
        close.clicked.connect(d.accept)
        lay.addWidget(close)
        return d

    def build_preferences(self):
        """Construct the dialog without showing it.

        Split from `open_preferences` so the controls can be built and inspected without entering a
        modal event loop: exec() blocks until a human closes the window, so a test that called it hung
        forever rather than failing.
        """
        d = QtWidgets.QDialog(self)
        d.setWindowTitle("Preferences")
        outer = QtWidgets.QVBoxLayout(d)
        self.pref_tabs = QtWidgets.QTabWidget()
        outer.addWidget(self.pref_tabs)
        appearance = QtWidgets.QWidget()
        form = QtWidgets.QFormLayout(appearance)
        self.pref_tabs.addTab(appearance, "Appearance")

        self.theme_box = QtWidgets.QComboBox()
        self.theme_box.addItems(TH.THEMES)
        self.theme_box.setCurrentText(self.theme)
        self.theme_box.setToolTip("dark and light are the general pair; slate is a low-contrast dark "
                                  "for long sessions, paper a high-contrast light for figures.")
        self.theme_box.currentTextChanged.connect(self.apply_theme)

        self.cmap_box = QtWidgets.QComboBox()
        self.cmap_box.addItem("auto (match the data)")
        self.cmap_box.addItems(list(TH.CMAPS))
        if self.cmap_name:
            self.cmap_box.setCurrentText(self.cmap_name)
        self.cmap_box.setToolTip(
            "Sequential for ordered quantities, diverging only where values straddle a midpoint, "
            "categorical for classes. A map of the wrong kind is ignored in favor of the right "
            "default, because a diverging ramp on a positive quantity invents a midpoint.")
        self.cmap_box.currentTextChanged.connect(self._on_cmap)

        self.point_box = QtWidgets.QComboBox()
        self.point_box.addItems(list(TH.POINT_STYLES))
        self.point_box.setCurrentText(self.point_style)
        self.point_box.setToolTip("pinpoint suits 8,000 genes at once; halo shows density; large is "
                                  "for a filtered subset.")
        self.point_box.currentTextChanged.connect(self._on_point_style)

        self.mode_box = QtWidgets.QComboBox()
        self.mode_box.addItems(TH.POINT_MODES)
        self.mode_box.setCurrentText(self.point_mode)
        self.mode_box.setToolTip(
            "occlude: nearer points hide farther ones — the correct default.\n"
            "additive: overlaps sum, which reads as density but saturates dense regions to white "
            "and destroys the color encoding.")
        self.mode_box.currentTextChanged.connect(self._on_point_mode)

        # GPU acceleration, as a switch rather than a checkbox: the same control this user has in
        # spacr, so a setting looks like a setting in both programs.
        from . import gpu
        self.gpu_switch = TH.Switch("", checked=self._gpu_wanted())
        self.gpu_note = QtWidgets.QLabel(gpu.describe())
        self.gpu_note.setWordWrap(True)
        self.gpu_switch.setToolTip(
            "cuml does UMAP and HDBSCAN themselves; cupy or torch do the array work -- scaling, "
            "ranking and the distance matrix the walk recomputes for every configuration. Nothing "
            "here is a dependency: with no backend installed the switch has nothing to turn on -- "
            'install them with `pip install starplast[gpu]`, or `pip install -e ".[gpu]"` in a '
            "checkout.\n\n"
            "A map built by cuml's UMAP is NOT the map umap-learn builds -- it is a different map "
            "of the same data -- so a walk whose rows came from both would compare the libraries "
            "rather than the settings. The arithmetic paths are checked against the CPU to 1e-5.")
        self.gpu_switch.toggled.connect(self._on_gpu)
        self.gpu_test = QtWidgets.QPushButton("compare CPU and GPU…")
        self.gpu_test.setToolTip(
            "Builds the same map twice, once each way, and shows both. The clock is the smaller "
            "half of the answer: cuml's UMAP is a different implementation, so turning this on does "
            "not speed a map up -- it produces a DIFFERENT map of the same data. The two are shown "
            "side by side with the share of each gene's nearest neighbours they agree on, which is "
            "the comparison that survives rotation, reflection and scale.")
        self.gpu_test.clicked.connect(self.compare_backends)

        self.spin_speed = QtWidgets.QDoubleSpinBox()
        self.spin_speed.setRange(0.05, 3.0)
        self.spin_speed.setSingleStep(0.05)
        self.spin_speed.setValue(self._spin_speed)
        self.spin_speed.setSuffix("  °/frame")
        self.spin_speed.setToolTip("Degrees of rotation per animation frame. Higher values rotate faster; this does not change the embedding or analysis.")
        self.spin_speed.valueChanged.connect(lambda v: setattr(self, "_spin_speed", v))

        self.depth_box = QtWidgets.QCheckBox("fade and shrink with distance")
        self.depth_box.setChecked(self.depth_cue)
        self.depth_box.setToolTip(
            "Without it, near and far points are equally bright and the map reads as a flat disc "
            "however far it is rotated. Turn it off to compare two points' colors exactly, since the "
            "fade changes apparent color with position.")
        self.depth_box.toggled.connect(self.set_depth_cue)

        self.ground_box = QtWidgets.QCheckBox("horizon grid")
        self.ground_box.setChecked(self.show_ground)
        self.ground_box.setToolTip(
            "A reference plane under the cloud. Spinning a bare point cloud, the eye cannot separate "
            "rotation from the points rearranging themselves.")
        self.ground_box.toggled.connect(self.set_show_ground)

        cfg = self.log_settings()
        self.log_box = QtWidgets.QCheckBox("keep a log file")
        self.log_box.setChecked(cfg["enabled"])
        self.log_box.setToolTip(
            "Write what happens to a rotating file under the cache directory: every job with how "
            "long it took, every fetch with its URL and outcome, every embedding with its recipe. "
            "Off by default because a tool that writes to your disk without being asked is one "
            "people stop trusting -- and on, it is what turns 'the search finished' into a record "
            "of which of its 288 configurations actually ran.")
        self.log_box.toggled.connect(lambda on: self.apply_log_settings(enabled=on))

        self.log_file_level = QtWidgets.QComboBox()
        self.log_file_level.addItems(logging_util.LEVELS)
        self.log_file_level.setCurrentText(cfg["file_level"])
        self.log_file_level.setToolTip(
            "How much detail the FILE keeps. DEBUG is the useful setting here: the file is read "
            "after something went wrong, and the line that explains it is usually the one nobody "
            "would have chosen to keep.")
        self.log_file_level.currentTextChanged.connect(
            lambda v: self.apply_log_settings(file_level=v))

        self.log_console_level = QtWidgets.QComboBox()
        self.log_console_level.addItems(logging_util.LEVELS)
        self.log_console_level.setCurrentText(cfg["console_level"])
        self.log_console_level.setToolTip(
            "How much reaches the console pane. Set independently of the file, because DEBUG is "
            "exactly what you want kept during a half-hour walk and exactly what you do not want "
            "scrolling past while you watch it.")
        self.log_console_level.currentTextChanged.connect(
            lambda v: self.apply_log_settings(console_level=v))

        self.log_path = QtWidgets.QLabel(logging_util.log_file() or "no file is being written")
        self.log_path.setToolTip("Where the log is being written. Rotating, so a long session "
                                 "cannot fill a disk.")
        self.log_path.setWordWrap(True)

        form.addRow("theme", self.theme_box)
        form.addRow("color map", self.cmap_box)
        form.addRow("points", self.point_box)
        form.addRow("rendering", self.mode_box)
        form.addRow("GPU acceleration", self.gpu_switch)
        form.addRow("", self.gpu_note)
        form.addRow("", self.gpu_test)
        form.addRow("spin speed", self.spin_speed)
        form.addRow("depth", self.depth_box)
        form.addRow("reference", self.ground_box)
        form.addRow("logging", self.log_box)
        form.addRow("keep at level", self.log_file_level)
        form.addRow("show at level", self.log_console_level)
        form.addRow("log file", self.log_path)
        self.pref_tabs.addTab(self._display_tab(), "Display")
        self.pref_tabs.addTab(self._window_tab(), "Window")
        close = QtWidgets.QPushButton("close")
        close.clicked.connect(d.accept)
        outer.addWidget(close)
        # The dialog is built after the window's own pass, so its tooltips are wrapped here.
        TH.wrap_tooltips(d)
        return d

    # ------------------------------------------------------------------ docks, jobs, progress
    def _tool_docks(self):
        """Console, jobs and assistant, tabbed together at the bottom.

        Hidden on startup. They are for when something is running or has gone wrong, and a browser
        that opens with three empty panels across the bottom has spent its screen space on nothing.
        """
        self.console = ConsolePanel()
        self.console.install()
        self.console_dock = QtWidgets.QDockWidget("console", self)
        self.console_dock.setWidget(self.console)

        jobs_w = QtWidgets.QWidget()
        jl = QtWidgets.QVBoxLayout(jobs_w)
        jl.setContentsMargins(6, 6, 6, 6)
        bar = QtWidgets.QHBoxLayout()
        self.stop_btn = QtWidgets.QPushButton("stop selected")
        self.stop_btn.setToolTip(
            "Ask the selected job to stop. Cooperative: it ends at the next point it reports "
            "progress, which for a hyperparameter walk is once per configuration -- seconds, not "
            "the rest of the run.")
        self.stop_btn.clicked.connect(self.stop_selected_job)
        self.stop_all_btn = QtWidgets.QPushButton("stop all")
        self.stop_all_btn.setToolTip("Ask every running job to stop.")
        self.stop_all_btn.clicked.connect(self.stop_all_jobs)
        self.free_btn = QtWidgets.QPushButton("free memory")
        self.free_btn.setToolTip(
            "Drop cached intermediates, collect garbage, and release the GPU cache if torch is "
            "loaded. Running jobs are untouched -- this reclaims what finished work left behind.")
        self.free_btn.clicked.connect(self.free_memory)
        for b in (self.stop_btn, self.stop_all_btn, self.free_btn):
            bar.addWidget(b)
        bar.addStretch(1)
        self.resources = QtWidgets.QLabel("")
        self.resources.setToolTip("Resident memory of this process, system memory, and GPU memory "
                                  "where a GPU is visible.")
        bar.addWidget(self.resources)
        jl.addLayout(bar)

        self.jobs_view = QtWidgets.QTreeWidget()
        self.jobs_view.setHeaderLabels(["job", "state", "detail"])
        self.jobs_view.setToolTip("Everything started this session, including what has finished. A "
                                  "failed job keeps its error here rather than vanishing with it.")
        self.jobs_view.setRootIsDecorated(False)
        self.jobs_view.setContextMenuPolicy(QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
        self.jobs_view.customContextMenuRequested.connect(self._job_menu)
        jl.addWidget(self.jobs_view, 1)

        self.jobs_dock = QtWidgets.QDockWidget("jobs", self)
        self.jobs_dock.setWidget(jobs_w)
        # A resource read costs a couple of file reads, so it is polled slowly and only while the
        # panel is visible -- a meter that samples every second is itself a background job.
        self._res_timer = QtCore.QTimer(self)
        self._res_timer.timeout.connect(self.refresh_resources)
        self._res_timer.start(3000)
        self.refresh_resources()

        self.chat = ChatPanel(context_provider=self.describe_state)
        self.chat_dock = QtWidgets.QDockWidget("assistant", self)
        self.chat_dock.setWidget(self.chat)

        area = QtCore.Qt.DockWidgetArea.BottomDockWidgetArea
        for d in (self.console_dock, self.jobs_dock, self.chat_dock):
            self.addDockWidget(area, d)
            d.hide()
        self.tabifyDockWidget(self.console_dock, self.jobs_dock)
        self.tabifyDockWidget(self.jobs_dock, self.chat_dock)

        self.jobs.started.connect(self._refresh_jobs)
        self.jobs.progress.connect(lambda *_: self._refresh_jobs())
        self.jobs.finished.connect(lambda *_: self._refresh_jobs())
        self.jobs.busy_changed.connect(self._on_busy)

    def _progress(self):
        """An indeterminate bar in the status bar, shown only while something is running."""
        self.progress_bar = QtWidgets.QProgressBar()
        self.progress_bar.setRange(0, 0)            # indeterminate until a job reports a percentage
        self.progress_bar.setMaximumWidth(180)
        self.progress_bar.setTextVisible(False)
        self.progress_bar.hide()
        self.status.addPermanentWidget(self.progress_bar)

    def _on_busy(self, busy: bool):
        """Show the bar, and optionally spin, while work is outstanding."""
        self.progress_bar.setVisible(busy)
        if self.spin_busy_act.isChecked():
            if busy and not self.spin_act.isChecked():
                self._spin_home = self._camera_state()
                self._spin_for_job = True
                self.spin_act.setChecked(True)
            elif not busy and getattr(self, "_spin_for_job", False):
                self._spin_for_job = False
                self.spin_act.setChecked(False)
                # Put the camera back where it was. A spinner that leaves the map at a random angle
                # has destroyed the view the user set up, which costs more than the reassurance.
                self._restore_camera(self._spin_home)
                self._spin_home = None

    def _camera_state(self):
        """Return camera azimuth, elevation, and distance for later restoration."""
        p = self.view.opts
        return (float(p.get("azimuth", 0.0)), float(p.get("elevation", 0.0)),
                float(p.get("distance", 40.0)))

    def _restore_camera(self, state):
        """Restore an (azimuth, elevation, distance) tuple; ignore an empty state."""
        if not state:
            return
        az, el, dist = state
        self.view.setCameraPosition(distance=dist, elevation=el, azimuth=az)

    def _refresh_jobs(self, *_):
        """Refresh the jobs list with status, errors, and stored traceback details."""
        self.jobs_view.clear()
        for j in sorted(self.jobs.jobs.values(), key=lambda x: -x.id):
            detail = j.error or j.note or ""
            it = QtWidgets.QTreeWidgetItem([j.name, j.state, detail])
            it.setData(0, QtCore.Qt.ItemDataRole.UserRole, j.id)
            if j.state == FAILED:
                it.setToolTip(2, j.traceback or j.error)
            self.jobs_view.addTopLevelItem(it)
        running = len(self.jobs.active())
        self.stop_btn.setEnabled(bool(running))
        self.stop_all_btn.setEnabled(bool(running))

    def selected_job(self):
        """The job selected in the Jobs panel, or None."""
        items = self.jobs_view.selectedItems()
        if not items:
            return None
        jid = items[0].data(0, QtCore.Qt.ItemDataRole.UserRole)
        return self.jobs.jobs.get(jid)

    def stop_selected_job(self):
        """Ask the selected job to stop, reporting why if it cannot."""
        job = self.selected_job()
        if job is None:
            self.status.showMessage("select a job in the list first")
            return
        if not job.active:
            self.status.showMessage(f"{job.name} has already finished ({job.state})")
            return
        job.cancel()
        self.status.showMessage(f"asked {job.name} to stop")
        self._refresh_jobs()

    def stop_all_jobs(self):
        """Ask every running job to stop."""
        n = len(self.jobs.active())
        if not n:
            self.status.showMessage("nothing is running")
            return
        self.jobs.cancel_all()
        self.status.showMessage(f"asked {n} job(s) to stop")
        self._refresh_jobs()

    def _job_menu(self, pos):
        """Right-click a job: stop it, or copy its error."""
        m = self.build_job_menu()
        m.exec(self.jobs_view.mapToGlobal(pos))
        return m

    def build_job_menu(self):
        """The Jobs right-click menu. Built rather than shown, because exec() blocks."""
        m = QtWidgets.QMenu(self)
        job = self.selected_job()
        stop = m.addAction("Stop this job")
        stop.setEnabled(bool(job and job.active))
        stop.triggered.connect(self.stop_selected_job)
        m.addAction("Stop all running jobs").triggered.connect(self.stop_all_jobs)
        m.addSeparator()
        err = m.addAction("Copy error / traceback")
        err.setEnabled(bool(job and job.traceback))
        err.triggered.connect(self.copy_job_error)
        return m

    def copy_job_error(self):
        """Copy a failed job's traceback to the clipboard."""
        job = self.selected_job()
        cb = QtWidgets.QApplication.clipboard()
        if job is not None and cb is not None:
            cb.setText(job.traceback or job.error or "")
            self.status.showMessage(f"copied the traceback from {job.name}")

    # ------------------------------------------------------------------ resources
    def resource_summary(self) -> str:
        """Process memory, system memory and GPU memory, as one line.

        Read from /proc and from torch where they exist, with no hard dependency on either: this is
        a status line, and a status line that raises on a machine without a GPU is worse than one
        that says nothing about GPUs.
        """
        parts = []
        try:
            with open("/proc/self/statm") as fh:
                rss_pages = int(fh.read().split()[1])
            parts.append(f"this process {rss_pages * os.sysconf('SC_PAGE_SIZE') / 2**30:.1f} GB")
        except Exception:
            pass
        try:
            info = {}
            with open("/proc/meminfo") as fh:
                for line in fh:
                    k, _, v = line.partition(":")
                    info[k] = float(v.strip().split()[0]) / 2**20      # kB -> GiB
            used = info["MemTotal"] - info["MemAvailable"]
            parts.append(f"system {used:.0f}/{info['MemTotal']:.0f} GB")
        except Exception:
            pass
        # Only if torch is ALREADY loaded, and never imported here. Importing it initialises CUDA,
        # and initialising CUDA inside a running OpenGL application segfaults -- which this did, on a
        # three-second timer, taking the next grabFramebuffer down with it. A status line must not be
        # able to crash the program it reports on.
        torch = sys.modules.get("torch")
        if torch is not None:
            try:
                if torch.cuda.is_initialized():
                    parts.append(f"GPU {torch.cuda.memory_reserved() / 2**30:.1f} GB reserved")
            except Exception:
                pass
        parts.append(f"{os.cpu_count()} CPUs")
        return "  ·  ".join(parts)

    def refresh_resources(self):
        """Update the memory and CPU line in the Jobs panel."""
        self.resources.setText(self.resource_summary())

    def free_memory(self) -> str:
        """Release what finished work left behind, and say what was actually released.

        Deliberately does NOT touch running jobs. "Clear RAM" that silently killed a search would be
        a data-loss button wearing a housekeeping label.
        """
        import gc
        freed = []
        n = gc.collect()
        freed.append(f"{n} unreachable objects")
        # Same rule: only if it is already loaded and CUDA already up. Freeing a cache that does not
        # exist is not worth initialising CUDA for.
        torch = sys.modules.get("torch")
        if torch is not None:
            try:
                if torch.cuda.is_initialized():
                    torch.cuda.empty_cache()
                    freed.append("GPU cache")
            except Exception:
                pass
        self.refresh_resources()
        msg = "freed: " + ", ".join(freed) + f"  ·  now {self.resource_summary()}"
        self.status.showMessage(msg)
        return msg

    def run_job(self, fn, name: str):
        """Submit background work. The one entry point, so everything slow is visible in one place."""
        job = self.jobs.submit(fn, name)
        self.status.showMessage(f"{name}…")
        return job

    def toggle_spin(self, on: bool):
        """Rotate the camera continuously. Depth in a 3D scatter only reads when it moves."""
        if not hasattr(self, "_spin_timer"):
            self._spin_timer = QtCore.QTimer(self)
            self._spin_timer.timeout.connect(self._spin_step)
        if on:
            self._spin_timer.start(33)          # ~30 fps; smooth without burning the GPU
        else:
            self._spin_timer.stop()

    def _spin_step(self):
        # orbit() moves the camera, not the data, so nothing is recomputed per frame.
        """Advance the camera orbit by the configured degrees per frame."""
        self.view.orbit(self._spin_speed, 0)

    def _on_cmap(self, name):
        """Apply a color map. Categorical colors are rebuilt, not just the continuous ramp.

        Only the ramp honoured this before, so choosing a map while coloring by compartment -- the
        default mode -- appeared to do nothing at all.
        """
        self.cmap_name = None if name.startswith("auto") else name
        kind = TH.CMAPS.get(self.cmap_name, (None,))[0] if self.cmap_name else None
        if kind in (None, "categorical"):
            self.color_of = dict(zip(self.comps, TH.categorical_colors(
                len(self.comps), self.theme, self.cmap_name if kind == "categorical" else None)))
            self.color_of["unassigned"] = TH.unknown_color(self.theme)[:3]
        self.redraw()

    def _on_point_style(self, name):
        """Apply the selected point size and opacity style, then redraw."""
        self.point_style = name
        self.apply_point_style()

    def _on_point_mode(self, name):
        """Apply the selected overlap mode without changing gene coordinates."""
        self.point_mode = name
        self.apply_point_style()

    # ------------------------------------------------------------------ panels
    def _left(self):
        """Search and the category filter. Everything else moved into the menus.

        This panel used to carry the level, the coloring, twelve edge checkboxes, spin, preferences
        and the compartment list all at once, which is most of the application's settings stacked in a
        column beside a 3D view. Settings that are chosen once belong in a menu; the two things used
        continuously -- finding a gene and narrowing the field -- stay on screen.
        """
        d = QtWidgets.QDockWidget("find")
        d.setFeatures(QtWidgets.QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
        w = QtWidgets.QWidget()
        # Kept, because the drifting background is parented to it and has to follow its size.
        self.left_panel = w
        L = QtWidgets.QVBoxLayout(w)
        L.setContentsMargins(8, 8, 8, 8)

        self.search = QtWidgets.QLineEdit(placeholderText="gene id or product…")
        self.search.setToolTip(
            "An exact accession first, then a partial one, then the product description. Most of this "
            "proteome has no symbol -- 6,500 genes are 'hypothetical protein' -- so product text is "
            "often the only handle you have on a gene.")
        self.search.returnPressed.connect(self.do_search)
        L.addWidget(self.search)

        L.addWidget(QtWidgets.QLabel("<b>color by</b>"))
        self.category_box = QtWidgets.QComboBox()
        self.category_box.setToolTip(
            "What color means right now, and what the list below filters and flies by. Three kinds "
            "of thing, because they answer the same question: any column with a manageable number "
            "of repeated values; any clustering you have kept, by name; and any quantity cut into "
            "bins, which is what makes a measurement comparable with a clustering. Localization is "
            "the worst-recovered property in this map, so it is a poor thing to be the only way in.")
        self.category_box.addItems(self.color_sources())
        if "compartment" in self.categories:
            self.category_box.setCurrentText("compartment")
        self.category_box.currentTextChanged.connect(self.on_category_changed)
        L.addWidget(self.category_box)

        self.bins_box = QtWidgets.QSpinBox()
        self.bins_box.setRange(2, 20)
        self.bins_box.setValue(self.bins)
        self.bins_box.setPrefix("bins: ")
        self.bins_box.setToolTip(
            "How many bins a quantity is cut into. Quantile bins, not equal-width: nearly every "
            "quantity in this table is heavy-tailed -- the fitness screens span 64x within "
            "themselves -- and equal-width bins put 95% of the genes in one color and call that a "
            "coloring. Only used when the source above is a binned quantity.")
        self.bins_box.valueChanged.connect(self.set_bins)
        self.bins_box.hide()
        L.addWidget(self.bins_box)

        self.comp_list = TH.CheckList()
        self.comp_list.setToolTip(
            "TICK to show only those classes; tick none to show everything. Drag across several rows "
            "and right-click to tick them all at once -- or click a single row's box. Space toggles "
            "whatever is selected.\n\n"
            "Ticks are what the map reads, and selection is only how you choose them, so an ordinary "
            "click can no longer destroy a set you built up over several ctrl-clicks.\n\n"
            "Double-click flies to a class's centroid. Counts are of genes with a value, so a class is "
            "at least this big and possibly bigger.")
        self.comp_list.checkedChanged.connect(self.redraw)
        self.comp_list.checkedChanged.connect(self._refresh_diagram)
        self.comp_list.itemDoubleClicked.connect(self.fly_to_compartment)
        L.addWidget(self.comp_list, 1)
        self._fill_category_list()

        # The cell, under the list, filled from the same palette the points are. Shown only for a
        # localization category: there is no sensible mapping from cell-cycle phase onto organelles.
        from .celldiagram import CellDiagram, available as diagram_available
        self.diagram = CellDiagram() if diagram_available() else None
        self.diagram_note = QtWidgets.QLabel("")
        self.diagram_note.setWordWrap(True)
        self.diagram_note.setStyleSheet("color: #888")
        # In a scroll area of FIXED height, because the note's own height changes with its text --
        # one line for a shared shape, three for the list of compartments the drawing has no
        # organelle for -- and a note that grows takes its height from the widget above it. The
        # parasite then jumps up and down as compartments are clicked, which reads as the drawing
        # being redrawn differently rather than as a caption reflowing. Scrolls when it overflows.
        self.diagram_note_area = QtWidgets.QScrollArea()
        self.diagram_note_area.setWidget(self.diagram_note)
        self.diagram_note_area.setWidgetResizable(True)
        self.diagram_note_area.setFixedHeight(DIAGRAM_NOTE_HEIGHT)
        self.diagram_note_area.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
        # Transparent, viewport included: a scroll area paints its own background, and a black card
        # under the caption covers whatever the theme is drawing behind the panel.
        self.diagram_note_area.setStyleSheet("background: transparent;")
        self.diagram_note_area.viewport().setAutoFillBackground(False)
        self.diagram_note_area.viewport().setStyleSheet("background: transparent;")
        self.diagram_note.setStyleSheet("color: #888; background: transparent;")
        self.diagram_note_area.setHorizontalScrollBarPolicy(
            QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        if self.diagram is not None:
            self.diagram.compartment_clicked.connect(self.select_compartment)
            L.addWidget(self.diagram, 1)
            L.addWidget(self.diagram_note_area)

        rename = QtWidgets.QHBoxLayout()
        self.run_name = QtWidgets.QLineEdit()
        self.run_name.setPlaceholderText("name this run…")
        self.run_name.setToolTip(
            "Rename the clustering selected above. Runs are named by the clock to the second so two "
            "a minute apart are distinguishable without anyone typing anything, but "
            "'hdbscan_60 on expression' is what you will look for a week later.")
        self.rename_btn = QtWidgets.QPushButton("save name")
        self.rename_btn.setToolTip(
            "Rename and re-save the run with its full recipe -- the embedding spec and the "
            "clustering parameters. Without the recipe a name is a label on nothing: the run cannot "
            "be rebuilt and two runs cannot be told apart except by their numbers.")
        self.rename_btn.clicked.connect(self._rename_current_run)
        rename.addWidget(self.run_name, 1)
        rename.addWidget(self.rename_btn)
        self.rename_row = QtWidgets.QWidget()
        self.rename_row.setLayout(rename)
        self.rename_row.hide()
        L.addWidget(self.rename_row)

        b = QtWidgets.QPushButton("reset view / clear filters")
        b.setToolTip("Clear the selection and every class filter, and frame the whole map again. "
                     "The way back when a filter has left you looking at forty genes and it is no "
                     "longer obvious which one.")
        b.clicked.connect(self.reset)
        L.addWidget(b)
        d.setWidget(w)
        w.setMinimumWidth(260)
        return d

    def _fill_category_list(self):
        """Rebuild the value list for the current source, ordered by size with absences last."""
        self.comp_list.blockSignals(True)
        self.comp_list.clear()
        vals = self.category_values()
        counts = vals.value_counts()
        absent = {str(x).lower() for x in ABSENCE}
        named = [v for v in counts.index if str(v).lower() not in absent]
        # Absence values are kept but sunk to the bottom: "unassigned" is usually the largest class in
        # this proteome, and sorting by size alone would put "we do not know" at the top of every list.
        tail = [v for v in counts.index if str(v).lower() in absent]
        for v in named + tail:
            it = self.comp_list.add(f"{v}  ({int(counts[v]):,})", v)
            col = self.color_of.get(v)
            if col is not None:
                it.setForeground(QtGui.QColor.fromRgbF(*col))
            if str(v).lower() in absent:
                it.setToolTip("Absence, not a class: these genes have no measurement, which is not "
                              "the same as measuring zero.")
        self.comp_list.blockSignals(False)

    def set_bins(self, n: int):
        """Change how finely a binned quantity is cut, and redraw if that is what is showing."""
        self.bins = int(n)
        if self.category.startswith(BIN_PREFIX):
            self.on_category_changed(self.category)

    def _rename_current_run(self):
        """Rename the run currently selected in the color-by list."""
        if not self.category.startswith(RUN_PREFIX):
            return
        self.rename_run(self.category[len(RUN_PREFIX):], self.run_name.text().strip())

    def on_category_changed(self, name: str):
        """Switch what color means, rebuild its palette and its list of values."""
        self.category = name
        # The controls that only apply to one kind of source appear only for that kind. A bins spin
        # box beside a compartment list is a control that does nothing, which teaches people that
        # controls here might do nothing.
        is_run = name.startswith(RUN_PREFIX)
        self.bins_box.setVisible(name.startswith(BIN_PREFIX))
        self.rename_row.setVisible(is_run)
        if is_run:
            self.run_name.setText(name[len(RUN_PREFIX):])
        # Coloring follows the panel: this IS the color choice, not a filter beside one.
        self.color_mode = "compartment" if name == "compartment" else self.color_mode
        vals = sorted(self.category_values().unique())
        absent = {str(x).lower() for x in ABSENCE}
        self.comps = ([v for v in vals if str(v).lower() not in absent]
                      + [v for v in vals if str(v).lower() in absent])
        self.color_of = dict(zip(self.comps,
                                  TH.categorical_colors(len(self.comps), self.theme, self.cmap_name)))
        for v in self.comps:
            if str(v).lower() in absent:
                self.color_of[v] = TH.unknown_color(self.theme)[:3]
        self._fill_category_list()
        self._refresh_diagram()
        self.redraw()

    def select_compartment(self, name: str):
        """Tick a compartment in the list, from the diagram. The other half of both directions.

        Ticks it rather than merely highlighting it, because the list now filters on ticks: clicking an
        organelle in the drawing has to do the same thing as ticking its row, or the two halves of "both
        directions" would no longer be the same operation.
        """
        for item in self.comp_list.items():
            if item.data(QtCore.Qt.ItemDataRole.UserRole) == name:
                self.comp_list.set_checked([name])
                self.comp_list.clearSelection()
                item.setSelected(True)
                self.comp_list.setCurrentItem(item)
                self.status.showMessage(f"{name} ticked from the diagram")
                return
        self.status.showMessage(f"{name} is not in this list")

    def _refresh_diagram(self):
        """Show the cell for a localization category, filled from the map's palette.

        Hidden for anything else: coloring organelles by cell-cycle phase would be a picture of a
        relationship that does not exist. The fills come from `color_of`, the same dict the points
        are drawn from, so the diagram and the map cannot disagree.
        """
        if self.diagram is None:
            return
        from .celldiagram import UNMAPPED_NOTE, missing_from_drawing
        show = self.category in ("compartment", "compartment_best")
        self.diagram.setVisible(show)
        # The AREA, not the label: hiding the label inside a fixed-height area leaves the area's
        # height behind, which is the gap this was supposed to remove.
        self.diagram_note_area.setVisible(show)
        if not show:
            return
        sel = self.comp_list.checked()
        self.diagram.set_palette(self.color_of, sel[0] if len(sel) == 1 else "")
        absent = [c for c in missing_from_drawing(self.comps, self.diagram.svg)
                  if str(c).lower() not in {str(x).lower() for x in ABSENCE}]
        note = self.diagram.showing()
        if not note and absent:
            # Named rather than dropped: between them these are a large part of the proteome, and a
            # compartment that vanishes from the legend reads as one that does not exist.
            note = f"{', '.join(absent)}: {UNMAPPED_NOTE}"
        self.diagram_note.setText(note)

    def _right(self):
        """Create the dock containing gene evidence and links to source records."""
        d = QtWidgets.QDockWidget("evidence")
        d.setFeatures(QtWidgets.QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
        self.detail = QtWidgets.QTextBrowser()
        self.detail.setOpenExternalLinks(True)
        self.detail.setHtml("<p style='color:#888'>Click a gene.</p>")
        self.detail.setMinimumWidth(400)
        d.setWidget(self.detail)
        self.right_dock = d          # the analysis dock tabs against this one
        return d

    # ------------------------------------------------------------------ drawing
    def color_sources(self) -> list:
        """Everything the map can be colored by, in one list: columns, runs, binned quantities.

        One list rather than three controls, because they answer the same question -- what should
        color mean right now -- and having to know which of three places to look for an answer is
        the state this panel replaced.
        """
        return (list(self.categories)
                + [RUN_PREFIX + r.name for r in self.runs.runs]
                + [BIN_PREFIX + c for c in self.numerics])

    def category_values(self) -> pd.Series:
        """The current color source as a string column over every gene.

        Absence is "" in all three cases and means the same thing each time: no value measured, no
        position in that run, no number to bin. It is drawn gray, never as a category.
        """
        src = self.category
        if src.startswith(RUN_PREFIX):
            run = self.runs.get(src[len(RUN_PREFIX):])
            return (run.values(self.n) if run is not None
                    else pd.Series([""] * self.n, dtype=object))
        if src.startswith(BIN_PREFIX):
            from .runs import bin_column
            col = src[len(BIN_PREFIX):]
            if col not in self.nodes.columns:
                return pd.Series([""] * self.n, dtype=object)
            return bin_column(self.nodes[col], self.bins,
                              log=lambda m: self.status.showMessage(f"{col}: {m}"))
        if src not in self.nodes.columns:
            return pd.Series([""] * self.n, dtype=object)
        return as_text(self.nodes[src])

    def keep_run(self, labels, recipe=None, name: str = ""):
        """Keep a clustering, list it in the color-by panel, and color the map by it.

        Kept rather than drawn and forgotten: the second run used to replace the first with no way
        back, which makes the one comparison this application is for -- does this structure survive
        different settings -- impossible to make by looking.
        """
        labels = np.asarray(labels)
        genes = getattr(self, "placed", None)
        genes = np.ones(self.n, bool) if genes is None else np.asarray(genes, bool)
        if len(labels) == self.n:
            genes = np.ones(self.n, bool)
        run = self.runs.add(labels, genes, recipe=recipe or {}, name=name)
        self.runs.save(run)
        self._refresh_sources()
        self.category_box.setCurrentText(RUN_PREFIX + run.name)
        self.status.showMessage(f"kept as {run.describe()}")
        return run

    def _refresh_sources(self):
        """Rebuild the color-by list, keeping the current choice if it still exists."""
        want = self.category_box.currentText()
        self.category_box.blockSignals(True)
        self.category_box.clear()
        self.category_box.addItems(self.color_sources())
        i = self.category_box.findText(want)
        self.category_box.setCurrentIndex(max(i, 0))
        self.category_box.blockSignals(False)

    def rename_run(self, old: str, new: str) -> bool:
        """Rename a kept run and keep the panel pointing at it."""
        if not self.runs.rename(old, new):
            self.status.showMessage(f"cannot rename to {new!r} -- that name is taken")
            return False
        self._refresh_sources()
        self.category_box.setCurrentText(RUN_PREFIX + new)
        self.status.showMessage(f"renamed to {new}, saved with its recipe")
        return True

    def visible_mask(self):
        """Boolean mask of the genes passing the current class filter AND having a position.

        The two are combined here rather than at each drawing site because everything downstream --
        points, edges, the gene count in the status bar -- reads this one mask. A gene the
        displayed embedding does not cover has no position to draw, and drawing it anyway would put
        absence on the map as though it were a measurement.
        """
        sel = self.comp_list.checked()
        vis = (np.ones(self.n, bool) if not sel
               else self.category_values().isin(sel).to_numpy())
        placed = getattr(self, "placed", None)
        return vis if placed is None else (vis & placed)

    def colors(self, vis):
        """An RGBA color per gene under the current color mode. Gray always means unknown."""
        mode = self.color_mode
        c = np.zeros((self.n, 4), dtype=np.float32)
        if mode.startswith("compartment"):
            # The default colors whatever the color-by panel names -- a column, a kept clustering
            # or a binned quantity -- because that panel IS the choice of what color means. The
            # "incl. transferred" variant is the one exception, a preset that names a column of its
            # own: ortholog-transferred labels are inferences from another species, offered because
            # coverage matters and kept separate because provenance matters more.
            if "transferred" in mode and "compartment_best" in self.nodes.columns:
                vals = as_text(self.nodes["compartment_best"])
            else:
                vals = self.category_values()
            for comp, col in self.color_of.items():
                c[(vals == comp).to_numpy(), :3] = col
            for absent in ("unassigned", ""):
                c[(vals == absent).to_numpy(), :3] = TH.unknown_color(self.theme)[:3]
        elif mode == "clusters":
            # Noise stays grey, with everything else that is unknown. HDBSCAN calling a gene
            # unclustered is a finding about that gene, not a gap in the drawing.
            if self.cluster_labels is None or len(self.cluster_labels) != self.n:
                c[:, :3] = TH.unknown_color(self.theme)[:3]
            else:
                lab = self.cluster_labels
                ids = sorted(set(lab[lab >= 0]))
                palette = TH.categorical_colors(max(len(ids), 1), self.theme, self.cmap_name)
                for col, k in zip(palette, ids):
                    c[lab == k, :3] = col
                c[lab < 0, :3] = TH.unknown_color(self.theme)[:3]
        elif mode == "annotations":
            # The fourth color, used for nothing else. Everything unannotated is grey -- not a
            # category, not zero: "nobody has proposed anything for this gene".
            from .annotations import ANNOTATION_COLOR
            c[:, :3] = TH.unknown_color(self.theme)[:3]
            if self.annotated is not None and len(self.annotated) == self.n:
                c[self.annotated, :3] = ANNOTATION_COLOR
        elif mode == "depth of attention":
            # Categorical, not a scale: these tiers are read off document structure (title / abstract /
            # body-only), so shading them along a gradient would imply a quantity that does not exist.
            # Never-named genes stay grey with everything else that is unknown rather than absent.
            if "attention_depth" not in self.nodes.columns:
                c[:, :3] = TH.unknown_color(self.theme)[:3]
            else:
                d = as_text(self.nodes.attention_depth).to_numpy()
                for tier, col in DEPTH_COLOR.items():
                    c[d == tier, :3] = col
                c[d == "", :3] = TH.unknown_color(self.theme)[:3]
        else:
            col = {"in vitro fitness": "fit_invitro_hff", "publications": "n_publications",
                   "structure confidence (pLDDT)": "mean_plddt"}.get(mode)
            if mode == "cyst / tachyzoite expression":
                v = (self.nodes.expr_cyst - self.nodes.expr_tachy).to_numpy(dtype=float)
            else:
                v = self.nodes[col].to_numpy(dtype=float)
                if mode == "publications":
                    v = np.log10(v + 1.0)
            ok = np.isfinite(v)
            if ok.sum():
                lo, hi = np.nanpercentile(v[ok], [2, 98])
                t = np.clip((v - lo) / max(hi - lo, 1e-9), 0, 1)
                cm = TH.resolve_cmap(self._cmap_for(v))
                c[:, :3] = cm.map(np.nan_to_num(t, nan=0.0), mode="float")[:, :3]
            c[~ok, :3] = TH.unknown_color(self.theme)[:3]   # missingness stays grey, never a value
        alpha = float(TH.POINT_STYLES.get(self.point_style, {}).get('alpha', 0.95))
        c[:, 3] = np.where(vis, alpha, 0.06)
        if self.sel is not None:
            c[self.sel] = (1.0, 1.0, 1.0, 1.0)
        if getattr(self, "gated", None) is not None and len(self.gated):
            # A gate is a SELECTION, not a claim about the data, so it must not recolor anything --
            # measurement, inference, absence and annotation are what color means here. The gated
            # genes keep their own color and everything else recedes.
            m = np.zeros(self.n, dtype=bool)
            m[self.gated] = True
            c[~m, 3] *= 0.12
        if getattr(self, "placed", None) is not None:
            # Fully transparent, not dimmed. A gene the displayed embedding has no coordinates for is
            # absent from it, and a faint point is still a point -- at 0.06 alpha, several thousand of
            # them overlapping read as a real feature of the map.
            c[~self.placed, 3] = 0.0
        return c

    def redraw(self):
        """Redraw everything from the current state: gene points, ground, selection halo and edges."""
        vis = self.visible_mask()

        for it in self.edge_items:
            self.view.removeItem(it)
        self.edge_items = []
        st = TH.POINT_STYLES.get(self.point_style, TH.POINT_STYLES[TH.DEFAULT_POINT_STYLE])
        # An explicit point size overrides the style's. The style still supplies the rest of its
        # look, so "large" plus an explicit 4 px is a large-style point drawn at 4 px.
        base = float(self.point_size if self.point_size is not None else st["size"])
        # redraw() sets sizes every time, so the chosen point style has to be applied here rather than
        # only in apply_point_style -- otherwise picking a style changed nothing the moment anything
        # else triggered a redraw.
        sizes = np.where(vis, base, max(base * 0.4, 1.5)).astype(np.float32)
        if getattr(self, "gated", None) is not None and len(self.gated):
            sizes[self.gated] = max(base * 1.7, 7.0)
        if self.sel is not None:
            sizes[self.sel] = max(base * 3.0, 12.0)
        colors = self.colors(vis)
        self._depth_ctx = None
        if self.depth_cue:
            # Fade and shrink with distance so the cloud has a front and a back. Without it the map is
            # a flat disc of color however much it is rotated.
            cam = self.view.cameraPosition()
            cam = (cam.x(), cam.y(), cam.z())
            d = np.linalg.norm(self.xyz - np.asarray(cam, np.float32), axis=1)
            # One range shared by points and edges. Normalizing each set against its own extent makes
            # an edge and the point it touches fade by different amounts, which reads as flicker.
            self._depth_ctx = (cam, (float(d.min()), float(d.max())))
            a, sizes = TH.depth_cue(self.xyz, cam, colors[:, 3], sizes, rng=self._depth_ctx[1])
            colors = colors.copy()
            colors[:, 3] = a
            sizes = sizes.astype(np.float32)

        # The flat colours are kept as they are: lighting modulates a COPY of them every frame, and
        # shading an already-shaded array would darken the map a little more on each pass until the
        # whole thing went black.
        self._base_colors, self._base_sizes = colors, sizes
        lit = None
        if self._lighting["mode"] != "off":
            from .lighting import shade
            lit = self.frame_lights()
            if self.scatter.gpu_material_enabled(self._lighting["point_mode"]):
                colors = np.array(colors, copy=True)
            else:
                colors = shade(self.xyz, colors, lit,
                               point_mode=self._lighting["point_mode"], eye=self._eye())
        # The ball each gene is drawn as. A display setting rather than a lighting one, so it is
        # built on every redraw and shows whether the moving lights are running or not.
        self._refresh_sprite(lit)
        self._draw_emitter(lit)
        from .lighting import POINT_MODES
        render_sizes = sizes * POINT_MODES[self._lighting["point_mode"]]["size"]
        self.scatter.setData(pos=self.xyz, color=colors, size=render_sizes)
        self._draw_ground()
        if self._lighting["mode"] != "off":
            # After `_draw_ground`, which builds a new grid item with the flat theme colour: lighting
            # it before that is lighting an object that is about to be replaced, which is why the
            # grid never appeared to be lit.
            self._light_ground(lit)
        self._draw_selection_halo()
        self.draw_edges(vis)

        act = [k for k, _ in EDGE_TYPES if self.edge_on.get(k)]
        note = "" if self.attn_on else "  ·  RAW co-mention (attention-biased)"
        if not self.edge_on.get("comention", False) and not act:
            note += "  ·  select a gene to see its edges, or Edges ▸ Draw all active edges"
        self.status.showMessage(
            f"{int(vis.sum()):,} / {self.n:,} genes shown  ·  edges: "
            f"{', '.join(act) if act else 'none'}{note}")

    def _draw_ground(self):
        """A faint grid under the cloud, so rotation has something to rotate against.

        A point cloud alone in black has no reference: spinning it, the eye cannot tell rotation from
        the points rearranging themselves, and depth cueing alone does not fix that because it gives
        near/far but not orientation. A horizon does.
        """
        if self.grid_item is not None:
            self.view.removeItem(self.grid_item)
            self.grid_item = None
        if not self.show_ground:
            return
        r = self.view.data_radius()
        centre = self.xyz.mean(0) if len(self.xyz) else np.zeros(3)
        g = gl.GLGridItem()
        g.setSize(x=r * 2.4, y=r * 2.4)
        g.setSpacing(x=r / 4.0, y=r / 4.0)
        # Sit it just below the lowest point rather than at the origin: the embedding is not centred on
        # zero and a grid cutting through the cloud reads as an artefact of the data.
        low = float(self.xyz[:, 2].min()) if len(self.xyz) else 0.0
        g.translate(float(centre[0]), float(centre[1]), low - r * 0.08)
        p = TH.palette_for(self.theme)
        rgb = TH.rgbf(p["border"])[:3]
        g.setColor((int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255),
                    70 if TH.is_light(p) else 40))
        g.setGLOptions("translucent")
        self.view.addItem(g)
        self.grid_item = g

    def _draw_selection_halo(self):
        """A soft ring behind the selected gene.

        A size bump alone is invisible in a dense region -- the neighbours are the same color and the
        selected point simply becomes a slightly bigger dot in a crowd. A translucent halo in the
        accent color separates it from its neighborhood at any density.
        """
        if self.halo_item is not None:
            self.view.removeItem(self.halo_item)
            self.halo_item = None
        if self.sel is None:
            return
        p = TH.palette_for(self.theme)
        base = float(TH.POINT_STYLES.get(self.point_style, {}).get("size", 5.0))
        rings = [(base * 6.0, 0.16), (base * 3.6, 0.30), (base * 2.2, 0.55)]
        if TH.is_light(p):
            # A translucent ring on a near-white ground has almost no contrast, so the same halo that
            # is unmistakable on dark is nearly invisible on light. Contrast has to go the other way:
            # more opacity, and the darker accent rather than the brighter one.
            rings = [(sz, min(a * 1.8, 0.92)) for sz, a in rings]
        pos = np.repeat(self.xyz[self.sel][None, :], len(rings), axis=0).astype(np.float32)
        key = "accent_lo" if TH.is_light(p) and "accent_lo" in p else "accent"
        col = np.array([(*TH.rgbf(p[key])[:3], a) for _, a in rings], dtype=np.float32)
        self.halo_item = gl.GLScatterPlotItem(
            pos=pos, color=col, size=np.array([s for s, _ in rings], np.float32), pxMode=True)
        self.halo_item.setGLOptions("translucent")
        self.view.addItem(self.halo_item)

    def draw_edges(self, vis):
        """Draw the active edge layers, capped and alpha-budgeted so they do not hide the map."""
        active = [k for k, _ in EDGE_TYPES if self.edge_on.get(k) and k in self.edges]
        if not active:
            return
        for k in active:
            e = self.edges[k]
            a, b = e["a"], e["b"]
            w = e["r"] if (k in COMENTION and self.attn_on) else e["w"]
            keep = vis[a] & vis[b]
            if k in COMENTION and self.attn_on:
                keep &= w > 0          # corrected mode shows only more-than-expected pairs
            if self.sel is not None and not self.all_edges_on:
                keep &= (a == self.sel) | (b == self.sel)
            elif not self.all_edges_on:
                continue               # nothing selected, whole map not requested: draw nothing
            idx = np.where(keep)[0]
            if idx.size == 0:
                continue
            if idx.size > 20000:       # cap is on drawing only; the cap is stated, not silent
                idx = idx[np.argsort(-w[idx])[:20000]]
                self.status.showMessage(f"{k}: showing strongest 20,000 of {int(keep.sum()):,} edges")
            seg = np.empty((idx.size * 2, 3), np.float32)
            seg[0::2] = self.xyz[a[idx]]
            seg[1::2] = self.xyz[b[idx]]
            col = {"comention": (0.95, 0.85, 0.35, 0.5), "comention_ft": (0.95, 0.62, 0.25, 0.45),
                   "orthogroup": (0.35, 0.85, 0.55, 0.5),
                   "coexpression": (0.40, 0.65, 0.95, 0.45), "compartment": (0.75, 0.75, 0.80, 0.25),
                   "cofitness": (0.95, 0.45, 0.75, 0.5), "domain": (0.60, 0.55, 0.45, 0.3),
                   # measured physical evidence gets its own cool, high-contrast family
                   "xlms": (0.30, 0.95, 0.90, 0.65), "ip_ms": (0.20, 1.00, 0.55, 0.90),
                   "struct": (0.70, 0.80, 0.30, 0.45),
                   # deliberately the loudest colors in the palette: these are the things to look at
                   "structural_hole": (1.00, 0.25, 0.25, 0.85),
                   "unwritten_interaction": (1.00, 0.55, 0.00, 0.90)}[k]
            # Fade each edge by its own weight. Drawn at one flat alpha, 7,733 full-text edges are an
            # opaque hairball in which the strongest and the weakest look identical -- which also made
            # the attention toggle almost invisible, though it reorders exactly this quantity. Scaling
            # alpha by weight is what lets the corrected view read differently from the raw one.
            cols = np.empty((idx.size * 2, 4), np.float32)
            cols[:, :3] = col[:3]
            ww = w[idx].astype(float)
            if idx.size > 20 and np.ptp(ww) > 0:
                lo, hi = np.percentile(ww, [10, 95])
                t = np.clip((ww - lo) / max(hi - lo, 1e-9), 0.0, 1.0)
            else:
                t = np.ones(idx.size)
            # Ink budget. Alpha per edge is not enough on its own: 20,000 translucent lines over the
            # same region sum to an opaque sheet, and the map underneath disappears entirely. (Rendered
            # with cofitness on, the whole cloud was one flat pink shape.) Total ink is what has to be
            # bounded, so alpha falls as the count rises -- sqrt, because coverage grows roughly with
            # line count over a fixed area. A small selected neighbourhood is untouched.
            ink = float(np.clip(np.sqrt(EDGE_INK_TARGET / max(idx.size, 1)), 0.12, 1.0))
            alpha = (col[3] * ink * (0.12 + 0.88 * t)).astype(np.float32)
            cols[0::2, 3] = alpha
            cols[1::2, 3] = alpha
            if self._depth_ctx is not None:
                # Fade by camera distance as well as by weight. In a dense layer every edge crosses
                # every other one; without the distance term the near neighbourhood -- the only part
                # anyone can actually trace -- is buried under edges from the far side of the cloud.
                cam, rng = self._depth_ctx
                t = TH.depth_t(seg, cam, rng)
                cols[:, 3] *= (1.0 - t * (1.0 - TH.EDGE_DEPTH_FADE)).astype(np.float32)
            it = gl.GLLinePlotItem(pos=seg, color=cols, width=1.0, mode="lines", antialias=True)
            self.view.addItem(it)
            self.edge_items.append(it)

    # ------------------------------------------------------------------ interaction
    def on_pick(self, i):
        """Select a gene and show its evidence."""
        self.sel = int(i)
        self.show_detail(self.sel)
        self.redraw()

    def do_search(self):
        """Find a gene by accession or product text and fly to it."""
        q = self.search.text().strip().lower()
        if not q:
            return
        gid = as_text(self.nodes.gene_id).str.lower()
        hit = np.where(gid == q)[0]
        if hit.size == 0:
            hit = np.where(gid.str.contains(q, regex=False))[0]
        if hit.size == 0:
            hit = np.where(as_text(self.nodes["product"]).str.lower()
                           .str.contains(q, regex=False))[0]
        if hit.size == 0:
            self.status.showMessage(f"no match for {q!r}")
            return
        self.on_pick(int(hit[0]))
        p = self.xyz[self.sel]
        self.view.setCameraPosition(pos=pg.Vector(*p), distance=45)
        self.status.showMessage(f"{hit.size} match(es); showing {self.nodes.gene_id.iloc[hit[0]]}")

    def fly_to_compartment(self, item):
        """Centre the camera on the genes in the selected category."""
        c = item.data(QtCore.Qt.ItemDataRole.UserRole)
        m = (self.category_values() == c).to_numpy()
        if m.sum():
            self.view.setCameraPosition(pos=pg.Vector(*self.xyz[m].mean(0)), distance=60)
            self.redraw()

    def reset(self):
        """Clear the selection and every filter, and frame the whole map again."""
        self.sel = None
        # Ticks are the filter now, so clearing the highlight alone would leave the map filtered while
        # the button said it had reset it.
        self.comp_list.set_checked([], emit=False)
        self.comp_list.clearSelection()
        self.view.fit_view()
        self.detail.setHtml("<p style='color:#888'>Click a gene.</p>")
        self.redraw()

    # ------------------------------------------------------------------ export
    def _ask_path(self, caption: str, filt: str, default: str) -> str:
        """Return the chosen export filename, or an empty string if cancelled."""
        path, _ = QtWidgets.QFileDialog.getSaveFileName(self, caption, default, filt)
        return path

    def export_image(self, path: str = ""):
        """The current frame as a PNG, at the size it is on screen."""
        path = path or self._ask_path("Export image", "PNG image (*.png)", "starplast.png")
        if not path:
            return None
        # Guarded, because grabbing a framebuffer needs a real GL context and there is not always
        # one: under the offscreen platform pyqtgraph itself warns that QOpenGLWidget is
        # unsupported, and grabbing there is undefined -- it segfaulted rather than failing. An
        # export that cannot happen must say so, not take the application down.
        img = None
        try:
            if self.view.isValid():
                img = self.view.grabFramebuffer()
        except Exception as exc:
            self.status.showMessage(f"could not capture the view: {type(exc).__name__}: {exc}")
            return None
        if img is None or img.isNull():
            self.status.showMessage("no OpenGL context to capture — export needs a real display")
            return None
        ok = img.save(path)
        self.status.showMessage(f"wrote {path}" if ok else f"could not write {path}")
        return path if ok else None

    def default_export_columns(self) -> list[str]:
        """What to tick when the picker opens: identity, the active filter, the active coloring.

        Taken from the current view rather than fixed, because the columns worth exporting are
        usually the ones being looked at.
        """
        want = ["gene_id", "product", self.category]
        want += {"in vitro fitness": ["fit_invitro_hff"],
                 "publications": ["n_publications"],
                 "structure confidence (pLDDT)": ["mean_plddt"],
                 "depth of attention": ["attention_depth"],
                 "cyst / tachyzoite expression": ["expr_cyst", "expr_tachy"],
                 }.get(self.color_mode, ["compartment"])
        seen, out = set(), []
        for c in want:
            if c in self.nodes.columns and c not in seen:
                seen.add(c)
                out.append(c)
        return out

    def choose_export_columns(self, preselect=None):
        """A checkable list of every column, for picking what an export carries.

        Built rather than shown, for the same reason as the context menu: exec() enters a modal loop
        and does not return until a human closes it.
        """
        d = QtWidgets.QDialog(self)
        d.setWindowTitle("Columns to export")
        L = QtWidgets.QVBoxLayout(d)
        L.addWidget(QtWidgets.QLabel(
            "The gene id and coordinates are always written. Tick anything else to carry along."))
        filt = QtWidgets.QLineEdit(placeholderText="filter columns…")
        filt.setToolTip("Narrow the list by name. Ticks are kept underneath, so filtering never "
                        "silently drops a column you had already chosen.")
        L.addWidget(filt)
        lst = QtWidgets.QListWidget()
        lst.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.NoSelection)
        chosen = set(preselect if preselect is not None else self.default_export_columns())
        for c in self.nodes.columns:
            it = QtWidgets.QListWidgetItem(str(c))
            it.setFlags(it.flags() | QtCore.Qt.ItemFlag.ItemIsUserCheckable)
            it.setCheckState(QtCore.Qt.CheckState.Checked if c in chosen
                             else QtCore.Qt.CheckState.Unchecked)
            lst.addItem(it)
        L.addWidget(lst, 1)

        def apply_filter(text):
            """Hide export columns whose names do not match the text; preserve their checks."""
            t = text.strip().lower()
            for i in range(lst.count()):
                lst.item(i).setHidden(bool(t) and t not in lst.item(i).text().lower())
        filt.textChanged.connect(apply_filter)

        row = QtWidgets.QHBoxLayout()
        for label, state in (("all", QtCore.Qt.CheckState.Checked),
                             ("none", QtCore.Qt.CheckState.Unchecked)):
            btn = QtWidgets.QPushButton(label)
            btn.setToolTip(f"Tick or untick every column the filter is currently showing. Acts on "
                           f"the visible ones only, so '{label}' after filtering does not disturb "
                           f"choices you cannot see.")
            # Only what the filter is showing, so "none" after filtering clears that group rather
            # than silently discarding ticks the user cannot currently see.
            btn.clicked.connect(lambda _c, s=state: [lst.item(i).setCheckState(s)
                                                     for i in range(lst.count())
                                                     if not lst.item(i).isHidden()])
            row.addWidget(btn)
        L.addLayout(row)
        bb = QtWidgets.QDialogButtonBox(
            QtWidgets.QDialogButtonBox.StandardButton.Ok
            | QtWidgets.QDialogButtonBox.StandardButton.Cancel)
        bb.accepted.connect(d.accept)
        bb.rejected.connect(d.reject)
        L.addWidget(bb)
        d.column_list = lst          # so a caller can read the ticks without entering exec()
        return d

    @staticmethod
    def _ticked(dialog) -> list[str]:
        """Return checked export column names, including rows hidden by the filter."""
        lst = dialog.column_list
        return [lst.item(i).text() for i in range(lst.count())
                if lst.item(i).checkState() == QtCore.Qt.CheckState.Checked]

    def export_visible(self, path: str = "", columns=None):
        """The genes currently passing the filter, with their coordinates and chosen columns.

        Coordinates travel with the rows because a gene list without them cannot reproduce what was
        on screen, and reproducing what was on screen is the point of exporting it.
        """
        path = path or self._ask_path("Export visible genes", "CSV (*.csv)", "starplast_genes.csv")
        if not path:
            return None
        if columns is None:
            d = self.choose_export_columns()
            if not d.exec():
                return None
            columns = self._ticked(d)
        vis = self.visible_mask()
        cols = ["gene_id"] + [c for c in columns if c != "gene_id" and c in self.nodes.columns]
        out = self.nodes.loc[vis, cols].copy()
        out["x"], out["y"], out["z"] = self.xyz[vis, 0], self.xyz[vis, 1], self.xyz[vis, 2]
        out.to_csv(path, index=False)
        self.status.showMessage(f"wrote {len(out):,} genes and {len(cols)} columns to {path}")
        return path

    def export_relationships(self, path: str = "", columns=None):
        """Every active edge between two visible genes, one row per edge, as CSV.

        The GraphML export is for Cytoscape; this is for a spreadsheet, and it is what "export the
        relationships I can see" actually means. Both endpoints carry the chosen columns, suffixed
        `_a` and `_b`, so the file can be read without joining it back against anything.

        The edge type is a column and is never collapsed. A row saying only that two genes are
        related at weight 0.8 has lost whether they were measured touching or merely mentioned
        together, which is the distinction the whole application is built around.
        """
        path = path or self._ask_path("Export relationships", "CSV (*.csv)", "starplast_edges.csv")
        if not path:
            return None
        if columns is None:
            d = self.choose_export_columns()
            if not d.exec():
                return None
            columns = self._ticked(d)
        vis = self.visible_mask()
        cols = [c for c in columns if c in self.nodes.columns and c != "gene_id"]
        active = [k for k, _ in EDGE_TYPES if self.edge_on.get(k) and k in self.edges]
        gid = self.nodes.gene_id.to_numpy()
        frames = []
        for k in active:
            e = self.edges[k]
            a, b = e["a"], e["b"]
            corrected = k in COMENTION and self.attn_on
            w = e["r"] if corrected else e["w"]
            keep = vis[a] & vis[b]
            if corrected:
                keep &= w > 0
            if self.sel is not None and not self.all_edges_on:
                # Match what is drawn. Exporting a whole layer while the screen shows one gene's
                # neighbourhood hands back a different graph from the one being looked at.
                keep &= (a == self.sel) | (b == self.sel)
            idx = np.flatnonzero(keep)
            if not idx.size:
                continue
            f = pd.DataFrame({
                "gene_a": gid[a[idx]], "gene_b": gid[b[idx]], "edge_type": k,
                "weight": w[idx],
                "weight_kind": "attention-corrected log2 obs/exp" if corrected else "raw",
            })
            for side, ends in (("a", a[idx]), ("b", b[idx])):
                for c in cols:
                    f[f"{c}_{side}"] = self.nodes[c].to_numpy()[ends]
            frames.append(f)
        if not frames:
            self.status.showMessage("no edges visible to export — turn on an edge type, then select "
                                    "a gene or use Edges ▸ Draw all active edges")
            return None
        out = pd.concat(frames, ignore_index=True)
        out.to_csv(path, index=False)
        self.status.showMessage(f"wrote {len(out):,} relationships over {len(active)} edge "
                                f"type(s) to {path}")
        return path

    def export_graphml(self, path: str = ""):
        """The active edge types over the visible genes, as GraphML for Cytoscape or networkx.

        Each edge keeps its type as an attribute rather than being merged into a generic edge, for
        the same reason the layers are never combined on screen -- see Edges ▸ Why are these never
        combined? An exported file outlives the session, so it is exactly where a merged score would
        do the most damage.
        """
        path = path or self._ask_path("Export graph", "GraphML (*.graphml)", "starplast.graphml")
        if not path:
            return None
        vis = self.visible_mask()
        active = [k for k, _ in EDGE_TYPES if self.edge_on.get(k) and k in self.edges]
        idx = np.flatnonzero(vis)
        gid = as_text(self.nodes.gene_id).to_numpy()
        from xml.sax.saxutils import escape as _esc
        lines = ['<?xml version="1.0" encoding="UTF-8"?>',
                 '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">',
                 '<key id="label" for="node" attr.name="label" attr.type="string"/>',
                 '<key id="etype" for="edge" attr.name="type" attr.type="string"/>',
                 '<key id="w" for="edge" attr.name="weight" attr.type="double"/>',
                 '<graph edgedefault="undirected">']
        for i in idx:
            lines.append(f'<node id="n{i}"><data key="label">{_esc(gid[i])}</data></node>')
        n_edges = 0
        for k in active:
            e = self.edges[k]
            a, b = e["a"], e["b"]
            w = e["r"] if (k in COMENTION and self.attn_on) else e["w"]
            keep = vis[a] & vis[b]
            for j in np.flatnonzero(keep):
                lines.append(f'<edge source="n{a[j]}" target="n{b[j]}">'
                             f'<data key="etype">{k}</data>'
                             f'<data key="w">{float(w[j]):.6g}</data></edge>')
                n_edges += 1
        lines += ["</graph>", "</graphml>"]
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(lines))
        self.status.showMessage(f"wrote {len(idx):,} nodes and {n_edges:,} edges to {path}")
        return path

    def describe_state(self) -> str:
        """A briefing for the assistant: what is on screen right now."""
        vis = self.visible_mask()
        bits = ["View: individual genes.",
                f"Coloring: {self.color_mode}.",
                f"Filter column: {self.category}.",
                f"{int(vis.sum()):,} of {self.n:,} genes visible."]
        active = [k for k, _ in EDGE_TYPES if self.edge_on.get(k)]
        bits.append(f"Edge layers on: {', '.join(active) if active else 'none'}"
                    + ("" if self.attn_on else " (co-mention shown RAW, attention-biased)") + ".")
        if self.sel is not None:
            r = self.nodes.iloc[self.sel]
            bits.append(f"\nSelected gene: {r.gene_id} — {r.get('product', '?')}")
            for c in ("compartment", "cellcycle_phase", "n_publications", "mean_plddt"):
                if c in self.nodes.columns:
                    bits.append(f"  {c}: {r[c]}")
        else:
            bits.append("No gene is selected.")
        return "\n".join(bits)

    # ------------------------------------------------------------------ detail
    def show_detail(self, i):
        """Fill the evidence panel for one gene, distinguishing absence from zero throughout."""
        r = self.nodes.iloc[i]
        gid = str(r.gene_id)

        def num(v, f="{:.2f}"):
            """Format a measured value, or show a dash for missing and non-finite values."""
            return "—" if v is None or pd.isna(v) or not np.isfinite(v) else f.format(v)

        compartment = r.get("compartment", "unassigned")
        if pd.isna(compartment):
            compartment = "unassigned"
        localization_label = "compartment (hyperLOPIT)" if "compartment" in r else "localization"
        rows = [(localization_label,
                 f"{compartment}" + (" <i>— unknown, not absent</i>"
                                     if compartment == "unassigned" else "")),
                ("orthogroup", str(r.get("orthogroup", "—"))),
                ("paralogs", num(r.get("paralog_number"), "{:.0f}")),
                ("InterPro domains", num(r.get("n_interpro"), "{:.0f}")),
                ("phosphosites", num(r.get("n_phosphosites"), "{:.0f}")),
                ("mean pLDDT", num(r.get("mean_plddt"))),
                ("abstracts naming it", num(r.get("n_publications"), "{:.0f}")),
                ("open-access full texts naming it", num(r.get("n_fulltext"), "{:.0f}")
                 + ("" if not r.get("lit_tier") else
                    f" <i>— by {r.get('lit_tier')}</i>")),
                ("papers with it in the title", num(r.get("n_papers_focal"), "{:.0f}")),
                ("papers with it in the abstract", num(r.get("n_papers_substantive"), "{:.0f}")),
                ("papers naming it only in passing",
                 num(r.get("n_papers_incidental"), "{:.0f}")),
                ("structural holes", num(r.get("n_holes"), "{:.0f}")
                 + (" <i>— genes it behaves like but is never discussed with</i>"
                    if r.get("n_holes", 0) else "")),
                ("crosslinked partners (XL-MS)", num(r.get("n_xlink_partners"), "{:.0f}")),
                ("IP-MS partners", num(r.get("n_ipms_partners"), "{:.0f}")),
                ("structurally similar (TM ≥ 0.7)", num(r.get("n_struct_similar"), "{:.0f}")),
                ("curated host targets", num(r.get("n_host_targets"), "{:.0f}"))]
        for column, label in (("expr_tachy", "log2 FPKM tachyzoite"),
                              ("expr_cyst", "log2 FPKM tissue cyst")):
            if column in r:
                rows.append((label, num(r[column])))
        tbl = "".join(f"<tr><td style='color:#888;padding-right:10px'>{k}</td>"
                      f"<td>{v}</td></tr>" for k, v in rows)

        fit = "".join(
            f"<tr><td style='color:#888;padding-right:10px'>{c.replace('fit_', '')}</td>"
            f"<td>{num(r.get(c))}</td></tr>" for c in FIT if c in self.nodes.columns)

        # Published screens, each with its own scope. Targeted libraries leave most genes untested, and
        # untested is shown as "—", never as zero effect.
        pub = [("GRA17 synthetic-lethal Δ", "crispr_gra17_synthlethal_delta"),
               ("GRA17 candidate", "crispr_gra17_candidate"),
               ("GRA12 screen 1 in vivo L2FC", "crispr_gra12s1_l2fc_invivo"),
               ("GRA12 screen 2 in vivo L2FC", "crispr_gra12s2_l2fc_invivo"),
               ("in vivo platform mean lfc", "crispr_invivo_platform_lfc"),
               ("host-transcription T²", "hosttx_T2"),
               ("protein abundance log2 iBAQ", "protein_ibaq_log2"),
               ("log2 FPKM sporulated oocyst", "expr_sporulated")]
        pub = "".join(
            f"<tr><td style='color:#888;padding-right:10px'>{lab}</td>"
            f"<td>{num(r.get(c))}</td></tr>" for lab, c in pub if c in self.nodes.columns)

        nb = []
        for k, label in EDGE_TYPES:
            if k not in self.edges:
                continue
            e = self.edges[k]
            m = (e["a"] == i) | (e["b"] == i)
            if not m.any():
                continue
            w = e["r"] if (k in COMENTION and self.attn_on) else e["w"]
            part = np.where(e["a"][m] == i, e["b"][m], e["a"][m])
            ww = w[m]
            o = np.argsort(-ww)[:8]
            names = ", ".join(f"{self.nodes.gene_id.iloc[int(part[j])]} ({ww[j]:.2f})" for j in o)
            nb.append(f"<p><b>{label}</b> — {int(m.sum())} edges<br>"
                      f"<span style='color:#aaa;font-size:11px'>{names}</span></p>")

        # How the measured binding is thought to happen: the crosslinked residues, the predicted
        # complexes, and whether those complexes place the crosslinks within reach. A model that does not
        # satisfy them is reported as such rather than dropped -- it says the model fails to explain the
        # measurement, not that the measurement is wrong.
        xl = ""
        if len(self.models):
            m = self.models[(self.models.gene_a == gid) | (self.models.gene_b == gid)]
            if len(m):
                items = []
                for row in m.sort_values("n_crosslinks", ascending=False).head(6).itertuples():
                    other = row.gene_b if row.gene_a == gid else row.gene_a
                    pos = json.loads(row.crosslink_positions or "[]")[:3]
                    res = ", ".join(f"{p[0]}–{p[1]}" for p in pos if p and p[0] is not None)
                    agree = ("model does not place them in contact"
                             if isinstance(row.frac_satisfied, float) and row.frac_satisfied == 0
                             else f"{row.frac_satisfied:.0%} of crosslinks satisfied"
                             if np.isfinite(row.frac_satisfied) else "not modeled")
                    nm = int(row.n_models) if np.isfinite(row.n_models) else 0
                    ok = bool(getattr(row, "model_trustworthy", False))
                    items.append(
                        f"<li>{other} — {row.n_crosslinks} crosslink(s)"
                        + (f" at residues {res}" if res else "")
                        + ("  <b style='color:#6c6'>model usable</b>" if ok else "")
                        + f"<br><span style='color:#888'>{agree}"
                        + (f" · {nm} Chai-1 models" if nm else "")
                        + (f" · ipTM {row.chai_iptm:.2f}" if np.isfinite(row.chai_iptm) else "")
                        + "</span></li>")
                where = m.model_dir.dropna().iloc[0] if m.model_dir.notna().any() else None
                xl = ("<p><b>How the binding is modeled</b> <span style='color:#888;"
                      "font-weight:normal;font-size:11px'>— the crosslink is the measurement; the "
                      "model is a guess at the pose, and 60% of them explain no crosslink at all"
                      "</span></p><ul style='margin-top:2px'>"
                      + "".join(items) + "</ul>"
                      + (f"<p style='color:#666;font-size:11px'>structures: {where}/"
                         f"&lt;id&gt;_model_&lt;0-3&gt;.cif</p>" if where else ""))

        # Being named is not being studied. A gene reached only through a screen's hit table would
        # otherwise read as attended-to simply because coverage counts every tier alike.
        if r.get("attention_depth", "") == "incidental":
            att = ("<p style='color:#c9a227'><b>Listed, not studied</b> — named in "
                   f"{int(r.get('n_papers_incidental', 0))} paper(s), never in a title or abstract. "
                   "Most such mentions are entries in a screen's hit table.</p>")
        elif r.get("n_publications", 0) > 6:
            att = ""
        else:
            att = ("<p style='color:#c9a227'><b>Effectively uncharacterised</b> — named in ≤6 abstracts"
                   + (f", {int(r.get('n_fulltext', 0))} open-access full texts"
                      if r.get("n_fulltext", 0) else "")
                   + ". Absence of evidence here is absence of attention, not absence of function.</p>")

        record_database = "PlasmoDB" if self.species == "Plasmodium falciparum" else "ToxoDB"
        record_base = ("https://plasmodb.org/plasmo" if record_database == "PlasmoDB"
                       else "https://toxodb.org/toxo")
        self.detail.setHtml(f"""
        <h2 style="margin-bottom:2px">{gid}</h2>
        <p style="color:#bbb;margin-top:0">{r.get('product', 'unannotated')}</p>
        {att}
        <table>{tbl}</table>
        <h4>CRISPR screens <span style="color:#888;font-weight:normal">— competitive growth,
        not essentiality; the screens do not agree with each other</span></h4>
        <table>{fit}</table>
        <h4>published screens &amp; abundance <span style="color:#888;font-weight:normal">— targeted
        libraries leave most genes untested; "—" means not measured, not no effect</span></h4>
        <table>{pub}</table>
        {xl}
        <h4>neighbours by edge type</h4>
        {''.join(nb) or '<p style="color:#888">no edges</p>'}
        <p><a href="{record_base}/app/record/gene/{gid}">{record_database} record</a> ·
           <a href="https://pubmed.ncbi.nlm.nih.gov/?term={gid}">PubMed</a></p>
        """)


def main():
    """Entry point for the `starplast` command: build the window and run the event loop."""
    # Check before building a window. A missing cache otherwise surfaces as a pandas error from inside
    # a constructor, which tells the user nothing about what to do next.
    ok, msg = paths.check()
    if not ok:
        print(f"starplast: {msg}", file=sys.stderr)
        print(paths.describe(), file=sys.stderr)
        sys.exit(1)
    pg.setConfigOptions(antialias=True)
    app = QtWidgets.QApplication(sys.argv)
    app.setApplicationName("starplast")
    w = Window()
    w.show()
    sys.exit(app.exec())


# No `if __name__ == "__main__"` guard here. There are already two ways in and this would be a third
# that nothing documents: the console script `starplast` points at `starplast.app:main`, and
# `starplast/__main__.py` makes `python -m starplast` work. Both are tested. A guard reachable only by
# `python -m starplast.app` is a line no test can execute without opening a real window, and an
# unreachable line is evidence the line should not exist.
