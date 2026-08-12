#!/usr/bin/env python3
"""starplast — a 3D browser for the Toxoplasma knowledge map.

Nodes are genes, positioned by a precomputed UMAP embedding so that proximity means biological similarity.
Six edge types are kept separate and toggled independently. Co-mention edges default to their
attention-corrected form, because the raw form reproduces the literature's popularity contest rather than
biology (see ../HANDOFF.md, decision 3).

Level of detail follows the data, not invented tiers:
    compartment (26 hyperLOPIT classes) -> orthogroup / module -> gene -> that gene's evidence
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

from . import lod  # noqa: E402
from . import paths  # noqa: E402
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

COLOUR_MODES = ["compartment", "compartment (incl. transferred)", "clusters", "in vitro fitness",
                "publications", "depth of attention", "structure confidence (pLDDT)",
                "cyst / tachyzoite expression"]

# None means "follow the point style". The rest are absolute pixel sizes.
POINT_SIZES = [("Automatic", None), ("Tiny (2 px)", 2.0), ("Small (4 px)", 4.0),
               ("Medium (7 px)", 7.0), ("Large (11 px)", 11.0), ("Huge (16 px)", 16.0)]

EDGE_EXPLANATION = (
    "They are kept as separate layers rather than added together into one 'interaction' edge.\n\n"
    "The twelve types are not twelve measurements of the same thing. A crosslink-MS edge is a "
    "measured physical contact. A co-mention edge is two genes appearing in one abstract, which "
    "happens to popular genes far more than to related ones. A co-expression edge is a correlation "
    "across a stage series.\n\n"
    "Combining them would let the weakest inference borrow the credibility of the strongest "
    "measurement: a merged score of 0.8 cannot tell you whether two proteins were measured touching "
    "or merely mentioned together, and there is no honest weighting that recovers the difference "
    "afterwards. So they are drawn in separate colours and toggled independently, and the answer to "
    "'are these two related' is 'by which evidence'.")

MAP_EXPLANATION = (
    "Position is similarity in expression, fitness screens, protein features and literature "
    "co-mention, reduced to three dimensions by UMAP. That is all it is.\n\n"
    "Held-out testing on the full proteome found that no target is reliably recovered: cell cycle "
    "recovers 0 of 5 phases and localisation 1 of 24 compartments. Localisation scores below the "
    "negative control — the map reflects how much a gene has been studied better than it reflects "
    "where the protein is.\n\n"
    "So proximity here is a hypothesis to check, never evidence on its own. No predictions are "
    "issued from it. Grey means unknown; it never means zero and never means a category.")

# How many distinct values a column may have and still be offered as a filter category. Above this it
# is an identifier rather than a class -- orthogroup has 7,331 values, and a list that long is not a
# filter, it is a scrolling exercise.
MAX_CATEGORY_VALUES = 60
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
# Kept as a fallback: the categorical colour map chosen in the UI supersedes it.
PALETTE = [
    (0.90, 0.24, 0.24), (0.20, 0.55, 0.90), (0.25, 0.75, 0.35), (0.95, 0.65, 0.15),
    (0.65, 0.35, 0.85), (0.15, 0.80, 0.78), (0.95, 0.45, 0.70), (0.55, 0.75, 0.20),
    (0.85, 0.35, 0.10), (0.35, 0.45, 0.85), (0.10, 0.65, 0.50), (0.80, 0.80, 0.20),
    (0.60, 0.20, 0.45), (0.30, 0.70, 0.95), (0.75, 0.55, 0.30), (0.45, 0.35, 0.70),
    (0.95, 0.80, 0.45), (0.20, 0.40, 0.35), (0.85, 0.55, 0.55), (0.40, 0.60, 0.50),
    (0.70, 0.70, 0.90), (0.55, 0.45, 0.20), (0.30, 0.85, 0.60), (0.90, 0.40, 0.45),
    (0.50, 0.50, 0.95), (0.65, 0.85, 0.35),
]
GREY = (0.45, 0.45, 0.48)   # fallback; the live value comes from TH.unknown_colour(theme)

# An orthogroup needs this many visible members before it earns a marker at system level. Below it the
# view fills with thousands of singleton markers, which is the gene level with extra steps.
MIN_ORTHOGROUP_FOR_SYSTEM = 4
# Roughly how many edges can be drawn at full alpha before they stop being separable lines and become a
# filled region. Beyond it, alpha is scaled down rather than edges being dropped -- density stays
# visible as brightness instead of being silently truncated.
EDGE_INK_TARGET = 1500

# Depth of attention is categorical (see literature.DEPTH_OF): named in a title / in an abstract / only in
# a body or caption. Distinct hues rather than a ramp, because the tiers are not a measured quantity.
DEPTH_COLOUR = {"focal": (0.98, 0.86, 0.30),          # the paper is about this gene
                "substantive": (0.35, 0.70, 0.95),    # a stated part of the paper's claims
                "incidental": (0.55, 0.35, 0.60)}     # mentioned in passing, or listed in a table


def load():
    """Load the built cache: node table, coordinates, edge layers and stored models.

    Raises with the command that builds the cache if it is absent, rather than letting pandas raise
    a file-not-found from inside a constructor, which tells the user nothing about what to do.
    """
    npz, pq = os.path.join(DATA, "graph.npz"), os.path.join(DATA, "nodes.parquet")
    if not (os.path.exists(npz) and os.path.exists(pq)):
        raise SystemExit("No cached graph. Run:  python -m starplast.build_graph")
    nodes = pd.read_parquet(pq)
    z = np.load(npz)
    xyz = z["xyz"].astype(np.float32)
    edges = {}
    for k, _ in EDGE_TYPES:
        if f"{k}__a" in z.files:
            edges[k] = {"a": z[f"{k}__a"], "b": z[f"{k}__b"],
                        "w": z[f"{k}__w"], "r": z[f"{k}__r"]}
    # Optional: the crosslink model table, which says how a measured interaction is thought to happen.
    mp = os.path.join(DATA, "crosslink_models.parquet")
    models = pd.read_parquet(mp) if os.path.exists(mp) else pd.DataFrame()
    return nodes, xyz, edges, models


class Map3D(gl.GLViewWidget):
    """GLViewWidget plus click-picking, which pyqtgraph does not provide."""

    picked = QtCore.pyqtSignal(int)

    def __init__(self, xyz):
        super().__init__()
        self.xyz = xyz
        self._proj_kind = None      # resolved once, see _projection_matrix
        # Which points may be picked. None means all of them; a mask is set when the displayed
        # embedding covers only some of the genes, so a hidden gene cannot be selected by clicking
        # where it would have been.
        self.pickable = None
        self.fit_view()

    def fit_view(self, margin=1.35):
        """Frame the data rather than assuming a fixed distance.

        The embedding is renormalised on every build and its extent changes with the feature set, so a
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
        """Radius of the point cloud about its centre, for framing and grid sizing."""
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

    def mouseReleaseEvent(self, ev):
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

    Holds the view state that the menus mutate -- level of detail, colouring, point size, which edge
    layers are active, the filter category -- and redraws from it. State lives here as plain
    attributes rather than being read back out of widgets, so a headless caller can drive the whole
    interface without a window manager, which is how the tests exercise it.
    """
    def __init__(self):
        super().__init__()
        self.theme = 'dark'
        self.point_style = TH.DEFAULT_POINT_STYLE
        self.point_mode = 'occlude'
        self._spin_speed = 0.35
        self.depth_cue = True
        self.show_ground = True
        self.cmap_name = None
        self.apply_log_settings()
        self.nodes, self.xyz, self.edges, self.models = load()
        self.n = len(self.nodes)
        self.sel = None
        self.setWindowTitle("starplast — Toxoplasma knowledge map")
        self.resize(1580, 950)

        self.categories = category_columns(self.nodes)
        self.category = "compartment" if "compartment" in self.categories else self.categories[0]
        comps = sorted(as_text(self.nodes[self.category]).unique())
        absent = {str(x).lower() for x in ABSENCE}
        self.comps = ([c for c in comps if c.lower() not in absent]
                      + [c for c in comps if c.lower() in absent])
        self.colour_of = dict(zip(self.comps,
                                  TH.categorical_colours(len(self.comps), self.theme)))
        for c in self.comps:
            if c.lower() in absent:
                self.colour_of[c] = GREY

        # Level of detail, colouring and edge state now live in the menus, so they are plain
        # attributes with menu actions over them rather than widgets read out of a side panel.
        self.level_idx = 2                    # gene tier: the map as it actually is
        self.colour_mode = "compartment"
        self.point_size = None                # None means "whatever the point style says"
        self.edge_on = {k: (k in ("comention", "cofitness") and k in self.edges)
                        for k, _ in EDGE_TYPES}
        self.attn_on = True                   # a correctness default, not a preference
        self.all_edges_on = False
        self._galaxies = None                 # computed lazily; the grid pass is not free
        # A clustering from the analysis panel, so the map can be coloured by it. Looking at
        # structure beside a held-out variable is what this application is for, and until this
        # existed the Clusters tab computed labels and discarded them.
        self.cluster_labels = None
        # Which genes the displayed embedding has coordinates for. None means "the built cache",
        # which covers all of them; a mask arrives with any map built over a subsample.
        self.placed = None
        self._spin_home = None                # where spin started, so it can be put back

        self.view = Map3D(self.xyz)
        self.view.picked.connect(self.on_pick)
        self.scatter = gl.GLScatterPlotItem(pos=self.xyz, size=5.0, pxMode=True)
        # GLScatterPlotItem blends additively by default, which sums the colours of overlapping points.
        # With 8,140 genes in dense UMAP clusters every mode rendered as one white blob and the colour
        # encoding -- the thing the map is for -- was invisible. Translucent blending with depth testing
        # makes nearer points occlude farther ones instead of adding to them.
        self.scatter.setGLOptions("translucent")
        self.view.addItem(self.scatter)
        self.centroid_item = None
        self.halo_item = None
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
        self.redraw()

    # ------------------------------------------------------------------ logging
    def settings(self):
        """Where preferences persist. One place, so a new setting cannot invent its own file."""
        return QtCore.QSettings("starplast", "starplast")

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
        """The colour map for a column: the user's choice if it suits the data, else the right default.

        A diverging map on a strictly positive quantity invents a midpoint, and a sequential map on a
        residual hides its sign, so the column's kind has the final say over an inappropriate choice.
        """
        kind = TH.kind_for_column(values)
        if self.cmap_name and TH.CMAPS.get(self.cmap_name, (None,))[0] == kind:
            return TH.CMAPS[self.cmap_name][1]
        return TH.DEFAULT_CMAP[kind]

    def apply_theme(self, name: str):
        """Repaint everything from the palette -- widgets, GL background, and the compartment colours."""
        self.theme = name
        app = QtWidgets.QApplication.instance()
        if app is not None:
            app.setStyleSheet(TH.stylesheet(name))
        self.view.setBackgroundColor(pg.mkColor(TH.palette_for(name)["bg"]))
        # Recolour the classes for this ground, then restore the deliberate grey for "unknown".
        self.colour_of = dict(zip(self.comps,
                                  TH.categorical_colours(len(self.comps), name, self.cmap_name)))
        self.colour_of["unassigned"] = TH.unknown_colour(name)[:3]
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
        self.scatter.setData(pos=self.xyz, color=self.colours(self.visible_mask()),
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
            from .tuning import EmbeddingStore
        except Exception as e:
            self.statusBar().showMessage(f"analysis panel unavailable: {e}")
            return
        d = QtWidgets.QDockWidget("analysis")
        d.setFeatures(QtWidgets.QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
        # Shares the window's runner, so an analysis appears in the Jobs panel and can be stopped
        # there like anything else. With its own private thread it was invisible and unstoppable,
        # and it also blocked every other tab for the several minutes a search takes.
        panel = AnalysisPanel(self.nodes, store=EmbeddingStore(os.path.join(DATA, "embeddings")),
                              runner=self.jobs)
        panel.status.connect(lambda m: self.statusBar().showMessage(m))
        panel.embedding_ready.connect(self.use_embedding)
        panel.clusters_ready.connect(self.use_clusters)
        d.setWidget(panel)
        self.addDockWidget(QtCore.Qt.DockWidgetArea.RightDockWidgetArea, d)
        self.tabifyDockWidget(self.right_dock, d)
        self.right_dock.raise_()
        self.analysis_dock = d
        self.panel = panel
        self._gallery()

    def _gallery(self):
        """The walk gallery, along the bottom where a wall of thumbnails has room to be a wall.

        Fed from the analysis panel's per-configuration signal, so a thumbnail appears as each map is
        computed rather than when the walk ends. Clicking one shows it in the central view, which is
        why this is a dock beside the map rather than a window over it: the gallery picks, the map
        displays, and what is displayed is the application's own 3D view with everything it can do.
        """
        from .gallery import GalleryPanel
        bg = TH.rgbf(TH.palette_for(self.theme)["bg"])[:3]
        self.gallery = GalleryPanel(colour_fn=self.colours_for_genes, background=bg)
        self.gallery.chosen.connect(self.show_walk_map)
        self.gallery_dock = QtWidgets.QDockWidget("gallery", self)
        self.gallery_dock.setWidget(self.gallery)
        self.addDockWidget(QtCore.Qt.DockWidgetArea.BottomDockWidgetArea, self.gallery_dock)
        self.gallery_dock.hide()
        self.panel.walk_started.connect(self._walk_started)
        self.panel.walk_step.connect(self.gallery.add)

    def _walk_started(self):
        """Clear the gallery and show it, so the first thumbnail lands somewhere visible."""
        self.gallery.clear()
        self.gallery_dock.show()
        self.gallery_dock.raise_()

    def colours_for_genes(self, mask):
        """The current colouring, restricted to a subset of genes -- what the gallery paints with.

        Taken from `colours` rather than reimplemented so a thumbnail is coloured by exactly what the
        map is coloured by, including the rule that grey means unknown. The visibility filter is
        deliberately not applied: a thumbnail showing only the filtered classes would look like a
        different embedding rather than the same one seen through a filter.
        """
        c = self.colours(np.ones(self.n, bool))
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
        # The galaxy tier is derived from the coordinates, so a new embedding invalidates it. Left
        # cached, the coarse tier would go on describing the map that was replaced.
        self._galaxies = None
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
        map", so genes are clickable, every colour mode applies and edges draw as usual.
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
        f.addAction("Export relationships (CSV)…").triggered.connect(self.export_relationships)
        f.addAction("Export graph (GraphML)…").triggered.connect(self.export_graphml)
        f.addSeparator()
        q = f.addAction("Quit")
        q.setShortcut("Ctrl+Q")
        q.triggered.connect(self.close)

        # ---- View
        v = mb.addMenu("&View")
        lvl = v.addMenu("Level of detail")
        lvl.setToolTipsVisible(True)
        self.level_group = QtGui.QActionGroup(self)
        for i, (name, tip) in enumerate([
                ("Galaxy — spatial structures", "Coarse structures found in the embedding itself."),
                ("System — orthogroup centroids", "One point per orthogroup with at least four members."),
                ("Planet — every gene", "All 8,140 genes. The map as it is.")]):
            act = lvl.addAction(name)
            act.setCheckable(True)
            act.setChecked(i == self.level_idx)
            act.setToolTip(tip)
            act.setData(i)
            self.level_group.addAction(act)
            act.triggered.connect(lambda _c, k=i: self.set_level(k))

        col = v.addMenu("Colour by")
        self.colour_group = QtGui.QActionGroup(self)
        for name in COLOUR_MODES:
            act = col.addAction(name)
            act.setCheckable(True)
            act.setChecked(name == self.colour_mode)
            self.colour_group.addAction(act)
            act.triggered.connect(lambda _c, n=name: self.set_colour_mode(n))

        ps = v.addMenu("Point size")
        self.size_group = QtGui.QActionGroup(self)
        for label, val in POINT_SIZES:
            act = ps.addAction(label)
            act.setCheckable(True)
            act.setChecked(val is None)
            act.setData(val)
            self.size_group.addAction(act)
            act.triggered.connect(lambda _c, s=val: self.set_point_size(s))

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
        v.addAction("Preferences…").triggered.connect(self.open_preferences)

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
        for dock in (self.console_dock, self.jobs_dock, self.chat_dock,
                     getattr(self, "analysis_dock", None), getattr(self, "gallery_dock", None),
                     self.right_dock):
            if dock is not None:
                t.addAction(dock.toggleViewAction())

        h = mb.addMenu("&Help")
        h.addAction("What this map does and does not show").triggered.connect(self.explain_map)
        h.addAction("Precision, recall, and how each can be gamed").triggered.connect(
            self.explain_scoring)

    def _context_menu(self, pos):
        """Right-click on the map. Shows the menu; `build_context_menu` makes it."""
        m = self.build_context_menu()
        m.exec(self.view.mapToGlobal(pos))
        return m

    def build_context_menu(self):
        """Construct the menu without showing it.

        Split from `_context_menu` for the same reason `build_preferences` is split from
        `open_preferences`: exec() enters a modal loop and does not return until a human closes the
        menu, so a test that called it hung forever instead of failing.
        """
        m = QtWidgets.QMenu(self)
        m.addAction(self.spin_act)
        m.addSeparator()
        lvl = m.addMenu("Level of detail")
        for act in self.level_group.actions():
            lvl.addAction(act)
        ps = m.addMenu("Point size")
        for act in self.size_group.actions():
            ps.addAction(act)
        m.addSeparator()
        m.addAction("Export image…").triggered.connect(self.export_image)
        m.addAction("Export visible genes (CSV)…").triggered.connect(self.export_visible)
        m.addAction("Export relationships (CSV)…").triggered.connect(self.export_relationships)
        m.addAction("Export graph (GraphML)…").triggered.connect(self.export_graphml)
        m.addSeparator()
        m.addAction("Reset view / clear filters").triggered.connect(self.reset)
        return m

    # ------------------------------------------------------------------ menu state
    def set_level(self, i: int):
        """Switch level of detail, easing the camera to suit the new tier."""
        self.level_idx = int(i)
        self.on_level_changed()

    def set_colour_mode(self, name: str):
        """Change what colour encodes. Each mode is a different claim about the data."""
        self.colour_mode = name
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
        """Explain, in words, why the twelve relation types are kept separate."""
        QtWidgets.QMessageBox.information(self, "Why edge types are kept separate", EDGE_EXPLANATION)

    def explain_map(self):
        """Explain what the map is and what held-out testing says it does not support."""
        QtWidgets.QMessageBox.information(self, "What this map shows", MAP_EXPLANATION)

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

    def use_clusters(self, labels):
        """Take a clustering from the analysis panel and colour the map by it.

        Switches the colouring automatically, because a user who has just pressed "cluster this map"
        wants to see the clusters -- and leaving it on compartment made the button look inert.
        """
        self.cluster_labels = np.asarray(labels)
        self.set_colour_mode("clusters")
        n = len(set(self.cluster_labels[self.cluster_labels >= 0]))
        self.status.showMessage(f"colouring by {n} clusters; grey is unclustered, which is a real "
                                f"answer and not a missing one")

    def open_preferences(self):
        """Appearance settings, gathered in one place rather than crowding the map panel."""
        self.build_preferences().exec()

    def build_preferences(self):
        """Construct the dialog without showing it.

        Split from `open_preferences` so the controls can be built and inspected without entering a
        modal event loop: exec() blocks until a human closes the window, so a test that called it hung
        forever rather than failing.
        """
        d = QtWidgets.QDialog(self)
        d.setWindowTitle("Preferences")
        form = QtWidgets.QFormLayout(d)

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
            "categorical for classes. A map of the wrong kind is ignored in favour of the right "
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
            "and destroys the colour encoding.")
        self.mode_box.currentTextChanged.connect(self._on_point_mode)

        self.spin_speed = QtWidgets.QDoubleSpinBox()
        self.spin_speed.setRange(0.05, 3.0)
        self.spin_speed.setSingleStep(0.05)
        self.spin_speed.setValue(self._spin_speed)
        self.spin_speed.setSuffix("  °/frame")
        self.spin_speed.valueChanged.connect(lambda v: setattr(self, "_spin_speed", v))

        self.depth_box = QtWidgets.QCheckBox("fade and shrink with distance")
        self.depth_box.setChecked(self.depth_cue)
        self.depth_box.setToolTip(
            "Without it, near and far points are equally bright and the map reads as a flat disc "
            "however far it is rotated. Turn it off to compare two points' colours exactly, since the "
            "fade changes apparent colour with position.")
        self.depth_box.toggled.connect(lambda v: (setattr(self, "depth_cue", v), self.redraw()))

        self.ground_box = QtWidgets.QCheckBox("horizon grid")
        self.ground_box.setChecked(self.show_ground)
        self.ground_box.setToolTip(
            "A reference plane under the cloud. Spinning a bare point cloud, the eye cannot separate "
            "rotation from the points rearranging themselves.")
        self.ground_box.toggled.connect(lambda v: (setattr(self, "show_ground", v), self.redraw()))

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
        form.addRow("colour map", self.cmap_box)
        form.addRow("points", self.point_box)
        form.addRow("rendering", self.mode_box)
        form.addRow("spin speed", self.spin_speed)
        form.addRow("depth", self.depth_box)
        form.addRow("reference", self.ground_box)
        form.addRow("logging", self.log_box)
        form.addRow("keep at level", self.log_file_level)
        form.addRow("show at level", self.log_console_level)
        form.addRow("log file", self.log_path)
        close = QtWidgets.QPushButton("close")
        close.clicked.connect(d.accept)
        form.addRow(close)
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
        p = self.view.opts
        return (float(p.get("azimuth", 0.0)), float(p.get("elevation", 0.0)),
                float(p.get("distance", 40.0)))

    def _restore_camera(self, state):
        if not state:
            return
        az, el, dist = state
        self.view.setCameraPosition(distance=dist, elevation=el, azimuth=az)

    def _refresh_jobs(self, *_):
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
        # The lazily-built level-of-detail cache is the largest thing this application holds that it
        # can rebuild for free.
        if self._galaxies is not None:
            self._galaxies = None
            freed.append("level-of-detail cache")
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
        self.view.orbit(self._spin_speed, 0)

    def _on_cmap(self, name):
        """Apply a colour map. Categorical colours are rebuilt, not just the continuous ramp.

        Only the ramp honoured this before, so choosing a map while colouring by compartment -- the
        default mode -- appeared to do nothing at all.
        """
        self.cmap_name = None if name.startswith("auto") else name
        kind = TH.CMAPS.get(self.cmap_name, (None,))[0] if self.cmap_name else None
        if kind in (None, "categorical"):
            self.colour_of = dict(zip(self.comps, TH.categorical_colours(
                len(self.comps), self.theme, self.cmap_name if kind == "categorical" else None)))
            self.colour_of["unassigned"] = TH.unknown_colour(self.theme)[:3]
        self.redraw()

    def _on_point_style(self, name):
        self.point_style = name
        self.apply_point_style()

    def _on_point_mode(self, name):
        self.point_mode = name
        self.apply_point_style()

    # ------------------------------------------------------------------ panels
    def _left(self):
        """Search and the category filter. Everything else moved into the menus.

        This panel used to carry the level, the colouring, twelve edge checkboxes, spin, preferences
        and the compartment list all at once, which is most of the application's settings stacked in a
        column beside a 3D view. Settings that are chosen once belong in a menu; the two things used
        continuously -- finding a gene and narrowing the field -- stay on screen.
        """
        d = QtWidgets.QDockWidget("find")
        d.setFeatures(QtWidgets.QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
        w = QtWidgets.QWidget()
        L = QtWidgets.QVBoxLayout(w)
        L.setContentsMargins(8, 8, 8, 8)

        self.search = QtWidgets.QLineEdit(placeholderText="gene id or product…")
        self.search.setToolTip(
            "An exact accession first, then a partial one, then the product description. Most of this "
            "proteome has no symbol -- 6,500 genes are 'hypothetical protein' -- so product text is "
            "often the only handle you have on a gene.")
        self.search.returnPressed.connect(self.do_search)
        L.addWidget(self.search)

        L.addWidget(QtWidgets.QLabel("<b>filter by</b>"))
        self.category_box = QtWidgets.QComboBox()
        self.category_box.setToolTip(
            "Which classification to filter and fly by. Any column with a manageable number of "
            "repeated values is offered, not compartment alone -- localisation is the worst-recovered "
            "property in this map, so it is a poor thing to be the only way in.")
        self.category_box.addItems(self.categories)
        if "compartment" in self.categories:
            self.category_box.setCurrentText("compartment")
        self.category_box.currentTextChanged.connect(self.on_category_changed)
        L.addWidget(self.category_box)

        self.comp_list = QtWidgets.QListWidget()
        self.comp_list.setToolTip(
            "Select to show only those classes; select none to show everything. Double-click flies to "
            "a class's centroid. Counts are of genes with a value, so a class is at least this big and "
            "possibly bigger.")
        self.comp_list.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.comp_list.itemSelectionChanged.connect(self.redraw)
        self.comp_list.itemDoubleClicked.connect(self.fly_to_compartment)
        L.addWidget(self.comp_list, 1)
        self._fill_category_list()

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
        """Rebuild the value list for the current category, ordered by size with absences last."""
        self.comp_list.blockSignals(True)
        self.comp_list.clear()
        vals = as_text(self.nodes[self.category])
        counts = vals.value_counts()
        absent = {str(x).lower() for x in ABSENCE}
        named = [v for v in counts.index if str(v).lower() not in absent]
        # Absence values are kept but sunk to the bottom: "unassigned" is usually the largest class in
        # this proteome, and sorting by size alone would put "we do not know" at the top of every list.
        tail = [v for v in counts.index if str(v).lower() in absent]
        for v in named + tail:
            it = QtWidgets.QListWidgetItem(f"{v}  ({int(counts[v]):,})")
            it.setData(QtCore.Qt.ItemDataRole.UserRole, v)
            col = self.colour_of.get(v)
            if col is not None:
                it.setForeground(QtGui.QColor.fromRgbF(*col))
            if str(v).lower() in absent:
                it.setToolTip("Absence, not a class: these genes have no measurement, which is not "
                              "the same as measuring zero.")
            self.comp_list.addItem(it)
        self.comp_list.blockSignals(False)

    def on_category_changed(self, name: str):
        """Switch the filter column, rebuild its palette and its list."""
        self.category = name
        vals = sorted(as_text(self.nodes[name]).unique())
        absent = {str(x).lower() for x in ABSENCE}
        self.comps = ([v for v in vals if str(v).lower() not in absent]
                      + [v for v in vals if str(v).lower() in absent])
        self.colour_of = dict(zip(self.comps,
                                  TH.categorical_colours(len(self.comps), self.theme, self.cmap_name)))
        for v in self.comps:
            if str(v).lower() in absent:
                self.colour_of[v] = TH.unknown_colour(self.theme)[:3]
        self._fill_category_list()
        self.redraw()

    def _right(self):
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
    def visible_mask(self):
        """Boolean mask of the genes passing the current class filter AND having a position.

        The two are combined here rather than at each drawing site because everything downstream --
        points, centroids, edges, the gene count in the status bar -- reads this one mask. A gene the
        displayed embedding does not cover has no position to draw, and drawing it anyway would put
        absence on the map as though it were a measurement.
        """
        sel = [i.data(QtCore.Qt.ItemDataRole.UserRole) for i in self.comp_list.selectedItems()]
        vis = (np.ones(self.n, bool) if not sel
               else as_text(self.nodes[self.category]).isin(sel).to_numpy())
        placed = getattr(self, "placed", None)
        return vis if placed is None else (vis & placed)

    def colours(self, vis):
        """An RGBA colour per gene under the current colour mode. Grey always means unknown."""
        mode = self.colour_mode
        c = np.zeros((self.n, 4), dtype=np.float32)
        if mode.startswith("compartment"):
            # The default colours the MEASURED hyperLOPIT call only. The second mode fills in
            # ortholog-transferred labels, which are inferences from another species -- offered because
            # coverage matters, kept separate because provenance matters more.
            col_name = ("compartment_best" if "transferred" in mode
                        and "compartment_best" in self.nodes.columns else "compartment")
            vals = as_text(self.nodes[col_name])
            for comp, col in self.colour_of.items():
                c[(vals == comp).to_numpy(), :3] = col
            c[(vals == "unassigned").to_numpy(), :3] = TH.unknown_colour(self.theme)[:3]
        elif mode == "clusters":
            # Noise stays grey, with everything else that is unknown. HDBSCAN calling a gene
            # unclustered is a finding about that gene, not a gap in the drawing.
            if self.cluster_labels is None or len(self.cluster_labels) != self.n:
                c[:, :3] = TH.unknown_colour(self.theme)[:3]
            else:
                lab = self.cluster_labels
                ids = sorted(set(lab[lab >= 0]))
                palette = TH.categorical_colours(max(len(ids), 1), self.theme, self.cmap_name)
                for col, k in zip(palette, ids):
                    c[lab == k, :3] = col
                c[lab < 0, :3] = TH.unknown_colour(self.theme)[:3]
        elif mode == "depth of attention":
            # Categorical, not a scale: these tiers are read off document structure (title / abstract /
            # body-only), so shading them along a gradient would imply a quantity that does not exist.
            # Never-named genes stay grey with everything else that is unknown rather than absent.
            if "attention_depth" not in self.nodes.columns:
                c[:, :3] = TH.unknown_colour(self.theme)[:3]
            else:
                d = as_text(self.nodes.attention_depth).to_numpy()
                for tier, col in DEPTH_COLOUR.items():
                    c[d == tier, :3] = col
                c[d == "", :3] = TH.unknown_colour(self.theme)[:3]
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
            c[~ok, :3] = TH.unknown_colour(self.theme)[:3]   # missingness stays grey, never a value
        alpha = float(TH.POINT_STYLES.get(self.point_style, {}).get('alpha', 0.95))
        c[:, 3] = np.where(vis, alpha, 0.06)
        if self.sel is not None:
            c[self.sel] = (1.0, 1.0, 1.0, 1.0)
        if getattr(self, "placed", None) is not None:
            # Fully transparent, not dimmed. A gene the displayed embedding has no coordinates for is
            # absent from it, and a faint point is still a point -- at 0.06 alpha, several thousand of
            # them overlapping read as a real feature of the map.
            c[~self.placed, 3] = 0.0
        return c

    def galaxy_labels(self):
        """Coarse spatial components of the embedding, computed once and cached.

        Cached because the grid pass walks every occupied cell and the tier is switched to more often
        than the embedding changes. `use_embedding` clears it.
        """
        if self._galaxies is None:
            self._galaxies = lod.galaxies(self.xyz)
        return self._galaxies

    def redraw(self):
        """Redraw everything from the current state: points, tiers, ground, halo and edges."""
        vis = self.visible_mask()
        lvl = self.level_idx

        for it in self.edge_items:
            self.view.removeItem(it)
        self.edge_items = []
        if self.centroid_item is not None:
            self.view.removeItem(self.centroid_item)
            self.centroid_item = None

        st = TH.POINT_STYLES.get(self.point_style, TH.POINT_STYLES[TH.DEFAULT_POINT_STYLE])
        # An explicit point size overrides the style's. The style still supplies the rest of its
        # look, so "large" plus an explicit 4 px is a large-style point drawn at 4 px.
        base = float(self.point_size if self.point_size is not None else st["size"])
        # redraw() sets sizes every time, so the chosen point style has to be applied here rather than
        # only in apply_point_style -- otherwise picking a style changed nothing the moment anything
        # else triggered a redraw.
        sizes = np.where(vis, base, max(base * 0.4, 1.5)).astype(np.float32)
        if self.sel is not None:
            sizes[self.sel] = max(base * 3.0, 12.0)
        if lvl == 0:
            # Galaxy level: genes recede and coarse SPATIAL structures carry the map.
            #
            # This used to draw one centroid per compartment, and it produced a cluster of dots in the
            # middle of the screen -- correctly, which was the problem. A compartment's genes are
            # scattered across the whole embedding, so their mean lands near the centre of it: the
            # median compartment centroid sits 0.17 of the map radius from the middle while the
            # median compartment's own members spread 0.37 of the radius. The centroids were an
            # honest average of something with no spatial meaning, drawn as though it had one, and
            # the held-out search says the same thing from the other side -- compartment is the
            # worst-recovered target in this map, below the study-effort control.
            #
            # So the tier is computed from the embedding itself. See lod.py.
            sizes = np.full(self.n, max(base * 0.4, 1.5), np.float32)
            lab = self.galaxy_labels()
            pos, num, spread = lod.centroids(self.xyz, np.where(vis, lab, -1))
            dom = lod.dominant(as_text(self.nodes[self.category]).to_numpy(),
                               np.where(vis, lab, -1), ignore=ABSENCE)
            self._galaxy_info = []
            col, ssz = [], []
            for i in range(len(pos)):
                name, frac = dom.get(i, ("mixed", 0.0))
                col.append((*self.colour_of.get(name, TH.unknown_colour(self.theme)[:3]), 0.95))
                ssz.append(float(10 + 30 * np.sqrt(num[i] / max(vis.sum(), 1))))
                self._galaxy_info.append((int(num[i]), name, frac))
            if len(pos):
                self.centroid_item = gl.GLScatterPlotItem(
                    pos=pos.astype(np.float32), color=np.array(col, np.float32),
                    size=np.array(ssz, np.float32), pxMode=True)
                self.centroid_item.setGLOptions("translucent")
                self.view.addItem(self.centroid_item)
        elif lvl == 1:
            # System level: orthogroups get their own centroids, the way compartments do at galaxy
            # level. Previously this tier did nothing at all unless a gene happened to be selected, so
            # it was indistinguishable from the gene level -- the middle of a three-tier hierarchy
            # silently missing.
            og = as_text(self.nodes.orthogroup).to_numpy()
            sizes = np.full(self.n, max(base * 0.5, 2.0), np.float32)
            pos, col, ssz = [], [], []
            groups = {}
            for i, g in enumerate(og):
                if g not in ("", "nan", "None") and vis[i]:
                    groups.setdefault(g, []).append(i)
            for g, idx in groups.items():
                if len(idx) < MIN_ORTHOGROUP_FOR_SYSTEM:
                    continue
                pos.append(self.xyz[idx].mean(0))
                comp = as_text(self.nodes[self.category]).iloc[idx[0]]
                col.append((*self.colour_of.get(comp, TH.unknown_colour(self.theme)[:3]), 0.9))
                ssz.append(float(6 + 18 * np.sqrt(len(idx) / 40.0)))
            if pos:
                self.centroid_item = gl.GLScatterPlotItem(
                    pos=np.array(pos, np.float32), color=np.array(col, np.float32),
                    size=np.array(ssz, np.float32), pxMode=True)
                self.centroid_item.setGLOptions("translucent")
                self.view.addItem(self.centroid_item)
            if self.sel is not None and og[self.sel] not in ("", "nan", "None"):
                sizes[og == og[self.sel]] = max(base * 2.0, 9.0)

        colours = self.colours(vis)
        self._depth_ctx = None
        if self.depth_cue:
            # Fade and shrink with distance so the cloud has a front and a back. Without it the map is
            # a flat disc of colour however much it is rotated.
            cam = self.view.cameraPosition()
            cam = (cam.x(), cam.y(), cam.z())
            d = np.linalg.norm(self.xyz - np.asarray(cam, np.float32), axis=1)
            # One range shared by points and edges. Normalising each set against its own extent makes
            # an edge and the point it touches fade by different amounts, which reads as flicker.
            self._depth_ctx = (cam, (float(d.min()), float(d.max())))
            a, sizes = TH.depth_cue(self.xyz, cam, colours[:, 3], sizes, rng=self._depth_ctx[1])
            colours = colours.copy()
            colours[:, 3] = a
            sizes = sizes.astype(np.float32)

        self.scatter.setData(pos=self.xyz, color=colours, size=sizes)
        self._draw_ground()
        self._draw_selection_halo()
        self.draw_edges(vis)

        act = [k for k, _ in EDGE_TYPES if self.edge_on.get(k)]
        note = "" if self.attn_on else "  ·  RAW co-mention (attention-biased)"
        extra = ""
        if lvl == 0 and getattr(self, "_galaxy_info", None):
            # Say what the blobs are. A tier that draws five unexplained spheres is not an
            # abstraction, it is a mystery.
            extra = ("  ·  " + ", ".join(f"{n:,} genes ({name} {frac*100:.0f}%)"
                                         for n, name, frac in self._galaxy_info[:4]))
        if not self.edge_on.get("comention", False) and not act:
            note += "  ·  select a gene to see its edges, or Edges ▸ Draw all active edges"
        self.status.showMessage(
            f"{int(vis.sum()):,} / {self.n:,} genes shown  ·  edges: "
            f"{', '.join(act) if act else 'none'}{note}{extra}")

    def on_level_changed(self):
        """Move the camera to suit the tier, easing rather than cutting.

        Pulling back for the galaxy tier and moving in for the gene tier makes the hierarchy legible as
        a change of scale, which is what the three tiers actually are.
        """
        r = self.view.data_radius()
        self.view.animate_distance(max(r * (2.0, 1.6, 1.25)[self.level_idx], 10.0))
        self.redraw()

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

        A size bump alone is invisible in a dense region -- the neighbours are the same colour and the
        selected point simply becomes a slightly bigger dot in a crowd. A translucent halo in the
        accent colour separates it from its neighbourhood at any density.
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
                   # deliberately the loudest colours in the palette: these are the things to look at
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
        """Move the camera to a class's centroid and drop to the gene tier."""
        c = item.data(QtCore.Qt.ItemDataRole.UserRole)
        m = (as_text(self.nodes[self.category]) == c).to_numpy()
        if m.sum():
            self.view.setCameraPosition(pos=pg.Vector(*self.xyz[m].mean(0)), distance=60)
            self.set_level(2)

    def reset(self):
        """Clear the selection and every filter, and frame the whole map again."""
        self.sel = None
        self.comp_list.clearSelection()
        self.view.fit_view()
        self.detail.setHtml("<p style='color:#888'>Click a gene.</p>")
        self.redraw()

    # ------------------------------------------------------------------ export
    def _ask_path(self, caption: str, filt: str, default: str) -> str:
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
        """What to tick when the picker opens: identity, the active filter, the active colouring.

        Taken from the current view rather than fixed, because the columns worth exporting are
        usually the ones being looked at.
        """
        want = ["gene_id", "product", self.category]
        want += {"in vitro fitness": ["fit_invitro_hff"],
                 "publications": ["n_publications"],
                 "structure confidence (pLDDT)": ["mean_plddt"],
                 "depth of attention": ["attention_depth"],
                 "cyst / tachyzoite expression": ["expr_cyst", "expr_tachy"],
                 }.get(self.colour_mode, ["compartment"])
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
        bits = [f"Level: {['galaxy','system','planet'][self.level_idx]}.",
                f"Colouring: {self.colour_mode}.",
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
            return "—" if v is None or (isinstance(v, float) and not np.isfinite(v)) else f.format(v)

        rows = [("compartment (hyperLOPIT)",
                 f"{r.compartment}" + (" <i>— unknown, not absent</i>"
                                       if r.compartment == "unassigned" else "")),
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
                ("curated host targets", num(r.get("n_host_targets"), "{:.0f}")),
                ("log2 FPKM tachyzoite", num(r.get("expr_tachy"))),
                ("log2 FPKM tissue cyst", num(r.get("expr_cyst")))]
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
                             if np.isfinite(row.frac_satisfied) else "not modelled")
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
                xl = ("<p><b>How the binding is modelled</b> <span style='color:#888;"
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
        <p><a href="https://toxodb.org/toxo/app/record/gene/{gid}">ToxoDB record</a> ·
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
