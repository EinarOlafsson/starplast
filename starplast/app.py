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
import pyqtgraph as pg  # noqa: E402
import pyqtgraph.opengl as gl  # noqa: E402

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")

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

# Qualitative palette; "unassigned" is deliberately grey, because a missing hyperLOPIT call means
# unknown (assignment tracks abundance) and must not read as a 27th compartment.
PALETTE = [
    (0.90, 0.24, 0.24), (0.20, 0.55, 0.90), (0.25, 0.75, 0.35), (0.95, 0.65, 0.15),
    (0.65, 0.35, 0.85), (0.15, 0.80, 0.78), (0.95, 0.45, 0.70), (0.55, 0.75, 0.20),
    (0.85, 0.35, 0.10), (0.35, 0.45, 0.85), (0.10, 0.65, 0.50), (0.80, 0.80, 0.20),
    (0.60, 0.20, 0.45), (0.30, 0.70, 0.95), (0.75, 0.55, 0.30), (0.45, 0.35, 0.70),
    (0.95, 0.80, 0.45), (0.20, 0.40, 0.35), (0.85, 0.55, 0.55), (0.40, 0.60, 0.50),
    (0.70, 0.70, 0.90), (0.55, 0.45, 0.20), (0.30, 0.85, 0.60), (0.90, 0.40, 0.45),
    (0.50, 0.50, 0.95), (0.65, 0.85, 0.35),
]
GREY = (0.45, 0.45, 0.48)

# Depth of attention is categorical (see literature.DEPTH_OF): named in a title / in an abstract / only in
# a body or caption. Distinct hues rather than a ramp, because the tiers are not a measured quantity.
DEPTH_COLOUR = {"focal": (0.98, 0.86, 0.30),          # the paper is about this gene
                "substantive": (0.35, 0.70, 0.95),    # a stated part of the paper's claims
                "incidental": (0.55, 0.35, 0.60)}     # mentioned in passing, or listed in a table


def load():
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
        self.setCameraPosition(distance=170)

    def _mvp(self):
        m = self.projectionMatrix() * self.viewMatrix()
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
        sx, sy = self.project()
        d = np.hypot(sx - p.x(), sy - p.y())
        if np.all(np.isnan(d)):
            return
        i = int(np.nanargmin(d))
        if d[i] < 14:
            self.picked.emit(i)


class Window(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.nodes, self.xyz, self.edges, self.models = load()
        self.n = len(self.nodes)
        self.sel = None
        self.setWindowTitle("starplast — Toxoplasma knowledge map")
        self.resize(1580, 950)

        comps = sorted(self.nodes.compartment.astype(str).unique())
        self.comps = [c for c in comps if c != "unassigned"] + \
                     (["unassigned"] if "unassigned" in comps else [])
        self.colour_of = {c: PALETTE[i % len(PALETTE)] for i, c in enumerate(self.comps)}
        self.colour_of["unassigned"] = GREY

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
        self.edge_items = []

        self.setCentralWidget(self.view)
        self.addDockWidget(QtCore.Qt.DockWidgetArea.LeftDockWidgetArea, self._left())
        self.addDockWidget(QtCore.Qt.DockWidgetArea.RightDockWidgetArea, self._right())
        self.status = self.statusBar()
        self.redraw()

    # ------------------------------------------------------------------ panels
    def _left(self):
        d = QtWidgets.QDockWidget("map")
        d.setFeatures(QtWidgets.QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
        w = QtWidgets.QWidget()
        L = QtWidgets.QVBoxLayout(w)
        L.setContentsMargins(8, 8, 8, 8)

        self.search = QtWidgets.QLineEdit(placeholderText="gene id or product…")
        self.search.returnPressed.connect(self.do_search)
        L.addWidget(self.search)

        L.addWidget(QtWidgets.QLabel("<b>level of detail</b>"))
        self.level = QtWidgets.QComboBox()
        self.level.addItems(["compartment (galaxy)", "orthogroup / module (system)",
                             "gene (planet)"])
        self.level.setCurrentIndex(2)
        self.level.currentIndexChanged.connect(self.redraw)
        L.addWidget(self.level)

        L.addWidget(QtWidgets.QLabel("<b>colour by</b>"))
        self.colour_by = QtWidgets.QComboBox()
        self.colour_by.addItems(["compartment", "compartment (incl. transferred)",
                                 "in vitro fitness", "publications",
                                 "depth of attention",
                                 "structure confidence (pLDDT)", "cyst / tachyzoite expression"])
        self.colour_by.currentIndexChanged.connect(self.redraw)
        L.addWidget(self.colour_by)

        L.addWidget(QtWidgets.QLabel("<b>edge types</b> (never merged)"))
        self.edge_cb = {}
        for k, label in EDGE_TYPES:
            cb = QtWidgets.QCheckBox(label)
            cb.setEnabled(k in self.edges)
            cb.setChecked(k in ("comention", "cofitness") and k in self.edges)
            cb.stateChanged.connect(self.redraw)
            n = len(self.edges[k]["a"]) if k in self.edges else 0
            cb.setToolTip(f"{n:,} edges")
            self.edge_cb[k] = cb
            L.addWidget(cb)

        self.attn = QtWidgets.QCheckBox("attention-corrected co-mention")
        self.attn.setChecked(True)     # a correctness default, not a preference
        self.attn.setToolTip("Raw co-mention counts track how often a gene is studied, not how "
                             "related two genes are. Corrected mode shows log2 observed/expected "
                             "given each gene's own publication count.")
        self.attn.stateChanged.connect(self.redraw)
        L.addWidget(self.attn)

        self.all_edges = QtWidgets.QCheckBox("draw all active edges (capped)")
        self.all_edges.stateChanged.connect(self.redraw)
        L.addWidget(self.all_edges)

        L.addWidget(QtWidgets.QLabel("<b>compartments</b>"))
        self.comp_list = QtWidgets.QListWidget()
        self.comp_list.setSelectionMode(
            QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        for c in self.comps:
            it = QtWidgets.QListWidgetItem(f"{c}  ({int((self.nodes.compartment == c).sum())})")
            it.setData(QtCore.Qt.ItemDataRole.UserRole, c)
            r, g, b = self.colour_of[c]
            it.setForeground(QtGui.QColor.fromRgbF(r, g, b))
            self.comp_list.addItem(it)
        self.comp_list.itemSelectionChanged.connect(self.redraw)
        self.comp_list.itemDoubleClicked.connect(self.fly_to_compartment)
        L.addWidget(self.comp_list, 1)

        b = QtWidgets.QPushButton("reset view / clear filters")
        b.clicked.connect(self.reset)
        L.addWidget(b)
        d.setWidget(w)
        w.setMinimumWidth(310)
        return d

    def _right(self):
        d = QtWidgets.QDockWidget("evidence")
        d.setFeatures(QtWidgets.QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
        self.detail = QtWidgets.QTextBrowser()
        self.detail.setOpenExternalLinks(True)
        self.detail.setHtml("<p style='color:#888'>Click a gene.</p>")
        self.detail.setMinimumWidth(400)
        d.setWidget(self.detail)
        return d

    # ------------------------------------------------------------------ drawing
    def visible_mask(self):
        sel = [i.data(QtCore.Qt.ItemDataRole.UserRole) for i in self.comp_list.selectedItems()]
        if not sel:
            return np.ones(self.n, bool)
        return self.nodes.compartment.astype(str).isin(sel).to_numpy()

    def colours(self, vis):
        mode = self.colour_by.currentText()
        c = np.zeros((self.n, 4), dtype=np.float32)
        if mode.startswith("compartment"):
            # The default colours the MEASURED hyperLOPIT call only. The second mode fills in
            # ortholog-transferred labels, which are inferences from another species -- offered because
            # coverage matters, kept separate because provenance matters more.
            col_name = ("compartment_best" if "transferred" in mode
                        and "compartment_best" in self.nodes.columns else "compartment")
            vals = self.nodes[col_name].astype(str)
            for comp, col in self.colour_of.items():
                c[(vals == comp).to_numpy(), :3] = col
            c[(vals == "unassigned").to_numpy(), :3] = GREY
        elif mode == "depth of attention":
            # Categorical, not a scale: these tiers are read off document structure (title / abstract /
            # body-only), so shading them along a gradient would imply a quantity that does not exist.
            # Never-named genes stay grey with everything else that is unknown rather than absent.
            if "attention_depth" not in self.nodes.columns:
                c[:, :3] = GREY
            else:
                d = self.nodes.attention_depth.astype(str).to_numpy()
                for tier, col in DEPTH_COLOUR.items():
                    c[d == tier, :3] = col
                c[d == "", :3] = GREY
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
                cm = pg.colormap.get("viridis")
                c[:, :3] = cm.map(np.nan_to_num(t, nan=0.0), mode="float")[:, :3]
            c[~ok, :3] = GREY          # informative missingness stays grey, never mapped to a value
        c[:, 3] = np.where(vis, 0.95, 0.06)
        if self.sel is not None:
            c[self.sel] = (1.0, 1.0, 1.0, 1.0)
        return c

    def redraw(self):
        vis = self.visible_mask()
        lvl = self.level.currentIndex()

        for it in self.edge_items:
            self.view.removeItem(it)
        self.edge_items = []
        if self.centroid_item is not None:
            self.view.removeItem(self.centroid_item)
            self.centroid_item = None

        sizes = np.where(vis, 5.0, 2.0).astype(np.float32)
        if self.sel is not None:
            sizes[self.sel] = 15.0
        if lvl == 0:
            # galaxy level: genes recede, compartment centroids carry the map
            sizes = np.full(self.n, 2.0, np.float32)
            pos, col, ssz = [], [], []
            for c in self.comps:
                m = (self.nodes.compartment.astype(str) == c).to_numpy() & vis
                if m.sum() < 3:
                    continue
                pos.append(self.xyz[m].mean(0))
                col.append((*self.colour_of[c], 0.95))
                ssz.append(float(8 + 26 * np.sqrt(m.sum() / max(vis.sum(), 1))))
            if pos:
                self.centroid_item = gl.GLScatterPlotItem(
                    pos=np.array(pos, np.float32), color=np.array(col, np.float32),
                    size=np.array(ssz, np.float32), pxMode=True)
                self.view.addItem(self.centroid_item)
        elif lvl == 1 and self.sel is not None:
            og = self.nodes.orthogroup.astype(str).to_numpy()
            same = og == og[self.sel]
            if og[self.sel] not in ("", "nan", "None"):
                sizes[same] = 10.0

        self.scatter.setData(pos=self.xyz, color=self.colours(vis), size=sizes)
        self.draw_edges(vis)

        act = [k for k, _ in EDGE_TYPES if self.edge_cb[k].isChecked()]
        note = "" if self.attn.isChecked() else "  ·  RAW co-mention (attention-biased)"
        self.status.showMessage(
            f"{int(vis.sum()):,} / {self.n:,} genes shown  ·  edges: "
            f"{', '.join(act) if act else 'none'}{note}")

    def draw_edges(self, vis):
        active = [k for k, _ in EDGE_TYPES if self.edge_cb[k].isChecked() and k in self.edges]
        if not active:
            return
        for k in active:
            e = self.edges[k]
            a, b = e["a"], e["b"]
            w = e["r"] if (k in COMENTION and self.attn.isChecked()) else e["w"]
            keep = vis[a] & vis[b]
            if k in COMENTION and self.attn.isChecked():
                keep &= w > 0          # corrected mode shows only more-than-expected pairs
            if self.sel is not None and not self.all_edges.isChecked():
                keep &= (a == self.sel) | (b == self.sel)
            elif not self.all_edges.isChecked():
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
            alpha = (col[3] * (0.12 + 0.88 * t)).astype(np.float32)
            cols[0::2, 3] = alpha
            cols[1::2, 3] = alpha
            it = gl.GLLinePlotItem(pos=seg, color=cols, width=1.0, mode="lines", antialias=True)
            self.view.addItem(it)
            self.edge_items.append(it)

    # ------------------------------------------------------------------ interaction
    def on_pick(self, i):
        self.sel = int(i)
        self.show_detail(self.sel)
        self.redraw()

    def do_search(self):
        q = self.search.text().strip().lower()
        if not q:
            return
        gid = self.nodes.gene_id.astype(str).str.lower()
        hit = np.where(gid == q)[0]
        if hit.size == 0:
            hit = np.where(gid.str.contains(q, regex=False))[0]
        if hit.size == 0:
            hit = np.where(self.nodes["product"].astype(str).str.lower()
                           .str.contains(q, regex=False))[0]
        if hit.size == 0:
            self.status.showMessage(f"no match for {q!r}")
            return
        self.on_pick(int(hit[0]))
        p = self.xyz[self.sel]
        self.view.setCameraPosition(pos=pg.Vector(*p), distance=45)
        self.status.showMessage(f"{hit.size} match(es); showing {self.nodes.gene_id.iloc[hit[0]]}")

    def fly_to_compartment(self, item):
        c = item.data(QtCore.Qt.ItemDataRole.UserRole)
        m = (self.nodes.compartment.astype(str) == c).to_numpy()
        if m.sum():
            self.view.setCameraPosition(pos=pg.Vector(*self.xyz[m].mean(0)), distance=60)
            self.level.setCurrentIndex(2)

    def reset(self):
        self.sel = None
        self.comp_list.clearSelection()
        self.view.setCameraPosition(pos=pg.Vector(0, 0, 0), distance=170)
        self.detail.setHtml("<p style='color:#888'>Click a gene.</p>")
        self.redraw()

    # ------------------------------------------------------------------ detail
    def show_detail(self, i):
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
            w = e["r"] if (k in COMENTION and self.attn.isChecked()) else e["w"]
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
    pg.setConfigOptions(antialias=True)
    app = QtWidgets.QApplication(sys.argv)
    app.setApplicationName("starplast")
    w = Window()
    w.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
