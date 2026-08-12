#!/usr/bin/env python3
"""The analysis tab: choose the data, tune the map, cluster it, and ask what the clusters mean.

This is the interface to `embedding`, `clustering`, `tuning` and `search`. It is deliberately arranged in
the order the work is actually done, top to bottom, because the order matters scientifically:

    1  data      which blocks feed the map, and how missing values and scales are handled
    2  map       UMAP hyperparameters, with a seeded walk to pick them
    3  clusters  DBSCAN / HDBSCAN, with a walk to pick those too
    4  inference what the clusters correspond to, held-out features only
    5  search    walk dataset combinations looking for structure that recovers a held-out label
    6  validation hide labels you already have and see whether the clustering puts them back

Long jobs run on a worker thread. The signal is relayed through a bound method rather than connected
directly, because a directly-connected `finished` handler runs on the worker thread and touching widgets
from there is a crash waiting for a slow machine.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from PyQt6 import QtCore, QtWidgets

from .embedding import BLOCKS, EmbeddingSpec, NA_POLICIES, SCALINGS, columns_for, variance_share
from .theme import CMAPS, POINT_MODES, POINT_STYLES, THEMES, cmaps_of, kind_for_column


class Worker(QtCore.QObject):
    """Runs one callable off the GUI thread and reports back."""
    done = QtCore.pyqtSignal(object, object)     # result, error
    progress = QtCore.pyqtSignal(str)

    def __init__(self, fn):
        super().__init__()
        self.fn = fn

    @QtCore.pyqtSlot()
    def run(self):
        """Run the callable and emit its result, or its exception, never raising into Qt."""
        try:
            self.done.emit(self.fn(self.progress.emit), None)
        except Exception as e:                                   # pragma: no cover - GUI path
            self.done.emit(None, e)



def wrap_tip(text: str, width: int = 64) -> str:
    """A tooltip as one block of even lines, wrapped once and left that way.

    Qt lays a plain tooltip out on a single line, so any real explanation becomes a strip wider than
    the screen. Inserting line breaks alone does not fix it: Qt re-wraps rich text at a width of its
    own choosing, so the manual breaks land inside Qt's lines and strand two words on a row -- which
    is exactly what the first version did.

    `white-space: pre` is what stops the second wrap. Measured on a typical tooltip it takes the
    laid-out width from 1188 pixels to 347 and keeps the breaks where they were put. A fixed table
    width and a styled div were both tried first and neither constrains a tooltip, because Qt only
    wraps when it is given an explicit text width and a tooltip sets its own.

    Lines are padded to equal length so the block is a rectangle rather than a ragged edge. True
    justification -- flush on both sides with stretched spaces -- is not available in Qt's rich-text
    subset, so this is the closest honest thing to it.
    """
    import textwrap
    from html import escape
    blocks = []
    for para in [" ".join(p.split()) for p in text.split("\n\n") if p.strip()]:
        lines = textwrap.wrap(para, width) or [""]
        longest = max(len(x) for x in lines)
        # Padded with non-breaking spaces, which Qt keeps; ordinary trailing spaces are dropped.
        blocks.append("\n".join(escape(x) + "&#160;" * (longest - len(x)) for x in lines))
    return '<div style="white-space:pre">' + "\n\n".join(blocks) + "</div>"


def range_note(widget) -> str:
    """Why a spin box stops where it does, read off the widget itself.

    A minimum and a maximum are decisions, and an unexplained one reads as arbitrary. Generated
    rather than written per control so it cannot fall out of step with the range actually set.
    """
    lo, hi = widget.minimum(), widget.maximum()
    fmt = (lambda v: f"{v:g}")
    return f"Range {fmt(lo)} to {fmt(hi)}."

#: Why each bounded control stops where it does. A minimum and a maximum are decisions, and an
#: unexplained bound reads as arbitrary -- or worse, as a limit of the method rather than a choice.
LIMITS = {
    "nn": "Below 2 there is no neighbourhood to embed from. The upper bound is well past "
          "useful: at 400 neighbours on 8,140 genes the map is nearly global structure only.",
    "md": "0 packs points as tightly as the layout allows; 1 spreads them as far as it can. "
          "Both extremes are legal and both are misleading, for opposite reasons.",
    "seed": "Any integer. It changes the answer, which is why it is recorded with every "
            "stored embedding.",
    "sample": "Below about 200 genes a UMAP is noise. The ceiling is the whole proteome, and "
              "the whole proteome is what a published number should be computed on.",
    "max_missing": "0 keeps only columns with no missing value at all; 1 keeps every column "
                   "however empty. Most of this proteome is unmeasured, so the top of this "
                   "range fills the map with absence indicators.",
    "mcs": "Below 3 a 'cluster' is a pair. The floor is what stops the degenerate answer, "
           "since any purity score is won outright by singletons.",
    "eps": "In embedding units, so the useful value depends entirely on min_dist and the "
           "scale of the coordinates. Ignored by HDBSCAN.",
    "search_sample": "Same range as the walk, same caveat: small subsamples score better "
                     "than the full proteome will reproduce.",
    "max_blocks": "One block per configuration up to four. Combinations grow factorially, so "
                  "four is already thousands of runs.",
    "val_folds": "Two folds is the least that gives a spread; twenty is well past the point "
                 "where the estimate stops moving.",
    "val_hold": "Hide too little and the estimate is noisy; hide more than half and the "
                "visible members no longer pick the cluster a user would pick.",
    "score_min_cluster": "Below 1 there is no cluster. The ceiling is arbitrary; what matters is "
                         "that this is a floor, and raising it is how you stop tiny clusters "
                         "winning a purity score.",
    "score_min_label": "Two is the least that can be a class. Raising it drops rare categories "
                       "from the score entirely, which is a trade rather than an improvement.",
    "min_recall": "0 accepts any cluster however small, which is the failure this floor "
                  "exists to prevent. 1 demands a cluster holding the entire label.",
}

#: Why each control exists, not what it is called. Applied after the tabs are built so every widget
#: is covered in one auditable place rather than scattered through five builders -- and so a new
#: control without an explanation is a visible omission rather than a silent one.
TOOLTIPS = {
    # 1 Data
    "cat_cb": "One-hot the measured hyperLOPIT compartment INTO the map. Leave it OFF when you "
              "intend to hold localisation out and test whether the map recovers it: a map built on "
              "a label separates that label by construction, and the result means nothing.",
    "nn_grid": "The n_neighbors values the WALK sweeps, as opposed to the single value above, "
               "which only affects 'build this map'. Comma-separated (5, 15, 25) or a range "
               "(5:100:5). Every extra value multiplies the length of the walk.",
    "md_grid": "The min_dist values the WALK sweeps. Comma-separated or min:max:step. Low values "
               "make dense clumps that look like clusters whether or not they are, so a grid that "
               "only samples low min_dist will flatter every clustering that follows.",
    "nn_grid2": "The n_neighbors values the search sweeps. Mirrors the Map tab.",
    "md_grid2": "The min_dist values the search sweeps. Mirrors the Map tab.",
    "mcs_grid": "The min_cluster_size values the SEARCH sweeps when clustering each embedding. "
                "This is the guard against the degenerate answer, so sweeping it low is sweeping "
                "toward solutions that win on purity by shattering into singletons.",
    "max_missing": "Drop a column missing in more than this fraction of genes. Most of this "
                   "proteome is unmeasured, so a permissive setting fills the map with columns that "
                   "are mostly absence indicators rather than measurements.",
    # 2 Map
    "nn": "UMAP n_neighbors: how much of the neighbourhood each point is placed by. Small values "
          "preserve local detail and fragment the map; large values preserve global shape and merge "
          "genuinely distinct groups. This is the single most consequential hyperparameter here.",
    "md": "UMAP min_dist: how tightly points may pack. Low values make dense, visually separated "
          "clumps that look like clusters whether or not they are; higher values spread points and "
          "make the density gradient honest. Low min_dist flatters every clustering that follows.",
    "seed": "Random seed. Recorded with every stored embedding, because a structure nobody can "
            "rebuild is not a result -- and UMAP moves noticeably between seeds at these sizes.",
    "sample": "How many genes the hyperparameter walk uses per configuration. Smaller is faster and "
              "FLATTERING: measured on this data, a 3,000-gene subsample scored compartment at "
              "0.228 where the full proteome scored 0.193. Confirm any winner at full size.",
    "emb_name": "Name this embedding so it can be reloaded and compared. Stored with its full "
                "recipe -- blocks, scaling, seed and hyperparameters -- so a map can be rebuilt "
                "exactly rather than approximately.",
    # 3 Clusters
    "algo": "HDBSCAN finds clusters of varying density and labels the rest noise, which suits a "
            "proteome where most genes belong to no tight group. DBSCAN needs one density for "
            "everything, so it either splits the sparse regions or merges the dense ones.",
    "mcs": "The smallest group that counts as a cluster. This is the guard against the degenerate "
           "answer: any purity objective is won outright by shattering the map into singletons, "
           "because a cluster of one is perfectly pure. Raise it if clusters look suspiciously tidy.",
    "eps": "DBSCAN neighbourhood radius, in embedding units. Ignored by HDBSCAN, which infers the "
           "equivalent per cluster instead of taking one value for the whole map.",
    # 5 Search
    "target": "The label to hold out and try to recover. It is excluded from the features along "
              "with anything that substantially restates it, so the map cannot see the answer it is "
              "being scored on -- which is the whole point of the exercise.",
    "search_sample": "Genes per configuration in the search. Same caveat as the walk sample: "
                     "smaller subsamples produce tighter, purer clusters and therefore better "
                     "scores than the full proteome will reproduce.",
    "max_blocks": "How many feature blocks may be combined in one configuration. The number of "
                  "combinations grows fast, so this bounds a walk that would otherwise run for "
                  "hours -- the count of configurations is reported when the search starts.",
    # 6 Validation
    "val_target": "The label to hide and try to recover. It must NOT be among the features the map "
                  "was built from: a cluster matching something the embedding already saw is "
                  "circular, and this refuses to score it rather than returning a flattering number.",
    "val_folds": "How many times to repeat the hide-and-recover test with a different random "
                 "selection. More folds give a steadier estimate; the spread across folds is what "
                 "tells you whether a single good result was luck.",
    "val_hold": "Fraction of each category's labelled genes hidden per fold. These are the genes "
                "the score is computed on, and they take no part in choosing which cluster to "
                "annotate from -- otherwise the test would be marking its own homework.",
    "val_refit": "Build a fresh map and clustering for every fold, with a different seed. The label "
                 "never feeds the embedding either way, so what this buys is that the answer stops "
                 "being a fact about one particular layout -- UMAP moves noticeably between seeds "
                 "at this size. It costs a full embedding per fold per category, so it is off by "
                 "default, and every row records which way it was obtained.",
}

#: Buttons, keyed by the method they call.
BUTTON_TOOLTIPS = {
    "show_variance": "Report how much of the feature matrix each block actually contributes. The "
                     "check that catches a block being named as an input while carrying almost "
                     "nothing -- hyperLOPIT came to 1.1%.",
    "run_umap_walk": "Score a grid of UMAP hyperparameters and rank them. Each configuration appears "
                     "in this table and as a thumbnail in the gallery the moment it is computed, so "
                     "a long sweep can be read while it runs; click either to open that map.",
    "run_embed": "Build ONE embedding from the current settings and show it in the 3D view, "
                 "replacing what is there.",
    "save_embedding": "Store this embedding with its full recipe, so it can be reloaded and "
                      "compared rather than rebuilt from memory of what the settings were.",
    "run_cluster_walk": "Score a grid of clustering hyperparameters against the current map.",
    "run_cluster": "Cluster the current map with these settings. Needs a map built first: this "
                   "clusters an embedding, it does not make one.",
    "run_battery": "Test what the clusters correspond to, using ONLY features the map was not built "
                   "from. A feature that fed the map separates the clusters by construction, so it "
                   "is evidence of nothing.",
    "run_search": "Walk dataset combinations and hyperparameters, scoring each by how well the "
                  "structure recovers the held-out label. This is the long one -- hundreds of "
                  "configurations, minutes to hours. It appears in Jobs and can be stopped there.",
    "run_validation": "Put an error rate on an annotation. Hides some genes that already carry a "
                      "category, picks the cluster holding most of the rest, and scores against the "
                      "hidden ones. A candidate list without this number is a list of guesses.",
}


from .jobs import Stopped as Cancelled  # noqa: E402  -- shared so the runner can recognise it

#: The walk's job name. Named once because it is also how a second walk is recognised and
#: refused: two sweeps filling one table and one gallery would read as a single sweep.
WALK_JOB = "UMAP hyperparameter walk"


class _Progress:
    """The `log` callable handed to a long analysis: reports progress, and honours a stop.

    These functions all report by calling `log` once per configuration, which makes it the one place
    that is guaranteed to be reached repeatedly without threading a cancellation flag through
    search.py, tuning.py and embedding.py. Raising from here unwinds the worker at the next line it
    prints, so a 328-configuration search stops in seconds rather than in half an hour.

    Cooperative by construction, which is the point: a search killed mid-write would leave a
    half-written embedding in the store.
    """

    def __init__(self, panel, job):
        self.panel, self.job = panel, job

    def __call__(self, message):
        if self.job.cancelled:
            raise Cancelled(f"{self.job.name} stopped after {self.job.note or 'some work'}")
        text = str(message)
        self.job.note = text
        # Queued to the GUI thread by Qt, because this is called from the worker.
        self.panel.status.emit(text)
        print(text, flush=True)          # and into the console pane, which is where a walk is read


class AnalysisPanel(QtWidgets.QWidget):
    """Data selection, tuning, clustering and the hypothesis battery."""

    embedding_ready = QtCore.pyqtSignal(object, object)          # coords, gene mask
    #: A clustering, so the map can colour by it. Without this the Clusters tab computed labels,
    #: printed how many there were, and threw them away -- which made the one thing this application
    #: is for, looking at structure coloured by a held-out variable, impossible to actually do.
    clusters_ready = QtCore.pyqtSignal(object)                   # labels, -1 for noise
    #: One finished configuration of a walk, as `tuning.WalkStep`, emitted from the worker thread the
    #: moment it is computed. Qt queues it to the GUI thread, which is what lets a row and a thumbnail
    #: appear while the walk is still running rather than all at once when it ends.
    walk_step = QtCore.pyqtSignal(object)
    #: A walk is starting: whatever the last one left on screen belongs to a different sweep.
    walk_started = QtCore.pyqtSignal()
    status = QtCore.pyqtSignal(str)

    def __init__(self, nodes: pd.DataFrame, store=None, parent=None, runner=None):
        super().__init__(parent)
        self.nodes = nodes
        self.store = store
        self.labels = None
        self.coords = None
        self.rows = None
        self._thread = None
        self._worker = None
        # The window's JobRunner, so this panel's work is visible and stoppable alongside everything
        # else. None is allowed: the panel is constructed without one in tests.
        self.runner = runner
        self._jobs = {}
        # What each results table is showing, what clicking one of its rows means, and what to call
        # the file if it is saved. Keyed by the widget so a table cannot be registered twice or
        # forgotten -- see `results_table`.
        self._frames, self._row_action, self._table_what = {}, {}, {}
        if runner is not None:
            runner.finished.connect(self._on_job_finished)
        # Connected to its own signal rather than called from the walk directly: `on_step` runs on
        # the worker thread, and filling a table from there is the crash this panel already documents
        # once. Going through the signal makes Qt queue it onto the GUI thread.
        self.walk_step.connect(self._walk_step_arrived)

        tabs = QtWidgets.QTabWidget()
        tabs.addTab(self._data_tab(), "1 · Data")
        tabs.addTab(self._map_tab(), "2 · Map")
        tabs.addTab(self._cluster_tab(), "3 · Clusters")
        tabs.addTab(self._meaning_tab(), "4 · Inference")
        tabs.addTab(self._search_tab(), "5 · Search")
        tabs.addTab(self._validation_tab(), "6 · Validation")
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(tabs)
        self._apply_tooltips()

    def _apply_tooltips(self):
        """Give every control an explanation, and put it on the label as well as the field.

        Three things this gets right that the first pass did not. The text is wrapped into a block
        rather than one line running off the screen. A bounded control says why its bounds are where
        they are, generated from the widget so it cannot drift from the range actually set. And the
        tooltip is attached to the form LABEL too -- which is the word people hover over, while the
        spin box beside it is the thing they click.
        """
        # Every attribute named in either table, so a control whose tooltip was set inline at
        # construction still gets wrapped and still gets its bounds explained.
        for attr in dict.fromkeys(list(TOOLTIPS) + list(LIMITS)):
            w = getattr(self, attr, None)
            if w is None:
                continue
            tip = TOOLTIPS.get(attr) or w.toolTip() or ""
            if attr in LIMITS and hasattr(w, "minimum"):
                tip = f"{tip}\n\n{range_note(w)} {LIMITS[attr]}"
            if tip:
                self._set_tip(w, wrap_tip(tip))
        for btn in self.findChildren(QtWidgets.QPushButton):
            tip = BUTTON_TOOLTIPS.get(self._handler_name(btn))
            if tip:
                btn.setToolTip(wrap_tip(tip))
        for name, cb in getattr(self, "block_cb", {}).items():
            cols = columns_for(self.nodes, EmbeddingSpec(blocks=(name,))).get(name, [])
            cb.setToolTip(wrap_tip(
                f"Feed the {name} block into the map: {len(cols)} columns"
                + (f", including {', '.join(cols[:4])}" if cols else "")
                + ". Anything fed in here cannot afterwards be used as evidence about the "
                  "clusters, because a feature the map was built on separates them by "
                  "construction."))

    def _set_tip(self, widget, html: str):
        """Set a tooltip on a widget and on its form label, so hovering the name works too."""
        widget.setToolTip(html)
        for form in self.findChildren(QtWidgets.QFormLayout):
            label = form.labelForField(widget)
            if label is not None:
                label.setToolTip(html)
                return

    @staticmethod
    def _handler_name(btn) -> str:
        """The name of the method a button is connected to, for keying its tooltip."""
        # Qt exposes no public way to read back a connection, so the button's own text is the key.
        # It is stable and unique across this panel, and a mismatch shows up as a missing tooltip,
        # which the audit test catches.
        return {
            "show what each block contributes": "show_variance",
            "walk hyperparameters": "run_umap_walk",
            "build this map": "run_embed",
            "save": "save_embedding",
            "walk clustering": "run_cluster_walk",
            "cluster this map": "run_cluster",
            "run the battery": "run_battery",
            "run the search": "run_search",
            "test the annotation": "run_validation",
        }.get(btn.text().strip(), "")

    # ------------------------------------------------------------------ 1 data
    def _data_tab(self):
        w = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(w)

        box = QtWidgets.QGroupBox("feature blocks")
        bl = QtWidgets.QVBoxLayout(box)
        self.block_cb = {}
        for b in BLOCKS:
            cols = columns_for(self.nodes, EmbeddingSpec(blocks=(b,))).get(b, [])
            cb = QtWidgets.QCheckBox(f"{b}  ({len(cols)} columns)")
            cb.setEnabled(bool(cols))
            cb.setChecked(b in ("expression_summary", "fitness_screens", "protein_features"))
            self.block_cb[b] = cb
            bl.addWidget(cb)
        self.cat_cb = QtWidgets.QCheckBox("compartment (hyperLOPIT, one-hot)")
        self.cat_cb.setChecked(True)
        bl.addWidget(self.cat_cb)
        v.addWidget(box)

        form = QtWidgets.QFormLayout()
        self.na_policy = QtWidgets.QComboBox()
        self.na_policy.addItems(NA_POLICIES)
        self.na_policy.setToolTip(
            "indicator adds a 0/1 column recording that a value was missing.\n"
            "drop_genes keeps only genes with no missing value at all -- on this table that is 3 of "
            "8,140, so it is offered rather than recommended.")
        self.scaling = QtWidgets.QComboBox()
        self.scaling.addItems(SCALINGS)
        self.scaling.setCurrentText("rank")
        self.scaling.setToolTip(
            "rank is the safe default: the published screens carry inverted sign conventions, ~64x "
            "differences in spread and heavy tails that z-scoring does not tame.")
        self.max_missing = QtWidgets.QDoubleSpinBox()
        self.max_missing.setRange(0.0, 1.0)
        self.max_missing.setSingleStep(0.05)
        self.max_missing.setValue(0.5)
        form.addRow("missing values", self.na_policy)
        form.addRow("scaling", self.scaling)
        form.addRow("max missing per column", self.max_missing)
        v.addLayout(form)

        self.variance_view = QtWidgets.QTextBrowser()
        self.variance_view.setMaximumHeight(150)
        btn = QtWidgets.QPushButton("show what each block contributes")
        btn.clicked.connect(self.show_variance)
        v.addWidget(btn)
        v.addWidget(self.variance_view)
        v.addStretch(1)
        return w

    def spec(self) -> EmbeddingSpec:
        """The EmbeddingSpec described by the current controls."""
        return EmbeddingSpec(
            blocks=tuple(b for b, cb in self.block_cb.items() if cb.isChecked()),
            categorical=("compartment",) if self.cat_cb.isChecked() else (),
            na_policy=self.na_policy.currentText(),
            scaling=self.scaling.currentText(),
            max_missing=self.max_missing.value(),
            n_neighbors=self.nn.value(), min_dist=self.md.value(),
            random_state=self.seed.value())

    def show_variance(self):
        """Report how much of the feature matrix each block actually carries."""
        try:
            v = variance_share(self.nodes, self.spec())
        except ValueError as e:
            self.variance_view.setPlainText(str(e))
            return
        rows = "".join(f"<tr><td>{i}</td><td align='right'>{s:.1%}</td></tr>"
                       for i, s in v["share"].items())
        self.variance_view.setHtml(
            "<p style='color:#888'>Share of the feature matrix each block carries. Equal weights give "
            "equal shares; if one block dominates, the map is mostly about that block.</p>"
            f"<table width='100%'>{rows}</table>")

    # ------------------------------------------------------------------ 2 map
    def _map_tab(self):
        w = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(w)
        form = QtWidgets.QFormLayout()
        self.nn = QtWidgets.QSpinBox(); self.nn.setRange(2, 400); self.nn.setValue(25)
        self.md = QtWidgets.QDoubleSpinBox(); self.md.setRange(0.0, 1.0)
        self.md.setSingleStep(0.05); self.md.setValue(0.25)
        self.seed = QtWidgets.QSpinBox(); self.seed.setRange(0, 10 ** 6); self.seed.setValue(42)
        self.sample = QtWidgets.QSpinBox(); self.sample.setRange(200, 20000)
        self.sample.setSingleStep(500); self.sample.setValue(2000)

        # The GRID the walk sweeps, as opposed to the single values above. These were hardcoded in
        # tuning.walk_umap and unreachable from the interface, so "walk hyperparameters" swept a set
        # nobody could see or change -- and the spin boxes above, which look like they control it,
        # only ever affected "build this map".
        self.nn_grid = QtWidgets.QLineEdit("5, 15, 25, 50, 100")
        self.md_grid = QtWidgets.QLineEdit("0.0, 0.1, 0.25, 0.5")
        self.mcs_grid = QtWidgets.QLineEdit("25, 60")
        form.addRow("n_neighbors", self.nn)
        form.addRow("min_dist", self.md)
        form.addRow("seed", self.seed)
        form.addRow("walk sample size", self.sample)
        form.addRow("walk: n_neighbors values", self.nn_grid)
        form.addRow("walk: min_dist values", self.md_grid)
        v.addLayout(form)

        row = QtWidgets.QHBoxLayout()
        b1 = QtWidgets.QPushButton("walk hyperparameters")
        b1.clicked.connect(self.run_umap_walk)
        b2 = QtWidgets.QPushButton("build this map")
        b2.setProperty("primary", True)
        b2.clicked.connect(self.run_embed)
        row.addWidget(b1); row.addWidget(b2)
        v.addLayout(row)

        self.walk_table = self.results_table(
            QtWidgets.QTableWidget(), self.show_walk_row, "umap_walk")
        self.walk_table.setToolTip(
            "Click a row to build that configuration and show it in the 3D view, clustered the way "
            "the walk clustered it. A table of scores is not a map, and the point of a walk is to "
            "look at the ones that scored well. Right-click to save the whole table as CSV.")
        v.addWidget(self.walk_table, 1)
        hint = QtWidgets.QLabel("<i>Click a row to build and show that map.</i>")
        hint.setWordWrap(True)
        v.addWidget(hint)

        save = QtWidgets.QHBoxLayout()
        self.emb_name = QtWidgets.QLineEdit(); self.emb_name.setPlaceholderText("name this embedding")
        bs = QtWidgets.QPushButton("save"); bs.clicked.connect(self.save_embedding)
        save.addWidget(self.emb_name, 1); save.addWidget(bs)
        v.addLayout(save)
        return w

    # ------------------------------------------------------------------ 3 clusters
    def _cluster_tab(self):
        w = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(w)
        form = QtWidgets.QFormLayout()
        self.algo = QtWidgets.QComboBox(); self.algo.addItems(["hdbscan", "dbscan"])
        self.mcs = QtWidgets.QSpinBox(); self.mcs.setRange(3, 2000); self.mcs.setValue(25)
        self.eps = QtWidgets.QDoubleSpinBox(); self.eps.setRange(0.001, 50.0)
        self.eps.setDecimals(3); self.eps.setValue(0.5)
        form.addRow("algorithm", self.algo)
        form.addRow("min_cluster_size / min_samples", self.mcs)
        form.addRow("eps (dbscan)", self.eps)
        v.addLayout(form)
        row = QtWidgets.QHBoxLayout()
        b1 = QtWidgets.QPushButton("walk clustering")
        b1.clicked.connect(self.run_cluster_walk)
        b2 = QtWidgets.QPushButton("cluster this map")
        b2.setProperty("primary", True)
        b2.clicked.connect(self.run_cluster)
        row.addWidget(b1); row.addWidget(b2)
        v.addLayout(row)
        self.cluster_table = self.results_table(
            QtWidgets.QTableWidget(), self.show_cluster_row, "clustering_walk")
        self.cluster_table.setToolTip(
            "Click a row to cluster the current map with those settings and colour it by the "
            "result. A silhouette for a clustering nobody can see is a number about nothing. "
            "Right-click to save the whole table as CSV.")
        v.addWidget(self.cluster_table, 1)
        return w

    # ------------------------------------------------------------------ 4 meaning
    def _meaning_tab(self):
        w = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(w)
        note = QtWidgets.QLabel(
            "Only <b>held-out</b> features are evidence. A feature that fed the map separates the "
            "clusters by construction, and anything that restates it does too.")
        note.setWordWrap(True)
        v.addWidget(note)
        b = QtWidgets.QPushButton("run the battery")
        b.setProperty("primary", True)
        b.clicked.connect(self.run_battery)
        v.addWidget(b)
        self.findings = QtWidgets.QTextBrowser()
        v.addWidget(self.findings, 1)
        self.battery_table = self.results_table(
            QtWidgets.QTableWidget(), None, "held_out_battery")
        self.battery_table.setToolTip(
            "One row per held-out feature. The score is over the whole feature at once, which is "
            "why the per-category table below it matters: the same Cramer's V describes 27 "
            "compartments weakly smeared everywhere and one compartment falling out cleanly.")
        v.addWidget(self.battery_table, 1)

        v.addWidget(QtWidgets.QLabel(
            "<i>Per category: the best single cluster for each value, which is what a feature-level "
            "score hides.</i>"))
        self.category_table = self.results_table(
            QtWidgets.QTableWidget(), self.show_inference_row, "per_category")
        self.category_table.setToolTip(
            "For each category of each held-out feature, the cluster that matches it best: "
            "precision over that cluster, recall over that category, and F1. Read both — a cluster "
            "that is 100% apicoplast holding 5% of apicoplast proteins is useless for inference, "
            "and one number cannot tell you which of the two you have.")
        v.addWidget(self.category_table, 1)
        return w

    # ------------------------------------------------------------------ 5 search
    def _search_tab(self):
        w = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(w)
        note = QtWidgets.QLabel(
            "Walk combinations of datasets and hyperparameters looking for a structure that recovers a "
            "label it was never given. The target and everything that restates it are excluded from "
            "every map, which is what makes the recovery meaningful — and what turns it into a "
            "prediction for the unlabelled genes in a pure cluster.")
        note.setWordWrap(True)
        v.addWidget(note)
        form = QtWidgets.QFormLayout()
        self.target = QtWidgets.QComboBox()
        # Read from search.TARGETS rather than listed here, so a target added to the registry is
        # reachable from the interface. Hardcoded, `cellcycle_phase` existed in the module and could
        # not be chosen -- half the project's stated purpose, present in the data and absent from the
        # menu. Filtered by what this table actually has, so an entry can never name a missing column.
        from .search import TARGETS
        self.target.addItems([c for c in TARGETS.values() if c in self.nodes.columns])
        self.target.setToolTip(
            "The label held out of every embedding and then scored. Measured labels "
            "(compartment, cellcycle_phase) are evidence; stage_enriched_derived is computed from "
            "expression and is a positive control, and attention_depth is a negative control -- if a "
            "structure recovers how much a gene has been STUDIED, the map is measuring the literature.")
        self.search_sample = QtWidgets.QSpinBox()
        self.search_sample.setRange(500, 20000); self.search_sample.setSingleStep(500)
        self.search_sample.setValue(3000)
        self.max_blocks = QtWidgets.QSpinBox(); self.max_blocks.setRange(1, 4); self.max_blocks.setValue(2)

        # Which objective a walk maximises is the most consequential choice in it, and it used to be
        # made for the user. See objectives.EXPLANATION, shown under Help.
        from .objectives import OBJECTIVES
        self.objective = QtWidgets.QComboBox()
        for name, what in OBJECTIVES.items():
            self.objective.addItem(f"{what}  [{name}]", name)
        self.objective.setCurrentIndex(list(OBJECTIVES).index("mean_f1"))
        self.objective.setToolTip(
            "What counts as good structure. Precision asks of a CLUSTER what fraction of its members "
            "share a label; recall asks of a LABEL what fraction of its genes share a cluster. Every "
            "one of these has a degenerate solution that wins it outright -- purity by shattering "
            "into singletons, recall by one giant cluster -- so read the coverage and cluster count "
            "beside the score. Help has the full table.")
        self.weighting = QtWidgets.QComboBox()
        self.weighting.addItems(["macro", "size"])
        self.weighting.setToolTip(
            "macro weights every label equally; size weights every gene equally. On this proteome "
            "nucleus-chromatin has 769 genes and dense granules 167, so size weighting means a "
            "search for dense granules is decided by the nucleus. macro is the default for that.")
        # A list rather than a dropdown: "precision for dense granules" and "recall for dense
        # granules and rhoptries" are both ordinary questions, and the second needs more than one.
        # Select none to score every category.
        self.focus = QtWidgets.QListWidget()
        self.focus.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.ExtendedSelection)
        self.focus.setMaximumHeight(120)
        self.focus.setToolTip(
            "Which categories to score. Select none for all of them; select one or several to "
            "optimise for exactly those. The objective and the label set are independent, so "
            "'precision for dense granules' and 'recall for dense granules and rhoptries' are both "
            "reachable. With several, mean objectives average over them and best objectives take "
            "the best among them.")
        # The SCORING floors, which are not the same thing as HDBSCAN's min_cluster_size. That one
        # decides what the algorithm FORMS as a cluster and is swept in the grid above; these decide
        # what is allowed to COUNT once clusters exist. Both matter and they were only half exposed:
        # a cluster too small to be a landmark could still win a purity objective.
        self.score_min_cluster = QtWidgets.QSpinBox()
        self.score_min_cluster.setRange(1, 2000); self.score_min_cluster.setValue(10)
        self.score_min_cluster.setToolTip(
            "The smallest cluster allowed to COUNT when scoring. Distinct from HDBSCAN's "
            "min_cluster_size above, which decides what gets formed in the first place. This is the "
            "guard against the singleton exploit: any purity objective is won outright by a cluster "
            "of one, so nothing below this is scored at all.")
        self.score_min_label = QtWidgets.QSpinBox()
        self.score_min_label.setRange(2, 2000); self.score_min_label.setValue(15)
        self.score_min_label.setToolTip(
            "The smallest label allowed to COUNT. A category with five genes can be captured "
            "perfectly by accident, and scoring it rewards luck. Raising this drops rare "
            "compartments from the score entirely, so it trades noise for coverage.")
        self.min_recall = QtWidgets.QDoubleSpinBox()
        self.min_recall.setRange(0.0, 1.0); self.min_recall.setSingleStep(0.05)
        self.min_recall.setValue(0.25)
        self.min_recall.setToolTip(
            "Used only by precision_at_recall: the least of a label a cluster must hold to count. "
            "Precision alone selects three co-located genes at 1.00 and gives you nothing to "
            "annotate from, so this is the floor that makes the answer usable.")
        # The same grid the Map tab shows, mirrored here because this is where it is swept. Kept in
        # step both ways so the two tabs cannot disagree about what a walk will do.
        self.nn_grid2 = QtWidgets.QLineEdit(self.nn_grid.text())
        self.md_grid2 = QtWidgets.QLineEdit(self.md_grid.text())
        # Mirrored with an explicit guard rather than by connecting each to the other's setText.
        # A plain two-way binding re-enters -- setText emits textChanged, which sets the first again
        # -- and the pair can still be firing at each other while Qt is deleting them.
        def mirror(src, dst):
            def on_change(text):
                if dst.text() != text:
                    dst.blockSignals(True)
                    dst.setText(text)
                    dst.blockSignals(False)
            src.textChanged.connect(on_change)
        for a, b in ((self.nn_grid, self.nn_grid2), (self.md_grid, self.md_grid2)):
            mirror(a, b)
            mirror(b, a)
        self.target.currentTextChanged.connect(self._refresh_focus)
        self._refresh_focus(self.target.currentText())

        form.addRow("target to recover", self.target)
        form.addRow("optimise for", self.objective)
        form.addRow("weighting", self.weighting)
        form.addRow("categories to score", self.focus)
        form.addRow("minimum recall", self.min_recall)
        form.addRow("min cluster size to score", self.score_min_cluster)
        form.addRow("min label size to score", self.score_min_label)
        form.addRow("walk: n_neighbors values", self.nn_grid2)
        form.addRow("walk: min_dist values", self.md_grid2)
        form.addRow("walk: min_cluster_size values", self.mcs_grid)
        form.addRow("sample size", self.search_sample)
        form.addRow("max blocks per combination", self.max_blocks)
        v.addLayout(form)
        b = QtWidgets.QPushButton("run the search")
        b.setProperty("primary", True)
        b.clicked.connect(self.run_search)
        v.addWidget(b)
        self.search_table = self.results_table(
            QtWidgets.QTableWidget(), self.show_search_row, "recovery_search")
        self.search_table.setToolTip(
            "Click a row to rebuild that exact configuration -- same blocks, same policy, same "
            "seed, same subsample, same excluded columns -- and show it with its clustering. A "
            "recovery score with no way to look at the structure it scored is a 'trust me'. "
            "Right-click to save the whole table as CSV.")
        v.addWidget(self.search_table, 1)
        return w

    @staticmethod
    def parse_grid(text: str, cast=float, fallback=()):
        """A comma-separated list of values to sweep, or a min:max:step range.

        Accepts "5, 15, 25" and "5:100:5" alike, because both are natural ways to say the same
        thing and refusing one of them is a papercut on the control used most. An unparsable entry
        falls back to the default rather than raising: a walk is expensive to start, and losing one
        to a stray comma is worse than sweeping the default set.
        """
        text = (text or "").strip()
        if not text:
            return tuple(fallback)
        try:
            if ":" in text:
                parts = [float(x) for x in text.split(":")]
                lo, hi = parts[0], parts[1]
                step = parts[2] if len(parts) > 2 else 1.0
                if step <= 0:
                    return tuple(fallback)
                out, v = [], lo
                while v <= hi + 1e-9:
                    out.append(cast(v))
                    v += step
                return tuple(out)
            return tuple(cast(x) for x in text.replace(",", " ").split())
        except (ValueError, IndexError):
            return tuple(fallback)

    def walk_grid(self) -> dict:
        """The hyperparameter grid the walk sweeps, parsed from the interface."""
        return {
            "n_neighbors_values": self.parse_grid(self.nn_grid.text(), int, (5, 15, 25, 50, 100)),
            "min_dist_values": self.parse_grid(self.md_grid.text(), float, (0.0, 0.1, 0.25, 0.5)),
            "min_cluster_sizes": self.parse_grid(self.mcs_grid.text(), int, (25, 60)),
        }

    def _refresh_focus(self, target: str):
        """Offer the categories of the chosen target, so 'focus on one' means something.

        Rebuilt whenever the target changes: the categories of `compartment` are not the categories
        of `cellcycle_phase`, and a stale list would let someone optimise for a label that no longer
        exists in the column being scored.
        """
        from .search import ABSENCE_LABELS
        self.focus.blockSignals(True)
        self.focus.clear()
        if target in self.nodes.columns:
            v = self.nodes[target].astype("object").where(self.nodes[target].notna(), "").astype(str)
            counts = v.value_counts()
            for name, n in counts.items():
                if str(name).lower() in ABSENCE_LABELS or n < 15:
                    continue
                it = QtWidgets.QListWidgetItem(f"{name}  ({n})")
                it.setData(QtCore.Qt.ItemDataRole.UserRole, str(name))
                self.focus.addItem(it)
        self.focus.blockSignals(False)

    def selected_categories(self) -> list:
        """The categories highlighted in the list, or an empty list meaning all of them."""
        return [i.data(QtCore.Qt.ItemDataRole.UserRole) for i in self.focus.selectedItems()]

    def objective_settings(self) -> dict:
        """The scoring choices, as the keyword arguments `objectives.score` takes.

        Returned as one dict so the settings travel together into the job and into whatever records
        the run -- a score whose objective is not recorded beside it cannot be compared with another.
        """
        return {"objective": self.objective.currentData(),
                "weighting": self.weighting.currentText(),
                # None means every category, which is what selecting nothing should mean.
                "category": self.selected_categories() or None,
                "min_recall": self.min_recall.value(),
                "min_cluster": self.score_min_cluster.value(),
                "min_label": self.score_min_label.value()}

    # ------------------------------------------------------------------ 6 validation
    def _validation_tab(self):
        """Test an annotation by hiding labels we already have and seeing whether they come back."""
        w = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(w)
        note = QtWidgets.QLabel(
            "A cluster that looks pure gives you a candidate list and <b>no error rate</b>. This "
            "hides a fraction of the genes that already carry a category, picks the cluster holding "
            "most of the rest \u2014 as you would by eye \u2014 and scores against the hidden ones, "
            "which had no say in that choice. Reported per category, because a method that recovers "
            "hidden dense granules but not hidden rhoptries is not one accuracy.")
        note.setWordWrap(True)
        v.addWidget(note)

        form = QtWidgets.QFormLayout()
        self.val_target = QtWidgets.QComboBox()
        self.val_target.addItems([c for c in ("compartment", "compartment_best", "cellcycle_phase",
                                              "stage_enriched_derived")
                                  if c in self.nodes.columns])
        self.val_target.setToolTip(
            "The label to hide. It must NOT be among the features the map was built from: a cluster "
            "matching something the embedding already saw is circular, and this refuses to score it.")
        self.val_folds = QtWidgets.QSpinBox(); self.val_folds.setRange(2, 20)
        self.val_folds.setValue(5)
        self.val_hold = QtWidgets.QDoubleSpinBox(); self.val_hold.setRange(0.05, 0.5)
        self.val_hold.setSingleStep(0.05); self.val_hold.setValue(0.2)
        self.val_hold.setToolTip("Fraction of each category's labelled genes hidden per fold.")
        self.val_refit = QtWidgets.QCheckBox("re-embed and re-cluster for every fold")
        form.addRow("hold out", self.val_target)
        form.addRow("folds", self.val_folds)
        form.addRow("fraction hidden", self.val_hold)
        form.addRow("stricter", self.val_refit)
        v.addLayout(form)

        b = QtWidgets.QPushButton("test the annotation")
        b.setProperty("primary", True)
        b.clicked.connect(self.run_validation)
        v.addWidget(b)

        self.val_note = QtWidgets.QLabel("")
        self.val_note.setWordWrap(True)
        self.val_note.hide()
        v.addWidget(self.val_note)

        self.val_table = self.results_table(
            QtWidgets.QTableWidget(), self.show_candidates, "validation")
        self.val_table.setToolTip(
            "One row per category, because a method that recovers hidden dense granules but not "
            "hidden rhoptries is not one accuracy. `refit` says whether each number came from a "
            "fresh map per fold or from one fixed map. Click a row for that category's candidates.")
        v.addWidget(self.val_table, 1)

        v.addWidget(QtWidgets.QLabel(
            "<i>Click a category above for the genes its cluster would have you annotate.</i>"))
        self.cand_table = self.results_table(
            QtWidgets.QTableWidget(), None, "candidates")
        self.cand_table.setToolTip(
            "The unlabelled members of that category's cluster: the list this whole tab exists to "
            "put a number on. Every row carries how much of its cluster already carries the "
            "category and how much carries something else, and how many genes already labelled it "
            "share an orthogroup, a Pfam or an InterPro domain with the candidate — agreement "
            "from evidence the map never saw. Zero support is the common answer.")
        v.addWidget(self.cand_table, 1)
        return w

    def used_columns(self) -> list:
        """The node-table columns the current embedding is built from, categoricals included.

        The real column names, not the block names. The block names were what the validation tab
        passed as "the columns the embedding used", which made the circularity guard compare a
        compartment against "expression_summary" -- so it never fired, and the tab would have
        scored a target the map was built on.
        """
        cols = [c for group in columns_for(self.nodes, self.spec()).values() for c in group]
        if self.cat_cb.isChecked():
            cols.append("compartment")
        return cols

    def _aligned_truth(self, target: str):
        """The label column restricted to the genes the clustering actually covers.

        `labels` comes from an embedding that may have dropped genes -- `drop_genes` keeps 3 of
        8,140 on this table -- and scoring a clustering of one set of genes against the labels of
        another compares gene i's cluster with gene j's compartment. The battery already aligned
        through `self.rows`; validation did not.
        """
        truth = self.nodes[target]
        if self.rows is not None and len(self.labels) != len(truth):
            truth = truth[self.rows]
        return truth

    def run_validation(self):
        """Put an error rate on an annotation by hiding labels that are already known."""
        from .validate import circularity_error, validate_all
        if self.labels is None:
            self.status.emit("cluster a map first -- validation scores a clustering, not a map")
            return
        target = self.val_target.currentText()
        used = self.used_columns()
        # Checked here as well as inside the job. `validate_all` refuses too -- that is the guard
        # that matters -- but going through the worker to find out means waiting for a job to fail
        # in order to be told the question cannot be asked, and a failed job is exactly the wrong
        # shape for that answer.
        problem = circularity_error(self.nodes[target], used, target)
        if problem:
            self._validation_refused(ValueError(problem))
            return
        truth, labels = self._aligned_truth(target), self.labels
        folds, frac = self.val_folds.value(), self.val_hold.value()
        refit = self.val_refit.isChecked()
        rebuild = self._refit_labels if refit else None
        if len(labels) != len(truth):
            self.status.emit(f"the clustering covers {len(labels):,} genes and {target} has "
                             f"{len(truth):,} -- rebuild the map, then cluster it again")
            return

        def job(p):
            p(f"hiding {frac:.0%} of each category in {target}, {folds} folds"
              + (" , re-embedding each one" if refit else ""))
            return validate_all(labels, truth, folds=folds, hold_frac=frac, used_columns=used,
                                target_column=target, refit=refit, rebuild=rebuild, log=p)

        self._validation_target = target
        self._run(job, self._validation_done, name=f"validate annotation ({target})",
                  on_error=self._validation_refused)

    def _refit_labels(self, fold: int):
        """One fold's clustering, from a map built again with a different seed.

        What re-fitting buys is not that the label is hidden from the embedding -- it never fed it,
        and the tab refuses when it did. It is that the estimate stops being conditional on one map:
        UMAP moves noticeably between seeds at this size, and a precision measured on a single
        layout is a fact about that layout. It costs a full embedding per fold, which is why it is
        off by default and says so.
        """
        from dataclasses import asdict
        from .clustering import cluster
        from .embedding import EmbeddingSpec, embed
        spec = EmbeddingSpec(**{**asdict(self.spec()),
                                "random_state": int(self.seed.value()) + 1 + int(fold)})
        coords, _, _ = embed(self.nodes, spec, log=lambda *a: None)
        return cluster(coords, algorithm=self.algo.currentText(),
                       min_cluster_size=self.mcs.value(), min_samples=self.mcs.value(),
                       eps=self.eps.value())

    def _validation_done(self, d):
        """Show the per-category scores, and clear whatever candidates the last run left."""
        self.val_note.hide()
        self._validation_scores = d
        self._fill(self.val_table, d)
        self.cand_table.clear()
        self.cand_table.setRowCount(0)
        self.cand_table.setColumnCount(0)
        self.status.emit(self._validation_verdict(d))

    def _validation_refused(self, error) -> bool:
        """A refusal is an explanation, not a failure. Returns whether it was handled here.

        `validate` raises when the label fed the embedding, which is the guard working. Reported as
        a red failed job with a traceback, it would read as the tab being broken -- and the next
        move would be to look for the bug rather than to rebuild the map without that column.
        """
        if not isinstance(error, ValueError):
            return False
        # Rendered and looked at: the note said "Not scored" while the previous run's table sat
        # underneath it, and a table of precisions under an explanation reads as the explanation's
        # result. A refusal must leave no number on screen that could be taken for this one's.
        self._validation_scores = None
        for table in (self.val_table, self.cand_table):
            table.clear()
            table.setRowCount(0)
            table.setColumnCount(0)
        self.val_note.setText(f"<b>Not scored.</b> {error}")
        self.val_note.show()
        self.status.emit(str(error))
        return True

    def show_candidates(self, row: int, _col: int = 0):
        """The genes one category's cluster would have you annotate, with the numbers attached.

        Never a bare list. Each candidate carries the composition of the cluster it comes from and
        whatever independent agreement exists, because the same list looks identical whether it is
        90% right or 6% right -- and on this proteome it has been 6%.
        """
        from .validate import best_cluster, candidates, orthogonal_support
        d = getattr(self, "_validation_scores", None)
        target = getattr(self, "_validation_target", None)
        if d is None or not len(d) or row >= len(d) or self.labels is None:
            return
        item = self.val_table.item(row, 0)
        category = item.text() if item is not None else str(d.iloc[row].category)
        truth = self._aligned_truth(target)
        cl = best_cluster(self.labels, truth, category)
        if cl is None:
            self.status.emit(f"no cluster holds any gene labelled {category!r}")
            return
        genes = self.nodes.gene_id[self.rows] if self.rows is not None else self.nodes.gene_id
        cand = candidates(self.labels, truth, category, cl, genes)
        if len(cand):
            support = orthogonal_support(self.nodes, cand.gene_id, truth, category,
                                         used_columns=self.used_columns(),
                                         log=lambda m: self.status.emit(m))
            cand = cand.merge(support, on="gene_id", how="left")
        self._fill(self.cand_table, cand)
        score = d.iloc[row]
        self.status.emit(
            f"{len(cand):,} candidates for {category} from cluster {cl} -- validated precision "
            f"{score.precision:.2f}, so roughly {score.precision:.0%} of them would be right")

    @staticmethod
    def _validation_verdict(d) -> str:
        """One line leading with the number a candidate list is meaningless without."""
        if d is None or not len(d):
            return "no category had enough labelled genes to hide any"
        best = d.iloc[0]
        return (f"best: {best.category} at precision {best.precision:.2f} -- annotate from that "
                f"cluster and roughly {best.precision:.0%} would be right")

    # ------------------------------------------------------------------ jobs
    def _run(self, fn, on_done, name: str = "analysis", on_error=None):
        """Run `fn` off the GUI thread, through the window's job runner when there is one.

        Routed through JobRunner rather than through a private QThread. The private one was the cause
        of three separate complaints at once: its work never appeared in the Jobs panel, there was no
        way to stop it, and -- worst -- it refused to start while anything else was running. A search
        takes minutes, so for those minutes every button on every other tab silently did nothing but
        write "a job is already running" into the status bar, which reads exactly like four broken
        tabs.

        `on_done` is invoked on the GUI thread. A handler connected straight to a worker signal runs
        on the worker, and touching widgets from there crashes on a slow machine.

        `on_error` is given the exception and returns whether it handled it. Some failures are not
        failures: validation raising because the target fed the embedding is the guard doing its
        job, and reporting it in red with a traceback teaches people to distrust the guard.
        """
        if self.runner is not None:
            job = self.runner.submit(lambda j: fn(_Progress(self, j)), name)
            self._jobs[job.id] = (on_done, on_error)
            return job
        # No runner (the panel used standalone, or in a test): fall back to a private thread.
        if self._thread is not None:
            self.status.emit("a job is already running")
            return None
        self._thread = QtCore.QThread(self)
        self._worker = Worker(fn)
        self._worker.moveToThread(self._thread)
        self._worker.progress.connect(self.status.emit)
        self._thread.started.connect(self._worker.run)

        def relay(result, error):
            self._thread.quit(); self._thread.wait()
            self._thread = None; self._worker = None
            if error is not None:
                if not (on_error is not None and on_error(error)):
                    self.status.emit(f"failed: {error}")
            else:
                on_done(result)

        self._worker.done.connect(relay)
        self._thread.start()
        return None

    def _on_job_finished(self, jid: int, ok: bool):
        """Deliver a finished job's result to its handler, on the GUI thread."""
        handlers = self._jobs.pop(jid, None)
        if handlers is None:
            return
        on_done, on_error = handlers
        job = self.runner.jobs.get(jid)
        if job is None:
            return
        if job.state == "cancelled":
            self.status.emit(f"{job.name}: stopped")
            return
        if not ok:
            # The exception itself where the runner kept it, so a handler can tell a deliberate
            # refusal from a crash by type rather than by reading the formatted message.
            exc = job.exception if job.exception is not None else RuntimeError(job.error)
            if on_error is not None and on_error(exc):
                return
            self.status.emit(f"{job.name} failed: {job.error}")
            return
        on_done(job.result)

    def results_table(self, table: QtWidgets.QTableWidget, on_row=None, what: str = "these results"):
        """Give a results table the two things every results table needs.

        Clicking a row shows the map that row is about -- rebuilt from the row's own configuration
        where the row names one, and the current map with its clustering where the row is about a
        category of it. A table of scores is not a result; the map it describes is, and a score
        nobody can look at is the thing this application exists not to produce.

        Right-clicking saves it. A table that can only be read on screen has to be re-derived
        anywhere else it is needed, and the run that produced it is minutes long.

        Wired here rather than per tab so a new table cannot arrive without either.
        """
        table.setAlternatingRowColors(True)
        table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        table.setContextMenuPolicy(QtCore.Qt.ContextMenuPolicy.CustomContextMenu)
        table.customContextMenuRequested.connect(
            lambda pos, t=table: self._table_menu(t, pos))
        self._row_action[table] = on_row
        self._table_what[table] = what
        if on_row is not None:
            table.cellClicked.connect(lambda row, _col, t=table: self._row_clicked(t, row))
        return table

    def _row_clicked(self, table, row: int):
        action = self._row_action.get(table)
        if action is not None:
            action(row)

    def _table_menu(self, table, pos):
        """Show the right-click menu for a results table. `build_table_menu` makes it."""
        menu = self.build_table_menu(table)
        menu.exec(table.viewport().mapToGlobal(pos))
        return menu

    def build_table_menu(self, table) -> QtWidgets.QMenu:
        """Construct a results table's menu without showing it.

        Split from `_table_menu` for the reason the window's menus are: `exec` enters a modal loop
        and does not return until a human closes the menu, so a test that called it would hang
        rather than fail.
        """
        m = QtWidgets.QMenu(self)
        df = self._frames.get(table)
        act = m.addAction("Save this table as CSV…")
        act.setEnabled(df is not None and len(df) > 0)
        act.triggered.connect(lambda: self.save_table(table))
        copy = m.addAction("Copy selected rows")
        copy.setEnabled(bool(table.selectedItems()))
        copy.triggered.connect(lambda: self.copy_rows(table))
        if self._row_action.get(table) is not None:
            m.addSeparator()
            show = m.addAction("Show this row's map")
            row = table.currentRow()
            show.setEnabled(row >= 0)
            show.triggered.connect(lambda: self._row_clicked(table, table.currentRow()))
        return m

    def save_table(self, table, path: str = "") -> str:
        """Write a results table to CSV, whole rather than as displayed.

        The full frame, not the 200 rows the widget shows: the truncation is there to keep the
        window responsive, and a file that silently stopped at row 200 would be a different result
        from the one that was computed. How many rows were written is reported for that reason.
        """
        df = self._frames.get(table)
        if df is None or not len(df):
            self.status.emit("nothing to save -- run something first")
            return ""
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(
                self, "Save table as CSV", f"{self._table_what.get(table, 'results')}.csv",
                "CSV (*.csv)")
        if not path:
            return ""
        df.to_csv(path, index=False)
        shown = min(len(df), 200)
        self.status.emit(f"wrote {len(df):,} rows to {path}"
                         + (f" (the table shows the first {shown})" if len(df) > shown else ""))
        return path

    def copy_rows(self, table) -> str:
        """Put the selected rows on the clipboard, tab-separated, with their headers."""
        rows = sorted({i.row() for i in table.selectedItems()})
        if not rows:
            self.status.emit("select a row first")
            return ""
        headers = [table.horizontalHeaderItem(c).text() for c in range(table.columnCount())]
        lines = ["\t".join(headers)]
        for r in rows:
            lines.append("\t".join((table.item(r, c).text() if table.item(r, c) else "")
                                   for c in range(table.columnCount())))
        text = "\n".join(lines)
        cb = QtWidgets.QApplication.clipboard()
        if cb is not None:
            cb.setText(text)
        self.status.emit(f"copied {len(rows)} row(s)")
        return text

    def row_values(self, table, row: int) -> dict:
        """One row as {column: text}, which is how a row is turned back into a configuration."""
        out = {}
        for c in range(table.columnCount()):
            head = table.horizontalHeaderItem(c)
            item = table.item(row, c)
            if head is not None and item is not None:
                out[head.text()] = item.text()
        return out

    def _publish_clusters(self, labels, genes=None):
        """Send a clustering to the map, over the whole node table.

        A clustering of a subsample covers only the genes that map covers, and the window colours
        8,140 points by it. Expanded here with -1 -- unclustered, which is drawn grey -- rather than
        left short, because a labels array of the wrong length made the window fall back to
        colouring everything grey, which reads as "this clustering found nothing".
        """
        labels = np.asarray(labels)
        rows = self.rows if genes is None else genes
        if rows is not None and np.size(rows) == len(self.nodes) and len(labels) != len(self.nodes):
            rows = np.asarray(rows).astype(bool)
            if int(rows.sum()) != len(labels):
                # Neither the map's length nor the table's, which means the clustering and the map
                # on screen are not of the same genes. Emitted unchanged so the window's own length
                # check draws it grey rather than putting cluster i's colour on gene j -- and said
                # out loud, because grey everywhere otherwise reads as "this found nothing".
                self.status.emit(f"the clustering covers {len(labels):,} genes and the map covers "
                                 f"{int(rows.sum()):,} -- cluster this map again")
                self.clusters_ready.emit(labels)
                return
            full = np.full(len(self.nodes), -1, dtype=int)
            full[rows] = labels
            labels = full
        self.clusters_ready.emit(labels)

    def _fill(self, table: QtWidgets.QTableWidget, df: pd.DataFrame, limit=200):
        """Fill a table, sortable by any column.

        Sorting is disabled while the rows go in and re-enabled afterwards: with it left on, Qt
        re-sorts after every insertion and the rows end up interleaved. Numbers are stored as
        numbers rather than as their formatted text, so a score column sorts 0.9 above 0.10 instead
        of lexically.
        """
        # The whole frame is remembered before the display is truncated, so saving writes the
        # result rather than the first screenful of it.
        self._frames[table] = df
        df = df.head(limit)
        table.setSortingEnabled(False)
        table.clear()
        table.setRowCount(len(df)); table.setColumnCount(len(df.columns))
        table.setHorizontalHeaderLabels([str(c) for c in df.columns])
        for i, (_, r) in enumerate(df.iterrows()):
            for j, v in enumerate(r):
                table.setItem(i, j, self._cell(v))
        table.resizeColumnsToContents()
        table.setSortingEnabled(True)

    @staticmethod
    def _cell(v) -> QtWidgets.QTableWidgetItem:
        """One table cell: a number stored as a number, so the column sorts 0.9 above 0.10."""
        item = QtWidgets.QTableWidgetItem()
        if isinstance(v, (int, float, np.integer, np.floating)) and np.isfinite(v):
            item.setData(QtCore.Qt.ItemDataRole.DisplayRole, float(v))
            item.setText(f"{v:.3f}" if isinstance(v, (float, np.floating)) else str(v))
        else:
            item.setText(str(v))
        return item

    def _start_table(self, table: QtWidgets.QTableWidget, columns):
        """Empty a table and give it headers, ready for rows to arrive one at a time."""
        table.setSortingEnabled(False)          # re-enabled by _fill when the run finishes
        # The last run's frame goes with it: saving a table that has been emptied on screen should
        # not write the previous walk's rows.
        self._frames.pop(table, None)
        table.clear()
        table.setRowCount(0)
        table.setColumnCount(len(columns))
        table.setHorizontalHeaderLabels([str(c) for c in columns])

    def _append(self, table: QtWidgets.QTableWidget, row: dict):
        """Add one result to the bottom of a table, in the order it was computed.

        Sorting stays off while a run streams: with it on, Qt re-sorts after every insertion and a
        row the user is reading moves under the pointer. The ranked table replaces this one when the
        run finishes, which is the point at which a ranking means anything -- ranking a walk that is
        one configuration in says only that one configuration has run.
        """
        if table.columnCount() == 0:
            self._start_table(table, list(row))
        headers = [table.horizontalHeaderItem(c).text() for c in range(table.columnCount())]
        i = table.rowCount()
        # Streamed rows are collected as they arrive, so a walk that is still running -- or one that
        # was stopped half way -- can be saved for what it found.
        prev = self._frames.get(table)
        one = pd.DataFrame([row])
        self._frames[table] = one if prev is None else pd.concat([prev, one], ignore_index=True)
        table.insertRow(i)
        for j, name in enumerate(headers):
            if name in row:
                table.setItem(i, j, self._cell(row[name]))
        table.resizeColumnsToContents()

    def _walk_step_arrived(self, step):
        """One configuration finished: put its scores in the table and say where the walk is.

        Runs on the GUI thread -- see the connection in `__init__`. The gallery is fed from the same
        signal by the window, so the row and the thumbnail appear together.
        """
        self._append(self.walk_table, step.row)
        self.status.emit(f"walk {step.index} of {step.total}: {step.label}")

    def run_umap_walk(self):
        """Sweep a grid of UMAP hyperparameters, reporting each configuration as it finishes.

        The walk emits per configuration and the table fills a row at a time, so a sweep of 288
        settings -- half an hour -- can be read while it runs instead of showing nothing until it
        ends. The ranking still arrives at the end, because ranking needs the whole sweep.

        Each embedding is written through the store as it is computed, so a walk that is stopped
        half way still leaves behind every configuration it finished, with the recipe to rebuild it.
        """
        from .tuning import walk_umap
        # One walk at a time. Two running together interleave their configurations into one table and
        # one gallery, which reads as a single sweep and invites a comparison between settings that
        # were never compared. The runner allows concurrent jobs deliberately -- that is what stops a
        # search blocking every other tab -- so the constraint belongs here, on the one job whose
        # output accumulates in a shared place.
        if any(j.active and j.name == WALK_JOB for j in getattr(self.runner, "jobs", {}).values()):
            self.status.emit("a walk is already running -- stop it in Jobs first")
            return
        spec, n, size, seed = self.spec(), self.nodes, self.sample.value(), self.seed.value()
        grid = {k: v for k, v in self.walk_grid().items() if k != 'min_cluster_sizes'}
        # Cleared before the signal, not after: the gallery and anything else listening should see a
        # panel that has already forgotten the last walk, rather than one still holding its rows.
        self._start_table(self.walk_table, [])
        self.walk_started.emit()
        self._run(lambda p: walk_umap(n, spec, sample_size=size, seed=seed, log=p,
                                      store=self.store, on_step=self.walk_step.emit, **grid),
                  lambda d: (self._fill(self.walk_table, d), self.status.emit("walk complete")),
                  name=WALK_JOB)

    def show_walk_row(self, row: int, _col: int = 0):
        """Build and display the configuration on one row of the walk table, and cluster it.

        The walk scores configurations and returns numbers; this is what turns a number back into
        something you can look at. It rebuilds rather than caching all of them, because a walk of
        288 embeddings over 8,140 genes is gigabytes and one rebuild is seconds -- and it goes
        through the same path as "build this map", so the result behaves identically.

        The clustering comes with it when the row reports one, at the same `min_cluster_size` the
        walk used: the row says "11 clusters", and a map shown without them leaves the reader to
        take that number on trust.
        """
        from .tuning import WALK_MIN_CLUSTER_SIZE
        values = self.row_values(self.walk_table, row)
        try:
            if "n_neighbors" in values:
                self.nn.setValue(int(float(values["n_neighbors"])))
            if "min_dist" in values:
                self.md.setValue(float(values["min_dist"]))
        except ValueError:
            self.status.emit("that row does not name a configuration this can rebuild")
            return
        cluster_it = "n_clusters_hdbscan" in values
        self.status.emit(f"building n_neighbors={self.nn.value()}, min_dist={self.md.value():g}"
                         + (", with the clustering the walk scored" if cluster_it else ""))
        self.run_embed(then_cluster=WALK_MIN_CLUSTER_SIZE if cluster_it else None)

    def show_cluster_row(self, row: int):
        """Re-cluster the current map with the settings on one row of the clustering walk.

        The walk scores clusterings of the map that is already on screen, so there is nothing to
        rebuild -- the row is a set of parameters, and this applies them and colours the map by the
        result. Without it the tab reported silhouettes for clusterings nobody could see.
        """
        from .clustering import cluster
        if self.coords is None:
            self.status.emit("build a map first")
            return
        v = self.row_values(self.cluster_table, row)
        algo = v.get("algorithm", self.algo.currentText())
        try:
            mcs = int(float(v.get("min_cluster_size", self.mcs.value())))
            eps = float(v.get("eps", self.eps.value()))
            ms = v.get("min_samples", "")
            ms = mcs if ms in ("", "None", "nan") else int(float(ms))
        except ValueError:
            self.status.emit("that row does not name a clustering this can rebuild")
            return
        # Put the controls where the row says, so the settings on screen describe the map on screen.
        self.algo.setCurrentText(algo)
        self.mcs.setValue(mcs)
        self.eps.setValue(eps)
        Y = self.coords
        self.status.emit(f"clustering: {algo}, min_cluster_size={mcs}"
                         + (f", eps={eps:g}" if algo == "dbscan" else ""))
        self._run(lambda p: cluster(Y, algorithm=algo, min_cluster_size=mcs, min_samples=ms,
                                    eps=eps),
                  self._clustered, name=f"cluster ({algo})")

    def show_search_row(self, row: int):
        """Rebuild the exact configuration on one row of the search table, and show it clustered.

        This is the table where a row is a whole recipe -- blocks, missing-value policy, scaling,
        both UMAP hyperparameters, the clustering size, the seed and the subsample -- and until now
        it was the one table whose rows could not be looked at. A recovery score with no way to see
        the structure it scored is exactly the "trust me" this project refuses elsewhere.

        Rebuilt on the SAME subsample, from the seed and sample size the row records, and with the
        same columns excluded. Rebuilding at full size, or over a different draw, would put a
        different map on screen from the one the row's numbers describe.
        """
        from .search import rebuild
        v = self.row_values(self.search_table, row)
        if not v.get("blocks"):
            self.status.emit("that row does not name a configuration this can rebuild")
            return
        nodes, blocks = self.nodes, v["blocks"]
        self.status.emit(f"rebuilding {blocks} nn={v.get('n_neighbors', '?')} "
                         f"md={v.get('min_dist', '?')} mcs={v.get('min_cluster_size', '?')} "
                         f"-- the same map, on the same genes, with its clustering")

        def job(p):
            coords, genes, labels, features = rebuild(nodes, v, log=p)
            return coords, features, genes, labels

        self._run(job, self._search_row_built, name=f"rebuild search row ({blocks})",
                  on_error=self._row_rebuild_failed)

    def _row_rebuild_failed(self, error) -> bool:
        """A row that cannot be rebuilt says why, rather than failing in red."""
        if not isinstance(error, ValueError):
            return False
        self.status.emit(f"cannot rebuild that row: {error}")
        return True

    def _search_row_built(self, result):
        from .clustering import NOISE
        coords, features, genes, labels = result
        self.coords, self.features, self.rows, self.labels = coords, features, genes, labels
        self.embedding_ready.emit(coords, genes)
        self._publish_clusters(labels, genes)
        k = len(set(labels[labels != NOISE]))
        self.status.emit(f"showing that configuration: {len(coords):,} genes, {k} clusters, "
                         f"{100 * (labels == NOISE).mean():.0f}% unclustered")

    def show_inference_row(self, row: int):
        """Colour the map by the clustering a battery row is about, and say which cluster it names.

        An Inference row is not a configuration -- it is a feature of the map already on screen --
        so the map does not change. What clicking it does is put the clustering the row was scored
        against back on the map, because reading "cluster 3 is 90% apicoplast" while looking at a
        map coloured by compartment is a needless act of translation.
        """
        if self.labels is None:
            self.status.emit("cluster a map first -- these rows describe a clustering")
            return
        v = self.row_values(self.category_table, row)
        self._publish_clusters(self.labels)
        cl, cat = v.get("cluster", "?"), v.get("category", v.get("feature", "that value"))
        lift = v.get("lift")
        self.status.emit(
            f"cluster {cl} is the best match for {cat}"
            + (f" at {float(lift):.1f}x its prevalence" if lift not in (None, "", "nan") else "")
            + " -- the map is coloured by that clustering")

    def run_embed(self, then_cluster=None):
        """Build one embedding from the current spec and show it in the 3D view.

        `then_cluster` clusters it in the same job at that `min_cluster_size`, which is how a walk
        row arrives with the clustering its score counted. One job rather than two, because the two
        belong together: a map that appears for a moment without the clusters the row promised
        reads as the clustering having failed.
        """
        from .clustering import cluster
        from .embedding import embed
        spec, n = self.spec(), self.nodes

        def job(p):
            coords, features, rows = embed(n, spec, log=p)
            labels = None
            if then_cluster:
                p(f"clustering at min_cluster_size={then_cluster}")
                labels = cluster(coords, algorithm="hdbscan", min_cluster_size=int(then_cluster))
            return coords, features, rows, labels

        self._run(job, self._embedded, name="build map")

    def _embedded(self, result):
        from .clustering import NOISE
        self.coords, self.features, self.rows, labels = result
        self.embedding_ready.emit(self.coords, self.rows)
        note = ""
        if labels is not None:
            self.labels = labels
            self._publish_clusters(labels)
            k = len(set(labels[labels != NOISE]))
            note = f", {k} clusters"
        self.status.emit(f"map built: {len(self.coords):,} genes, "
                         f"{len(self.features)} features{note}")

    def save_embedding(self):
        """Store the current embedding with its full recipe, so it can be rebuilt exactly."""
        if self.coords is None or self.store is None:
            self.status.emit("build a map first")
            return
        name = self.emb_name.text().strip() or "unnamed"
        self.store.save(name, self.coords, self.spec(),
                        gene_ids=self.nodes.gene_id[self.rows], features=self.features)
        self.status.emit(f"saved embedding {name!r} with its full recipe")

    def run_cluster_walk(self):
        """Score a grid of clustering hyperparameters against the current map."""
        from .clustering import walk_dbscan, walk_hdbscan
        if self.coords is None:
            self.status.emit("build a map first"); return
        Y, algo = self.coords, self.algo.currentText()
        fn = walk_hdbscan if algo == "hdbscan" else walk_dbscan
        self._run(lambda p: fn(Y, log=p),
                  lambda d: (self._fill(self.cluster_table, d), self.status.emit("walk complete")),
                  name=f"{algo} hyperparameter walk")

    def run_cluster(self):
        """Cluster the current map with the chosen algorithm and settings."""
        from .clustering import NOISE, cluster
        if self.coords is None:
            self.status.emit("build a map first"); return
        Y, algo, mcs, eps = self.coords, self.algo.currentText(), self.mcs.value(), self.eps.value()
        self._run(lambda p: cluster(Y, algorithm=algo, min_cluster_size=mcs,
                                    min_samples=mcs, eps=eps),
                  self._clustered, name=f"cluster ({algo})")

    def _clustered(self, labels):
        from .clustering import NOISE
        self.labels = labels
        k = len(set(labels[labels != NOISE]))
        # Expanded to the node table before it leaves: a clustering of a map built over a subsample
        # covers only that subsample, and the window colours all 8,140 points by it.
        self._publish_clusters(labels)
        self.status.emit(f"{k} clusters, {100 * (labels == NOISE).mean():.0f}% unassigned "
                         f"-- colour the map by 'clusters' to see them")

    def run_battery(self):
        """Test what the clusters correspond to, using only features the map never saw."""
        from .clustering import battery, describe
        if self.labels is None:
            self.status.emit("cluster the map first"); return
        used = [c for cols in columns_for(self.nodes, self.spec()).values() for c in cols]
        if self.cat_cb.isChecked():
            used.append("compartment")
        sub, lab = self.nodes.loc[self.rows], self.labels

        def job(p):
            S, D = battery(sub, lab, used_features=used, log=p)
            return S, D, describe(S, D)

        self._run(job, self._battery_done, name="held-out battery")

    def _battery_done(self, result):
        from .clustering import per_category
        S, D, lines = result
        held = S[S.evidence == "held_out"] if not S.empty else S
        # `+` binds tighter than `or`, so the header made the whole expression truthy and the
        # fallback was unreachable: a clustering that organised nothing showed the explanation alone,
        # which reads as a result that has not loaded rather than as an honest empty answer.
        body = ("".join(f"<p>• {ln}</p>" for ln in lines)
                or "<p>nothing separated the clusters</p>")
        self.findings.setHtml(
            "<p style='color:#888'>Held-out features only. <b>assoc</b> is how much each feature "
            "restates something the map already used — a high value means the finding is close to "
            "circular even when it is technically held out.</p>" + body)
        cols = [c for c in ("feature", "evidence", "score_type", "score", "assoc_with_input", "n", "q")
                if c in held.columns]
        self._fill(self.battery_table, held[cols])
        per = per_category(D)
        self._fill(self.category_table, per)
        recovered = int((per.f1 >= 0.5).sum()) if len(per) else 0
        enriched = int((per.lift >= 2).sum()) if len(per) else 0
        self.status.emit(
            f"battery: {len(held)} held-out features tested; of {len(per)} categories big enough to "
            f"score, {enriched} sit in a cluster at twice their own prevalence and {recovered} "
            f"reach F1 0.5")

    def run_search(self):
        """Walk dataset combinations, scoring each by how well it recovers the held-out label."""
        from .search import search
        import itertools
        n, target = self.nodes, self.target.currentText()
        size, seed = self.search_sample.value(), self.seed.value()
        obj, grid = self.objective_settings(), self.walk_grid()
        base = [b for b in BLOCKS if columns_for(n, EmbeddingSpec(blocks=(b,))).get(b)]
        r = self.max_blocks.value()
        sets = [tuple(c) for k in range(1, r + 1) for c in itertools.combinations(base, k)]

        def job(p):
            R, P = search(n, target=target, block_sets=sets, sample_size=size, seed=seed,
                          store=self.store, objective=obj, log=p, **grid)
            return R

        self._run(job, lambda R: (self._fill(self.search_table, R),
                                  self.status.emit(f"search complete: {len(R)} runs scored")),
                  name=f"recovery search ({target})")
