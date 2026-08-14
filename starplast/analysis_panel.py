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



from .theme import tip as wrap_tip  # noqa: E402  -- one implementation, in the theme


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
    "nn": "Below 2 there is no neighborhood to embed from. The upper bound is well past "
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
              "intend to hold localization out and test whether the map recovers it: a map built on "
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
    "nn": "UMAP n_neighbors: how much of the neighborhood each point is placed by. Small values "
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
    "eps": "DBSCAN neighborhood radius, in embedding units. Ignored by HDBSCAN, which infers the "
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
    "val_hold": "Fraction of each category's labeled genes hidden per fold. These are the genes "
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
from .logging_util import get_logger  # noqa: E402

_log = get_logger(__name__)

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

    def stopped(self) -> bool:
        """Whether the user has asked this job to stop, without raising.

        The raising form above lands at the next log line, and the search logs every fortieth run --
        four minutes of a full-proteome sweep, which reads as a button that does nothing. Passed to
        `search` and `walk_umap` as `should_stop`, they check it every configuration and return what
        they have.
        """
        return bool(self.job.cancelled)

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
    #: A clustering, so the map can color by it. Without this the Clusters tab computed labels,
    #: printed how many there were, and threw them away -- which made the one thing this application
    #: is for, looking at structure colored by a held-out variable, impossible to actually do.
    clusters_ready = QtCore.pyqtSignal(object)                   # labels, -1 for noise
    #: One finished configuration of a walk, as `tuning.WalkStep`, emitted from the worker thread the
    #: moment it is computed. Qt queues it to the GUI thread, which is what lets a row and a thumbnail
    #: appear while the walk is still running rather than all at once when it ends.
    walk_step = QtCore.pyqtSignal(object)
    #: A walk is starting: whatever the last one left on screen belongs to a different sweep.
    walk_started = QtCore.pyqtSignal()
    #: An annotation was saved or withdrawn, so the map's fourth color has changed.
    annotations_changed = QtCore.pyqtSignal()
    #: One scored configuration of a search, as `search.RunStep`, emitted from the worker thread as
    #: each embedding finishes. Carries its clustering, because here the clustering is half of what
    #: was scored and a map shown without it is a map shown without the result.
    search_step = QtCore.pyqtSignal(object)
    #: One configuration of a discovery climb, as a plain row. A row rather than a step object
    #: because the climb's artefacts -- the labels, the findings -- are held for the reading rather
    #: than drawn, and shipping them through a GUI signal every step would copy them for nothing.
    discovery_step = QtCore.pyqtSignal(object)
    status = QtCore.pyqtSignal(str)

    def __init__(self, nodes: pd.DataFrame, store=None, parent=None, runner=None,
                 annotations=None):
        super().__init__(parent)
        self.nodes = nodes
        self.store = store
        #: Where annotations are written. A separate file from everything else, always: the node
        #: table is measurement, and an inference stored beside it becomes indistinguishable from
        #: one the moment anybody reads the table without knowing which columns are which.
        self.store_annotations = annotations
        self.labels = None
        self.coords = None
        self.rows = None
        self._thread = None
        self._worker = None
        # The window's JobRunner, so this panel's work is visible and stoppable alongside everything
        # else. None is allowed: the panel is constructed without one in tests.
        self.runner = runner
        self._jobs = {}
        #: Buttons that stop whatever this panel is running, and the jobs they stop. A stop that
        #: lives only in the Jobs dock's right-click menu is a stop nobody finds: the first thing
        #: asked about it was "I can't see the stop button".
        self._stop_btns, self._running = [], []
        #: Rows are written to disk as they arrive, one file per table per run -- see `_row_log`.
        self._row_logs = {}
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
        self.search_step.connect(self._search_step_arrived)
        self.discovery_step.connect(self._discovery_step_arrived)

        tabs = QtWidgets.QTabWidget()
        tabs.addTab(self._data_tab(), "1 · Data")
        tabs.addTab(self._map_tab(), "2 · Map")
        tabs.addTab(self._cluster_tab(), "3 · Clusters")
        tabs.addTab(self._meaning_tab(), "4 · Inference")
        tabs.addTab(self._search_tab(), "5 · Search")
        tabs.addTab(self._validation_tab(), "6 · Validation")
        tabs.addTab(self._discover_tab(), "7 · Discover")
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(tabs)
        # `_apply_tooltips` already wraps every one of these through `theme.tip` and copies it onto
        # the form label. Walking the children again from here crashed the interpreter: this runs
        # inside __init__, and findChildren hands back wrappers around C++ objects that widgets
        # still under construction destroy as they go.
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
        # Imported columns get their own block, added when something is imported: a block that is
        # always there and always empty is a control that does nothing.
        self.imported_cb = QtWidgets.QCheckBox("imported  (0 columns)")
        self.imported_cb.setToolTip(
            "Columns from a table you imported yourself, with the preprocessing you chose recorded "
            "beside them. Ticked, they feed the map like any other block -- and like any other "
            "block, anything fed in here cannot afterwards be evidence about the clusters.")
        self.imported_cb.setEnabled(False)
        self.imported_columns = []
        bl.addWidget(self.imported_cb)
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

    def add_imported(self, columns):
        """Take newly imported columns, so the Data tab can offer them as a block."""
        for c in columns:
            if c not in self.imported_columns:
                self.imported_columns.append(c)
        self.imported_cb.setText(f"imported  ({len(self.imported_columns)} columns)")
        self.imported_cb.setEnabled(bool(self.imported_columns))
        self.imported_cb.setChecked(bool(self.imported_columns))
        return self.imported_columns

    def spec(self) -> EmbeddingSpec:
        """The EmbeddingSpec described by the current controls."""
        return EmbeddingSpec(
            extra_columns=(tuple(self.imported_columns)
                           if self.imported_cb.isChecked() else ()),
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
        row.addWidget(b1); row.addWidget(b2); row.addWidget(self.stop_button())
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
        row.addWidget(b1); row.addWidget(b2); row.addWidget(self.stop_button())
        v.addLayout(row)
        self.cluster_table = self.results_table(
            QtWidgets.QTableWidget(), self.show_cluster_row, "clustering_walk")
        self.cluster_table.setToolTip(
            "Click a row to cluster the current map with those settings and color it by the "
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
    def _discover_tab(self):
        """Search the space of maps for the two claims a map can make, and read the winner back.

        A separate tab from Search because it asks the opposite question. Search asks whether a map
        can be trusted -- hold out a label, see if it comes back. This asks whether a map is USEFUL:
        how much does it say about genes nobody has measured, and where does one layer split a
        category another layer calls uniform. A configuration can be excellent at one and useless at
        the other, and putting them in one tab under one "optimize for" box would hide that.
        """
        from . import discovery, optimize
        w = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(w)
        note = QtWidgets.QLabel(
            "Hill-climb the space of embeddings and clusterings looking for structure that PREDICTS "
            "rather than structure that agrees. Guilt by association finds clusters whose labelled "
            "members agree and whose unlabelled members inherit the claim; layer disagreement finds "
            "clusters that agree about one measurement and split on another. Then press "
            "\u201cread the results\u201d and the run is written back as ranked claims with their "
            "numbers and their caveats.")
        note.setWordWrap(True)
        v.addWidget(note)

        form = QtWidgets.QFormLayout()
        self.discover_mode = QtWidgets.QComboBox()
        self.discover_mode.addItems(list(optimize.MODES))
        self.discover_mode.setToolTip(
            "What the climb maximises.\n\n"
            "guilt -- how much the map predicts about genes nobody has measured.\n"
            "disagreement -- how many categories one layer splits along another.\n"
            "both -- the sum, for \u201cfind me anything\u201d.\n"
            "recovery -- the old objective: how well a held-out label comes back. Measures trust "
            "rather than yield, and a map tuned until it recovers a label perfectly has often "
            "found nothing new.\n"
            "auprc / auroc -- how good a shortlist the map gives per category, as a ranking rather "
            "than a partition. Unlike the partition scores these cannot be won by merging "
            "everything into one cluster. AUPRC is reported as lift over prevalence, because a raw "
            "AUPRC is uninterpretable without knowing how common the class is.\n"
            "knn -- scores the EMBEDDING alone, with no clustering in it. Optimise this first when "
            "the clustering is the thing that keeps going wrong.")

        columns = list(self.nodes.columns) if self.nodes is not None else []
        def _layers(numeric):
            out = []
            for c in columns:
                s = self.nodes[c]
                is_num = pd.api.types.is_numeric_dtype(s)
                if is_num == numeric and s.notna().sum() > 100 and (numeric or s.nunique() <= 40):
                    out.append(c)
            return out

        self.discover_layer = QtWidgets.QComboBox()
        self.discover_layer.addItems(_layers(False) + _layers(True))
        for i, name in enumerate(("compartment_best", "compartment")):
            if self.discover_layer.findText(name) >= 0:
                self.discover_layer.setCurrentText(name)
                break
        self.discover_layer.setToolTip(
            "The layer a claim is about. A category (localisation, cell-cycle phase) produces "
            "enrichment claims about the genes in a cluster that carry no label; a quantity "
            "(fitness, abundance) produces claims about the genes nobody measured.")

        self.discover_against = QtWidgets.QComboBox()
        self.discover_against.addItems(_layers(True) + _layers(False))
        # Fitness first where there is one: "same compartment, opposite fitness" is the disagreement
        # people actually come here for, and the alphabetical default was protein length.
        for name in ("fit_invitro_hff", "fitness", "cellcycle_phase"):
            if self.discover_against.findText(name) >= 0:
                self.discover_against.setCurrentText(name)
                break
        self.discover_against.setToolTip(
            "The second layer, for disagreement: the one a cluster is allowed to disagree about "
            "while agreeing about the first. Same compartment, opposite fitness.")

        self.discover_budget = QtWidgets.QSpinBox()
        self.discover_budget.setRange(4, 400)
        self.discover_budget.setValue(40)
        self.discover_budget.setToolTip(
            "How many configurations the climb may evaluate. Each one is an embedding and a "
            "clustering, so this is the run\u2019s cost in minutes as much as its thoroughness. "
            "Embeddings are cached across steps that change only the clustering, which is most of "
            "them, so the true cost is far below the count.")

        self.discover_restarts = QtWidgets.QSpinBox()
        self.discover_restarts.setRange(1, 10)
        self.discover_restarts.setValue(2)
        self.discover_restarts.setToolTip(
            "A hill climber finds the top of whatever hill it started on. Restarts are the cheapest "
            "defence against reporting a local optimum as the answer.")

        form.addRow("optimize for", self.discover_mode)
        form.addRow("layer", self.discover_layer)
        form.addRow("disagreeing with", self.discover_against)
        form.addRow("configurations to try", self.discover_budget)
        form.addRow("restarts", self.discover_restarts)
        v.addLayout(form)

        run = QtWidgets.QHBoxLayout()
        go = QtWidgets.QPushButton("search for structure")
        go.setProperty("primary", True)
        go.setToolTip(
            "Start the climb. It steps one coordinate at a time -- one hyperparameter, or one "
            "dataset added or removed -- keeps the step when the score improves, and restarts "
            "elsewhere when it can no longer improve.\n\n"
            "Every configuration appears in the table as it finishes, with its score and every "
            "other metric alongside, so a run optimised for one thing can be re-read against "
            "another afterwards. Stopping keeps everything already evaluated.")
        go.clicked.connect(self.run_discovery)
        self.read_button = QtWidgets.QPushButton("read the results")
        self.read_button.clicked.connect(self.read_discovery)
        self.read_button.setEnabled(False)
        self.read_button.setToolTip(
            "Rank what the winning configuration found and write it out as claims: which clusters "
            "matter, what they say, which genes they are about, and what is wrong with each one. "
            "Ranked by strength x reach x novelty -- a statistically overwhelming claim about two "
            "well-published genes is not the finding to read first.")
        run.addWidget(go)
        run.addWidget(self.read_button)
        run.addWidget(self.stop_button())
        v.addLayout(run)

        self.discover_table = self.results_table(
            QtWidgets.QTableWidget(), self.show_discovery_row, "discovery")
        v.addWidget(self.discover_table, 1)
        self.discover_report = QtWidgets.QTextBrowser()
        self.discover_report.setOpenExternalLinks(True)
        self.discover_report.setMinimumHeight(180)
        v.addWidget(self.discover_report, 1)
        self._discovered = None
        return w

    def run_discovery(self):
        """Climb, streaming each configuration into the table as it finishes."""
        from . import optimize
        if self.nodes is None:
            return
        mode = self.discover_mode.currentText()
        layer, against = self.discover_layer.currentText(), self.discover_against.currentText()
        budget, restarts = self.discover_budget.value(), self.discover_restarts.value()
        pool = optimize.block_pool(self.nodes)
        seed = self.seed.value()
        start = {**{k: v[len(v) // 2] for k, v in optimize.SPACE.items()},
                 "method": "umap", "algorithm": "hdbscan", "blocks": tuple(pool[:3]) or ("",)}
        self._start_table(self.discover_table, [])
        self.discover_report.setPlainText("")

        def job(p):
            evaluate = optimize.evaluator(self.nodes, mode=mode, layers=(layer,),
                                          against=(against,), seed=seed, log=p)
            return optimize.climb(evaluate, start, block_pool=pool, restarts=restarts,
                                  max_evaluations=budget, seed=seed,
                                  on_step=lambda row, cfg, extras: self.discovery_step.emit(row),
                                  should_stop=getattr(p, "stopped", None), log=p)

        self._run(job, self._discovery_done, name=f"discovery ({mode}, {layer})")

    def _discovery_step_arrived(self, row):
        """One configuration finished. Its private artefacts are not shown -- they are held for the
        reading -- so the table stays a table."""
        self._append(self.discover_table, {k: v for k, v in row.items() if not k.startswith("_")})

    def _discovery_done(self, result):
        self._discovered = result
        got = result is not None and not result.empty
        self.read_button.setEnabled(bool(got))
        if got:
            best = result.iloc[0]
            self.status.emit(f"best {best.score:.3f} from {len(result)} configurations "
                             f"({int(best.get('n_findings', 0))} findings)")

    def read_discovery(self):
        """Write the winning configuration's findings back as ranked claims."""
        from . import interpret
        if self._discovered is None or self._discovered.empty:
            return
        findings = self._discovered.iloc[0].get("_findings")
        self.discover_report.setMarkdown(interpret.report(findings, self.nodes, top=8))

    def show_discovery_row(self, row: int, table=None):
        """Rebuild the configuration on one row and put it on screen with its clustering.

        Through `search.rebuild`, which is what the recovery table uses: one implementation, so the
        two tables cannot come to disagree about what a row means. A row here carries the newer
        coordinates as well -- the method, the algorithm -- and rebuild reads what it recognises.
        """
        from .search import rebuild
        v = self.row_values(table if table is not None else self.discover_table, row)
        if not v.get("blocks"):
            self.status.emit("that row does not name a configuration this can rebuild")
            return
        nodes = self.nodes
        self.status.emit(f"rebuilding {v['blocks']} {v.get('method', 'umap')} / "
                         f"{v.get('algorithm', 'hdbscan')} -- the map that scored "
                         f"{v.get('score', '?')}")

        def job(p):
            coords, genes, labels, features = rebuild(nodes, v, log=p)
            return coords, features, genes, labels

        self._run(job, self._search_row_built, name=f"rebuild discovery row ({v['blocks']})")

    def _search_tab(self):
        w = QtWidgets.QWidget()
        v = QtWidgets.QVBoxLayout(w)
        note = QtWidgets.QLabel(
            "Walk combinations of datasets and hyperparameters looking for a structure that recovers a "
            "label it was never given. The target and everything that restates it are excluded from "
            "every map, which is what makes the recovery meaningful — and what turns it into a "
            "prediction for the unlabeled genes in a pure cluster.")
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
            "optimize for exactly those. The objective and the label set are independent, so "
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
        form.addRow("optimize for", self.objective)
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
        run = QtWidgets.QHBoxLayout()
        b = QtWidgets.QPushButton("run the search")
        b.setProperty("primary", True)
        b.clicked.connect(self.run_search)
        run.addWidget(b); run.addWidget(self.stop_button())
        v.addLayout(run)
        self.search_table = self.results_table(
            QtWidgets.QTableWidget(), self.show_search_row, "recovery_search")
        self.category_search_table = self.results_table(
            QtWidgets.QTableWidget(), self.show_search_category_row, "recovery_per_category")
        self.category_search_table.setToolTip(
            "One row per configuration and category: how well that map isolates that compartment, "
            "phase or class, rather than how it did on average. A mean hides the case this table "
            "exists for -- a map where the GRAs are clean and everything else is a mess is exactly "
            "what you want when you are looking for GRAs. Click a row to rebuild that map.")
        self.search_table.setToolTip(
            "Click a row to rebuild that exact configuration -- same blocks, same policy, same "
            "seed, same subsample, same excluded columns -- and show it with its clustering. A "
            "recovery score with no way to look at the structure it scored is a 'trust me'. "
            "Right-click to save the whole table as CSV.")
        v.addWidget(self.search_table, 1)
        v.addWidget(QtWidgets.QLabel(
            "<i>Per configuration and category — a mean over categories hides the map that nails "
            "one of them.</i>"))
        v.addWidget(self.category_search_table, 1)
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
        of `cellcycle_phase`, and a stale list would let someone optimize for a label that no longer
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
        self.val_hold.setToolTip("Fraction of each category's labeled genes hidden per fold.")
        self.val_refit = QtWidgets.QCheckBox("re-embed and re-cluster for every fold")
        form.addRow("hold out", self.val_target)
        form.addRow("folds", self.val_folds)
        form.addRow("fraction hidden", self.val_hold)
        form.addRow("stricter", self.val_refit)
        v.addLayout(form)

        run = QtWidgets.QHBoxLayout()
        b = QtWidgets.QPushButton("test the annotation")
        b.setProperty("primary", True)
        b.clicked.connect(self.run_validation)
        run.addWidget(b); run.addWidget(self.stop_button())
        v.addLayout(run)

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
            "The unlabeled members of that category's cluster: the list this whole tab exists to "
            "put a number on. Every row carries how much of its cluster already carries the "
            "category and how much carries something else, and how many genes already labeled it "
            "share an orthogroup, a Pfam or an InterPro domain with the candidate — agreement "
            "from evidence the map never saw. Zero support is the common answer.")
        v.addWidget(self.cand_table, 1)

        save = QtWidgets.QHBoxLayout()
        self.reasoning = QtWidgets.QLineEdit()
        self.reasoning.setPlaceholderText("why this cluster is worth annotating from…")
        self.reasoning.setToolTip(
            "Free text, saved with every annotation. The numbers say how often it would be right; "
            "this says why you thought it was worth proposing, which is the part a collaborator "
            "will want to argue with and the part no column can hold.")
        self.save_annotations = QtWidgets.QPushButton("save these as annotations")
        self.save_annotations.setToolTip(
            "Write the candidates above to the annotations file, each carrying the validated "
            "precision and recall for its category, the cluster's composition and the "
            "configuration. Refused outright when there is no validated precision, or when it is "
            "below 10% -- that is not an annotation, it is a coin toss with a record attached. "
            "Never written into the node table: measurement and inference do not share a file.")
        self.save_annotations.clicked.connect(self.save_candidates)
        save.addWidget(self.reasoning, 1)
        save.addWidget(self.save_annotations)
        v.addLayout(save)
        return w

    def save_candidates(self):
        """Save the candidate list as annotations, with the numbers that justify each one.

        Only from here, and only from a validated run: the refusal is the feature. The candidates
        and their error rate are joined in `annotations.from_candidates` so a precision from a
        different category or a different clustering cannot travel with a row.
        """
        import datetime
        from .annotations import from_candidates
        if self.store_annotations is None:
            self.status.emit("no annotations file is configured")
            return
        cand = self._frames.get(self.cand_table)
        d = getattr(self, "_validation_scores", None)
        row = getattr(self, "_candidate_row", None)
        if cand is None or not len(cand) or d is None or row is None:
            self.status.emit("no candidates -- validate a target, then click a category")
            return
        scores = d.iloc[row]
        ann = from_candidates(cand, getattr(self, "_validation_target", ""), scores,
                              configuration=self._configuration(),
                              reasoning=self.reasoning.text().strip(),
                              date=datetime.date.today().isoformat())
        try:
            out = self.store_annotations.save(ann)
        except ValueError as exc:
            # Shown as an explanation, like the circularity refusal: this is the store working.
            self.val_note.setText(f"<b>Not saved.</b> {exc}")
            self.val_note.show()
            self.status.emit(str(exc))
            return
        self.val_note.hide()
        self.annotations_changed.emit()
        self.status.emit(f"saved {len(ann):,} annotations at precision {scores.precision:.2f} "
                         f"({len(out):,} in the file)")

    def _configuration(self) -> dict:
        """The recipe behind the map the candidates came from, as the store records it."""
        spec = self.spec()
        return {"blocks": "+".join(spec.blocks), "na_policy": spec.na_policy,
                "scaling": spec.scaling, "n_neighbors": spec.n_neighbors,
                "min_dist": spec.min_dist, "min_cluster_size": self.mcs.value(),
                "seed": spec.random_state}

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
            self.status.emit(f"no cluster holds any gene labeled {category!r}")
            return
        genes = self.nodes.gene_id[self.rows] if self.rows is not None else self.nodes.gene_id
        self._candidate_row = row
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
            return "no category had enough labeled genes to hide any"
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
            self._running = [j for j in self._running if j.active] + [job]
            self._set_stoppable(True)
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
        self._running = [j for j in self._running if j.active and j.id != jid]
        if not self._running:
            self._set_stoppable(False)
        saved = self._close_row_logs()
        where = f" -- rows saved to {'; '.join(saved)}" if saved else ""
        if job.state == "cancelled":
            # A stopped search RETURNS what it computed, so there is usually a result to deliver.
            # Discarding it because the user pressed stop would make stopping cost the whole run,
            # which is the opposite of what the button is for.
            kept = job.result
            if kept is not None and on_done is not None:
                try:
                    on_done(kept)
                except Exception as exc:                      # a partial result may be shaped oddly
                    _log.warning("stopped job %d: %s: %s", jid, type(exc).__name__, exc)
            self.status.emit(f"{job.name}: stopped -- what it finished is kept{where}")
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
        if saved:
            self.status.emit(f"{job.name}: done{where}")

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

    def table_kind(self, table) -> str:
        """The key a saved file records for this table -- what makes a load refuse a mismatch."""
        return self._table_what.get(table, "")

    def results_tables(self) -> dict:
        """Every results table by its key. What "all tabs at once" means, in one place."""
        return {self._table_what[t]: t for t in self._table_what}

    def save_results(self, table, path: str = "") -> str:
        """Write one table so it can be loaded back and clicked, not merely read."""
        from .results import save_table
        df = self._frames.get(table)
        kind = self.table_kind(table)
        if df is None or not len(df):
            self.status.emit("nothing to save -- run something first")
            return ""
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(
                self, "Save results", f"{kind}.csv", "starplast results (*.csv)")
        if not path:
            return ""
        save_table(path, kind, df)
        self.status.emit(f"wrote {len(df):,} rows of {kind} to {path} -- load it back and its rows "
                         f"are still clickable")
        return path

    def load_results(self, table, path: str = "") -> bool:
        """Load a saved table back into this one, refusing a file that came from another.

        Refused rather than loaded anyway: the rows of one table are recipes for rebuilding a map and
        the rows of another are per-category scores, so a validation file dropped into the walk tab
        produces rows nobody can rebuild and an error message about the wrong thing.
        """
        from .results import load_table, table_kind
        kind = self.table_kind(table)
        if not path:
            path, _ = QtWidgets.QFileDialog.getOpenFileName(
                self, f"Load {kind}", "", "starplast results (*.csv);;All files (*)")
        if not path:
            return False
        found = table_kind(path)
        if found and found != kind:
            self.status.emit(f"that file holds {found}, not {kind} -- load it into the tab it came "
                             f"from, where its rows mean something")
            return False
        try:
            df = load_table(path)
        except Exception as exc:
            self.status.emit(f"could not read {path}: {type(exc).__name__}: {exc}")
            return False
        self._fill(table, df)
        self.status.emit(f"loaded {len(df):,} rows into {kind}"
                         + ("" if found else " -- the file did not say which table it came from"))
        return True

    def save_all_results(self, path: str = "") -> str:
        """Every tab's results in one file."""
        from .results import save_bundle
        tables = {kind: self._frames.get(t) for kind, t in self.results_tables().items()}
        if not any(df is not None and len(df) for df in tables.values()):
            self.status.emit("no results to save -- run something first")
            return ""
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(
                self, "Save all results", "starplast_results.starplast",
                "starplast results (*.starplast)")
        if not path:
            return ""
        save_bundle(path, tables, meta={"target": self.target.currentText(),
                                        "spec": self.spec().to_dict()})
        n = sum(1 for df in tables.values() if df is not None and len(df))
        self.status.emit(f"wrote {n} table(s) to {path}")
        return path

    def load_all_results(self, path: str = "") -> int:
        """Load a bundle back into every tab it names."""
        from .results import load_bundle
        if not path:
            path, _ = QtWidgets.QFileDialog.getOpenFileName(
                self, "Load results", "", "starplast results (*.starplast);;All files (*)")
        if not path:
            return 0
        try:
            tables, meta = load_bundle(path)
        except Exception as exc:
            self.status.emit(f"could not read {path}: {type(exc).__name__}: {exc}")
            return 0
        known = self.results_tables()
        loaded = 0
        for kind, df in tables.items():
            if kind in known and len(df):
                self._fill(known[kind], df)
                loaded += 1
        # Named rather than counted: a bundle whose tables this version does not have is a file from
        # a later one, and "3 of 5 loaded" with no names is not something anyone can act on.
        unknown = [k for k in tables if k not in known]
        self.status.emit(f"loaded {loaded} table(s) from {path}"
                         + (f"; this version has no tab for {', '.join(unknown)}" if unknown else ""))
        return loaded

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
        m.addSeparator()
        keep = m.addAction("Save these results (reloadable)…")
        keep.setEnabled(df is not None and len(df) > 0)
        keep.triggered.connect(lambda: self.save_results(table))
        back = m.addAction("Load results into this table…")
        back.triggered.connect(lambda: self.load_results(table))
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

        A clustering of a subsample covers only the genes that map covers, and the window colors
        8,140 points by it. Expanded here with -1 -- unclustered, which is drawn gray -- rather than
        left short, because a labels array of the wrong length made the window fall back to
        coloring everything gray, which reads as "this clustering found nothing".
        """
        labels = np.asarray(labels)
        rows = self.rows if genes is None else genes
        if rows is not None and np.size(rows) == len(self.nodes) and len(labels) != len(self.nodes):
            rows = np.asarray(rows).astype(bool)
            if int(rows.sum()) != len(labels):
                # Neither the map's length nor the table's, which means the clustering and the map
                # on screen are not of the same genes. Emitted unchanged so the window's own length
                # check draws it grey rather than putting cluster i's color on gene j -- and said
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

    def stop_button(self) -> QtWidgets.QPushButton:
        """A stop button for whatever this panel is running, disabled until something is."""
        b = QtWidgets.QPushButton("stop")
        b.setEnabled(False)
        b.setToolTip("Stop the running analysis after the configuration it is on. Everything "
                     "already computed is kept: the rows are in the table, each embedding is in "
                     "the store, and both have been written to disk as they were produced.")
        b.clicked.connect(self.stop_running)
        self._stop_btns.append(b)
        return b

    def stop_running(self) -> int:
        """Ask every job this panel started to stop. Returns how many were asked."""
        asked = [j for j in self._running if j is not None and j.active]
        for job in asked:
            job.cancel()
        if asked:
            self.status.emit(f"stopping after the current configuration -- what is already "
                             f"computed is kept and saved")
        else:
            self.status.emit("nothing is running")
        return len(asked)

    def _set_stoppable(self, on: bool) -> None:
        for b in self._stop_btns:
            b.setEnabled(on)

    def _row_log(self, table: QtWidgets.QTableWidget):
        """The file this table's rows are being written to as they arrive, opening one if needed."""
        from .results import RowLog, autosave_path
        log = self._row_logs.get(table)
        if log is None:
            import datetime
            from . import paths
            kind = self._table_what.get(table, "results")
            stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            log = RowLog(autosave_path(paths.user_cache_dir(), kind, stamp), kind)
            self._row_logs[table] = log
        return log

    def _close_row_logs(self) -> list:
        """Close every open autosave file and return the paths that got rows."""
        out = []
        for table, log in list(self._row_logs.items()):
            path = log.close()
            if path:
                out.append(path)
            self._row_logs.pop(table, None)
        return out

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
        # And onto disk, immediately. In memory is enough to survive a stop; it is not enough to
        # survive quitting, a crash or a power cut, and an hour of search is too much to hold in
        # a process nobody promised to keep alive.
        try:
            self._row_log(table).append(row)
        except OSError as exc:
            _log.warning("autosave: %s: %s", type(exc).__name__, exc)
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
        # `getattr`, because `log` is documented as any callable: the panel hands in a _Progress,
        # which carries the stop flag, and a caller with a plain function still gets a walk.
        self._run(lambda p: walk_umap(n, spec, sample_size=size, seed=seed, log=p,
                                      store=self.store, on_step=self.walk_step.emit,
                                      should_stop=getattr(p, "stopped", None), **grid),
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
        rebuild -- the row is a set of parameters, and this applies them and colors the map by the
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

    def show_search_row(self, row: int, table=None):
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
        # The table is a parameter because the per-category rows carry the same recipe and rebuild
        # the same way: one implementation, so the two cannot come to disagree about what a row means.
        v = self.row_values(table if table is not None else self.search_table, row)
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
        note = ""
        want = getattr(self, "_category_of_interest", None)
        if want and want[0]:
            note = f"  ·  {want[0]} was matched by cluster {want[1]}"
        self._category_of_interest = None
        self.status.emit(f"showing that configuration: {len(coords):,} genes, {k} clusters, "
                         f"{100 * (labels == NOISE).mean():.0f}% unclustered{note}")

    def show_inference_row(self, row: int):
        """Color the map by the clustering a battery row is about, and say which cluster it names.

        An Inference row is not a configuration -- it is a feature of the map already on screen --
        so the map does not change. What clicking it does is put the clustering the row was scored
        against back on the map, because reading "cluster 3 is 90% apicoplast" while looking at a
        map colored by compartment is a needless act of translation.
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
            + " -- the map is colored by that clustering")

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
        # covers only that subsample, and the window colors all 8,140 points by it.
        self._publish_clusters(labels)
        self.status.emit(f"{k} clusters, {100 * (labels == NOISE).mean():.0f}% unassigned "
                         f"-- color the map by 'clusters' to see them")

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
        """Walk dataset combinations, scoring each by how well it recovers the held-out label.

        Every configuration is embedded, clustered at each `min_cluster_size` in the grid, and
        scored per category -- the automated walk. Results stream: the configuration table and the
        per-category table fill as they are computed, and each embedding's best clustering appears
        in the gallery, colored by that clustering rather than by the map's color mode. A search is
        hundreds of runs and a table that arrives at the end is a table nobody watches.
        """
        from .search import search
        import itertools
        n, target = self.nodes, self.target.currentText()
        size, seed = self.search_sample.value(), self.seed.value()
        obj, grid = self.objective_settings(), self.walk_grid()
        base = [b for b in BLOCKS if columns_for(n, EmbeddingSpec(blocks=(b,))).get(b)]
        r = self.max_blocks.value()
        sets = [tuple(c) for k in range(1, r + 1) for c in itertools.combinations(base, k)]
        for table in (self.search_table, self.category_search_table):
            self._start_table(table, [])
        self.walk_started.emit()

        def job(p):
            R, P = search(n, target=target, block_sets=sets, sample_size=size, seed=seed,
                          store=self.store, objective=obj, on_run=self.search_step.emit,
                          should_stop=getattr(p, "stopped", None), log=p, **grid)
            return R, P

        self._run(job, self._search_done, name=f"recovery search ({target})")

    def _search_step_arrived(self, step):
        """One configuration finished: its rows into both tables, its map into the gallery.

        Runs on the GUI thread -- see the connection in `__init__`. The per-category rows go in as
        they arrive rather than being collected, so the question this walk exists to answer ("is
        there a map where THIS category comes out clean") can be asked while it is still running.
        """
        self._append(self.search_table, step.row)
        for _, row in step.per.iterrows():
            self._append(self.category_search_table, self._category_row(step, row))
        self.status.emit(f"search {step.index} of {step.total}: {step.label}")

    @staticmethod
    def _category_row(step, row) -> dict:
        """One (configuration x category) row, with enough of the recipe to rebuild it."""
        out = {"category": row.get("label"), "precision": row.get("precision"),
               "recall": row.get("recall"), "f1": row.get("f1"),
               "n_label": row.get("n_label"), "n_in_cluster": row.get("n_in_cluster"),
               "cluster": row.get("cluster")}
        # The recipe travels on every row: a per-category score whose configuration is only in
        # another table cannot be rebuilt from what is on screen, which is what clicking it needs.
        for k in ("blocks", "na_policy", "scaling", "n_neighbors", "min_dist", "min_cluster_size",
                  "seed", "sample_size", "excluded", "mean_f1", "best_f1"):
            if k in step.row:
                out[k] = step.row[k]
        return out

    def _search_done(self, result):
        """Replace the streamed rows with the ranked tables, and say what was found."""
        R, P = result
        self._fill(self.search_table, R)
        if not R.empty:
            per = P.copy()
            if not per.empty:
                per = per.rename(columns={"label": "category"}).sort_values("f1", ascending=False)
                self._fill(self.category_search_table, per)
            best = R.iloc[0]
            n_front = int(R.on_frontier.sum()) if "on_frontier" in R.columns else 0
            # Read defensively: this is the display path, and a results frame that is missing a
            # column should cost the sentence, not the whole run's output.
            score = f"{best['mean_f1']:.3f}" if "mean_f1" in R.columns else "-"
            blocks = best["blocks"] if "blocks" in R.columns else "?"
            self.status.emit(
                f"search complete: {len(R)} runs scored; best mean F1 {score} ({blocks}), "
                f"{n_front} on the frontier of mean and best F1")
        else:
            self.status.emit("search complete: nothing was scorable")

    def show_search_category_row(self, row: int):
        """Rebuild the configuration behind one per-category row, and say which cluster to look at.

        The same rebuild as the configuration table -- the recipe is on the row -- so a category
        score can be taken from a number to a map in one click. The cluster it names is the one the
        score is about, which is not obvious from a map colored by 40 clusters.
        """
        v = self.row_values(self.category_search_table, row)
        if not v.get("blocks"):
            self.status.emit("that row does not name a configuration this can rebuild")
            return
        self._category_of_interest = (v.get("category"), v.get("cluster"))
        self.show_search_row(row, table=self.category_search_table)
