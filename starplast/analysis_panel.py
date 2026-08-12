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


#: Why each control exists, not what it is called. Applied after the tabs are built so every widget
#: is covered in one auditable place rather than scattered through five builders -- and so a new
#: control without an explanation is a visible omission rather than a silent one.
TOOLTIPS = {
    # 1 Data
    "cat_cb": "One-hot the measured hyperLOPIT compartment INTO the map. Leave it OFF when you "
              "intend to hold localisation out and test whether the map recovers it: a map built on "
              "a label separates that label by construction, and the result means nothing.",
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
}

#: Buttons, keyed by the method they call.
BUTTON_TOOLTIPS = {
    "show_variance": "Report how much of the feature matrix each block actually contributes. The "
                     "check that catches a block being named as an input while carrying almost "
                     "nothing -- hyperLOPIT came to 1.1%.",
    "run_umap_walk": "Score a grid of UMAP hyperparameters and rank them. Produces a table, not "
                     "maps: seeing the embeddings themselves is the gallery, which is not built yet.",
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
        if runner is not None:
            runner.finished.connect(self._on_job_finished)

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
        """Give every control an explanation, from the tables at the top of this module.

        Applied in one place rather than at each widget's construction so that the set of
        explanations is auditable: a control with no entry is a visible omission, and a test asserts
        there are none. The convention throughout is that a tooltip says WHY the control exists and
        what choosing wrongly costs -- "n_neighbors" is a label, not an explanation.
        """
        for attr, tip in TOOLTIPS.items():
            w = getattr(self, attr, None)
            if w is not None and not w.toolTip():
                w.setToolTip(tip)
        # Buttons are keyed by the handler they call, since they are built inline and not stored.
        for btn in self.findChildren(QtWidgets.QPushButton):
            tip = BUTTON_TOOLTIPS.get(self._handler_name(btn))
            if tip and not btn.toolTip():
                btn.setToolTip(tip)
        # The feature blocks explain themselves from the registry rather than from a hand-written
        # list that would drift as blocks are added.
        for name, cb in getattr(self, "block_cb", {}).items():
            if not cb.toolTip():
                cols = columns_for(self.nodes, EmbeddingSpec(blocks=(name,))).get(name, [])
                cb.setToolTip(
                    f"Feed the {name} block into the map: {len(cols)} columns"
                    + (f", including {', '.join(cols[:4])}" if cols else "")
                    + ". Anything fed in here cannot afterwards be used as evidence about the "
                      "clusters, because a feature the map was built on separates them by "
                      "construction.")

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
        form.addRow("n_neighbors", self.nn)
        form.addRow("min_dist", self.md)
        form.addRow("seed", self.seed)
        form.addRow("walk sample size", self.sample)
        v.addLayout(form)

        row = QtWidgets.QHBoxLayout()
        b1 = QtWidgets.QPushButton("walk hyperparameters")
        b1.clicked.connect(self.run_umap_walk)
        b2 = QtWidgets.QPushButton("build this map")
        b2.setProperty("primary", True)
        b2.clicked.connect(self.run_embed)
        row.addWidget(b1); row.addWidget(b2)
        v.addLayout(row)

        self.walk_table = QtWidgets.QTableWidget()
        self.walk_table.setAlternatingRowColors(True)
        v.addWidget(self.walk_table, 1)

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
        self.cluster_table = QtWidgets.QTableWidget()
        self.cluster_table.setAlternatingRowColors(True)
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
        self.battery_table = QtWidgets.QTableWidget()
        self.battery_table.setAlternatingRowColors(True)
        v.addWidget(self.battery_table, 1)
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
        form.addRow("target to recover", self.target)
        form.addRow("sample size", self.search_sample)
        form.addRow("max blocks per combination", self.max_blocks)
        v.addLayout(form)
        b = QtWidgets.QPushButton("run the search")
        b.setProperty("primary", True)
        b.clicked.connect(self.run_search)
        v.addWidget(b)
        self.search_table = QtWidgets.QTableWidget()
        self.search_table.setAlternatingRowColors(True)
        v.addWidget(self.search_table, 1)
        return w

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
        form.addRow("hold out", self.val_target)
        form.addRow("folds", self.val_folds)
        form.addRow("fraction hidden", self.val_hold)
        v.addLayout(form)

        b = QtWidgets.QPushButton("test the annotation")
        b.setProperty("primary", True)
        b.clicked.connect(self.run_validation)
        v.addWidget(b)

        self.val_table = QtWidgets.QTableWidget()
        self.val_table.setAlternatingRowColors(True)
        v.addWidget(self.val_table, 1)
        return w

    def run_validation(self):
        """Put an error rate on an annotation by hiding labels that are already known."""
        from .validate import validate_all
        if self.labels is None:
            self.status.emit("cluster a map first -- validation scores a clustering, not a map")
            return
        target = self.val_target.currentText()
        used = list(columns_for(self.nodes, self.spec()).keys())
        truth, labels = self.nodes[target], self.labels
        folds, frac = self.val_folds.value(), self.val_hold.value()

        def job(p):
            p(f"hiding {frac:.0%} of each category in {target}, {folds} folds")
            return validate_all(labels, truth, folds=folds, hold_frac=frac, used_columns=used)

        self._run(job, lambda d: (self._fill(self.val_table, d),
                                  self.status.emit(self._validation_verdict(d))),
                  name=f"validate annotation ({target})")

    @staticmethod
    def _validation_verdict(d) -> str:
        """One line leading with the number a candidate list is meaningless without."""
        if d is None or not len(d):
            return "no category had enough labelled genes to hide any"
        best = d.iloc[0]
        return (f"best: {best.category} at precision {best.precision:.2f} -- annotate from that "
                f"cluster and roughly {best.precision:.0%} would be right")

    # ------------------------------------------------------------------ jobs
    def _run(self, fn, on_done, name: str = "analysis"):
        """Run `fn` off the GUI thread, through the window's job runner when there is one.

        Routed through JobRunner rather than through a private QThread. The private one was the cause
        of three separate complaints at once: its work never appeared in the Jobs panel, there was no
        way to stop it, and -- worst -- it refused to start while anything else was running. A search
        takes minutes, so for those minutes every button on every other tab silently did nothing but
        write "a job is already running" into the status bar, which reads exactly like four broken
        tabs.

        `on_done` is invoked on the GUI thread. A handler connected straight to a worker signal runs
        on the worker, and touching widgets from there crashes on a slow machine.
        """
        if self.runner is not None:
            job = self.runner.submit(lambda j: fn(_Progress(self, j)), name)
            self._jobs[job.id] = on_done
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
                self.status.emit(f"failed: {error}")
            else:
                on_done(result)

        self._worker.done.connect(relay)
        self._thread.start()
        return None

    def _on_job_finished(self, jid: int, ok: bool):
        """Deliver a finished job's result to its handler, on the GUI thread."""
        on_done = self._jobs.pop(jid, None)
        if on_done is None:
            return
        job = self.runner.jobs.get(jid)
        if job is None:
            return
        if job.state == "cancelled":
            self.status.emit(f"{job.name}: stopped")
            return
        if not ok:
            self.status.emit(f"{job.name} failed: {job.error}")
            return
        on_done(job.result)

    def _fill(self, table: QtWidgets.QTableWidget, df: pd.DataFrame, limit=200):
        df = df.head(limit)
        table.clear()
        table.setRowCount(len(df)); table.setColumnCount(len(df.columns))
        table.setHorizontalHeaderLabels([str(c) for c in df.columns])
        for i, (_, r) in enumerate(df.iterrows()):
            for j, v in enumerate(r):
                s = f"{v:.3f}" if isinstance(v, float) and np.isfinite(v) else str(v)
                table.setItem(i, j, QtWidgets.QTableWidgetItem(s))
        table.resizeColumnsToContents()

    def run_umap_walk(self):
        """Score a grid of UMAP hyperparameters and fill the table with the ranking."""
        from .tuning import walk_umap
        spec, n, size, seed = self.spec(), self.nodes, self.sample.value(), self.seed.value()
        self._run(lambda p: walk_umap(n, spec, sample_size=size, seed=seed, log=p),
                  lambda d: (self._fill(self.walk_table, d), self.status.emit("walk complete")),
                  name="UMAP hyperparameter walk")

    def run_embed(self):
        """Build one embedding from the current spec and show it in the 3D view."""
        from .embedding import embed
        spec, n = self.spec(), self.nodes
        self._run(lambda p: embed(n, spec, log=p), self._embedded, name="build map")

    def _embedded(self, result):
        self.coords, self.features, self.rows = result
        self.embedding_ready.emit(self.coords, self.rows)
        self.status.emit(f"map built: {len(self.coords):,} genes, {len(self.features)} features")

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
        self.status.emit(f"{k} clusters, {100 * (labels == NOISE).mean():.0f}% unassigned")

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
        self.status.emit(f"battery: {len(held)} held-out features tested")

    def run_search(self):
        """Walk dataset combinations, scoring each by how well it recovers the held-out label."""
        from .search import search
        import itertools
        n, target = self.nodes, self.target.currentText()
        size, seed = self.search_sample.value(), self.seed.value()
        base = [b for b in BLOCKS if columns_for(n, EmbeddingSpec(blocks=(b,))).get(b)]
        r = self.max_blocks.value()
        sets = [tuple(c) for k in range(1, r + 1) for c in itertools.combinations(base, k)]

        def job(p):
            R, P = search(n, target=target, block_sets=sets, sample_size=size, seed=seed,
                          store=self.store, log=p)
            return R

        self._run(job, lambda R: (self._fill(self.search_table, R),
                                  self.status.emit(f"search complete: {len(R)} runs scored")),
                  name=f"recovery search ({target})")
