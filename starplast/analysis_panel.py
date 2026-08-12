#!/usr/bin/env python3
"""The analysis tab: choose the data, tune the map, cluster it, and ask what the clusters mean.

This is the interface to `embedding`, `clustering`, `tuning` and `search`. It is deliberately arranged in
the order the work is actually done, top to bottom, because the order matters scientifically:

    1  data      which blocks feed the map, and how missing values and scales are handled
    2  map       UMAP hyperparameters, with a seeded walk to pick them
    3  clusters  DBSCAN / HDBSCAN, with a walk to pick those too
    4  meaning   what the clusters correspond to, held-out features only
    5  search    walk dataset combinations looking for structure that recovers a held-out label

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
        try:
            self.done.emit(self.fn(self.progress.emit), None)
        except Exception as e:                                   # pragma: no cover - GUI path
            self.done.emit(None, e)


class AnalysisPanel(QtWidgets.QWidget):
    """Data selection, tuning, clustering and the hypothesis battery."""

    embedding_ready = QtCore.pyqtSignal(object, object)          # coords, gene mask
    status = QtCore.pyqtSignal(str)

    def __init__(self, nodes: pd.DataFrame, store=None, parent=None):
        super().__init__(parent)
        self.nodes = nodes
        self.store = store
        self.labels = None
        self.coords = None
        self.rows = None
        self._thread = None
        self._worker = None

        tabs = QtWidgets.QTabWidget()
        tabs.addTab(self._data_tab(), "1 · Data")
        tabs.addTab(self._map_tab(), "2 · Map")
        tabs.addTab(self._cluster_tab(), "3 · Clusters")
        tabs.addTab(self._meaning_tab(), "4 · Meaning")
        tabs.addTab(self._search_tab(), "5 · Search")
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(tabs)

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
        return EmbeddingSpec(
            blocks=tuple(b for b, cb in self.block_cb.items() if cb.isChecked()),
            categorical=("compartment",) if self.cat_cb.isChecked() else (),
            na_policy=self.na_policy.currentText(),
            scaling=self.scaling.currentText(),
            max_missing=self.max_missing.value(),
            n_neighbors=self.nn.value(), min_dist=self.md.value(),
            random_state=self.seed.value())

    def show_variance(self):
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

    # ------------------------------------------------------------------ jobs
    def _run(self, fn, on_done):
        """Run `fn` off the GUI thread.

        `done` is relayed through this bound method rather than connected straight to `on_done`, so the
        handler runs on the GUI thread. A directly-connected handler executes on the worker and touching
        widgets from there crashes on a slow machine.
        """
        if self._thread is not None:
            self.status.emit("a job is already running")
            return
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
        from .tuning import walk_umap
        spec, n, size, seed = self.spec(), self.nodes, self.sample.value(), self.seed.value()
        self._run(lambda p: walk_umap(n, spec, sample_size=size, seed=seed, log=p),
                  lambda d: (self._fill(self.walk_table, d), self.status.emit("walk complete")))

    def run_embed(self):
        from .embedding import embed
        spec, n = self.spec(), self.nodes
        self._run(lambda p: embed(n, spec, log=p), self._embedded)

    def _embedded(self, result):
        self.coords, self.features, self.rows = result
        self.embedding_ready.emit(self.coords, self.rows)
        self.status.emit(f"map built: {len(self.coords):,} genes, {len(self.features)} features")

    def save_embedding(self):
        if self.coords is None or self.store is None:
            self.status.emit("build a map first")
            return
        name = self.emb_name.text().strip() or "unnamed"
        self.store.save(name, self.coords, self.spec(),
                        gene_ids=self.nodes.gene_id[self.rows], features=self.features)
        self.status.emit(f"saved embedding {name!r} with its full recipe")

    def run_cluster_walk(self):
        from .clustering import walk_dbscan, walk_hdbscan
        if self.coords is None:
            self.status.emit("build a map first"); return
        Y, algo = self.coords, self.algo.currentText()
        fn = walk_hdbscan if algo == "hdbscan" else walk_dbscan
        self._run(lambda p: fn(Y, log=p),
                  lambda d: (self._fill(self.cluster_table, d), self.status.emit("walk complete")))

    def run_cluster(self):
        from .clustering import NOISE, cluster
        if self.coords is None:
            self.status.emit("build a map first"); return
        Y, algo, mcs, eps = self.coords, self.algo.currentText(), self.mcs.value(), self.eps.value()
        self._run(lambda p: cluster(Y, algorithm=algo, min_cluster_size=mcs,
                                    min_samples=mcs, eps=eps),
                  self._clustered)

    def _clustered(self, labels):
        from .clustering import NOISE
        self.labels = labels
        k = len(set(labels[labels != NOISE]))
        self.status.emit(f"{k} clusters, {100 * (labels == NOISE).mean():.0f}% unassigned")

    def run_battery(self):
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

        self._run(job, self._battery_done)

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
                                  self.status.emit(f"search complete: {len(R)} runs scored")))
