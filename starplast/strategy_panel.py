"""The Strategies tab: thirty-two ways to infer something, each explained, runnable and self-testing.

Docked to the right of Evidence and Analysis. The top half lists the strategies by family, with the
verdict each one earned when its self-test was run on the shipped data; the bottom half has three
tabs for the selected one:

    Guide     what it infers, why that works, how it fails, a step-by-step walkthrough, and how
              it is tested -- the text a person needs before believing its output
    Settings  its parameters, each with a tooltip saying why it exists, and Run / Test / Stop
    Results   the summary, the tables it produced, and the verdict of the last self-test

A strategy runs on the window's job runner, so it appears in the Jobs panel and can be stopped
there like anything else. Results tables are registered through `AnalysisPanel.results_table`
when the analysis panel is present, so they carry the same row action and right-click save as
every other results table in the application.

Nothing here computes anything: `strategies` and `strategy_catalog` do, and this is their face.
"""
from __future__ import annotations

import html
import json
import os

import numpy as np
import pandas as pd
from PyQt6 import QtCore, QtGui, QtWidgets

from . import strategies as S
from . import theme as TH
from .jobs import Stopped

DATA = os.path.join(os.path.dirname(__file__), "data")
#: The verdicts each strategy's self-test earned on the shipped tables, written by
#: `scripts/strategy_selftests.py`. Shown beside each strategy so its trustworthiness is visible
#: before anything is run.
MEASURED = os.path.join(DATA, "strategy_selftests.json")
#: How many result tables a strategy can show at once. A fixed pool, registered once, so every
#: table has its row action and save menu from the moment it exists.
TABLES = 6
VERDICT_COLORS = {"PASS": "#2e8b57", "FAIL": "#c0392b", "INCONCLUSIVE": "#888888"}

BUTTON_TIPS = {
    "run": "Run the selected strategy with these settings on the whole table, in the background. "
           "Its tables appear under Results; a strategy that builds a map can then be shown on it.",
    "test": "Run the strategy's self-test: hide information that is already known, ask the strategy "
            "for it back, and compare the answer with the same procedure on shuffled data. PASS "
            "means it beat that null by the stated margin on these settings.",
    "stop": "Ask the running strategy to stop. It stops at its next checkpoint, usually within a "
            "few seconds; nothing is written anywhere by a stopped run.",
    "map": "Put the map this result was computed on into the central view, colored by its clusters "
           "or modules, so the result can be looked at rather than only read.",
    "genes_file": "Read gene identifiers from a text or CSV file: every token in the first column, "
                  "or the whole file if it is a plain list. Anything not found is reported.",
    "genes_example": "Fill the list with the members of one known category of the default label, "
                     "so a gene-list strategy can be tried before you have a list of your own.",
    "genes_gate": "Use the genes currently gated on the map, so a structure you drew around can be "
                  "handed straight to a gene-list strategy.",
    "filter": "Type to show only strategies whose title, question or family contains the text -- "
              "for example 'list', 'network', 'held-out' or 'Plasmodium'.",
}


def load_measured(organism: str) -> dict:
    """The shipped self-test verdicts for one organism, or nothing if they were never measured."""
    try:
        with open(MEASURED) as fh:
            return json.load(fh).get(organism, {})
    except (OSError, ValueError):
        return {}


class _Progress:
    """The log a strategy job reports through: notes the job, and unwinds it when a stop is asked."""

    def __init__(self, panel, job):
        """Tie progress messages to a panel and a job."""
        self.panel, self.job = panel, job

    def stopped(self) -> bool:
        """Whether this job has been asked to stop."""
        return bool(getattr(self.job, "cancelled", False))

    def __call__(self, message):
        """Record a progress message, or raise `jobs.Stopped` if a stop was asked for."""
        if self.stopped():
            raise Stopped(f"{self.job.name} stopped")
        self.job.note = str(message)
        self.panel.status.emit(str(message))


class StrategyPanel(QtWidgets.QWidget):
    """The strategy list, a guide to the selected strategy, its settings, and its results."""

    #: A map to show: coordinates and the genes they place, as `AnalysisPanel.embedding_ready`.
    embedding_ready = QtCore.pyqtSignal(object, object)
    #: A clustering or module assignment over every gene, -1 for none.
    clusters_ready = QtCore.pyqtSignal(object)
    #: A gene someone clicked in a results table.
    gene_selected = QtCore.pyqtSignal(str)
    status = QtCore.pyqtSignal(str)
    #: A finished run (`strategies.StrategyResult`) and a finished self-test (`TestResult`).
    result_ready = QtCore.pyqtSignal(object)
    test_ready = QtCore.pyqtSignal(object)

    def __init__(self, nodes: pd.DataFrame, runner=None, analysis=None, organism: str | None = None,
                 graph=None, other=None, gated=None, parent=None):
        """Build the panel over one organism's table. Nothing is computed until something is run."""
        super().__init__(parent)
        self.ctx = S.Context(nodes, graph=graph, other=other, organism=organism)
        self.runner, self.analysis, self.gated = runner, analysis, gated
        self.measured = load_measured(self.ctx.organism)
        self.current: S.Strategy | None = None
        self.inputs: dict = {}
        self.last_result: S.StrategyResult | None = None
        self.last_test: S.TestResult | None = None
        self._jobs: dict = {}
        self._frames: dict = {}
        if runner is not None:
            runner.finished.connect(self._on_job_finished)

        self.filter = QtWidgets.QLineEdit()
        self.filter.setPlaceholderText("filter strategies…")
        self.filter.setToolTip(TH.tip(BUTTON_TIPS["filter"]))
        self.filter.textChanged.connect(self._apply_filter)
        self.tree = QtWidgets.QTreeWidget()
        self.tree.setHeaderLabels(["strategy", "on the shipped data"])
        self.tree.setRootIsDecorated(True)
        self.tree.itemSelectionChanged.connect(self._tree_selected)
        self._populate()

        self.guide = QtWidgets.QTextBrowser()
        self.guide.setOpenExternalLinks(True)

        self.form_host = QtWidgets.QWidget()
        self.form = QtWidgets.QFormLayout(self.form_host)
        self.run_btn = QtWidgets.QPushButton("Run")
        self.test_btn = QtWidgets.QPushButton("Test (hold-out)")
        self.stop_btn = QtWidgets.QPushButton("Stop")
        self.map_btn = QtWidgets.QPushButton("Show on map")
        for b, key in ((self.run_btn, "run"), (self.test_btn, "test"), (self.stop_btn, "stop"),
                       (self.map_btn, "map")):
            b.setToolTip(TH.tip(BUTTON_TIPS[key]))
        self.run_btn.clicked.connect(self.run_current)
        self.test_btn.clicked.connect(self.test_current)
        self.stop_btn.clicked.connect(self.stop_running)
        self.map_btn.clicked.connect(self.show_on_map)
        self.stop_btn.setEnabled(False)
        self.map_btn.setEnabled(False)
        buttons = QtWidgets.QHBoxLayout()
        for b in (self.run_btn, self.test_btn, self.stop_btn):
            buttons.addWidget(b)
        settings = QtWidgets.QWidget()
        sl = QtWidgets.QVBoxLayout(settings)
        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(self.form_host)
        sl.addWidget(scroll)
        sl.addLayout(buttons)

        self.summary = QtWidgets.QLabel("Nothing run yet.")
        self.summary.setWordWrap(True)
        self.summary.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        self.verdict = QtWidgets.QLabel("")
        self.verdict.setWordWrap(True)
        self.verdict.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        self.result_tabs = QtWidgets.QTabWidget()
        self.tables = [QtWidgets.QTableWidget() for _ in range(TABLES)]
        for i, t in enumerate(self.tables):
            self._wire_table(t, i)
        results = QtWidgets.QWidget()
        rl = QtWidgets.QVBoxLayout(results)
        rl.addWidget(self.verdict)
        rl.addWidget(self.summary)
        rl.addWidget(self.map_btn)
        rl.addWidget(self.result_tabs, 1)

        self.tabs = QtWidgets.QTabWidget()
        self.tabs.addTab(self.guide, "Guide")
        self.tabs.addTab(settings, "Settings")
        self.tabs.addTab(results, "Results")

        top = QtWidgets.QWidget()
        tl = QtWidgets.QVBoxLayout(top)
        tl.setContentsMargins(0, 0, 0, 0)
        tl.addWidget(self.filter)
        tl.addWidget(self.tree)
        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        split.addWidget(top)
        split.addWidget(self.tabs)
        split.setStretchFactor(0, 2)
        split.setStretchFactor(1, 3)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.addWidget(split)
        first = self.tree.topLevelItem(0)
        if first is not None and first.childCount():
            self.tree.setCurrentItem(first.child(0))

    # ------------------------------------------------------------------ the list
    def _populate(self):
        """Fill the tree: one branch per family, one leaf per strategy, each with its verdict."""
        self.tree.clear()
        self.items: dict = {}
        families: dict = {}
        for s in S.catalog():
            fam = families.get(s.family)
            if fam is None:
                fam = QtWidgets.QTreeWidgetItem([s.family, ""])
                fam.setFlags(fam.flags() & ~QtCore.Qt.ItemFlag.ItemIsSelectable)
                font = fam.font(0)
                font.setBold(True)
                fam.setFont(0, font)
                self.tree.addTopLevelItem(fam)
                families[s.family] = fam
            m = self.measured.get(s.key, {})
            verdict = m.get("verdict", "") if isinstance(m, dict) else ""
            item = QtWidgets.QTreeWidgetItem([f"{s.number:02d} · {s.title}", verdict])
            item.setData(0, QtCore.Qt.ItemDataRole.UserRole, s.key)
            item.setToolTip(0, TH.tip(s.tooltip))
            item.setToolTip(1, TH.tip(self._measured_line(s.key) or
                                      "Not yet measured on the shipped data: press Test."))
            if verdict in VERDICT_COLORS:
                item.setForeground(1, QtGui.QBrush(QtGui.QColor(VERDICT_COLORS[verdict])))
            fam.addChild(item)
            self.items[s.key] = item
        self.tree.expandAll()
        header = self.tree.header()
        header.setStretchLastSection(False)
        header.setSectionResizeMode(0, QtWidgets.QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeMode.ResizeToContents)

    def _apply_filter(self, text: str):
        """Hide the strategies that do not match, and the families left empty by that."""
        needle = text.strip().lower()
        for i in range(self.tree.topLevelItemCount()):
            fam = self.tree.topLevelItem(i)
            shown = 0
            for j in range(fam.childCount()):
                item = fam.child(j)
                s = S.get(item.data(0, QtCore.Qt.ItemDataRole.UserRole))
                hay = " ".join([s.title, s.question, s.family, s.tooltip, s.key]).lower()
                hide = bool(needle) and needle not in hay
                item.setHidden(hide)
                shown += not hide
            fam.setHidden(shown == 0)

    def _tree_selected(self):
        items = self.tree.selectedItems()
        key = items[0].data(0, QtCore.Qt.ItemDataRole.UserRole) if items else None
        if key:
            self.select(key)

    def select(self, key: str):
        """Make `key` the current strategy: its guide, and a fresh form of its parameters."""
        self.current = S.get(key)
        self.guide.setHtml(self.guide_html(self.current))
        self._build_form(self.current)
        item = self.items.get(key)
        if item is not None and self.tree.currentItem() is not item:
            self.tree.setCurrentItem(item)

    def _measured_line(self, key: str) -> str:
        m = self.measured.get(key)
        if not isinstance(m, dict) or "verdict" not in m:
            return ""
        if m["verdict"] == "INCONCLUSIVE":
            return f"INCONCLUSIVE on the shipped data -- {m.get('note', '')}"
        return (f"{m['verdict']} on the shipped data: {m['metric']} {m['observed']:.3f} against "
                f"{m['null_mean']:.3f} under {m['null_kind']} ({m['n_hidden']:,} hidden).")

    def guide_html(self, s: S.Strategy) -> str:
        """The Guide tab for one strategy, as HTML."""
        e = html.escape
        paras = "".join(f"<p>{e(p)}</p>" for p in s.explanation.split("\n\n"))
        steps = "".join(f"<li>{e(x)}</li>" for x in s.walkthrough)
        params = "".join(f"<li><b>{e(p.label)}</b> -- {e(p.tip)}</li>" for p in s.params)
        needs = f"<p><i>Needs {e(', '.join(s.needs))}.</i></p>" if s.needs else ""
        measured = self._measured_line(s.key)
        verdict = (self.measured.get(s.key) or {}).get("verdict", "")
        colour = VERDICT_COLORS.get(verdict, "#888888")
        mline = (f"<p style='color:{colour}'><b>{e(measured)}</b></p>" if measured else
                 "<p><i>Not yet measured on the shipped data -- press Test to measure it.</i></p>")
        return (f"<h3>{s.number:02d} · {e(s.title)}</h3><p><i>{e(s.question)}</i></p>{mline}"
                f"{paras}<h4>Walkthrough</h4><ol>{steps}</ol>"
                f"<h4>How it is tested</h4><p>{e(s.test_description)}</p>"
                f"<h4>Settings</h4><ul>{params}</ul>{needs}"
                f"<p><i>Typical cost: {e(s.cost)}.</i></p>")

    # ------------------------------------------------------------------ the settings form
    def _options(self, p: S.Param) -> list:
        ctx = self.ctx.other() if p.space == "other" else self.ctx
        if ctx is None:
            return []
        if p.kind == "category":
            return ctx.categorical_columns()
        if p.kind == "number":
            return ctx.numeric_targets()
        if p.kind == "column":
            return ctx.categorical_columns() + [c for c in ctx.numeric_targets()
                                                if c not in ctx.categorical_columns()]
        if p.kind == "layer":
            return ctx.layers()
        return list(p.options(self.ctx))

    def _build_form(self, s: S.Strategy):
        """One row per parameter, with its tooltip on both the label and the control."""
        while self.form.rowCount():
            self.form.removeRow(0)
        self.inputs = {}
        defaults = s.defaults(self.ctx)
        for p in s.params:
            w = self._widget(p, defaults.get(p.name))
            w.setToolTip(TH.tip(p.tip))
            label = QtWidgets.QLabel(p.label)
            label.setToolTip(TH.tip(p.tip))
            self.form.addRow(label, w)
            self.inputs[p.name] = (p, w)

    def _widget(self, p: S.Param, value):
        if p.kind in ("category", "number", "column", "layer", "choice"):
            w = QtWidgets.QComboBox()
            if p.optional:
                w.addItem("(none)", None)
            for o in self._options(p):
                w.addItem(str(o), o)
            i = w.findData(value)
            if i >= 0:
                w.setCurrentIndex(i)
            return w
        if p.kind == "int":
            w = QtWidgets.QSpinBox()
            w.setRange(int(p.lo), int(p.hi))
            w.setSingleStep(max(1, int(p.step)))
            w.setValue(int(value if value is not None else p.lo))
            return w
        if p.kind == "float":
            w = QtWidgets.QDoubleSpinBox()
            w.setDecimals(3)
            w.setRange(float(p.lo), float(p.hi))
            w.setSingleStep(float(p.step))
            w.setValue(float(value if value is not None else p.lo))
            return w
        if p.kind == "grid":
            return QtWidgets.QLineEdit(str(value or ""))
        return self._genes_widget(p, value)

    def _genes_widget(self, p: S.Param, value):
        host = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(host)
        lay.setContentsMargins(0, 0, 0, 0)
        box = QtWidgets.QPlainTextEdit(str(value or ""))
        box.setPlaceholderText("one accession per line, or separated by commas or spaces")
        box.setToolTip(TH.tip(p.tip))
        row = QtWidgets.QHBoxLayout()
        load = QtWidgets.QPushButton("Load file…")
        example = QtWidgets.QPushButton("Example set")
        gate = QtWidgets.QPushButton("From gate")
        for b, key in ((load, "genes_file"), (example, "genes_example"), (gate, "genes_gate")):
            b.setToolTip(TH.tip(BUTTON_TIPS[key]))
            row.addWidget(b)
        load.clicked.connect(lambda: self.load_genes_file(box))
        example.clicked.connect(lambda: self.fill_example(box))
        gate.clicked.connect(lambda: self.fill_from_gate(box))
        gate.setEnabled(self.gated is not None)
        lay.addWidget(box)
        lay.addLayout(row)
        host.box = box
        return host

    def load_genes_file(self, box, path: str = "") -> int:
        """Fill a gene box from a file. Returns how many identifiers were read."""
        if not path:
            path, _ = QtWidgets.QFileDialog.getOpenFileName(self, "Gene list", "",
                                                            "Text or CSV (*.txt *.csv *.tsv);;All (*)")
        if not path:
            return 0
        with open(path, errors="replace") as fh:
            lines = [ln.strip() for ln in fh if ln.strip()]
        tokens = [ln.replace("\t", ",").split(",")[0].strip().strip('"') for ln in lines]
        box.setPlainText("\n".join(t for t in tokens if t))
        found, missing = self.ctx.resolve_genes(tokens)
        self.status.emit(f"{len(found)} genes found in this table, {len(missing)} not found")
        return len(tokens)

    def fill_example(self, box) -> int:
        """Fill a gene box with a known category's members, so a list strategy can be tried."""
        cat, members = S.example_set(self.ctx, S.default_category(self.ctx))
        box.setPlainText("\n".join(self.ctx.gene_ids[members]))
        self.status.emit(f"example set: the {len(members)} genes of {cat!r}" if cat else
                         "no category of a usable size to take an example from")
        return len(members)

    def fill_from_gate(self, box) -> int:
        """Fill a gene box with the genes gated on the map."""
        genes = list(self.gated() if callable(self.gated) else (self.gated or []))
        box.setPlainText("\n".join(map(str, genes)))
        return len(genes)

    def settings(self) -> dict:
        """The current form as the keyword arguments the strategy takes."""
        out = {}
        for name, (p, w) in self.inputs.items():
            if isinstance(w, QtWidgets.QComboBox):
                out[name] = w.currentData()
            elif isinstance(w, (QtWidgets.QSpinBox, QtWidgets.QDoubleSpinBox)):
                out[name] = w.value()
            elif isinstance(w, QtWidgets.QLineEdit):
                out[name] = w.text()
            else:
                out[name] = w.box.toPlainText()
        return out

    def set_setting(self, name: str, value):
        """Set one parameter's control -- for scripting the panel and for tests."""
        p, w = self.inputs[name]
        if isinstance(w, QtWidgets.QComboBox):
            i = w.findData(value)
            if i < 0:
                raise ValueError(f"{value!r} is not an option for {name}")
            w.setCurrentIndex(i)
        elif isinstance(w, (QtWidgets.QSpinBox, QtWidgets.QDoubleSpinBox)):
            w.setValue(value)
        elif isinstance(w, QtWidgets.QLineEdit):
            w.setText(str(value))
        else:
            w.box.setPlainText(value if isinstance(value, str) else "\n".join(map(str, value)))

    # ------------------------------------------------------------------ running
    def run_current(self):
        """Run the selected strategy with the current settings."""
        return self._submit("run")

    def test_current(self):
        """Run the selected strategy's self-test with the current settings."""
        return self._submit("test")

    def _submit(self, kind: str):
        if self.current is None:
            return None
        s, params = self.current, self.settings()
        name = f"strategy {s.number:02d} {'test' if kind == 'test' else 'run'}: {s.title}"

        def work(log, should_stop):
            ctx = self.ctx.bound(log=log, should_stop=should_stop)
            return s.test(ctx, **params) if kind == "test" else s.run(ctx, **params)

        self.status.emit(f"{name} -- started")
        if self.runner is None:
            try:
                out = work(self.status.emit, None)
            except Exception as exc:                  # shown, never raised into Qt
                self._failed(exc)
                return None
            self._deliver(kind, out)
            return out
        job = self.runner.submit(lambda j: work(_Progress(self, j), lambda: j.cancelled), name)
        self._jobs[job.id] = (kind, job)
        self.stop_btn.setEnabled(True)
        return job

    def _on_job_finished(self, jid: int, ok: bool):
        entry = self._jobs.pop(jid, None)
        if entry is None:
            return
        kind, job = entry
        self.stop_btn.setEnabled(bool(self._jobs))
        if job.state == "cancelled":
            self.status.emit(f"{job.name}: stopped")
            return
        if not ok:
            self._failed(job.exception or RuntimeError(job.error))
            return
        self._deliver(kind, job.result)

    def _failed(self, exc):
        text = f"{type(exc).__name__}: {exc}"
        self.summary.setText(f"Could not run: {text}")
        self.status.emit(f"strategy failed -- {text}")
        self.tabs.setCurrentIndex(2)

    def stop_running(self) -> int:
        """Ask every job this panel started to stop. Returns how many were asked."""
        n = 0
        for _kind, job in list(self._jobs.values()):
            if job.active:
                job.cancel()
                n += 1
        return n

    def _deliver(self, kind: str, out):
        if kind == "test":
            self.show_test(out)
        else:
            self.show_result(out)

    # ------------------------------------------------------------------ results
    def _wire_table(self, table, index: int):
        table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        if self.analysis is not None:
            self.analysis.results_table(table, on_row=lambda row, t=table: self._row(t, row),
                                        what=f"strategy result {index + 1}")
        else:
            table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
            table.cellClicked.connect(lambda row, _col, t=table: self._row(t, row))

    def _fill(self, table, frame: pd.DataFrame):
        if self.analysis is not None:
            self.analysis._fill(table, frame)
            return
        self._frames[table] = frame
        shown = frame.head(200)
        table.setSortingEnabled(False)
        table.clear()
        table.setRowCount(len(shown))
        table.setColumnCount(len(shown.columns))
        table.setHorizontalHeaderLabels([str(c) for c in shown.columns])
        for i, row in enumerate(shown.itertuples(index=False)):
            for j, v in enumerate(row):
                item = QtWidgets.QTableWidgetItem()
                if isinstance(v, (int, float, np.integer, np.floating)) and np.isfinite(v):
                    item.setData(QtCore.Qt.ItemDataRole.DisplayRole, float(v))
                    item.setText(f"{v:.3f}" if isinstance(v, (float, np.floating)) else str(v))
                else:
                    item.setText(str(v))
                table.setItem(i, j, item)
        table.resizeColumnsToContents()
        table.setSortingEnabled(True)

    def frame(self, index: int) -> pd.DataFrame:
        """The whole frame behind result table `index`, not only the rows on screen."""
        t = self.tables[index]
        frames = self.analysis._frames if self.analysis is not None else self._frames
        return frames.get(t, pd.DataFrame())

    def _row(self, table, row: int):
        """A clicked row: select its gene on the map if it names one."""
        header = [table.horizontalHeaderItem(j).text() if table.horizontalHeaderItem(j) else ""
                  for j in range(table.columnCount())]
        for col in ("gene_id", "gene_a"):
            if col in header:
                item = table.item(row, header.index(col))
                if item is not None and item.text():
                    self.gene_selected.emit(item.text())
                return

    def _show_tables(self, tables: dict):
        self.result_tabs.clear()
        for i, (name, frame) in enumerate(list(tables.items())[:TABLES]):
            frame = frame if isinstance(frame, pd.DataFrame) else pd.DataFrame(frame)
            self._fill(self.tables[i], frame)
            self.result_tabs.addTab(self.tables[i], f"{name} ({len(frame):,})")

    def show_result(self, result: S.StrategyResult):
        """Put a finished run on screen: its summary, its tables, and a way to see its map."""
        self.last_result = result
        self.summary.setText(result.summary)
        self._show_tables(result.tables)
        self.map_btn.setEnabled(result.coords is not None or result.labels is not None)
        self.tabs.setCurrentIndex(2)
        self.status.emit(f"{result.strategy}: done in {result.seconds:.1f}s")
        self.result_ready.emit(result)

    def show_test(self, test: S.TestResult):
        """Put a finished self-test on screen, leading with its verdict."""
        self.last_test = test
        colour = VERDICT_COLORS.get(test.verdict, "#888888")
        # `summary()` leads with the verdict; it is drawn once, in color, and the rest follows.
        rest = test.summary().split(" -- ", 1)[-1]
        self.verdict.setText(f"<span style='color:{colour}'><b>{html.escape(test.verdict)}</b></span>"
                             f" -- {html.escape(rest)}<br><i>Hidden: "
                             f"{html.escape(test.hidden)}. Null: {html.escape(test.null_kind)}.</i>")
        tables = {"self-test": test.details} if len(test.details) else {}
        if test.numbers:
            tables["numbers"] = pd.DataFrame({"quantity": list(test.numbers),
                                              "value": list(test.numbers.values())})
        if tables:
            self._show_tables(tables)
        self.tabs.setCurrentIndex(2)
        self.status.emit(f"{test.strategy}: {test.verdict}")
        self.test_ready.emit(test)

    def show_on_map(self) -> bool:
        """Send the last result's map and clusters to the central view."""
        r = self.last_result
        if r is None:
            return False
        if r.coords is not None and r.positions is not None:
            order = np.argsort(np.asarray(r.positions))
            mask = np.zeros(self.ctx.n, bool)
            mask[np.asarray(r.positions, dtype=int)] = True
            self.embedding_ready.emit(np.asarray(r.coords)[order], mask)
        if r.labels is not None:
            self.clusters_ready.emit(np.asarray(r.labels))
        return True
