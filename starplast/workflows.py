"""Guided gene exploration, held-out prediction and measured-screen comparison.

The map remains a browser. Prediction jobs run in the shared background runner,
and exports retain the full evaluation contract rather than just a ranked score.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PyQt6 import QtCore, QtWidgets

from . import paths
from .jobs import JobRunner, Stopped


def _combo(items, tip):
    box = QtWidgets.QComboBox()
    for label, value in items:
        box.addItem(label, value)
    box.setToolTip(tip)
    return box


def _table():
    table = QtWidgets.QTableWidget()
    table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
    table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
    table.setAlternatingRowColors(True)
    table.horizontalHeader().setStretchLastSection(True)
    return table


def _fill(table, frame, limit=1000):
    frame = frame.head(limit)
    table.setSortingEnabled(False)
    table.clear()
    table.setColumnCount(len(frame.columns)); table.setRowCount(len(frame))
    table.setHorizontalHeaderLabels(list(map(str, frame.columns)))
    for i, row in enumerate(frame.itertuples(index=False, name=None)):
        for j, value in enumerate(row):
            text = '' if value is None or (np.isscalar(value) and pd.isna(value)) else str(value)
            item = QtWidgets.QTableWidgetItem(text)
            item.setToolTip(text)
            table.setItem(i,j,item)
    table.resizeColumnsToContents()
    table.setSortingEnabled(True)


class WorkflowDialog(QtWidgets.QDialog):
    """Three task-oriented views sharing the application's gene table and job runner."""
    gene_selected = QtCore.pyqtSignal(str)

    def __init__(self, nodes, runner=None, parent=None):
        super().__init__(parent)
        self.setWindowTitle('Starplast · Guided workflows')
        self.resize(1060,740)
        self.nodes = nodes.copy()
        self.runner = runner or JobRunner(self)
        self.result = None; self.comparison = None; self.screen = None; self.job = None
        self.runner.finished.connect(self._finished)
        self.runner.progress.connect(self._progress)
        self.tabs = QtWidgets.QTabWidget()
        self.status = QtWidgets.QLabel('Choose a task. Double-click a gene in a table to show it on the map.')
        self.status.setWordWrap(True)
        layout = QtWidgets.QVBoxLayout(self)
        layout.addWidget(self.tabs); layout.addWidget(self.status)
        self._explore_tab(); self._predict_tab(); self._screen_tab()

    def _explore_tab(self):
        page = QtWidgets.QWidget(); layout = QtWidgets.QVBoxLayout(page)
        line = QtWidgets.QHBoxLayout()
        self.gene = QtWidgets.QLineEdit()
        self.gene.setPlaceholderText('Gene accession, symbol or product')
        self.gene.setToolTip('Search this organism’s identifiers, symbols and product descriptions.')
        find = QtWidgets.QPushButton('Explore gene'); find.clicked.connect(self._explore)
        self.gene.returnPressed.connect(self._explore)
        line.addWidget(self.gene); line.addWidget(find); layout.addLayout(line)
        self.matches = _table(); self.matches.setMaximumHeight(160)
        self.matches.cellDoubleClicked.connect(self._open_match)
        layout.addWidget(self.matches)
        self.evidence = _table(); layout.addWidget(self.evidence)
        note = QtWidgets.QLabel('Each row is an existing gene summary. Missing measurements remain unknown; '
                               'source links describe where a value came from. Replicates require the original source table.')
        note.setWordWrap(True); layout.addWidget(note)
        self.tabs.addTab(page,'Explore a gene')

    def _explore(self):
        query = self.gene.text().strip()
        if not query:
            self.status.setText('Enter an accession, symbol or product name.'); return
        mask = np.zeros(len(self.nodes),dtype=bool)
        for col in ('gene_id','symbol','product'):
            if col in self.nodes:
                mask |= self.nodes[col].astype(str).str.contains(query,case=False,regex=False).to_numpy()
        rows = self.nodes.loc[mask]
        columns = [c for c in ('gene_id','symbol','product') if c in rows]
        _fill(self.matches,rows[columns],100)
        self.status.setText(f'{len(rows):,} matches. Showing at most 100; double-click to inspect another gene.')
        if len(rows): self._show_evidence(str(rows.gene_id.iloc[0]))
        else: self.evidence.setRowCount(0)

    def _open_match(self, row, _column):
        if self.matches.item(row,0): self._show_evidence(self.matches.item(row,0).text())

    def _show_evidence(self, gene):
        from .evidence import from_gene_table
        from .slots import table_organism
        subset = self.nodes[self.nodes.gene_id.astype(str).eq(gene)]
        records = list(from_gene_table(subset,table_organism(self.nodes) or 'Tg','current_session',include_missing=True))
        _fill(self.evidence,pd.DataFrame([{'trait':r.trait,'value':r.value,'missingness':r.missing_state,
                                         'source':r.source_id,'evidence status':r.evidence_status,
                                         'source URL':r.source_url} for r in records]))
        self.gene_selected.emit(gene)

    def _predict_tab(self):
        page = QtWidgets.QWidget(); layout = QtWidgets.QVBoxLayout(page)
        form = QtWidgets.QFormLayout()
        targets = [c for c in self.nodes if c!='gene_id']
        self.target = _combo([(c,c) for c in targets],
                            'The observed outcome to hold out. Unknown entries never become negative examples.')
        self.target.setEditable(True)
        self.target.setInsertPolicy(QtWidgets.QComboBox.InsertPolicy.NoInsert)
        if 'compartment' in targets: self.target.setCurrentText('compartment')
        self.kind = _combo([('Category / class','classification'),('Continuous measurement','regression')],
                           'Classification predicts labels. Regression preserves the numerical outcome and reports residual intervals.')
        self.method = _combo([('Linear baseline','linear'),('Boosted trees','boosted'),
                              ('Feature-space neighbours','neighbors'),('PCA neighbours','pca'),
                              ('UMAP neighbours','umap'),('Masked factors','multiview'),('Prior baseline','prior')],
                             'Compare methods on identical held-out families. The display map is not used for prediction.')
        groups = [(c,c) for c in self.nodes if c!='gene_id' and not pd.api.types.is_numeric_dtype(self.nodes[c])]
        self.group = _combo([('Separate genes (random holdout)',None)]+groups,
                            'Keep related proteins or study groups together. Choose orthogroup for family holdout; '
                            'random gene folds can overestimate transfer to new families.')
        if 'orthogroup' in self.nodes: self.group.setCurrentText('orthogroup')
        self.folds = QtWidgets.QSpinBox(); self.folds.setRange(2,10); self.folds.setValue(3)
        self.folds.setToolTip('Number of outer evaluation folds. Calibration uses separate groups inside each training fold.')
        self.sequence = QtWidgets.QCheckBox('Include frozen protein sequence features')
        self.sequence.setToolTip('Load the bundled ESM-2 table by gene ID. This adds sequence evidence without folding proteins or downloading a model.')
        self.sequence.setEnabled(bool(self.nodes.gene_id.astype(str).str.startswith('TGME49_').any()))
        self.threshold = QtWidgets.QDoubleSpinBox(); self.threshold.setRange(0,1); self.threshold.setSingleStep(.05)
        self.threshold.setToolTip('Abstain below this model probability. A value of zero applies only the measured-feature coverage rule. '
                                  'Check calibration status before interpreting probabilities.')
        for label, widget in [('Outcome',self.target),('Outcome type',self.kind),('Method',self.method),
                              ('Hold out groups by',self.group),('Evaluation folds',self.folds),
                              ('Minimum model probability',self.threshold),('Sequence evidence',self.sequence)]:
            form.addRow(label,widget)
        layout.addLayout(form)
        line = QtWidgets.QHBoxLayout()
        self.run_button = QtWidgets.QPushButton('Evaluate and predict'); self.run_button.clicked.connect(self._run)
        self.cancel_button = QtWidgets.QPushButton('Stop'); self.cancel_button.setEnabled(False)
        self.cancel_button.clicked.connect(self._cancel)
        self.export_button = QtWidgets.QPushButton('Export full run…'); self.export_button.setEnabled(False)
        self.export_button.clicked.connect(self._export)
        for widget in (self.run_button,self.cancel_button,self.export_button): line.addWidget(widget)
        layout.addLayout(line)
        self.summary = QtWidgets.QLabel('Results report held-out performance separately from hypotheses for unlabelled genes.')
        self.summary.setWordWrap(True); self.summary.setTextInteractionFlags(QtCore.Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(self.summary)
        self.predictions = _table(); self.predictions.cellDoubleClicked.connect(self._select_prediction)
        layout.addWidget(self.predictions)
        self.tabs.addTab(page,'Predict a trait')

    def _run(self):
        from .prediction import TaskSpec, run
        if self.job is not None and self.job.active: return
        target = self.target.currentText()
        if target not in self.nodes:
            self.status.setText('Choose an outcome column present in this table.'); return
        spec = TaskSpec(target,kind=self.kind.currentData(),method=self.method.currentData(),
                        group_column=self.group.currentData(),folds=self.folds.value(),
                        min_probability=self.threshold.value())
        nodes = self.nodes.copy(); add_sequence = self.sequence.isChecked() and self.sequence.isEnabled()
        self.result = None; self.export_button.setEnabled(False); self.predictions.setRowCount(0)
        self.summary.setText('Running held-out evaluation…')
        self.run_button.setEnabled(False); self.cancel_button.setEnabled(True)
        def work(job):
            def log(message):
                if job.cancelled: raise Stopped('Prediction cancelled between model fits')
                job.note=message; self.runner._sig.progress.emit(job.id,-1,message)
            frame=nodes
            if add_sequence:
                log('Loading frozen sequence evidence')
                sequence=pd.read_parquet(paths.cache_file('esm_features.parquet'))
                frame=frame.merge(sequence,on='gene_id',how='left',validate='one_to_one')
            log('Preparing evaluation')
            from threadpoolctl import threadpool_limits
            with threadpool_limits(limits=2):
                return run(frame,spec,log=log)
        self.job = self.runner.submit(work,f'Predict {target} ({spec.method})')

    def _progress(self, job_id, _percent, note):
        if self.job is not None and job_id==self.job.id: self.status.setText(note)

    def _cancel(self):
        if self.job is not None: self.job.cancel()
        self.status.setText('Stop requested. The current model fit will finish before cancellation.')

    def _finished(self, job_id, ok):
        if self.job is None or job_id!=self.job.id: return
        self.run_button.setEnabled(True); self.cancel_button.setEnabled(False)
        if not ok:
            self.summary.setText('No result: '+(self.job.error or self.job.note or self.job.state))
            self.status.setText('The job record retains the error or cancellation details.'); return
        self.result = self.job.result; self.export_button.setEnabled(True)
        self.show_result(self.result)

    def show_result(self, result):
        """Present evaluated performance and candidate calls with their support status."""
        metrics = ', '.join(f'{key}: {value:.3f}' for key,value in result.metrics.items()
                            if key in {'accuracy','macro_f1','macro_average_precision','rmse','r2','coverage'})
        candidates=result.predictions[result.predictions.role.eq('unlabelled_candidate')]
        self.summary.setText(f'Held-out evaluation — {metrics}. '
                             f'{int(candidates.supported.sum()):,} supported hypotheses among {len(candidates):,} unlabelled genes. '
                             'These are model predictions, not measured annotations.')
        columns=[c for c in ('gene_id','prediction','score','feature_coverage','supported','abstention_reason',
                             'calibration_status','interval_lower','interval_upper') if c in candidates]
        _fill(self.predictions,candidates[columns])
        self.status.setText('Showing up to 1,000 unlabelled rows. Export includes every row, per-class results, split IDs and provenance.')

    def _select_prediction(self, row, _column):
        if self.predictions.item(row,0): self.gene_selected.emit(self.predictions.item(row,0).text())

    def _export(self):
        if self.result is None: return
        directory=QtWidgets.QFileDialog.getExistingDirectory(self,'Export evaluated run')
        if directory:
            self.result.save(directory); self.status.setText(f'Exported complete run to {directory}')

    def _screen_tab(self):
        page=QtWidgets.QWidget(); layout=QtWidgets.QVBoxLayout(page)
        note=QtWidgets.QLabel('Import a gene-level spaCR or other measured table. Choose the identifier and effect columns. '
                              'Unmapped identifiers and missing measurements remain visible; duplicate genes need explicit aggregation first.')
        note.setWordWrap(True); layout.addWidget(note)
        load=QtWidgets.QPushButton('Open CSV or TSV…'); load.clicked.connect(self._load_screen); layout.addWidget(load)
        form=QtWidgets.QFormLayout()
        self.screen_gene=_combo([], 'The column containing one gene identifier per row.')
        self.screen_value=_combo([], 'A measured effect or other numerical result; missing entries stay missing.')
        form.addRow('Gene identifier',self.screen_gene); form.addRow('Measured value',self.screen_value); layout.addLayout(form)
        compare=QtWidgets.QPushButton('Compare with gene evidence'); compare.clicked.connect(self._compare)
        layout.addWidget(compare)
        self.screen_table=_table(); layout.addWidget(self.screen_table)
        export=QtWidgets.QPushButton('Export comparison…'); export.clicked.connect(self._export_screen); layout.addWidget(export)
        self.tabs.addTab(page,'Compare a screen')

    def _load_screen(self):
        path,_=QtWidgets.QFileDialog.getOpenFileName(self,'Measured gene table','', 'Tables (*.csv *.tsv *.txt)')
        if path:
            try: self.load_screen(path)
            except Exception as exc: self.status.setText(f'Cannot read table: {exc}')

    def load_screen(self, path):
        """Read a measured table and preserve its original content hash for export."""
        path=Path(path)
        self.screen=pd.read_csv(path,sep='\t' if path.suffix.lower() in {'.tsv','.txt'} else ',')
        self.screen_source={'filename':path.name,'sha256':hashlib.sha256(path.read_bytes()).hexdigest()}
        self.comparison=None
        for box in (self.screen_gene,self.screen_value):
            box.clear(); box.addItems(list(self.screen.columns))
        if 'gene_id' in self.screen: self.screen_gene.setCurrentText('gene_id')
        numeric=list(self.screen.select_dtypes(include=np.number))
        if numeric: self.screen_value.setCurrentText(numeric[0])
        _fill(self.screen_table,self.screen)
        self.status.setText(f'Loaded {len(self.screen):,} rows. Select columns and compare.')

    def _compare(self):
        from .prioritization import compare_screen
        if self.screen is None: self.status.setText('Open a measured table first.'); return
        try:
            self.comparison=compare_screen(self.nodes,self.screen,self.screen_value.currentText(),self.screen_gene.currentText())
            self.screen_source.update(gene_column=self.screen_gene.currentText(),value_column=self.screen_value.currentText())
            columns=['input_gene_id','gene_id','mapping_status','screen_effect','measurement_status']
            if 'product' in self.comparison: columns.append('product')
            _fill(self.screen_table,self.comparison[columns])
            mapped=int(self.comparison.mapping_status.eq('mapped').sum())
            self.status.setText(f'{mapped:,} of {len(self.comparison):,} rows mapped. Unresolved and unmeasured rows remain in the export.')
        except Exception as exc:
            self.comparison=None; self.status.setText(f'Comparison stopped: {exc}')

    def _export_screen(self):
        if self.comparison is None: self.status.setText('Compare a measured table first.'); return
        path,_=QtWidgets.QFileDialog.getSaveFileName(self,'Export screen comparison','','CSV (*.csv)')
        if path:
            self.comparison.to_csv(path,index=False)
            Path(path).with_suffix('.json').write_text(json.dumps(self.screen_source,indent=2)+'\n')
            self.status.setText(f'Exported comparison and source record to {path}')
