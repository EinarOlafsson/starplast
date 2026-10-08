"""Browse recorded label memberships and separate mechanism evaluation scopes.

Labels begin with function and retain the backend's native class keys, unknown
counts and annotation interpretations. Every gene membership remains reachable
through bounded pages. Shared evidence and performance cards stay separate, and
multiple evaluations retain their original setting, seed and partition. This
widget routes recorded sources and genes without fitting strategies, pooling
cohorts, estimating confidence or admitting annotation recovery as biological
accuracy.
"""
import json

import pandas as pd
from PyQt6 import QtCore, QtWidgets

from .scorecard_browser import ScorecardBrowser


def _text(value):
    if value is None or value is pd.NA or value is pd.NaT:
        return 'Unavailable'
    missing = pd.isna(value)
    return 'Unavailable' if pd.api.types.is_scalar(missing) and bool(missing) else str(value)


def _table(columns):
    table = QtWidgets.QTableWidget(0, len(columns))
    table.setHorizontalHeaderLabels(columns)
    table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
    table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
    table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
    table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Interactive)
    table.horizontalHeader().setStretchLastSection(True)
    for i, label in enumerate(columns):
        table.horizontalHeaderItem(i).setToolTip(label)
    return table


def _fill(table, frame):
    table.setColumnCount(len(frame.columns))
    table.setHorizontalHeaderLabels([str(column) for column in frame.columns])
    table.setRowCount(len(frame))
    for j, column in enumerate(frame):
        table.horizontalHeaderItem(j).setToolTip(str(column))
        for i, value in enumerate(frame[column]):
            item = QtWidgets.QTableWidgetItem(_text(value))
            item.setToolTip(repr(value))
            table.setItem(i, j, item)


class LabelSpaceBrowser(QtWidgets.QWidget):
    """Explore one organism's known labels, native classes and recorded strategies.

    Label and family filters preserve the source inventory. Class balances and
    unknown counts remain distinct from evaluation metrics; selecting a class
    never transforms its native membership into a new major class. Gene records
    are paged without excluding overlap, while each strategy evaluation retains
    its supplied cohort and settings in a separate immutable scorecard. Missing
    tests show explicit gaps, and navigation stays within the selected organism.
    """

    gene_requested = QtCore.pyqtSignal(str, str)
    source_requested = QtCore.pyqtSignal(str)

    def __init__(self, model, parent=None):
        super().__init__(parent)
        self.model, self.organism = model, model.organism
        self.target, self.value = None, None
        self.page_start, self.page_size = 0, 200
        self.member_frame = pd.DataFrame()
        self.class_frame = pd.DataFrame()
        self.mechanism_records = []
        self.evaluations = []
        self._updating = False
        layout = QtWidgets.QVBoxLayout(self)
        controls = QtWidgets.QHBoxLayout()
        self.search = QtWidgets.QLineEdit(placeholderText='Search labels and recorded sources…')
        self.search.setToolTip('Filter available label names, families and recorded source identifiers without changing their memberships.')
        controls.addWidget(self.search, 1)
        self.family = QtWidgets.QComboBox()
        self.family.setToolTip('Filter the supplied labels by biological question family. Function labels appear first.')
        self.family.addItem('All families', None)
        for family in sorted(set(model.labels.family), key=lambda value: (value != 'Function', value)):
            self.family.addItem(family, family)
        controls.addWidget(self.family)
        self.target_choice = QtWidgets.QComboBox()
        self.target_choice.setToolTip('Choose a supplied label; known annotation membership and independent prediction accuracy remain separate.')
        controls.addWidget(self.target_choice, 1)
        self.class_choice = QtWidgets.QComboBox()
        self.class_choice.setToolTip('Choose an original annotation class or all classes. Overlapping memberships and known False values retain their meaning.')
        controls.addWidget(self.class_choice, 1)
        layout.addLayout(controls)
        self.summary = QtWidgets.QLabel()
        self.summary.setWordWrap(True)
        self.summary.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        layout.addWidget(self.summary)
        tabs = QtWidgets.QTabWidget()
        membership = QtWidgets.QWidget()
        members_layout = QtWidgets.QVBoxLayout(membership)
        self.class_table = _table(['Class', 'Description', 'Known genes', 'Fraction', 'Unknown genes'])
        for i, width in enumerate((160, 300, 120, 120, 120)):
            self.class_table.setColumnWidth(i, width)
        self.class_table.cellClicked.connect(self._class_clicked)
        members_layout.addWidget(self.class_table)
        self.evidence_card = ScorecardBrowser()
        members_layout.addWidget(self.evidence_card, 1)
        source_controls = QtWidgets.QHBoxLayout()
        self.source_choice = QtWidgets.QComboBox()
        self.source_choice.setToolTip('Choose an explicitly recorded source identifier for this label; an unresolved identifier supplies no source route.')
        source_controls.addWidget(self.source_choice, 1)
        self.source_open = QtWidgets.QPushButton('Open source')
        self.source_open.setToolTip('Open the selected recorded source identifier in the dataset browser, retaining the current organism.')
        self.source_open.clicked.connect(self._open_source)
        source_controls.addWidget(self.source_open)
        members_layout.addLayout(source_controls)
        self.member_table = _table([])
        self.member_table.cellDoubleClicked.connect(self._open_gene)
        members_layout.addWidget(self.member_table, 1)
        paging = QtWidgets.QHBoxLayout()
        self.page_previous, self.page_next = QtWidgets.QPushButton('Previous'), QtWidgets.QPushButton('Next')
        self.page_previous.setToolTip('Show the previous page of original membership records, including overlapping class memberships.')
        self.page_next.setToolTip('Show the next page; every recorded gene membership remains reachable.')
        self.page_note = QtWidgets.QLabel()
        self.page_previous.clicked.connect(lambda: self._change_page(-1))
        self.page_next.clicked.connect(lambda: self._change_page(1))
        paging.addWidget(self.page_previous)
        paging.addWidget(self.page_note, 1)
        paging.addWidget(self.page_next)
        members_layout.addLayout(paging)
        tabs.addTab(membership, 'Membership, definitions and sources')
        mechanisms = QtWidgets.QWidget()
        mechanisms_layout = QtWidgets.QVBoxLayout(mechanisms)
        self.strategy_choice = QtWidgets.QComboBox()
        self.strategy_choice.setToolTip('Inspect every registered strategy and its applicability or recorded evaluation gaps for the selected label/class.')
        self.evaluation_choice = QtWidgets.QComboBox()
        self.evaluation_choice.setToolTip('Choose one exact recorded evaluation scope. Setting, seed, partition and mode are preserved; unrelated cohorts are never combined.')
        mechanisms_layout.addWidget(self.strategy_choice)
        mechanisms_layout.addWidget(self.evaluation_choice)
        self.strategy_note = QtWidgets.QLabel()
        self.strategy_note.setWordWrap(True)
        self.strategy_note.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        mechanisms_layout.addWidget(self.strategy_note)
        self.strategy_card = ScorecardBrowser()
        mechanisms_layout.addWidget(self.strategy_card, 1)
        self.evaluation_details = QtWidgets.QPlainTextEdit()
        self.evaluation_details.setReadOnly(True)
        self.evaluation_details.setMaximumHeight(180)
        self.evaluation_details.setToolTip('Exact supplied class precision, recall, F1 and confusion records for this one evaluation scope. Missing fields remain unavailable; these records do not establish independent biological accuracy.')
        mechanisms_layout.addWidget(self.evaluation_details)
        tabs.addTab(mechanisms, 'Strategies and recorded tests')
        layout.addWidget(tabs, 1)
        self.search.textChanged.connect(self._filter_labels)
        self.family.currentIndexChanged.connect(self._filter_labels)
        self.target_choice.currentIndexChanged.connect(self._target_selected)
        self.class_choice.currentIndexChanged.connect(self._class_selected)
        self.strategy_choice.currentIndexChanged.connect(self._strategy_selected)
        self.evaluation_choice.currentIndexChanged.connect(self._evaluation_selected)
        self._filter_labels()

    def select(self, target, value=None):
        """Select an exact backend label and original class key for application routing.

        A valid route clears only filters that would hide its target. Native class
        keys come from the backend, so False-like string keys are never treated
        as absence. Unknown labels or classes raise an explicit error instead of
        falling back to another organism, a related class or a different source.
        """
        if target not in set(self.model.labels.target):
            raise ValueError('Label is unavailable in this organism')
        classes = self.model.classes(target)
        if value is not None and value not in classes.value.tolist():
            raise ValueError('Native class is unavailable for this label')
        self._updating = True
        self.search.clear()
        self.family.setCurrentIndex(0)
        self._updating = False
        self._filter_labels()
        if self.target_choice.currentData() != target:
            self.target_choice.setCurrentIndex(self.target_choice.findData(target))
        if self.class_choice.currentData() != value:
            self.class_choice.setCurrentIndex(self.class_choice.findData(value))

    def _filter_labels(self, *_):
        if self._updating:
            return
        current = self.target_choice.currentData()
        frame = self.model.labels.copy()
        if self.family.currentData() is not None:
            frame = frame[frame.family.eq(self.family.currentData())]
        text = self.search.text().strip().casefold()
        if text:
            frame = frame[frame.astype(str).apply(lambda column: column.str.casefold().str.contains(text, regex=False)).any(axis=1)]
        rows = sorted(frame.to_dict('records'), key=lambda row: (row['family'] != 'Function', row['family'], row['target']))
        self.target_choice.blockSignals(True)
        self.target_choice.clear()
        for row in rows:
            self.target_choice.addItem(row.get('title') or row['target'], row['target'])
            self.target_choice.setItemData(self.target_choice.count()-1, str(row), QtCore.Qt.ItemDataRole.ToolTipRole)
        index = self.target_choice.findData(current)
        self.target_choice.setCurrentIndex(max(0, index))
        self.target_choice.blockSignals(False)
        self._target_selected()

    def _target_selected(self, *_):
        if self._updating:
            return
        self.target = self.target_choice.currentData()
        self.class_choice.blockSignals(True)
        self.class_choice.clear()
        self.class_choice.addItem('All native classes', None)
        self.class_frame = self.model.classes(self.target) if self.target is not None else pd.DataFrame()
        for row in self.class_frame.to_dict('records'):
            self.class_choice.addItem(f"{row['value']} · {_text(row.get('description'))} ({row['annotated_genes']} genes)", row['value'])
            self.class_choice.setItemData(self.class_choice.count()-1, str(row), QtCore.Qt.ItemDataRole.ToolTipRole)
        self.class_choice.setEnabled(self.target is not None)
        self.class_choice.blockSignals(False)
        columns = [key for key in ('value', 'description', 'annotated_genes', 'class_fraction', 'unknown_genes') if key in self.class_frame]
        _fill(self.class_table, self.class_frame[columns])
        self._class_selected()

    def _class_selected(self, *_):
        if self._updating:
            return
        self.value = self.class_choice.currentData()
        self.page_start = 0
        if self.target is None:
            self.summary.setText('No label matches the selected filters; no membership or evaluation substituted.')
            self.evidence_card.clear()
            self.member_frame = pd.DataFrame()
            self.mechanism_records = []
            sources = []
        else:
            label = self.model.labels[self.model.labels.target.eq(self.target)].iloc[0]
            self.summary.setText(f"{self.organism} · {self.target}: {label['annotated_genes']} known annotated genes / "
                f"{label['genes']} genes; {label['unannotated_genes']} unknown. Native class memberships may overlap. "
                + str(label.get('truth_interpretation', 'Annotation coverage is not prediction accuracy.')))
            self.evidence_card.set_scorecard(self.model.evidence_card(self.target, self.value))
            self.member_frame = self.model.members(self.target, self.value)
            self.mechanism_records = self.model.mechanisms(self.target, self.value)
            declared = label.get('source_ids', '')
            sources = declared.split(';') if isinstance(declared, str) else list(declared or [])
            sources = [str(source).strip() for source in sources if str(source).strip().casefold() not in {'', 'unknown', 'unresolved', 'unavailable'}]
        self.source_choice.clear()
        for source in sources:
            self.source_choice.addItem(source, source)
        self.source_choice.setEnabled(bool(sources))
        self.source_open.setEnabled(bool(sources))
        self._page()
        self.strategy_choice.blockSignals(True)
        self.strategy_choice.clear()
        for record in self.mechanism_records:
            self.strategy_choice.addItem(f"{record.get('title') or record['strategy']} · {record['status']}", record)
            self.strategy_choice.setItemData(self.strategy_choice.count()-1,
                str({key: record.get(key) for key in ('strategy', 'applicable', 'status', 'gaps')}), QtCore.Qt.ItemDataRole.ToolTipRole)
        self.strategy_choice.blockSignals(False)
        self._strategy_selected()

    def _class_clicked(self, row, _column):
        if 0 <= row < len(self.class_frame):
            self.class_choice.setCurrentIndex(self.class_choice.findData(self.class_frame.iloc[row]['value']))

    def _page(self):
        self.page_frame = self.member_frame.iloc[self.page_start:self.page_start + self.page_size]
        _fill(self.member_table, self.page_frame)
        if 'gene_id' in self.page_frame:
            self.member_table.setColumnWidth(self.page_frame.columns.get_loc('gene_id'), 170)
        total = len(self.member_frame)
        self.page_previous.setEnabled(self.page_start > 0)
        self.page_next.setEnabled(self.page_start + len(self.page_frame) < total)
        self.page_note.setText(f'{self.page_start + 1 if len(self.page_frame) else 0}–{self.page_start + len(self.page_frame)} of {total} membership records')

    def _change_page(self, direction):
        start = self.page_start + direction * self.page_size
        if 0 <= start < len(self.member_frame):
            self.page_start = start
            self._page()

    def _open_gene(self, row, _column):
        if not 0 <= row < len(self.page_frame) or 'gene_id' not in self.page_frame:
            return
        record = self.page_frame.iloc[row]
        if record.get('organism', self.organism) == self.organism and isinstance(record.gene_id, str) and record.gene_id:
            self.gene_requested.emit(self.organism, record.gene_id)

    def _open_source(self):
        source = self.source_choice.currentData()
        if source is not None:
            self.source_requested.emit(source)

    def _strategy_selected(self, *_):
        record = self.strategy_choice.currentData()
        self.evaluations = list(record.get('evaluations', [])) if record else []
        self.evaluation_choice.blockSignals(True)
        self.evaluation_choice.clear()
        for index, evaluation in enumerate(self.evaluations):
            self.evaluation_choice.addItem(str(evaluation['scope']), index)
            self.evaluation_choice.setItemData(index, str(evaluation['scope']), QtCore.Qt.ItemDataRole.ToolTipRole)
        self.evaluation_choice.setEnabled(bool(self.evaluations))
        self.evaluation_choice.blockSignals(False)
        self._evaluation_selected()

    def _evaluation_selected(self, *_):
        record = self.strategy_choice.currentData()
        index = self.evaluation_choice.currentData()
        evaluation = self.evaluations[index] if index is not None else None
        card = evaluation.get('card') if evaluation else record.get('card') if record else None
        gaps = evaluation.get('gaps', []) if evaluation else record.get('gaps', []) if record else []
        status = record['status'] if record else 'Unavailable'
        self.strategy_note.setText(status + (' · ' + ' / '.join(map(str, gaps)) if gaps else '')
                                   + (' · No recorded evaluation; metrics unavailable.' if card is None else ''))
        details = {'Evaluation scope': evaluation['scope'],
                   'Class precision, recall and F1': evaluation.get('class_metrics', 'Unavailable'),
                   'Confusion records': evaluation.get('confusion', 'Unavailable')} if evaluation else None
        self.evaluation_details.setPlainText(json.dumps(details, ensure_ascii=False, indent=2, allow_nan=False)
                                            if details is not None else 'Class metrics and confusions unavailable: no selected recorded evaluation.')
        if card is None:
            self.strategy_card.clear()
        else:
            self.strategy_card.set_scorecard(card)
