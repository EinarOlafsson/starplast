"""Present original evidence for one resolved gene in an explicitly qualified space.

The condensed view begins with twelve nonmissing quantities and can reveal every
stored column, including missing values. Source choices retain their full inventory
scope and open shared evidence cards. Navigation signals refer only to supplied
sources and suitable categorical labels; the widget does not fit models, download
data, infer biological accuracy or replace missing measurements with negatives.
"""
from copy import deepcopy

import pandas as pd
from PyQt6 import QtCore, QtWidgets

from .gene_evidence import evidence_rows
from .scorecard_browser import ScorecardBrowser


def _missing(value):
    if value is None or value is pd.NA or value is pd.NaT:
        return True
    missing = pd.isna(value)
    return bool(missing) if pd.api.types.is_scalar(missing) else False


def _text(value):
    return 'Unknown / missing' if _missing(value) else str(value)


def _categories(nodes):
    denied = {'gene_id', 'product', 'sequence', 'symbol', 'gene_name', 'name', 'description',
              'orthogroup', 'alphafold_accession'}
    result = set()
    for column in nodes:
        if column in denied or str(column).endswith(('_desc', '_description', '_accession')):
            continue
        values = nodes[column].dropna().tolist()
        if not values or not all(isinstance(value, str) for value in values):
            continue
        classes = set(values)
        if (2 <= len(classes) <= 60 and max(map(len, classes)) <= 120
                and (len(classes) < len(values) or isinstance(nodes[column].dtype, pd.CategoricalDtype))):
            result.add(column)
    return result


class GeneEvidenceBrowser(QtWidgets.QWidget):
    """Browse exact evidence and declared source context for one organism-qualified gene.

    The initial table keeps a concise set of nonmissing values, with an explicit
    expansion showing every original column. Selecting an evidence row exposes
    each attributable source through the shared scorecard renderer. Optional
    source, label and class navigation is emitted to the parent application only
    when a corresponding recorded route is available; unknown measurements stay
    distinct from stored zero and False values.
    """

    source_requested = QtCore.pyqtSignal(object)
    label_requested = QtCore.pyqtSignal(str)
    class_requested = QtCore.pyqtSignal(str, str)

    def __init__(self, nodes, organism, gene_id, space, parent=None):
        super().__init__(parent)
        self.organism, self.gene_id, self.space = organism, gene_id, space
        self.unavailable = ''
        try:
            supplied_rows = evidence_rows(nodes, organism, gene_id, space)
        except ValueError as error:
            if (str(error) != 'Select one exact canonical gene; missing or ambiguous aliases are not canonical IDs'
                    or not isinstance(nodes, pd.DataFrame) or 'gene_id' not in nodes
                    or gene_id in nodes.gene_id.tolist()):
                raise
            supplied_rows = []
            self.unavailable = 'Requested target is unavailable in this organism’s current gene table; no evidence values were substituted.'
        self.rows = sorted(supplied_rows,
                           key=lambda row: (str(row['question']), row['column']))
        self.label_columns = _categories(nodes)
        self.selected_row = None
        self.shown_rows = []
        layout = QtWidgets.QVBoxLayout(self)
        self.summary = QtWidgets.QLabel()
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        self.show_all = QtWidgets.QCheckBox('Show all evidence (including missing)')
        self.show_all.setToolTip('Reveal every original evidence column, including unknown or missing values; zero and False are retained values.')
        self.show_all.toggled.connect(self._refresh)
        layout.addWidget(self.show_all)
        self.table = QtWidgets.QTableWidget(0, 7)
        self.table.setHorizontalHeaderLabels(['Question', 'Column', 'Original value', 'Context', 'Unit', 'Origin', 'Recorded gaps'])
        for column, width in enumerate((240, 220, 180, 240, 100, 180, 320)):
            self.table.setColumnWidth(column, width)
            self.table.horizontalHeaderItem(column).setToolTip(self.table.horizontalHeaderItem(column).text())
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.itemSelectionChanged.connect(self._selected)
        layout.addWidget(self.table, 1)
        controls = QtWidgets.QHBoxLayout()
        self.source_choice = QtWidgets.QComboBox()
        self.source_choice.setToolTip('Choose an attributable source for this exact organism, storage unit and evidence column to inspect its recorded provenance and coverage.')
        self.source_choice.currentIndexChanged.connect(self._source_selected)
        controls.addWidget(self.source_choice, 1)
        self.source_open = QtWidgets.QPushButton('Open source in datasets')
        self.source_open.setToolTip('Navigate to the selected source scope in the dataset browser; the source card preserves its recorded storage coverage.')
        self.source_open.clicked.connect(self._open_source)
        controls.addWidget(self.source_open)
        self.label_open = QtWidgets.QPushButton('Open label')
        self.label_open.setToolTip('Browse the selected categorical annotation column. Numeric quantities, identifiers and free text do not offer a label route.')
        self.label_open.clicked.connect(lambda: self.label_requested.emit(self.selected_row['column']))
        controls.addWidget(self.label_open)
        self.class_open = QtWidgets.QPushButton('Open class')
        self.class_open.setToolTip('Browse the recorded string class for the selected categorical label; missing values do not supply a class.')
        self.class_open.clicked.connect(lambda: self.class_requested.emit(self.selected_row['column'], self.selected_row['value']))
        controls.addWidget(self.class_open)
        layout.addLayout(controls)
        self.card = ScorecardBrowser()
        layout.addWidget(self.card, 1)
        self._refresh()

    def _refresh(self, *_):
        selected = self.selected_row['column'] if self.selected_row else None
        self.shown_rows = self.rows if self.show_all.isChecked() else [row for row in self.rows if not _missing(row['value'])][:12]
        self.table.blockSignals(True)
        self.table.clearContents()
        self.table.setRowCount(len(self.shown_rows))
        for i, row in enumerate(self.shown_rows):
            for j, key in enumerate(('question', 'column', 'value', 'contexts', 'quantity_unit', 'origins', 'gaps')):
                value = row[key]
                item = QtWidgets.QTableWidgetItem(_text(value))
                item.setToolTip(repr(value) if key == 'value' else _text(value))
                self.table.setItem(i, j, item)
        self.table.blockSignals(False)
        self.summary.setText(f'{self.organism} · {self.gene_id}: {len(self.shown_rows)} of {len(self.rows)} evidence columns shown. '
                             + (self.unavailable or 'Values are original records; missing is unknown, not a verified negative.'))
        if self.shown_rows:
            index = next((i for i, row in enumerate(self.shown_rows) if row['column'] == selected), 0)
            self.table.setCurrentCell(index, 0)
        else:
            self.table.setCurrentCell(-1, -1)
        self._selected()

    def _selected(self):
        index = self.table.currentRow()
        self.selected_row = self.shown_rows[index] if 0 <= index < len(self.shown_rows) else None
        self.source_choice.blockSignals(True)
        self.source_choice.clear()
        for row in self.selected_row['sources'] if self.selected_row else ():
            title = row.get('title') or row.get('name') or row['source_id']
            self.source_choice.addItem(title, deepcopy(row))
        self.source_choice.blockSignals(False)
        label = self.selected_row is not None and self.selected_row['column'] in self.label_columns
        self.label_open.setEnabled(label)
        self.class_open.setEnabled(label and isinstance(self.selected_row['value'], str) and bool(self.selected_row['value']))
        self.source_choice.setEnabled(self.source_choice.count() > 0)
        self._source_selected()

    def _source_selected(self, *_):
        row = self.source_choice.currentData()
        self.source_open.setEnabled(row is not None)
        if row is None:
            self.card.clear()
        else:
            self.card.set_scorecard(self.space.source_card(row))

    def _open_source(self):
        row = self.source_choice.currentData()
        if row is not None:
            self.source_requested.emit(deepcopy(row))
