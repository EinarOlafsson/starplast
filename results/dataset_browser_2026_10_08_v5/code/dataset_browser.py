"""Browse supplied dataset inventory, provenance and original stored values.

The widget filters an existing DatasetSpace without acquiring sources, fitting
models or projecting measurements. Coverage remains storage coverage with its
declared population. Source cards use the shared evidence renderer, while exact
stored entity identifiers drive organism-qualified navigation. Refused questions
and unavailable data keep their own recorded reasons and cannot borrow another
source's values or biological accuracy.
"""
from __future__ import annotations

import pandas as pd
from PyQt6 import QtCore, QtWidgets

from .scorecard_browser import ScorecardBrowser


def _records(value):
    return value.to_dict('records') if isinstance(value, pd.DataFrame) else list(value)


def _scope(row):
    return tuple(row.get(key, '') for key in ('source_id', 'organism', 'unit', 'question'))


def _display(value):
    if value is None or value is pd.NA or value is pd.NaT:
        return 'Unavailable'
    missing = pd.isna(value)
    if pd.api.types.is_scalar(missing) and bool(missing):
        return 'Unavailable'
    if isinstance(value, (list, tuple)):
        return ' / '.join(str(item) for item in value) or 'Unavailable'
    return str(value)


class _ValueModel(QtCore.QAbstractTableModel):
    def __init__(self, frame, parent=None):
        super().__init__(parent)
        self.frame = frame

    def rowCount(self, parent=QtCore.QModelIndex()):
        return 0 if parent.isValid() else len(self.frame)

    def columnCount(self, parent=QtCore.QModelIndex()):
        return 0 if parent.isValid() else len(self.frame.columns)

    def data(self, index, role=QtCore.Qt.ItemDataRole.DisplayRole):
        if index.isValid() and role in (QtCore.Qt.ItemDataRole.DisplayRole, QtCore.Qt.ItemDataRole.ToolTipRole):
            return _display(self.frame.iat[index.row(), index.column()])
        return None

    def headerData(self, section, orientation, role=QtCore.Qt.ItemDataRole.DisplayRole):
        if role not in (QtCore.Qt.ItemDataRole.DisplayRole, QtCore.Qt.ItemDataRole.ToolTipRole):
            return None
        if orientation == QtCore.Qt.Orientation.Horizontal:
            return str(self.frame.columns[section]) if 0 <= section < len(self.frame.columns) else None
        return str(section + 1) if 0 <= section < len(self.frame) else None


class DatasetBrowser(QtWidgets.QWidget):
    """Filter one supplied evidence inventory and inspect each source's original values.

    Source summaries retain scope, storage counts, context and origin. Selecting a
    source opens the shared immutable evidence card and its entity table; missing
    or refused records show their recorded reasons. Navigation emits explicitly
    qualified stored identifiers, leaving application routing to the parent host.
    No data download, model fit or biological-performance estimate occurs here.
    """

    entity_requested = QtCore.pyqtSignal(str, str, str)

    def __init__(self, space, *, organism=None, parent=None):
        super().__init__(parent)
        self.space = None
        self.shown_rows = []
        self.selected_row = None
        self.entity_frame = pd.DataFrame()
        self.page_start = 0
        self.page_size = 200
        self._updating = False
        layout = QtWidgets.QVBoxLayout(self)
        controls = QtWidgets.QGridLayout()
        self.filters = {}
        for column, (key, label) in enumerate((('organism', 'Organism / host'), ('unit', 'Unit'),
                ('family', 'Evidence family'), ('context', 'Measured context'), ('status', 'Availability'))):
            combo = QtWidgets.QComboBox()
            combo.currentIndexChanged.connect(self._refresh)
            self.filters[key] = combo
            controls.addWidget(QtWidgets.QLabel(label), 0, column)
            controls.addWidget(combo, 1, column)
        self.search = QtWidgets.QLineEdit(placeholderText='Filter sources, labels, columns and recorded gaps…')
        self.search.textChanged.connect(self._refresh)
        controls.addWidget(self.search, 2, 0, 1, 5)
        layout.addLayout(controls)
        self.summary = QtWidgets.QLabel('')
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)
        split = QtWidgets.QSplitter(QtCore.Qt.Orientation.Vertical)
        self.source_table = QtWidgets.QTableWidget(0, 9)
        self.source_table.setHorizontalHeaderLabels(['Source', 'Organism', 'Unit', 'Origin',
            'Family', 'Context', 'Availability', 'Stored / table rows', 'Question / gap'])
        source_header = self.source_table.horizontalHeader()
        source_header.setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Interactive)
        source_header.setStretchLastSection(True)
        for column, width in enumerate((220, 60, 90, 150, 200, 220, 100, 150, 160)):
            self.source_table.setColumnWidth(column, width)
            header = self.source_table.horizontalHeaderItem(column)
            header.setToolTip(header.text())
        self.source_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.source_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.source_table.setSelectionMode(QtWidgets.QAbstractItemView.SelectionMode.SingleSelection)
        self.source_table.verticalHeader().setVisible(False)
        self.source_table.itemSelectionChanged.connect(self._show_selected)
        split.addWidget(self.source_table)
        details = QtWidgets.QSplitter(QtCore.Qt.Orientation.Horizontal)
        self.card = ScorecardBrowser()
        details.addWidget(self.card)
        value_page = QtWidgets.QWidget()
        value_layout = QtWidgets.QVBoxLayout(value_page)
        self.value_note = QtWidgets.QLabel('')
        self.value_note.setWordWrap(True)
        value_layout.addWidget(self.value_note)
        self.entity_table = QtWidgets.QTableView()
        self.entity_table.horizontalHeader().setSectionResizeMode(QtWidgets.QHeaderView.ResizeMode.Interactive)
        self.entity_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.entity_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.entity_table.doubleClicked.connect(self._open_entity)
        value_layout.addWidget(self.entity_table, 1)
        paging = QtWidgets.QHBoxLayout()
        self.page_previous = QtWidgets.QPushButton('Previous')
        self.page_next = QtWidgets.QPushButton('Next')
        self.page_note = QtWidgets.QLabel('')
        self.page_previous.clicked.connect(lambda: self._change_page(-1))
        self.page_next.clicked.connect(lambda: self._change_page(1))
        paging.addWidget(self.page_previous)
        paging.addWidget(self.page_note, 1)
        paging.addWidget(self.page_next)
        value_layout.addLayout(paging)
        details.addWidget(value_page)
        split.addWidget(details)
        layout.addWidget(split, 1)
        self.set_space(space, organism=organism)

    def set_space(self, space, *, organism=None):
        """Refresh the supplied DatasetSpace and preserve an applicable scoped selection.

        Facets and counts come from the new inventory. A requested organism takes
        precedence over the previous organism filter; other surviving filters and
        the complete source/organism/unit/question address are retained. Removed
        sources clear their previous card and values, preventing a stale result
        from being exported or used for entity navigation after a refresh.
        """
        if not all(callable(getattr(space, name, None)) for name in ('filter_rows', 'facets', 'source_card', 'entities')):
            raise TypeError('A DatasetSpace interface is required')
        previous = {key: combo.currentData() for key, combo in self.filters.items()}
        self.space = space
        facets = space.facets()
        self._updating = True
        for key, combo in self.filters.items():
            combo.clear()
            combo.addItem('All', None)
            for value in facets.get(key, ()):
                combo.addItem(str(value), value)
            desired = organism if key == 'organism' and organism is not None else previous[key]
            index = combo.findData(desired)
            if key == 'organism' and organism is not None and index < 0:
                combo.addItem(str(organism) + ' (unavailable)', organism)
                index = combo.count() - 1
            combo.setCurrentIndex(index if index >= 0 else 0)
        self._updating = False
        self._refresh()

    def _refresh(self, *_):
        if self._updating or self.space is None:
            return
        selected = _scope(self.selected_row) if self.selected_row else None
        filters = {key: combo.currentData() for key, combo in self.filters.items()}
        self.shown_rows = _records(self.space.filter_rows(**filters, text=self.search.text()))
        self._updating = True
        self.source_table.setRowCount(len(self.shown_rows))
        self.source_table.clearContents()
        for i, row in enumerate(self.shown_rows):
            counts = ' / '.join(_display(row.get(key)) for key in ('stored_any_rows', 'table_rows'))
            if row.get('unit') == 'pair':
                counts = _display(row.get('pair_records')) + ' pair records; denominator unavailable'
            values = (row.get('title') or row.get('name') or row['source_id'], row['organism'], row['unit'],
                row.get('origin'), row.get('evidence_families'), row.get('contexts'), row['status'], counts,
                ' · '.join(str(row.get(key) or '') for key in ('question', 'gap')).strip(' ·'))
            for j, value in enumerate(values):
                text = _display(value)
                label = 'Unverified registry' if j == 3 and text == 'registry_asserted_unverified' else text
                item = QtWidgets.QTableWidgetItem(label)
                item.setToolTip(text)
                self.source_table.setItem(i, j, item)
        self._updating = False
        self.summary.setText(f'{len(self.shown_rows):,} source scopes shown. Coverage describes stored values '
            'in the declared table; missing values remain unknown. Counts are not biological accuracy.')
        if self.shown_rows:
            index = next((i for i, row in enumerate(self.shown_rows) if _scope(row) == selected), 0)
            self.source_table.setCurrentCell(index, 0)
            self._show_selected()
        else:
            self.selected_row = None
            self.card.clear()
            self._set_entity_frame(pd.DataFrame())
            self.page_start = 0
            self.page_previous.setEnabled(False)
            self.page_next.setEnabled(False)
            self.page_note.setText('0 stored records')
            self.value_note.setText('No source matches the selected filters.')

    def _show_selected(self):
        if self._updating:
            return
        index = self.source_table.currentRow()
        if not 0 <= index < len(self.shown_rows):
            return
        self.selected_row = self.shown_rows[index]
        self.card.set_scorecard(self.space.source_card(self.selected_row))
        self.page_start = 0
        self._load_page()

    def _set_entity_frame(self, frame):
        previous = self.entity_table.model()
        self.entity_frame = frame.copy(deep=True)
        self.entity_table.setModel(_ValueModel(self.entity_frame, self.entity_table))
        if len(self.entity_frame.columns):
            self.entity_table.setColumnWidth(0, 170)
        if previous is not None:
            previous.deleteLater()

    def _load_page(self):
        frame = self.space.entities(self.selected_row, start=self.page_start, limit=self.page_size)
        self._set_entity_frame(frame)
        total = self.entity_frame.attrs['total_rows']
        self.page_previous.setEnabled(self.page_start > 0)
        self.page_next.setEnabled(total is not None and self.page_start + len(self.entity_frame) < total)
        self.page_note.setText('0 shown; stored record total unavailable' if total is None else
            f'{self.page_start + 1 if len(self.entity_frame) else 0:,}–'
            f'{self.page_start + len(self.entity_frame):,} of {total:,} stored records')
        gap = self.entity_frame.attrs.get('gap') or self.selected_row.get('gap')
        self.value_note.setText(str(gap) if gap else f'{total:,} original stored rows. '
            'Unavailable cells retain their original missingness; double-click a stored entity to open it.')

    def _change_page(self, direction):
        if self.selected_row is None:
            return
        total = self.entity_frame.attrs['total_rows']
        start = self.page_start + direction * self.page_size
        if total is not None and 0 <= start < total:
            self.page_start = start
            self._load_page()

    def _open_entity(self, index):
        if not index.isValid() or self.selected_row is None or not 0 <= index.row() < len(self.entity_frame):
            return
        attrs = self.entity_frame.attrs
        column = attrs.get('entity_column')
        if (not column or column not in self.entity_frame or attrs.get('status') == 'unavailable'
                or attrs.get('organism') != self.selected_row['organism']
                or attrs.get('unit') != self.selected_row['unit']):
            return
        entity = self.entity_frame.iloc[index.row()][column]
        if isinstance(entity, str) and entity:
            self.entity_requested.emit(attrs['organism'], attrs['unit'], entity)
