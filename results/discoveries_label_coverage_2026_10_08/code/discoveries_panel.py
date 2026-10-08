"""Browse functional and other labels, their known members and precomputed claims.

Annotation membership is separate from prediction. Missing inference and validation
remain visible. The claims view preserves its legacy calibration and visible filters;
those legacy tests do not substitute for independently admitted biological benchmarks.
Functional classes retain overlapping memberships and organism-specific source
descriptions. Selecting a class opens its known genes; only recorded evaluations
enable scorecard navigation, so an annotation never becomes an accuracy estimate.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from PyQt6 import QtCore, QtWidgets

from . import claims as C
from . import discovery_labels as DL
from . import strategies as S
from . import theme as TH
from . import track_record as T

TIPS = {
    "label": "Browse functional domains, enzyme classes, phenotypes, stages, localization and other "
             "available labels. Known annotation and inferred claims have separate views; an absent "
             "claim or test remains visible.",
    "status": "Tested: a legacy verifier passed the recipe's independence screen and agreed or disagreed. Untested: nothing independent "
              "reaches the gene, so only the generator's certainty is known. Outside the tested "
              "range: unlike every gene certainty was measured on, so no number is attached.",
    "confidence": "Of held-out genes given claims this confident, at least this share were right. "
                  "Measured, not estimated: certainty is fit on genes whose answer was hidden.",
    "lift": "Confidence over how common the claimed class is among labelled genes. A claim of a class "
            "that is 80% of genes at 85% confidence says little; lift 2 or more says something.",
    "table": "Click a claim to open it in the evidence panel: what was claimed, by what, how it was "
             "tested and what this recipe's claims were worth on held-out genes.",
    "export": "Save the claims shown, with every number and verdict, as a table.",
    "colour": "Colour the map by this label: annotated genes in full colour, claimed genes in their "
              "claimed class faded by how uncertain the claim is, unclaimed genes grey.",
    "recipes": "Why each label has the claims it has: the recipe chosen, whether it is calibrated end "
               "to end, and what its claims were worth on held-out genes.",
}

STATUSES = ("tested", "untested", "outside tested range")
COLUMNS = ("gene_id", "product", "label", "claim", "confidence", "lift", "checks", "status")


class DiscoveriesPanel(QtWidgets.QWidget):
    """All available labels and classes, with a separate filtered claim view."""

    gene_chosen = QtCore.pyqtSignal(str)        # a gene id, for the evidence panel
    claim_chosen = QtCore.pyqtSignal(str, str)  # (gene id, label), for the claim's reasoning
    colour_by = QtCore.pyqtSignal(str)          # a label to colour the map by
    record_chosen = QtCore.pyqtSignal(str,str)  # target, optional class in the legacy scorecard

    def __init__(self, organism: str, products=None, parent=None, context=None):
        super().__init__(parent)
        self.organism = organism
        self.products = products or {}
        self.frame = C.shipped(organism)
        self.context = context or S.Context.shipped(organism)
        self.inventory, self.memberships = DL.catalogue(self.context,self.frame,C.recipes(organism),T.shipped(organism))
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(8, 8, 8, 8)
        self.summary = QtWidgets.QLabel("")
        self.summary.setWordWrap(True)

        row = QtWidgets.QHBoxLayout()
        self.label = QtWidgets.QComboBox()
        self.label.addItem("every label", "")
        for r in self.inventory.itertuples(index=False):
            self.label.addItem(f'{r.family}: {r.title}',r.target)
        self.label.setToolTip(TH.tip(TIPS["label"]))
        row.addWidget(self.label, 1)
        _strip(lay, row)
        status_row = QtWidgets.QHBoxLayout()
        self.status = {}
        for s in STATUSES:
            box = QtWidgets.QCheckBox(s)
            box.setChecked(s == "tested")
            box.setToolTip(TH.tip(TIPS["status"]))
            self.status[s] = box
            status_row.addWidget(box)

        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("confidence ≥"))
        self.confidence = QtWidgets.QDoubleSpinBox()
        self.confidence.setRange(0.0, 1.0)
        self.confidence.setSingleStep(0.05)
        self.confidence.setValue(0.8)
        self.confidence.setToolTip(TH.tip(TIPS["confidence"]))
        row.addWidget(self.confidence)
        row.addWidget(QtWidgets.QLabel("lift ≥"))
        self.lift = QtWidgets.QDoubleSpinBox()
        self.lift.setRange(0.0, 50.0)
        self.lift.setSingleStep(0.5)
        self.lift.setValue(2.0)
        self.lift.setToolTip(TH.tip(TIPS["lift"]))
        row.addWidget(self.lift)
        row.addStretch(1)
        self.colour = QtWidgets.QPushButton("Colour map")
        self.colour.setToolTip(TH.tip(TIPS["colour"]))
        row.addWidget(self.colour)
        self.export = QtWidgets.QPushButton("Save…")
        self.export.setToolTip(TH.tip(TIPS["export"]))
        row.addWidget(self.export)
        claim_controls = row

        self.tabs = QtWidgets.QTabWidget()
        self.annotations_page = QtWidgets.QWidget()
        annotations_lay = QtWidgets.QVBoxLayout(self.annotations_page)
        self.annotation_search = QtWidgets.QLineEdit()
        self.annotation_search.setPlaceholderText('Search functions, domains, enzyme classes or labels…')
        annotations_lay.addWidget(self.annotation_search)
        self.label_table = QtWidgets.QTableWidget(0,6)
        self.label_table.setHorizontalHeaderLabels(['label / function','kind','annotated genes','classes','inferred claims','evaluation'])
        self.label_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.label_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.label_table.verticalHeader().setVisible(False)
        annotations_lay.addWidget(self.label_table,1)
        self.annotation_note = QtWidgets.QLabel('')
        self.annotation_note.setWordWrap(True)
        annotations_lay.addWidget(self.annotation_note)
        self.annotation_class = QtWidgets.QComboBox()
        self.annotation_class.addItem('Choose a label above to browse its classes','')
        annotations_lay.addWidget(self.annotation_class)
        self.annotation_record = QtWidgets.QPushButton('Open held-out scorecard')
        self.annotation_record.setToolTip(TH.tip('Open existing held-out label/class results. Legacy source-recovery tests have their own scope; they are not newly admitted biological accuracy.'))
        annotations_lay.addWidget(self.annotation_record)
        self.member_table = QtWidgets.QTableWidget(0,3)
        self.member_table.setHorizontalHeaderLabels(['gene','product','known annotation'])
        self.member_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.member_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.member_table.verticalHeader().setVisible(False)
        annotations_lay.addWidget(self.member_table,1)
        self.tabs.addTab(self.annotations_page,'Labels and functions')
        self.claims_page = QtWidgets.QWidget()
        claims_lay = QtWidgets.QVBoxLayout(self.claims_page)
        claims_lay.addWidget(self.summary)
        _strip(claims_lay,status_row)
        _strip(claims_lay,claim_controls)
        self.table = QtWidgets.QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(list(COLUMNS))
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setToolTip(TH.tip(TIPS["table"]))
        self.table.setMinimumWidth(120)
        claims_lay.addWidget(self.table,1)
        self.tabs.addTab(self.claims_page,'Inferred claims')
        lay.addWidget(self.tabs,1)

        self.recipes = QtWidgets.QLabel("")
        self.recipes.setWordWrap(True)
        self.recipes.setToolTip(TH.tip(TIPS["recipes"]))
        claims_lay.addWidget(self.recipes)

        self.label.currentIndexChanged.connect(self.refresh)
        for box in self.status.values():
            box.toggled.connect(self.refresh)
        self.confidence.valueChanged.connect(self.refresh)
        self.lift.valueChanged.connect(self.refresh)
        self.table.cellClicked.connect(self._clicked)
        self.colour.clicked.connect(self._colour)
        self.export.clicked.connect(lambda: self.save())
        self.annotation_search.textChanged.connect(self._browse_labels)
        self.annotation_search.textChanged.connect(self._annotation_classes)
        self.label_table.cellClicked.connect(self._select_annotation_label)
        self.annotation_class.currentIndexChanged.connect(self._browse_members)
        self.member_table.cellClicked.connect(self._member_clicked)
        self.annotation_record.clicked.connect(self._open_annotation_record)
        self.shown = pd.DataFrame()
        self.refresh()
        self._browse_labels()

    # ------------------------------------------------------------------ state
    def selection(self) -> pd.DataFrame:
        """The claims the filters leave, most confident first."""
        f = self.frame
        if not len(f):
            return f
        status = [s for s, box in self.status.items() if box.isChecked()]
        keep = f["status"].astype(str).isin(status)
        target = self.label.currentData()
        if target:
            keep &= f["target"].astype(str) == target
        conf = f["confidence"]
        # "Outside tested range" has no number by design; a confidence floor must not silently drop it.
        numeric = conf.notna()
        keep &= (~numeric) | ((conf >= self.confidence.value()) & (f["lift"] >= self.lift.value()))
        return f[keep].sort_values("confidence", ascending=False, na_position="last")

    def refresh(self):
        """Re-apply the filters: the summary line, the table and the recipe note."""
        self.shown = self.selection()
        total = len(self.frame)
        tested = int((self.frame["status"].astype(str) == "tested").sum()) if total else 0
        self.summary.setText(
            f"<b>{len(self.shown):,}</b> claims shown, of {total:,} about genes with no label "
            f"({tested:,} tested by legacy verifiers)." if total else
            "No claims have been built for this organism.")
        self._fill()
        self._describe_recipes()
        self._annotation_classes()

    def _browse_labels(self):
        """Search variable names and individual functional terms without requiring claims."""
        query = self.annotation_search.text().strip().casefold()
        view = self.inventory
        if query:
            hits = self.memberships[self.memberships[['value','description']].astype(str).apply(
                lambda c:c.str.casefold().str.contains(query,regex=False)).any(axis=1)].target
            keep = view[['title','family','target']].astype(str).apply(
                lambda c:c.str.casefold().str.contains(query,regex=False)).any(axis=1)|view.target.isin(hits)
            view = view[keep]
        self.browsed_labels = view.reset_index(drop=True)
        self.label_table.setRowCount(len(view))
        for i,r in enumerate(view.itertuples(index=False)):
            values = [r.title,r.family,f'{r.annotated_genes:,} / {r.genes:,}',str(r.classes),str(r.claims),r.evaluation_status]
            for j,value in enumerate(values):self.label_table.setItem(i,j,QtWidgets.QTableWidgetItem(value))
        self.label_table.resizeColumnsToContents()

    def _select_annotation_label(self, row, _column):
        """Keep the selected organism/label address while opening its classes and claims."""
        if 0<=row<len(self.browsed_labels):
            target = self.browsed_labels.iloc[row].target
            self.label.setCurrentIndex(self.label.findData(target))

    def _annotation_classes(self):
        """Show measured membership and explicit evaluation gaps for the selected variable."""
        target = self.label.currentData()
        self.annotation_class.blockSignals(True)
        self.annotation_class.clear()
        self.annotation_class.addItem('All known classes','')
        classes = DL.class_summary(self.memberships,target) if target else pd.DataFrame()
        search = self.annotation_search.text().strip().casefold()
        for r in classes.itertuples(index=False):
            text = f'{r.value} — {r.description}' if r.description else str(r.value)
            if not search or search in text.casefold() or search in str(target).replace('_',' ').casefold():
                self.annotation_class.addItem(f'{text} ({r.annotated_genes:,} genes)',r.value)
        self.annotation_class.blockSignals(False)
        selected = self.inventory[self.inventory.target.eq(target)] if target else pd.DataFrame()
        self.annotation_record.setEnabled(bool(len(selected) and selected.held_out_rows.iloc[0]>0))
        if target and len(self.inventory[self.inventory.target.eq(target)]):
            r = self.inventory[self.inventory.target.eq(target)].iloc[0]
            self._annotation_note_base = (f'{r.annotated_genes:,} genes with known annotation; {r.unannotated_genes:,} without annotation. '
                f'{r.claims:,} precomputed inferred claims. {r.evaluation_status}. '
                f'{r.held_out_rows:,} legacy held-out rows across {r.evaluated_strategies} strategies. '
                f'Source: {r.source_ids}. Independent biological class precision/recall: not evaluated. '
                'Known annotation is separate from inference; missing annotation is unknown membership. '
                'False/0 annotation flags do not establish biological absence.')
        else:
            self._annotation_note_base = (f'{len(self.inventory)} available labels across function, phenotype, stage, localization and other categories. '
                'Select a label to inspect its known classes and genes. Inferred claims and their filters are in the next tab.')
        self._browse_members()

    def _browse_members(self):
        """Show gene membership in one class, retaining multi-valued functional annotations."""
        target,value = self.label.currentData(),self.annotation_class.currentData()
        rows = self.memberships[self.memberships.target.eq(target)].copy()
        if value:rows=rows[rows.value.eq(value)]
        search = self.annotation_search.text().strip().casefold()
        if search and target and search not in str(target).replace('_',' ').casefold():
            rows = rows[rows[['value','description']].astype(str).apply(lambda c:c.str.casefold().str.contains(search,regex=False)).any(axis=1)]
        self.shown_members = rows.reset_index(drop=True)
        shown = rows.head(2000)
        self.member_table.setRowCount(len(shown))
        for i,r in enumerate(shown.itertuples(index=False)):
            for j,text in enumerate((r.gene_id,self.products.get(r.gene_id,''),r.value)):
                self.member_table.setItem(i,j,QtWidgets.QTableWidgetItem(text))
        self.member_table.resizeColumnsToContents()
        self.annotation_note.setText(self._annotation_note_base)
        if len(rows)>2000:
            self.annotation_note.setText(self.annotation_note.text()+' First 2,000 membership rows shown; select a class or refine search.')

    def _member_clicked(self, row, _column):
        """Open a known member's gene without inventing an inferred claim for it."""
        if 0<=row<min(2000,len(self.shown_members)):
            self.gene_chosen.emit(str(self.shown_members.iloc[row].gene_id))

    def _open_annotation_record(self):
        """Link the selected label/class to existing held-out results when they exist."""
        if self.annotation_record.isEnabled():
            self.record_chosen.emit(str(self.label.currentData()),str(self.annotation_class.currentData() or ''))

    def _fill(self):
        rows = self.shown.head(2000)
        verdicts = [c for c in rows.columns if str(c).endswith(" verdict")]
        self.table.setRowCount(len(rows))
        for i, (_,rr) in enumerate(rows.iterrows()):
            said = [str(rr[c]) for c in verdicts if pd.notna(rr.get(c)) and str(rr[c]) != "silent"]
            checks = (f"{said.count('agrees')} agree, {said.count('disagrees')} disagree"
                      if said else "—")
            conf = rr["confidence"]
            values = [str(rr["gene_id"]), self.products.get(str(rr["gene_id"]), ""),
                      str(rr["target"]).replace("_", " "), str(rr["claim"]),
                      "not measured" if not np.isfinite(conf) else f"{conf:.0%}",
                      "" if not np.isfinite(rr["lift"]) else f"{rr['lift']:.1f}×",
                      checks, str(rr["status"])]
            for j, v in enumerate(values):
                self.table.setItem(i, j, QtWidgets.QTableWidgetItem(v))
        self.table.resizeColumnsToContents()

    def _describe_recipes(self):
        rec = C.recipes(self.organism)
        target = self.label.currentData()
        if target and len(rec):
            rec = rec[rec["target"].astype(str) == target]
        if not len(rec):
            self.recipes.setText("")
            return
        parts = []
        for q in rec.itertuples(index=False):
            name = str(q.target).replace("_", " ")
            if not q.proven:
                parts.append(f"{name}: no claims -- its best recipe is not calibrated end to end")
            elif q.claims == 0:
                parts.append(f"{name}: no claims -- every gene already has this label")
            else:
                checks = q.verifiers or "no independent check"
                parts.append(f"{name}: {str(q.generator).replace('_', ' ')} checked by "
                             f"{str(checks).replace('_', ' ')}")
        self.recipes.setText("<span style='color:#888'>" + "; ".join(parts) + "</span>")

    # ------------------------------------------------------------------ actions
    def _clicked(self, i: int, _j: int):
        if 0 <= i < len(self.shown):
            r = self.shown.iloc[i]
            self.gene_chosen.emit(str(r["gene_id"]))
            self.claim_chosen.emit(str(r["gene_id"]), str(r["target"]))

    def _colour(self):
        target = self.label.currentData() or (
            str(self.shown["target"].iloc[0]) if len(self.shown) else "")
        if target:
            self.colour_by.emit(target)

    def save(self, path: str = "") -> int:
        """Write the claims shown to CSV or TSV. Returns how many rows were written."""
        if not path:
            path, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Save claims", "claims.csv",
                                                            "CSV (*.csv);;TSV (*.tsv)")
        if not path:
            return 0
        out = self.shown.assign(product=self.shown["gene_id"].astype(str).map(self.products))
        out.to_csv(path, index=False, sep="\t" if path.endswith(".tsv") else ",")
        return len(out)


def _strip(lay: QtWidgets.QVBoxLayout, row: QtWidgets.QHBoxLayout) -> None:
    """Add a control row in a horizontal scroll strip, so the panel can be narrower than the row.

    A row of filters laid out directly sets the panel's minimum width, and through the dock the
    window's: on an 800-pixel screen that was 804, and the window could not fit.
    """
    holder = QtWidgets.QWidget()
    holder.setLayout(row)
    strip = QtWidgets.QScrollArea()
    strip.setWidget(holder)
    strip.setWidgetResizable(True)
    strip.setFrameShape(QtWidgets.QFrame.Shape.NoFrame)
    strip.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
    strip.setFixedHeight(holder.sizeHint().height() + 2)
    lay.addWidget(strip)

def install(window) -> DiscoveriesPanel | None:
    """Add the Discoveries tab beside Strategies and wire it to the evidence panel and the map."""
    from . import organisms
    code = organisms.by_species(window.species).code
    products = dict(zip(window.nodes["gene_id"].astype(str),
                        window.nodes["product"].astype(str) if "product" in window.nodes
                        else [""] * len(window.nodes)))
    panel = DiscoveriesPanel(code,products,context=S.Context(window.nodes,graph={},organism=code))

    def show_gene(gene_id):
        """Open the selected annotation member in the existing organism gene view."""
        hits = window.nodes.index[window.nodes.gene_id.astype(str).eq(gene_id)]
        if len(hits):window.show_detail(int(hits[0]))
        window.right_dock.raise_()

    def show_claim(gene_id: str, target: str):
        ids = window.nodes["gene_id"].astype(str)
        hit = ids[ids == gene_id]
        if len(hit):
            window.show_detail(int(hit.index[0]))
        window._detail_link(QtCore.QUrl(f"starplast://claim/{gene_id}/{target}"))
        window.right_dock.raise_()

    panel.claim_chosen.connect(show_claim)
    panel.gene_chosen.connect(show_gene)
    def show_record(target,label):
        """Open the existing scorecard without reinterpreting its truth source or test scope."""
        from urllib.parse import quote
        address = (f"starplast://class/{quote(target,safe='')}/{quote(label,safe='')}" if label
            else f"starplast://target/{quote(target,safe='')}")
        window._detail_link(QtCore.QUrl(address))
        window.right_dock.raise_()
    panel.record_chosen.connect(show_record)
    panel.colour_by.connect(window.colour_by_claims)
    dock = QtWidgets.QDockWidget("discoveries", window)
    dock.setFeatures(QtWidgets.QDockWidget.DockWidgetFeature.NoDockWidgetFeatures)
    dock.setWidget(panel)
    window.addDockWidget(QtCore.Qt.DockWidgetArea.RightDockWidgetArea, dock)
    anchor = getattr(window, "guided_dock", None) or getattr(window, "strategies_dock", None)
    if anchor is not None:
        window.tabifyDockWidget(anchor, dock)
    if getattr(window, "right_dock", None) is not None:
        window.right_dock.raise_()
    window.discoveries, window.discoveries_dock = panel, dock
    return panel
