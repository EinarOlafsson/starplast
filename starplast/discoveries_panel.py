"""The Discoveries tab: what Starplast infers about genes nobody has labelled, and how each claim was
tested -- condensed to one line per claim, with every claim a click from its reasoning.

Read from the shipped claims (`claims.shipped`), so opening it costs a lookup. The defaults show only
claims worth a person's time: tested by evidence measured to be independent of what generated them,
confident, and well above how common the claimed class is anyway. Every filter is visible and can be
loosened; nothing is hidden that the person cannot bring back.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from PyQt6 import QtCore, QtWidgets

from . import claims as C
from . import theme as TH

TIPS = {
    "label": "Which kind of knowledge: a localization, a phenotype, a stage. Each label has its own "
             "recipe -- the generator, and the checks measured to be independent of it.",
    "status": "Tested: an independent check agreed or disagreed. Untested: nothing independent "
              "reaches the gene, so only the generator's certainty is known. Outside the tested "
              "range: unlike every gene certainty was measured on, so no number is attached.",
    "confidence": "Of held-out genes given claims this confident, at least this share were right. "
                  "Measured, not estimated: certainty is fit on genes whose answer was hidden.",
    "lift": "Confidence over how common the claimed class is among labelled genes. A claim of a class "
            "that is 80% of genes at 85% confidence says little; lift 2 or more says something.",
    "table": "Click a claim to open it in the evidence panel: what was claimed, by what, how it was "
             "tested and what this recipe's claims were worth on held-out genes.",
    "export": "Save the claims shown, with every number and verdict, as a table.",
    "colour": "Colour the map by this label: measured genes in full colour, claimed genes in their "
              "claimed class faded by how uncertain the claim is, unclaimed genes grey.",
    "recipes": "Why each label has the claims it has: the recipe chosen, whether it is calibrated end "
               "to end, and what its claims were worth on held-out genes.",
}

STATUSES = ("tested", "untested", "outside tested range")
COLUMNS = ("gene_id", "product", "label", "claim", "confidence", "lift", "checks", "status")


class DiscoveriesPanel(QtWidgets.QWidget):
    """Filters on top, one condensed line of what they leave, then the claims."""

    gene_chosen = QtCore.pyqtSignal(str)        # a gene id, for the evidence panel
    claim_chosen = QtCore.pyqtSignal(str, str)  # (gene id, label), for the claim's reasoning
    colour_by = QtCore.pyqtSignal(str)          # a label to colour the map by

    def __init__(self, organism: str, products=None, parent=None):
        super().__init__(parent)
        self.organism = organism
        self.products = products or {}
        self.frame = C.shipped(organism)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(8, 8, 8, 8)
        self.summary = QtWidgets.QLabel("")
        self.summary.setWordWrap(True)
        lay.addWidget(self.summary)

        row = QtWidgets.QHBoxLayout()
        self.label = QtWidgets.QComboBox()
        self.label.addItem("every label", "")
        targets = (sorted(set(self.frame["target"].astype(str))) if len(self.frame) else [])
        for t in targets:
            self.label.addItem(t.replace("_", " "), t)
        self.label.setToolTip(TH.tip(TIPS["label"]))
        row.addWidget(self.label, 1)
        self.status = {}
        for s in STATUSES:
            box = QtWidgets.QCheckBox(s)
            box.setChecked(s == "tested")
            box.setToolTip(TH.tip(TIPS["status"]))
            self.status[s] = box
            row.addWidget(box)
        lay.addLayout(row)

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
        lay.addLayout(row)

        self.table = QtWidgets.QTableWidget(0, len(COLUMNS))
        self.table.setHorizontalHeaderLabels(list(COLUMNS))
        self.table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setToolTip(TH.tip(TIPS["table"]))
        lay.addWidget(self.table, 1)

        self.recipes = QtWidgets.QLabel("")
        self.recipes.setWordWrap(True)
        self.recipes.setToolTip(TH.tip(TIPS["recipes"]))
        lay.addWidget(self.recipes)

        self.label.currentIndexChanged.connect(self.refresh)
        for box in self.status.values():
            box.toggled.connect(self.refresh)
        self.confidence.valueChanged.connect(self.refresh)
        self.lift.valueChanged.connect(self.refresh)
        self.table.cellClicked.connect(self._clicked)
        self.colour.clicked.connect(self._colour)
        self.export.clicked.connect(lambda: self.save())
        self.shown = pd.DataFrame()
        self.refresh()

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
        self.shown = self.selection()
        total = len(self.frame)
        tested = int((self.frame["status"].astype(str) == "tested").sum()) if total else 0
        self.summary.setText(
            f"<b>{len(self.shown):,}</b> claims shown, of {total:,} about genes with no label "
            f"({tested:,} tested by independent evidence)." if total else
            "No claims have been built for this organism.")
        self._fill()
        self._describe_recipes()

    def _fill(self):
        rows = self.shown.head(2000)
        verdicts = [c for c in rows.columns if str(c).endswith(" verdict")]
        self.table.setRowCount(len(rows))
        for i, r in enumerate(rows.itertuples(index=False)):
            rr = pd.Series(r._asdict()) if hasattr(r, "_asdict") else r
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
        if target:
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


def install(window) -> DiscoveriesPanel | None:
    """Add the Discoveries tab beside Strategies and wire it to the evidence panel and the map."""
    from . import organisms
    code = organisms.by_species(window.species).code
    products = dict(zip(window.nodes["gene_id"].astype(str),
                        window.nodes["product"].astype(str) if "product" in window.nodes
                        else [""] * len(window.nodes)))
    panel = DiscoveriesPanel(code, products)

    def show_claim(gene_id: str, target: str):
        ids = window.nodes["gene_id"].astype(str)
        hit = ids[ids == gene_id]
        if len(hit):
            window.show_detail(int(hit.index[0]))
        window._detail_link(QtCore.QUrl(f"starplast://claim/{gene_id}/{target}"))
        window.right_dock.raise_()

    panel.claim_chosen.connect(show_claim)
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
