"""Browse functional and other labels, their known members and precomputed claims.

Annotation membership is separate from prediction. Missing inference and validation
remain visible. The claims view preserves its legacy calibration and visible filters;
those legacy tests do not substitute for independently admitted biological benchmarks.
Functional classes retain overlapping memberships and organism-specific source
descriptions. Selecting a class opens its known genes; only recorded evaluations
enable scorecard navigation, so an annotation never becomes an accuracy estimate.
"""
from __future__ import annotations

from html import escape
import numpy as np
import pandas as pd
from PyQt6 import QtCore, QtWidgets

from . import claims as C
from . import discovery_labels as DL
from . import functional_results as FR
from . import functional_coverage as FC
from . import scorecard_view as SV
from . import host_foundation_scorecard as HF
from .scorecard_browser import ScorecardBrowser
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


def _domain_provenance(record, *, detailed=False):
    """Explain displayed domain nomenclature while preserving original assignment scope."""
    status=str(record.get('ontology_status','') or '')
    if not status:return ''
    lines=[f"Domain {record.get('value','')}. Displayed name source: {record.get('description_source','') or 'not available'}. "
        f"Current nomenclature status: {status.replace('_',' ')}; release: {record.get('ontology_release','') or 'unknown'}. "
        f"Accession version: {str(record.get('accession_version_status','') or 'unknown').replace('_',' ')}."]
    if 'withdrawn' in status:lines.append('This identifier is withdrawn; original gene membership is retained. Forwarding identifiers were not assigned to the gene.')
    elif status!='current_metadata':lines.append('Missing current metadata does not establish domain absence or a replacement identifier.')
    lines.append('Current nomenclature names an identifier; it does not verify gene function or establish the original assignment release.')
    if detailed:
        lines.insert(1,'Original source description: '+(str(record.get('source_description','') or '') or 'not recorded')+'.')
        for field,title in (('metadata_source_url','Nomenclature source'),('metadata_license_url','Metadata license')):
            if record.get(field):lines.append(title+': '+str(record[field]))
    return '\n\n'.join(lines)


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
        if self.context.organism!=organism:
            raise ValueError('Discoveries context must match the selected organism')
        self.inventory, self.memberships = DL.catalogue(self.context,self.frame,C.recipes(organism),T.shipped(organism))
        self.functional_benchmarks,self.functional_unavailable = FR.shipped(organism)
        self.coverage_warning = ''
        try:
            if organism not in FC.DEFAULT_ORGANISMS:
                raise LookupError('Functional coverage adapters are not declared for this organism')
            for benchmark in self.functional_benchmarks:
                FR.require_context(benchmark,self.context)
            self.functional_coverage = FC.build({organism:self.inventory},
                benchmarks=self.functional_benchmarks,organisms=(organism,))
        except LookupError as exc:
            self.coverage_warning=str(exc)
            self.functional_coverage={'rows':[]}
        except (ValueError, TypeError) as exc:
            # A custom/imported context cannot inherit a frozen source population.
            self.coverage_warning = 'Frozen recovery does not match this context: '+str(exc)
            self.functional_benchmarks=[]
            self.functional_unavailable=self.coverage_warning
            self.functional_coverage = FC.build({organism:self.inventory},organisms=(organism,))
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
        self.annotation_note.setTextFormat(QtCore.Qt.TextFormat.PlainText)
        annotations_lay.addWidget(self.annotation_note)
        self.annotation_class = QtWidgets.QComboBox()
        self.annotation_class.addItem('Choose a label above to browse its classes','')
        annotations_lay.addWidget(self.annotation_class)
        self.annotation_record = QtWidgets.QPushButton('Open held-out scorecard')
        self.annotation_record.setToolTip(TH.tip('Open existing held-out label/class results. Legacy source-recovery tests have their own scope; they are not newly admitted biological accuracy.'))
        annotations_lay.addWidget(self.annotation_record)
        self.annotation_functional = QtWidgets.QPushButton('Open functional recovery tests')
        self.annotation_functional.setToolTip(TH.tip('Inspect frozen functional annotation recovery: label, class and strategy scorecards, matched controls and every held-out gene. Independent biological activity accuracy remains unknown.'))
        annotations_lay.addWidget(self.annotation_functional)
        self.annotation_coverage = QtWidgets.QPushButton('Inspect functional strategy coverage')
        self.annotation_coverage.clicked.connect(self._open_coverage)
        annotations_lay.addWidget(self.annotation_coverage)
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
        self._build_functional_page()
        self._build_coverage_page()
        self._build_host_page()
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
        self.annotation_functional.clicked.connect(self._open_functional_record)
        self.shown = pd.DataFrame()
        self.refresh()
        self._browse_labels()

    def _build_functional_page(self):
        """Expose all frozen recovery scorecards and retained outcomes without runtime fitting."""
        self.functional_page=QtWidgets.QWidget()
        layout=QtWidgets.QVBoxLayout(self.functional_page)
        self.functional_choice=QtWidgets.QComboBox()
        for benchmark in self.functional_benchmarks:
            meta=benchmark.metadata
            self.functional_choice.addItem(f"{meta['target'].replace('_',' ')} · {meta['strategy'].replace('_',' ')}",benchmark)
        layout.addWidget(self.functional_choice)
        self.functional_note=QtWidgets.QLabel('')
        self.functional_note.setWordWrap(True)
        layout.addWidget(self.functional_note)
        self.functional_view=QtWidgets.QComboBox()
        layout.addWidget(self.functional_view)
        self.functional_card=ScorecardBrowser()
        self.functional_card.setOpenExternalLinks(False)
        layout.addWidget(self.functional_card,1)
        self.functional_export=QtWidgets.QPushButton('Save selected scorecard…')
        self.functional_export.clicked.connect(self._save_functional_card)
        layout.addWidget(self.functional_export)
        self.functional_rows=QtWidgets.QTableWidget(0,6)
        self.functional_rows.setHorizontalHeaderLabels(['held-out gene','product','recorded profile','recovered profile','outcome','vote share (uncalibrated)'])
        self.functional_rows.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.functional_rows.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.functional_rows.verticalHeader().setVisible(False)
        layout.addWidget(self.functional_rows,1)
        self.functional_choice.currentIndexChanged.connect(self._functional_selected)
        self.functional_view.currentIndexChanged.connect(self._functional_view_selected)
        self.functional_rows.cellClicked.connect(self._functional_gene_clicked)
        self.functional_card.rows_requested.connect(lambda _key:self.functional_rows.setFocus())
        self.tabs.addTab(self.functional_page,'Functional tests')
        self._functional_selected()

    def _functional_selected(self):
        """Keep tested reference recovery and its source context explicit at every level."""
        benchmark=self.functional_choice.currentData()
        self.functional_view.blockSignals(True)
        self.functional_view.clear()
        if benchmark is None:
            self.functional_note.setText(self.functional_unavailable)
            self.functional_choice.setEnabled(False)
            self.functional_view.setEnabled(False)
        else:
            self.functional_choice.setEnabled(True)
            self.functional_view.setEnabled(True)
            self.functional_view.addItem('Label: complete annotation profiles','label:')
            self.functional_view.addItem('Strategy: '+benchmark.metadata['strategy'].replace('_',' '),'strategy:')
            for name in benchmark.baseline_cards:
                self.functional_view.addItem('Training-only control: '+name.replace('_',' '),'baseline:'+name)
            for card in benchmark.major_class_cards:
                self.functional_view.addItem('Class: '+FR.class_title(card['class']),'major:'+card['class'])
            for card in benchmark.profile_class_cards:
                self.functional_view.addItem('Complete profile: '+card['class'],'profile:'+card['class'])
            self.functional_note.setText(f"{len(benchmark.rows):,} held-out genes; frozen reference-annotation recovery. "
                "Independent biological activity accuracy: not evaluated. No calibrated probabilities or new verified unknown-gene claims. "
                "Unannotated/unresolved genes were excluded as unknown, not treated as biological negatives. "
                "Source context: "+benchmark.metadata['source_context'])
        self.functional_view.blockSignals(False)
        self._functional_view_selected()

    def _build_coverage_page(self):
        """Expose exact source/target/strategy/task availability without fitting."""
        self.coverage_page=QtWidgets.QWidget()
        layout=QtWidgets.QVBoxLayout(self.coverage_page)
        self.coverage_note=QtWidgets.QLabel('')
        self.coverage_note.setWordWrap(True)
        layout.addWidget(self.coverage_note)
        self.coverage_table=QtWidgets.QTableWidget(0,9)
        self.coverage_table.setHorizontalHeaderLabels(['source label','target','strategy','task',
            'adapter status','reference recovery','biological test','calibration','deployment'])
        self.coverage_table.setEditTriggers(QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        self.coverage_table.setSelectionBehavior(QtWidgets.QAbstractItemView.SelectionBehavior.SelectRows)
        self.coverage_table.verticalHeader().setVisible(False)
        layout.addWidget(self.coverage_table,1)
        self.coverage_card=ScorecardBrowser()
        layout.addWidget(self.coverage_card,1)
        self.coverage_open=QtWidgets.QPushButton('Open recorded recovery scorecard')
        self.coverage_open.setEnabled(False)
        self.coverage_open.clicked.connect(self._open_coverage_result)
        layout.addWidget(self.coverage_open)
        self.coverage_table.cellClicked.connect(self._coverage_selected)
        self.tabs.addTab(self.coverage_page,'Strategy coverage')

    def _build_host_page(self):
        """Show the reviewed human foundation while host gene inference remains unavailable."""
        self.host_page=QtWidgets.QWidget()
        layout=QtWidgets.QVBoxLayout(self.host_page)
        note=QtWidgets.QLabel('Human host source foundation: source preservation and mapping only. '
            'This separate human scope does not inherit parasite accuracy. Host gene spaces and inference packs are unfinished; '
            'mouse source-gene foundation: unavailable.')
        note.setWordWrap(True)
        layout.addWidget(note)
        self.host_card=ScorecardBrowser()
        layout.addWidget(self.host_card,1)
        view,reason=HF.shipped()
        if view is None:self.host_card.setPlainText(reason)
        else:self.host_card.set_scorecard(view)
        self.tabs.addTab(self.host_page,'Host source status')

    def _refresh_coverage(self):
        """Filter by the selected original source label and keep missing states visible."""
        target=self.label.currentData()
        self.coverage_shown=[row for row in self.functional_coverage['rows']
            if not target or row['source_label']==target]
        self.coverage_table.setRowCount(len(self.coverage_shown))
        for i,row in enumerate(self.coverage_shown):
            recovery=row['reference_recovery']
            recorded=(f"{recovery['rows']} test genes; {recovery['answered']} answered"
                if recovery['status']=='verified_reference_recovery' else 'unavailable')
            values=[row['source_label'],row['derived_target'],row['strategy'],row['task'],
                row['applicability_status'].replace('_',' '),recorded,
                row['independent_biological_test']['status'],row['calibration']['status'],row['deployment']['status']]
            for j,value in enumerate(values):
                self.coverage_table.setItem(i,j,QtWidgets.QTableWidgetItem(value))
        self.coverage_table.resizeColumnsToContents()
        self.coverage_open.setEnabled(False)
        self.coverage_card.clear()
        self.coverage_selected_row=None
        self.annotation_coverage.setEnabled(bool(self.coverage_shown))
        self.coverage_note.setText(f'{len(self.coverage_shown):,} source/target/strategy/task addresses for {self.organism}. '
            'Click a row for its source counts, applicability and exact evidence scope. '
            'Declarations do not establish runnable inputs. Reference recovery measures recorded annotations; '
            'independent biological testing, calibration and unknown-gene deployment remain unavailable. '
            'Counts are not pooled accuracy or independent studies. '+self.coverage_warning)

    def _open_coverage(self):
        """Open coverage for the current functional source without changing claim filters."""
        self.tabs.setCurrentWidget(self.coverage_page)

    def _coverage_selected(self,row,_column):
        """Display one scoped evidence-quality card and enable only its own result."""
        if not 0<=row<len(self.coverage_shown):return
        record=self.coverage_shown[row]
        self.coverage_selected_row=record
        self.coverage_card.set_scorecard(SV.build_scorecard_view({
            'scope':{'organism':record['organism'],'target':record['derived_target'],
                'strategy':record['strategy'],'unit':record['benchmark_unit']},
            'evidence':record},kind='evidence',title='Functional strategy coverage',
            status=record['reference_recovery']['status'].replace('_',' '),
            details={'source':{'name':record['source_label'],'grade':'Unresolved source admission',
                'negative_semantics':record['negative_semantics']},'gaps':record['applicability_status']}))
        self.coverage_open.setEnabled(record['reference_recovery']['status']=='verified_reference_recovery')

    def _open_coverage_result(self):
        """Navigate by the complete validated address, never by a borrowed accuracy."""
        row=self.coverage_selected_row
        if not row or not self.coverage_open.isEnabled():return
        recovery=row['reference_recovery']
        for i,benchmark in enumerate(self.functional_benchmarks):
            if (benchmark.metadata['source_artifact_identity']==recovery['source_artifact_identity']
                    and benchmark.source_matches(row['source_label'])
                    and benchmark.metadata['organism']==row['organism']
                    and benchmark.metadata['target']==row['derived_target']
                    and benchmark.metadata['strategy']==row['strategy']
                    and benchmark.card['scope']['task']==row['task']):
                self.functional_choice.setCurrentIndex(i)
                self.functional_view.setCurrentIndex(0)
                self.tabs.setCurrentWidget(self.functional_page)
                return

    def _save_functional_card(self):
        """Export the exact selected shared card including its supplied test record."""
        if self.functional_card.view is None:return
        path,_=QtWidgets.QFileDialog.getSaveFileName(self,'Save scorecard','','JSON (*.json)')
        if path:
            from pathlib import Path
            Path(path).write_text(self.functional_card.export_json(),encoding='utf-8')

    def _functional_view_selected(self):
        """Display one exact scoped card and include missed calls and abstentions in outcomes."""
        benchmark=self.functional_choice.currentData()
        address=self.functional_view.currentData()
        if benchmark is None or address is None:
            self.functional_card.clear();self.functional_rows.setRowCount(0);self.functional_shown_rows=[]
            self.functional_export.setEnabled(False)
            return
        kind,value=address.split(':',1)
        rows=benchmark.rows
        if kind=='major':
            card=next(card for card in benchmark.major_class_cards if card['class']==value)
            rows=[row for row in rows if value in FR.profile_classes(row['truth'])
                or value in (FR.profile_classes(row['prediction']) or ())]
            metrics={name:card.get(name) for name in ('precision','recall','f1','coverage','reference_prevalence')}
            counts={name:card[name] for name in ('eligible','known_positive_genes','true_positive','false_positive','false_negative')}
            title='Class: '+FR.class_title(value)
        else:
            card=(benchmark.baseline_cards[value] if kind=='baseline' else
                next(card for card in benchmark.profile_class_cards if card['class']==value) if kind=='profile' else benchmark.card)
            metrics=card['metrics'];counts=card['counts']
            if kind=='baseline':rows=[]
            if kind=='profile':
                # Keep false calls from other reference profiles visible too.
                rows=[row for row in rows if row['truth']==value or row['prediction']==value]
                metrics={**metrics,**card['class_metrics']}
            title=self.functional_view.currentText()
        def display(value):
            """Render null metrics as unavailable rather than as zero accuracy."""
            if value is None:return 'unavailable'
            return f'{value:.6g}' if isinstance(value,float) else str(value)
        values={**counts,**metrics}
        uncertainty=card.get('extra',{}).get('uncertainty',{})
        for name in ('accuracy_all_hidden','accuracy_among_calls'):
            if uncertainty.get(name) is not None:
                bounds=uncertainty[name]
                values[name+' descriptive group interval']=' – '.join(display(bound) for bound in bounds)
        table=''.join('<tr><td>'+escape(str(key).replace('_',' '))+'</td><td>'+escape(display(number))+'</td></tr>' for key,number in values.items())
        meta=benchmark.metadata;scope=benchmark.card['scope']
        summary=meta['summary']
        context=(f"Benchmark: {meta['benchmark_id']} · truth grade: {meta['truth_grade']} · "
            f"strategy: {meta['strategy']} · frozen seed: {scope['seed']} · partition: {scope['partition']}. "
            f"Frozen roles: {summary['train']} training, {summary['tune']} tuning, "
            f"{summary['calibration']} calibration, {summary['test']} test genes. "
            "This fixed adapter fitted training genes only; tuning and calibration roles were unused. "
            'Settings: '+scope['settings'])
        row_note=('Control card shown; native gene outcomes are hidden because they are not control predictions.' if kind=='baseline' else
            'Outcomes include reference members, missed calls and calls from other profiles; class precision counts those false calls.' if kind in {'major','profile'} else
            'Every frozen test gene is retained below, including abstentions and unsupported reference profiles.')
        metric_note=('Class precision is the share of class calls matching the recorded reference; recall is the share of recorded members recovered, including misses and abstentions in its denominator.' if kind in {'major','profile'} else
            'Accuracy uses all eligible test genes; coverage is the share given a profile call; precision of calls uses only answered genes. These measure annotation recovery in this cohort.')
        if kind=='profile':metric_note+=' Eligible, answered and accuracy on this card concern genes whose recorded profile is this class. Class precision also includes false calls from the other reference profiles; those rows are retained below.'
        limits=''.join('<li>'+escape(limit)+'</li>' for limit in scope.get('gaps',[]))
        self.functional_card.setHtml('<h3>'+escape(title)+'</h3><p>Recovery of recorded annotations; independent biological precision/recall unavailable.</p>'
            '<p>'+escape(context)+'</p><table>'+table+'</table><p>'+escape(metric_note)+'</p><p>'+escape(row_note)+'</p><p>Biological precision: unavailable. Biological recall: unavailable. Calibrated confidence: unavailable.</p>'
            '<p>Source artifact: '+escape(meta['source_artifact_identity'])+'</p><ul>'+limits+'</ul>')
        self.functional_card.view=None
        if kind in {'label','strategy','baseline'}:
            self.functional_card.set_scorecard(SV.build_scorecard_view(card,title=title,
                status='Reference-annotation recovery; independent biological activity accuracy unavailable. Calibrated confidence: unavailable.',
                details={'source':{'name':', '.join(meta['source_targets']),'grade':meta['truth_grade'],
                    'sha256':meta['source_artifact_identity'],'context':meta['source_context'],
                    'lineage':'Frozen recorded reference profiles; independent biological activity not admitted',
                    'negative_semantics':'Missing or unresolved annotation is unknown, not biological absence'},
                    'split':{'roles':{key:summary[key] for key in ('train','tune','calibration','test')},
                        'scope':scope,'fit_entities':meta['source_manifest']['spec']['fit_entities']},
                    'baseline':{name:{'counts':control['counts'],'metrics':control['metrics']}
                        for name,control in benchmark.baseline_cards.items()},
                    'controls':list(benchmark.baseline_cards),'rows':rows,
                    'failures':[row for row in rows if row['truth']!=row['prediction']],
                    'freshness':{'artifact_identity':meta['source_artifact_identity']},
                    'gaps':scope.get('gaps',[])+['Biological precision: unavailable. Biological recall: unavailable. Calibrated confidence: unavailable.',row_note]},
                links=(SV.ScorecardLink('Retained gene outcomes','scorecard:rows/selected'),)))
        self.functional_export.setEnabled(self.functional_card.view is not None)
        self.functional_shown_rows=rows
        self.functional_rows.setRowCount(len(rows))
        self.functional_rows.setHorizontalHeaderItem(4,QtWidgets.QTableWidgetItem('class recovery outcome' if kind=='major' else 'complete profile outcome'))
        for i,row in enumerate(rows):
            outcome='abstained' if row['prediction'] is None else 'recovered reference profile' if row['truth']==row['prediction'] else 'different reference profile'
            if kind=='major':
                actual=value in FR.profile_classes(row['truth'])
                recovered=value in (FR.profile_classes(row['prediction']) or ())
                outcome=('recovered reference class' if actual and recovered else 'unexpected reference-class call' if recovered else
                    'abstained; reference class missed' if row['prediction'] is None else 'reference class missed')
            values=(row['entity'],self.products.get(row['entity'],''),row['truth'],row['prediction'] or 'abstained',outcome,display(row.get('support')))
            for j,text in enumerate(values):self.functional_rows.setItem(i,j,QtWidgets.QTableWidgetItem(str(text)))
            for j,key in ((2,'truth'),(3,'prediction')):
                members=FR.profile_classes(row[key])
                if members:self.functional_rows.item(i,j).setToolTip('; '.join(FR.class_title(member) for member in members))
        self.functional_rows.resizeColumnsToContents()

    def _functional_gene_clicked(self,row,_column):
        """Open one retained held-out gene through the existing organism gene navigation."""
        if 0<=row<len(self.functional_shown_rows):self.gene_chosen.emit(self.functional_shown_rows[row]['entity'])

    def _open_functional_record(self):
        """Route an installed EC source label directly to its distinct frozen recovery target."""
        target=self.label.currentData()
        for i,benchmark in enumerate(self.functional_benchmarks):
            if benchmark.source_matches(target):
                self.functional_choice.setCurrentIndex(i)
                self.functional_view.setCurrentIndex(0)
                self.tabs.setCurrentWidget(self.functional_page)
                break

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
        self._refresh_coverage()

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
            status=r.evaluation_status
            if any(b.source_matches(r.target) for b in self.functional_benchmarks):status+='; frozen functional recovery available'
            values = [r.title,r.family,f'{r.annotated_genes:,} / {r.genes:,}',str(r.classes),str(r.claims),status]
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
        self._annotation_class_metadata={}
        search = self.annotation_search.text().strip().casefold()
        for r in classes.itertuples(index=False):
            text = f'{r.value} — {r.description}' if r.description else str(r.value)
            if not search or search in text.casefold() or search in str(target).replace('_',' ').casefold():
                self.annotation_class.addItem(f'{text} ({r.annotated_genes:,} genes)',r.value)
                metadata=r._asdict()
                self._annotation_class_metadata[str(r.value)]=metadata
                provenance=_domain_provenance(metadata,detailed=True)
                if provenance:self.annotation_class.setItemData(self.annotation_class.count()-1,TH.tip(provenance),QtCore.Qt.ItemDataRole.ToolTipRole)
        self.annotation_class.blockSignals(False)
        selected = self.inventory[self.inventory.target.eq(target)] if target else pd.DataFrame()
        self.annotation_record.setEnabled(bool(len(selected) and selected.held_out_rows.iloc[0]>0))
        matched=[b for b in self.functional_benchmarks if b.source_matches(target)]
        self.annotation_functional.setEnabled(bool(matched))
        if target and len(self.inventory[self.inventory.target.eq(target)]):
            r = self.inventory[self.inventory.target.eq(target)].iloc[0]
            self._annotation_note_base = (f'{r.annotated_genes:,} genes with known annotation; {r.unannotated_genes:,} without annotation. '
                f'{r.claims:,} precomputed inferred claims. {r.evaluation_status}. '
                f'{r.held_out_rows:,} legacy held-out rows across {r.evaluated_strategies} strategies. '
                f'Source: {r.source_ids}. Independent biological class precision/recall: not evaluated. '
                'Known annotation is separate from inference; missing annotation is unknown membership. '
                'False/0 annotation flags do not establish biological absence. '+r.annotation_warning)
            if matched:
                self._annotation_note_base+=f' {sum(len(b.rows) for b in matched):,} additional frozen functional recovery rows; open Functional tests for exact profile/class/strategy scorecards. These are reference recovery, not independent biological activity tests.'
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
            provenance=_domain_provenance(r._asdict(),detailed=True)
            if provenance:self.member_table.item(i,2).setToolTip(TH.tip(provenance))
        self.member_table.resizeColumnsToContents()
        self.annotation_note.setText(self._annotation_note_base)
        metadata=self._annotation_class_metadata.get(str(value),{})
        provenance=_domain_provenance(metadata)
        self.annotation_class.setToolTip(TH.tip(_domain_provenance(metadata,detailed=True)) if provenance else '')
        if provenance:self.annotation_note.setText(self.annotation_note.text()+'\n\n'+provenance)
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
