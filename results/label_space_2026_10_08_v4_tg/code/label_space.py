"""Scoped label and class entry points over recorded native annotations.

Membership, descriptions, unknown states and source coverage stay separate from
performance. Existing categorical ledger cohorts can be presented only against
their explicitly supplied unchanged table and unspecified legacy context. Each
setting, seed and hold-out mode retains its own shared scorecard. No classifier
fits, source acquisition, inferred protein functions, biological negatives or
cross-organism pooling occur. Every strategy remains visible with explicit gaps.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import asdict

import numpy as np
import pandas as pd

from . import capabilities as C, discovery_labels as DL, functional_domain_profiles as DP
from . import query as Q, record_scorecards as R, slots, strategies as S
from .scorecard_view import build_scorecard_view


def _items(value):
    if isinstance(value, (list, tuple, np.ndarray)):
        return list(value)
    return [value]


def _known(value):
    if value is None or value is pd.NA or value is pd.NaT:
        return False
    if pd.api.types.is_scalar(value) and pd.isna(value):
        return False
    return isinstance(value, (str, bool, int, float, np.generic)) and str(value).strip().lower() not in S.ABSENT


class LabelSpace:
    """One exact organism/table/context, its native classes and frozen ledger views."""

    def __init__(self, nodes, organism, *, ledger=None, ledger_nodes=None, context=None):
        """Copy supplied annotations and refuse changed or context-unbound ledger accuracy."""
        if not isinstance(nodes, pd.DataFrame) or nodes.columns.has_duplicates or 'gene_id' not in nodes:
            raise ValueError('An explicit unique-column canonical gene table is required')
        genes = nodes.gene_id.tolist()
        if len(set(genes)) != len(genes):
            raise ValueError('Canonical genes must be unique')
        for gene in genes:
            Q.EntityRef(organism, 'gene', gene)
        self.organism = organism
        self.context = Q.BiologicalContext() if context is None else context
        if not isinstance(self.context, Q.BiologicalContext):
            raise TypeError('An explicit BiologicalContext is required')
        self._nodes = nodes.copy(deep=True)
        for column in self._nodes:
            if pd.api.types.is_object_dtype(self._nodes[column].dtype):
                self._nodes[column] = pd.Series([deepcopy(value) for value in nodes[column]], index=nodes.index, dtype=object)
        self._ctx = S.Context(self._nodes, organism=organism, graph={})
        # Explicit collections in function fields represent complete assigned
        # source entries; no substring or delimiter guessing applies elsewhere.
        source_nodes = self._nodes.copy(deep=True)
        for target in DL.FUNCTION_FIELDS & set(source_nodes):
            source_nodes[target] = source_nodes[target].map(lambda value:
                '; '.join(str(item) for item in _items(value) if _known(item))
                if isinstance(value, (list, tuple, np.ndarray)) else value)
        source_ctx = S.Context(source_nodes, organism=organism, graph={})
        catalog, observed = DL.catalogue(source_ctx, pd.DataFrame(), pd.DataFrame())
        self._catalog, pieces = catalog.copy(deep=True), []
        self._native = {}
        for target in self._catalog.target:
            if target not in self._nodes:
                continue
            if target in DL.FUNCTION_FIELDS:
                membership = observed[observed.target.eq(target)].copy()
                if target in DP.TARGET_PATTERNS:
                    complete = {str(source_ctx.gene_ids[i]): set(DP.parse_source_profile(
                        str(value).strip() if _known(value) else None, target)['recorded_identifiers'])
                        for i, value in enumerate(source_nodes[target])
                        if DP.parse_source_profile(str(value).strip() if _known(value) else None, target)['eligible']}
                    membership = membership[membership.apply(lambda row:
                        row.gene_id in complete and row.value in complete[row.gene_id], axis=1)]
                native = {value: value for value in membership.value.unique()}
            else:
                rows, native = [], {}
                for i, cell in enumerate(self._nodes[target]):
                    for value in _items(cell):
                        if not _known(value):
                            continue
                        key = str(value).strip()
                        raw = value.item() if isinstance(value, np.generic) else value
                        if key in native and type(native[key]) is not type(raw):
                            raise ValueError('Distinct native classes share the same displayed class key')
                        native[key] = deepcopy(raw)
                        rows.append({'organism': organism, 'target': target, 'value': key,
                                     'description': '', 'gene_id': str(self._ctx.gene_ids[i])})
                membership = pd.DataFrame(rows, columns=['organism', 'target', 'value', 'description', 'gene_id'])
            membership = membership.drop_duplicates(['organism', 'target', 'value', 'gene_id'])
            pieces.append(membership)
            self._native[target] = native
            known = int(membership.gene_id.nunique())
            source_known = sum(any(_known(value) for value in _items(cell)) for cell in self._nodes[target])
            mask = self._catalog.target.eq(target)
            self._catalog.loc[mask, ['annotated_genes', 'unannotated_genes', 'classes']] = [known, len(genes) - known, len(native)]
            self._catalog.loc[mask, 'annotation_coverage'] = known / len(genes) if genes else None
            self._catalog.loc[mask, ['source_annotated_genes', 'unparsed_annotation_genes']] = [source_known, source_known - known]
        self._members = pd.concat(pieces, ignore_index=True) if pieces else observed.iloc[:0].copy()
        self._ledger = pd.DataFrame() if ledger is None else ledger.copy(deep=True)
        self.ledger_gap = 'Frozen ledger not supplied'
        self._ledger_bound = False
        if len(self._ledger):
            if not isinstance(ledger_nodes, pd.DataFrame) or not nodes.equals(ledger_nodes):
                self.ledger_gap = 'Frozen ledger table binding unavailable or changed; archived accuracy withheld'
            elif self.context != Q.BiologicalContext():
                self.ledger_gap = 'Legacy ledger context unspecified; requested context accuracy withheld'
            else:
                self._ledger_bound, self.ledger_gap = True, ''
        self._evaluations = {}
        if self._ledger_bound:
            for target in self._catalog.target:
                selected = self._ledger[self._ledger.organism.eq(organism) & self._ledger.target.eq(target)]
                truth = dict(zip(self._ctx.gene_ids, self._ctx.truth(target)))
                if len(selected) and all(gene in truth and truth[gene] == observed and observed in self._native[target]
                                         for gene, observed in zip(selected.gene_id, selected.truth)):
                    mask = self._catalog.target.eq(target)
                    self._catalog.loc[mask, ['held_out_rows', 'evaluated_strategies']] = [len(selected), selected.strategy.nunique()]
                    self._catalog.loc[mask, 'evaluation_status'] = 'Legacy held-out arithmetic; source and truth admission unresolved'

    @property
    def labels(self):
        """Copy the native label catalog and its exact membership/unknown counts."""
        return deepcopy(self._catalog)

    def _target(self, target):
        if target not in set(self._catalog.target):
            raise ValueError('Label is outside this exact organism/table scope')
        return self._catalog[self._catalog.target.eq(target)].iloc[0].to_dict()

    def _value(self, target, value):
        self._target(target)
        key = str(value).strip()
        if key not in self._native.get(target, {}):
            raise ValueError('Class is outside this exact label scope')
        return key

    def members(self, target, value=None):
        """Copy every recorded membership or one exact class, retaining overlaps."""
        self._target(target)
        selected = self._members[self._members.target.eq(target)]
        if value is not None:
            selected = selected[selected.value.eq(self._value(target, value))]
        return deepcopy(selected.reset_index(drop=True))

    def classes(self, target):
        """Return native class balance, observed definitions and unavailable-truth gaps."""
        label = self._target(target)
        frame = DL.class_summary(self._members, target)
        frame['native_value'] = pd.Series([deepcopy(self._native[target][value]) for value in frame.value], dtype=object)
        frame['class_fraction'] = frame.annotated_genes / label['genes'] if label['genes'] else None
        frame['unknown_genes'] = label['unannotated_genes']
        frame['gaps'] = [['Class omissions are reference annotations, not verified biological negatives'] for _ in range(len(frame))]
        return frame

    def query(self, target, value=None):
        """Build a typed label/class address without losing organism or requested context."""
        self._target(target)
        return Q.Query(self.organism, 'label' if value is None else 'class', target=target,
                       values=() if value is None else (self._value(target, value),),
                       context=self.context, output='scorecard')

    def evidence_card(self, target, value=None):
        """Present source memberships and unknown denominators without inferred performance."""
        label = self._target(target)
        members = self.members(target, value)
        classes = self.classes(target)
        contexts = sorted({' / '.join(slot.context_path) or slot.context
                           for slot in slots.all_slots(self.organism)
                           if slot.unit == 'gene' and target in slots.declared_columns(self._nodes, slot)})
        if value is not None:
            classes = classes[classes.value.eq(self._value(target, value))]
        card = {'scope': {'organism': self.organism, 'target': target, 'unit': 'gene'},
                'counts': {'whole_universe': int(label['genes']), 'annotated_genes': int(label['annotated_genes']),
                           'unknown_genes': int(label['unannotated_genes']), 'class_members': int(members.gene_id.nunique())},
                'source': {'name': label['source_ids'], 'grade': 'unresolved',
                           'context': asdict(self.context), 'negative_semantics': label['truth_interpretation']},
                'evidence': {'label': label, 'classes': classes.drop(columns=['native_value']).to_dict('records'),
                             'native_values': {key: deepcopy(raw) for key, raw in self._native[target].items()},
                             'query': self.query(target, value).to_dict(),
                             'recorded_contexts': contexts,
                             'hierarchy': {'status': 'unavailable', 'reason': 'Native source does not supply an admitted class-parent hierarchy'},
                             'unmeasured_classes': {'status': 'unavailable', 'reason': 'Class universe beyond observed memberships is not supplied'},
                             'balance_semantics': 'Classes may overlap; fractions use the whole original gene universe'}}
        return build_scorecard_view(card, kind='evidence', title=target if value is None else target + ': ' + self._value(target, value),
                                    status='Recorded annotation membership', details={'gaps': [label['annotation_warning'],
                                    'Gene assignment lineage and independent biology remain unresolved']})

    def _cohorts(self, target):
        if target in self._evaluations:
            return self._evaluations[target]
        evaluations = []
        if self._ledger_bound:
            selected = self._ledger[self._ledger.organism.eq(self.organism) & self._ledger.target.eq(target)]
            if len(selected):
                truth = dict(zip(self._ctx.gene_ids, self._ctx.truth(target)))
                for scope, rows in R.legacy_cohorts(selected):
                    # A single-label ledger cannot validate overlapping member
                    # classes by silently choosing one domain or enzyme code.
                    valid = all(gene in truth and truth[gene] == observed and observed in self._native[target]
                                for gene, observed in zip(rows.entity, rows.truth))
                    if not valid:
                        evaluations.append({'strategy': scope.strategy, 'scope': asdict(scope), 'gap':
                            'Frozen ledger truth/member identity differs from current native annotations'})
                        continue
                    views = R.views(rows, scope)
                    evaluations.append({'strategy': scope.strategy, 'scope': asdict(scope), 'views': views})
        self._evaluations[target] = evaluations
        return evaluations

    def mechanisms(self, target, value=None):
        """Expose all 39 strategies and separate frozen evaluations or explicit unavailable metrics."""
        self._target(target)
        key = None if value is None else self._value(target, value)
        query = self.query(target, value)
        results = []
        capabilities = {capability.strategy: capability for capability in C.catalog()}
        for strategy in S.catalog():
            capability = capabilities[strategy.key]
            supported = query.kind in capability.query_kinds and self.organism in capability.organisms
            evaluations, gaps = [], []
            for cohort in self._cohorts(target):
                if cohort['strategy'] != strategy.key:
                    continue
                if 'gap' in cohort:
                    gaps.append(cohort['gap'])
                    continue
                views = cohort['views']
                raw = views['target'] if key is None else next((card for card in views['class'] if card['class'] == key), None)
                if raw is None:
                    gaps.append('Requested class has no truth rows in this frozen cohort')
                    continue
                full = views['target']
                raw = deepcopy(raw)
                if key is not None:
                    raw['extra'] = {**raw['extra'], 'class': key, 'class_metrics': deepcopy(raw['class_metrics']),
                                    'full_cohort_confusion': deepcopy(full['extra'].get('confusion', [])),
                                    'class_interpretation': 'Precision includes calls from other reference classes; recall includes abstentions'}
                evaluations.append({'scope': deepcopy(cohort['scope']), 'raw_card': deepcopy(raw),
                    'card': build_scorecard_view(raw, details={'gaps': list(cohort['scope']['gaps'])}),
                    'class_metrics': deepcopy(raw.get('class_metrics', {})),
                    'confusion': deepcopy(full['extra'].get('confusion', [])), 'gaps': list(cohort['scope']['gaps'])})
            if not evaluations:
                gaps.append(self.ledger_gap or 'No matching frozen cohort for this strategy/label/class')
            status = 'Legacy held-out arithmetic; independent biology unadmitted' if evaluations else 'Metrics unavailable'
            results.append({'strategy': strategy.key, 'title': strategy.title, 'applicable': supported,
                            'status': status, 'metrics': deepcopy(evaluations[0]['raw_card']['metrics']) if len(evaluations) == 1 else None,
                            'card': evaluations[0]['card'] if len(evaluations) == 1 else None,
                            'evaluations': evaluations, 'gaps': list(dict.fromkeys(gaps)),
                            'applicability_interpretation': 'Declared query kind only; source-excluded input availability not certified'})
        return results
