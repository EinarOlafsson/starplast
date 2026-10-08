"""Browse explicitly qualified inventory records and their stored evidence.

This backend preserves the inventory's source, organism, storage unit and question
boundaries. It exposes original stored values and immutable evidence cards without
projecting identifiers across species, fitting graphs, admitting biological truth
or interpreting storage coverage as accuracy. Missing lineage, unavailable tables,
unattributed pairs and question-specific refusals remain visible as separate gaps.
"""
from __future__ import annotations

from collections.abc import Mapping
from copy import deepcopy
from dataclasses import asdict, is_dataclass
import math

import numpy as np
import pandas as pd

from . import datasets, inventory as I, organisms, slots
from .provenance import EVIDENCE_GRADES, MeasurementTrace
from .scorecard_view import ScorecardLink, build_scorecard_view

_ADDRESS = ('source_id', 'organism', 'unit', 'question')
_IDENTIFIERS = {'gene': 'gene_id', 'protein': 'host_id', 'metabolite': 'metabolite'}
_COUNTS = ('table_rows', 'stored_any_rows', 'missing_unknown_rows', 'false_cells', 'pair_records',
           'assayed_pair_denominator', 'unattributed_pair_records', 'unmeasured_records',
           'unmapped_records', 'below_detection_records', 'failed_qc_records',
           'unknown_missing_records', 'observed_records', 'explicit_negative_records')
_GRADES = {'direct_experiment': 'A', 'orthology_transfer': 'B', 'prediction': 'C', 'derived_quantity': 'C'}


def _plain(value):
    if is_dataclass(value):
        return _plain(asdict(value))
    if isinstance(value, Mapping):
        return {key: _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [_plain(item) for item in value]
    if value is pd.NA or value is pd.NaT:
        return None
    if isinstance(value, np.generic):
        return _plain(value.item())
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return deepcopy(value)


def _address(row):
    return tuple(row.get(key, '') or '' for key in _ADDRESS)


def _table_copy(table):
    copied = table.copy(deep=True)
    for column in copied:
        if pd.api.types.is_object_dtype(copied[column].dtype):
            copied[column] = pd.Series([deepcopy(value) for value in copied[column]],
                                       index=copied.index, dtype=object)
    copied.attrs = deepcopy(table.attrs)
    return copied


class DatasetSpace:
    """Copied scoped inventory, facets, provenance cards and original entity records.

    Tables are keyed by ``(organism, unit)``; graphs and bridges by organism.
    Origins are mappings keyed by the full four-part inventory address, or by a
    qualified table address for its origin only. Exact-address metadata may supply
    typed ``MeasurementTrace`` objects in ``traces``. Evidence grades never follow
    automatically from an installed table or a registry description.
    """

    def __init__(self, inventory, *, tables=None, graphs=None, bridges=None, sources=None, catalog=None, origins=None):
        """Snapshot supplied storage and validate source/question addressing without I/O."""
        self._tables = {}
        for key, table in (tables or {}).items():
            if not isinstance(key, tuple) or len(key) != 2 or key[1] not in _IDENTIFIERS:
                raise ValueError('Tables require explicit (organism, gene/protein/metabolite) keys')
            if not isinstance(table, pd.DataFrame) or table.columns.has_duplicates:
                raise ValueError('Tables require unique-column dataframes')
            self._tables[key] = _table_copy(table)
        self._bridges = {key: _table_copy(table) for key, table in (bridges or {}).items()}
        for key, table in self._bridges.items():
            if key not in organisms.SPACES and key not in organisms.HOST_TABLES:
                raise ValueError('Bridges require an explicit registered organism')
            if table.columns.has_duplicates:
                raise ValueError('Bridge columns must be unique')
            if 'gene_id' in table and key in organisms.SPACES:
                if not table.gene_id.map(organisms.get(key).matches).all():
                    raise ValueError('Bridge accessions do not match the explicit organism')
        # Reuse organism/accession and denominator validation, even when an
        # already-built inventory is supplied. No inventory counts are replaced.
        I.build_inventory(self._tables, sources=[], catalog=[])
        supplied_sources = list(datasets.registry() if sources is None else sources)
        self._sources = {source.key: source for source in supplied_sources}
        if len(self._sources) != len(supplied_sources):
            raise ValueError('Source keys must be unique')
        self._catalog = tuple(slots.all_slots() if catalog is None else catalog)
        self._origins = deepcopy(origins or {})
        for key, origin in self._origins.items():
            if not isinstance(key, tuple) or len(key) not in (2, 4) or not isinstance(origin, Mapping):
                raise ValueError('Origins require qualified table or full source/question keys')
        if inventory is None:
            inventory = I.build_inventory(self._tables, graphs=graphs, bridges=self._bridges,
                                          sources=list(self._sources.values()), catalog=self._catalog)
        if not isinstance(inventory, pd.DataFrame) or inventory.columns.has_duplicates:
            raise ValueError('A unique-column inventory dataframe is required')
        if not inventory.empty and not set(_ADDRESS) <= set(inventory):
            raise ValueError('Inventory requires source, organism, unit and question scope')
        self._rows = []
        for supplied in inventory.to_dict('records'):
            row = _plain(supplied)
            for field in _COUNTS:
                value = row.get(field)
                if value is not None:
                    if isinstance(value, bool) or not isinstance(value, (int, float)) or value < 0 or int(value) != value:
                        raise ValueError('Inventory counts must be nonnegative whole records or unknown')
                    row[field] = int(value)
            if any(not isinstance(row.get(key), str) or not row[key] for key in _ADDRESS[:-1]):
                raise ValueError('Inventory source/organism/unit addresses must be explicit')
            if row['organism'] not in organisms.SPACES and row['organism'] not in organisms.HOST_TABLES:
                raise ValueError('Inventory organism must be registered')
            row['question'] = row.get('question') or ''
            if not isinstance(row['question'], str) or not isinstance(row.get('status'), str):
                raise ValueError('Question and availability status must be explicit text')
            source = self._sources.get(row['source_id'])
            origin = self._metadata(row)
            traces = self._traces(row, origin)
            grades = sorted({trace.evidence_grade for trace in traces}) or [origin.get('evidence_grade', 'unresolved')]
            if any(grade not in EVIDENCE_GRADES for grade in grades):
                raise ValueError('Unknown supplied evidence grade')
            row.update(name=source.name if source else row.get('name', row['source_id']),
                       origin=origin.get('kind', 'registry_asserted_unverified' if source else 'unknown'),
                       declared_derivation=list(source.derived_from) if source else [], evidence_grade=' / '.join(grades),
                       slot_grades=sorted({_GRADES[grade] for grade in grades if grade in _GRADES}))
            row['title'] = row['name'] + (' — ' + row['question'] if row['question'] else '')
            row['gap'] = self._gap(row)
            self._rows.append(row)
        if len({_address(row) for row in self._rows}) != len(self._rows):
            raise ValueError('Inventory contains duplicate source/organism/unit/question scopes')
        # Qualified graph handles remain caller-owned and open. Read only a
        # requested layer during paging; eager snapshots could load every edge.
        self._graphs = dict(graphs or {})

    def _metadata(self, row):
        table_origin = self._origins.get((row['organism'], row['unit']), {})
        # A table origin cannot assign all its publications one evidence grade.
        common = {key: value for key, value in table_origin.items() if key in {'kind', 'path', 'sha256', 'version'}}
        return {**common, **self._origins.get(_address(row), {})}

    def _traces(self, row, metadata):
        traces = metadata.get('traces', ())
        if not isinstance(traces, (tuple, list)) or any(not isinstance(trace, MeasurementTrace) for trace in traces):
            raise ValueError('Source traces require typed MeasurementTrace records')
        for trace in traces:
            if (trace.source_id, trace.output_organism, trace.storage_unit) != _address(row)[:3] or trace.column not in row.get('columns', []):
                raise ValueError('Trace belongs to another source/organism/unit/column')
        return tuple(traces)

    def _gap(self, row):
        if row.get('rejection_reason'):
            return row['rejection_reason']
        missing = row.get('missing_columns') or []
        if row.get('status') == 'unattributed':
            return 'Pair records lack source attribution; they cannot be assigned to this publication'
        if row.get('status') in {'unavailable', 'source_only'}:
            return 'Stored entity values unavailable' + (': ' + ', '.join(missing) if missing else '')
        return 'Missing declared columns: ' + ', '.join(missing) if missing else ''

    @property
    def rows(self):
        """Return fresh plain inventory records, preserving every scoped refusal."""
        return deepcopy(self._rows)

    def _row(self, row):
        if not isinstance(row, Mapping):
            raise TypeError('Select a scoped inventory record')
        matches = [stored for stored in self._rows if _address(stored) == _address(row)]
        if len(matches) != 1:
            raise ValueError('Source/organism/unit/question is outside this space')
        return matches[0]

    def filter_rows(self, organism=None, unit=None, family=None, context=None, status=None, text=''):
        """Filter exact facets and case-insensitive text without merging source rows."""
        def chosen(row):
            if any(value is not None and row.get(key) != value for key, value in
                   (('organism', organism), ('unit', unit), ('status', status))):
                return False
            if family is not None and family not in row.get('evidence_families', []):
                return False
            if context is not None and context not in row.get('contexts', []):
                return False
            return text.casefold() in str(row).casefold()
        return deepcopy([row for row in self._rows if chosen(row)])

    def facets(self):
        """List exact available organism, unit, family, context and status filters."""
        return {key: sorted({value for row in self._rows for value in
                (row.get(field, []) if multiple else [row[field]])})
                for key, field, multiple in (('organism', 'organism', False), ('unit', 'unit', False),
                    ('family', 'evidence_families', True), ('context', 'contexts', True), ('status', 'status', False))}

    def source_card(self, row):
        """Present one immutable evidence card; storage counts are never accuracy."""
        row = self._row(row)
        metadata = self._metadata(row)
        traces = self._traces(row, metadata)
        source = self._sources.get(row['source_id'])
        links = []
        url = row.get('url') or (source.url if source else '')
        unsupported_url = False
        if url:
            try:
                links.append(ScorecardLink('Source', url))
            except ValueError:
                unsupported_url = True
        if source and source.pmid and str(source.pmid).isdigit():
            links.append(ScorecardLink('Publication', 'https://pubmed.ncbi.nlm.nih.gov/' + str(source.pmid) + '/'))
        evidence = {'inventory': row, 'origin': _plain(metadata), 'traces': _plain(traces),
                    'registry': _plain(source) if source else None,
                    'quantity_units': {trace.column: trace.quantity_unit for trace in traces},
                    'slot_grades': {'A': 'Supplied direct experiment', 'B': 'Supplied orthology transfer',
                                    'C': 'Supplied prediction or derived quantity'},
                    'storage_interpretation': 'Nonmissing stored values include zero and False; coverage is not assay sampling or accuracy'}
        gaps = [row['gap']] if row['gap'] else []
        gaps.extend(metadata.get('gaps', []))
        for trace in traces:
            gaps.extend(trace.gaps)
        if not traces and not metadata.get('lineage'):
            gaps.append('Source processing, mapping and evidence lineage not supplied')
        if unsupported_url:
            gaps.append('Recorded source URL is not an allowed HTTPS navigation destination')
        if row['unit'] == 'pair' and row.get('assayed_pair_denominator') is None:
            gaps.append('Assayed pair universe unavailable; pair records are not biological negatives')
        card = {'scope': {'organism': row['organism'], 'unit': row['unit'],
                          'target': row['question'] or row['source_id']},
                'counts': {key: row.get(key) for key in _COUNTS}, 'evidence': evidence,
                'source': {'name': row['name'], 'grade': row['evidence_grade'],
                           'version': metadata.get('version', 'unresolved'), 'sha256': metadata.get('sha256', 'unresolved'),
                           'context': row.get('contexts', []), 'lineage': _plain(traces) or metadata.get('lineage', 'unresolved'),
                           'negative_semantics': 'Missing stored values are unknown; zero/False remain original values, not automatically verified negatives'}}
        return build_scorecard_view(card, kind='evidence', title=row['title'], status=row['status'],
                                    details={'gaps': gaps, 'rows': {'scope': _address(row),
                                             'quantity_unit': metadata.get('quantity_unit', 'unresolved')}}, links=tuple(links))

    def entities(self, row, *, start=0, limit=None):
        """Copy original stored entities/values; unavailable scopes return an explicit gap."""
        if type(start) is not int or start < 0 or (limit is not None and (type(limit) is not int or limit < 1)):
            raise ValueError('Page start must be nonnegative and limit positive or None')
        row = self._row(row)
        identifier = _IDENTIFIERS.get(row['unit'])
        attrs = {'organism': row['organism'], 'unit': row['unit'], 'entity_column': identifier,
                 'source_id': row['source_id'], 'question': row['question'],
                 'status': row['status'], 'gap': row['gap'], 'start': start,
                 'total_rows': None, 'displayed_rows': 0}
        def unavailable(gap):
            result = pd.DataFrame(columns=([identifier] if identifier else []) + row.get('columns', []))
            result.attrs.update(attrs, status='unavailable', gap=gap)
            return result
        if row['status'] in {'rejected', 'unavailable', 'source_only'}:
            return unavailable(row['gap'] or 'No attributable stored entity records')
        if identifier:
            table = self._tables.get((row['organism'], row['unit']))
            if table is None or identifier not in table:
                return unavailable('Explicit stored entity identifiers unavailable; row positions are not entity IDs')
            columns = [column for column in row.get('columns', []) if column in table and column != identifier]
            if not columns:
                return unavailable('Declared stored value columns unavailable')
            stop = None if limit is None else start + limit
            result = _table_copy(table.iloc[start:stop][[identifier, *columns]])
            attrs['total_rows'] = len(table)
        elif row['unit'] == 'pair':
            frames, schema = [], []
            found = False
            graph = self._graphs.get(row['organism'], {})
            available = getattr(graph, 'files', graph.keys() if isinstance(graph, Mapping) else ())
            offset = 0
            stop = None if limit is None else start + limit
            nodes = self._tables.get((row['organism'], 'gene'))
            mapped = False
            if nodes is not None and 'gene_id' in nodes and 'gene_ids' in available:
                graph_ids = np.asarray(graph['gene_ids']).astype(str)
                mapped = tuple(graph_ids) == tuple(nodes.gene_id.astype(str))
            for column in row.get('columns', []):
                if column in row.get('missing_columns', []):
                    continue
                if column.startswith(slots.EDGE_PREFIX):
                    name = column[len(slots.EDGE_PREFIX):]
                    keys = [name + '__a', name + '__b']
                    if any(key not in available for key in keys):
                        continue
                    if name + '__w' in available:
                        keys.append(name + '__w')
                    arrays = {key: np.asarray(graph[key]) for key in keys}
                    n = len(arrays[keys[0]])
                    if any(array.ndim != 1 or len(array) != n for array in arrays.values()):
                        raise ValueError('Graph endpoints/weights do not form complete stored pair records')
                    found = True
                    schema.extend(keys + (['endpoint_a_gene_id', 'endpoint_b_gene_id'] if mapped else []) + ['source_column'])
                    lo, hi = max(0, start - offset), n if stop is None else min(n, stop - offset)
                    if lo < hi:
                        selected = pd.DataFrame({key: array[lo:hi] for key, array in arrays.items()})
                        if mapped:
                            for key, label in zip(keys[:2], ('endpoint_a_gene_id', 'endpoint_b_gene_id')):
                                endpoint = arrays[key][lo:hi]
                                if endpoint.dtype.kind not in 'iu' or (endpoint < 0).any() or (endpoint >= len(nodes)).any():
                                    raise ValueError('Pair endpoints exceed the explicit matching gene universe')
                                selected[label] = nodes.gene_id.iloc[endpoint].to_numpy()
                        selected['source_column'] = column
                        frames.append(selected)
                    offset += n
                elif column.startswith('bridge:'):
                    table = self._bridges.get(row['organism'])
                    if table is not None and {'source', 'bridge'} <= set(table):
                        found = True
                        schema.extend(table.columns)
                        positions = np.flatnonzero(((table.source == row['source_id']) &
                                                   (table.bridge == column.split(':', 1)[1])).to_numpy())
                        n = len(positions)
                        lo, hi = max(0, start - offset), n if stop is None else min(n, stop - offset)
                        if lo < hi:
                            frames.append(table.iloc[positions[lo:hi]])
                        offset += n
            if not found:
                return unavailable('Attributable pair endpoints unavailable')
            result = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=list(dict.fromkeys(schema)))
            attrs.update(total_rows=offset, endpoint_mapping='verified_gene_order' if mapped else 'unresolved')
            if not mapped and any(column.startswith(slots.EDGE_PREFIX) for column in row.get('columns', [])):
                attrs['gap'] = '; '.join(filter(None, (attrs['gap'], 'Graph endpoint gene mapping unavailable without matching explicit gene_ids and qualified table')))
        else:
            return unavailable('This storage unit has no qualified entity table')
        attrs['displayed_rows'] = len(result)
        result.attrs.update(attrs)
        return _table_copy(result)
