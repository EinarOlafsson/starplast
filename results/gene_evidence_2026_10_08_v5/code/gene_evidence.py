"""Resolve scoped gene identities and expose every original evidence column.

Alias mappings come only from the supplied canonical table or an explicitly
identified existing index. This backend preserves ambiguous choices and original
stored values, including missing cells, zero, False and nested collections. Source
cards supply provenance and units; absent lineage remains an explicit gap. Slot
questions organize evidence without implying measured biological truth, inference
accuracy, cross-organism projection or automatic source admission.
"""
from __future__ import annotations

from copy import deepcopy
import hashlib
import json

import numpy as np
import pandas as pd

from . import organisms as O, query, slots
from .identity import GeneIndex


def _validate(nodes, organism):
    if organism not in O.SPACES and organism not in O.HOST_TABLES:
        raise ValueError('An explicit registered organism is required')
    if not isinstance(nodes, pd.DataFrame) or nodes.columns.has_duplicates or 'gene_id' not in nodes:
        raise ValueError('A unique-column canonical gene table is required')
    ids = nodes.gene_id.tolist()
    if any(not isinstance(gene, str) or not gene or gene != gene.strip() for gene in ids) or len(set(ids)) != len(ids):
        raise ValueError('Canonical gene IDs must be unique, nonmissing exact strings')
    if any(not isinstance(column, str) for column in nodes):
        raise ValueError('Stored column names must be explicit strings')
    return {gene: query.EntityRef(organism, 'gene', gene) for gene in ids}


def _aliases(value):
    if isinstance(value, (tuple, list, np.ndarray)):
        return [alias for item in value for alias in _aliases(item)]
    if value is None or value is pd.NA or value is pd.NaT or (pd.api.types.is_scalar(value) and pd.isna(value)):
        return []
    if not isinstance(value, str):
        raise ValueError('Stored symbols must be text or an explicit collection of text')
    # A supplied cell is an exact alias. No substring, delimiter, prefix or
    # accession-suffix heuristic creates additional symbol mappings.
    value = value.strip()
    return [value] if value and query.norm(value) else []


def build_resolver(nodes, organism, *, index=None, index_source=''):
    """Build exact table aliases plus a scoped index, preserving every ambiguous choice.

    Alias-record sources identify the ordered table mapping content separately
    from a caller-supplied index artifact. Valid same-organism index genes outside
    the current table remain choices, preserving collisions with unavailable
    targets. Selecting their evidence still requires current table membership.
    """
    entities = _validate(nodes, organism)
    columns = [column for column in ('symbol', 'gene_name') if column in nodes]
    mappings = [(gene, column, alias) for position, gene in enumerate(entities)
                for column in columns for alias in _aliases(nodes[column].iat[position])]
    payload = {'organism': organism, 'canonical_ids': list(entities), 'stored_alias_mappings': mappings}
    identity = hashlib.sha256(json.dumps(payload, ensure_ascii=False, separators=(',', ':')).encode()).hexdigest()
    source = 'table_alias_mapping_sha256:' + identity
    records = [query.AliasRecord(gene, entity, 'accession', source + ':gene_id') for gene, entity in entities.items()]
    records.extend(query.AliasRecord(alias, entities[gene], 'symbol', source + ':' + column)
                   for gene, column, alias in mappings)
    if index is not None:
        if not isinstance(index, GeneIndex):
            raise TypeError('A supplied index must be GeneIndex')
        if not isinstance(index_source, str) or not index_source.strip() or index_source != index_source.strip():
            raise ValueError('A supplied index requires an explicit index_source identity')
        index_entities = {gene: query.EntityRef(organism, 'gene', gene) for gene in index.canonical}
        referenced = {gene for gene, _kind in index.lookup.values()}
        referenced.update(gene for genes in index.ambiguous.values() for gene in genes)
        if referenced - set(index.canonical):
            raise ValueError('Index alias references lack canonical index membership')
        records.extend(query.AliasRecord(gene, entity, 'accession', index_source)
                       for gene, entity in index_entities.items())
        records.extend(query.AliasRecord(alias, index_entities[gene], kind, index_source)
                       for alias, (gene, kind) in index.lookup.items())
        records.extend(query.AliasRecord(alias, index_entities[gene], 'ambiguous_alias', index_source)
                       for alias, genes in index.ambiguous.items() for gene in genes)
    return query.GeneResolver(organism, records)


def evidence_rows(nodes, organism, gene_id, space):
    """Return every original column with scoped source cards, slot questions and gaps.

    The caller must select an exact canonical gene. Aliases and ambiguities are
    resolved separately before entering this function. Source scopes and values
    remain distinct; grouping does not invent lineage or biological accuracy.
    """
    entities = _validate(nodes, organism)
    if not isinstance(gene_id, str) or gene_id not in entities:
        raise ValueError('Select one exact canonical gene; missing or ambiguous aliases are not canonical IDs')
    position = nodes.gene_id.tolist().index(gene_id)
    matched_slots = {}
    for slot in slots.all_slots(organism):
        if slot.organism == organism and slot.unit == 'gene':
            for column in slots.declared_columns(nodes, slot):
                matched_slots.setdefault(column, []).append(slot)
    source_rows = [deepcopy(row) for row in space.rows if row['organism'] == organism and row['unit'] == 'gene']
    metadata = {}
    for row in source_rows:
        address = tuple(row.get(key, '') for key in ('source_id', 'organism', 'unit', 'question'))
        card = space.source_card(row)
        saved = json.loads(card.snapshot_json)
        metadata[address] = saved
    result = []
    for column in nodes:
        if column == 'gene_id':
            continue
        questions = list(dict.fromkeys(slot.name for slot in matched_slots.get(column, [])))
        contexts = {' / '.join(slot.context_path) or slot.context for slot in matched_slots.get(column, [])}
        sources = [deepcopy(row) for row in source_rows if column in row.get('columns', [])]
        units, origins, gaps = set(), set(), []
        for row in sources:
            address = tuple(row.get(key, '') for key in ('source_id', 'organism', 'unit', 'question'))
            saved = metadata[address]
            evidence = saved['card'].get('evidence', {})
            origin = evidence.get('origin', {})
            unit = evidence.get('quantity_units', {}).get(column, origin.get('quantity_unit', 'unresolved'))
            if isinstance(unit, str) and unit.strip() and unit.casefold() not in {'unresolved', 'unknown', 'unavailable'}:
                units.add(unit)
            origins.add(row.get('origin', 'unknown'))
            contexts.update(row.get('contexts', []))
            gaps.extend(saved.get('details', {}).get('gaps', []))
        contexts.discard('')
        if not sources:
            gaps.append('Source attribution unavailable for stored column')
            gaps.append('Source processing, mapping and evidence lineage not supplied')
        if not units:
            gaps.append('Quantity unit unresolved')
        if not contexts:
            gaps.append('Measured context unspecified')
        if not questions:
            gaps.append('Stored column is unassigned to a biological question')
        if len(units) > 1:
            gaps.append('Supplied source quantity units conflict; no conversion or preferred unit selected')
        quantity = next(iter(units)) if len(units) == 1 else 'unresolved' if not units else 'conflicting: ' + ' / '.join(sorted(units))
        result.append({'column': column, 'value': deepcopy(nodes[column].iat[position]),
                       'question': ' / '.join(questions) if questions else 'Unassigned', 'questions': questions,
                       'contexts': sorted(contexts), 'quantity_unit': quantity, 'quantity_units': sorted(units),
                       'origins': sorted(origins) or ['unknown'], 'gaps': list(dict.fromkeys(gaps)), 'sources': sources})
    return result
