"""Assemble explicitly addressed local inputs for the dataset browser."""
from __future__ import annotations

import json
from collections.abc import Mapping

from . import datasets, inventory, organisms, source_refusals


class _Graph:
    def __init__(self, mapping):
        self.mapping = mapping
        self.files = tuple(mapping)

    def __getitem__(self, key):
        return self.mapping[key]


def from_sources(source_inputs, *, imports=(), imported_organism=None):
    """Build a live inventory from supplied tables, without downloads or projections.

    ``source_inputs`` has the same organism/nodes/graph/tables structure as the
    existing slot viewer. Session imports become separate source rows; they do
    not authenticate publication provenance or change native source identities.
    """
    from .dataset_space import DatasetSpace

    tables, graphs, bridges = {}, {}, {}
    for code, source in source_inputs.items():
        if source.get('nodes') is not None:
            tables[(code, 'gene')] = source['nodes']
        if source.get('graph') is not None:
            graph = source['graph']
            graphs[code] = _Graph(graph) if isinstance(graph, Mapping) else graph
        side = source.get('tables') or {}
        if side.get('metabolite') is not None:
            tables[(code, 'metabolite')] = side['metabolite']
        if side.get('bridge') is not None:
            bridges[code] = side['bridge']
        for host, frame in (side.get('host_gene') or {}).items():
            if host not in organisms.HOST_TABLES:
                raise ValueError('Host protein table requires an explicit host organism')
            if (host, 'protein') in tables and not tables[(host, 'protein')].equals(frame):
                raise ValueError('Conflicting host protein tables')
            tables[(host, 'protein')] = frame

    registry = list(datasets.registry())
    origins = {}
    for index, record in enumerate(imports):
        columns = tuple(record.get('columns', ()))
        if not columns:
            continue
        if imported_organism not in source_inputs or (imported_organism, 'gene') not in tables:
            raise ValueError('Session import requires its explicit gene-table organism')
        key = 'session_import:' + str(index + 1)
        registry.append(datasets.Dataset(key=key, name='Session import ' + str(index + 1),
            level='gene', kind='user_imported', provides='Imported values',
            organism=imported_organism, columns=columns,
            note=json.dumps(record, sort_keys=True, allow_nan=False)))
        origins[(key, imported_organism, 'gene', '')] = {
            'kind': 'user_imported', 'evidence_grade': 'unresolved',
            'lineage': 'Session preprocessing record; original source-file identity not supplied',
            'quantity_unit': 'unresolved',
            'gaps': ['Import provenance is a session record, not independent biological truth',
                     'Original source-file identity and publication association are unavailable']}

    # A processed application cache is not a located raw assay file.
    raw_paths = {d.key: datasets.local_path(d.key) for d in datasets.registry()
                 if d.path and not d.path.startswith('starplast/data/')}
    report = inventory.build_inventory(tables=tables, graphs=graphs, bridges=bridges,
        sources=registry, raw_paths=raw_paths, refusals=source_refusals.records())
    return DatasetSpace(report, tables=tables, graphs=graphs, bridges=bridges,
        sources=registry, origins=origins)
