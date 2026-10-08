"""Inventory registered evidence without treating absence as a negative assay.

One row per source/organism/storage unit retains slot context and evidence facets.
Coverage is the fraction of installed rows having a stored nonmissing value in
any declared column, not an assay's sampling fraction or a ground-truth accuracy.
Null causes remain unknown unless source-level observations explicitly record
them. Pair counts have no coverage denominator without an assayed pair universe.
Rejected proposals are scoped to their question, never to the entire publication.
"""
from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd

from . import datasets, organisms, slots


def _present(series):
    present = series.notna()
    if pd.api.types.is_numeric_dtype(series) and not pd.api.types.is_bool_dtype(series):
        return present & np.isfinite(series.fillna(0).astype(float))
    if not pd.api.types.is_bool_dtype(series):
        # This is storage coverage: an 'unknown' annotation is stored text, not
        # biological truth. Empty cells are absent; 0 and False are retained.
        present &= series.map(lambda value: not isinstance(value, str) or bool(value.strip()))
    return present


def _base(source_id, organism, unit):
    return {"source_id": source_id, "organism": organism, "unit": unit,
            "registry_organism": organism, "status": "unavailable", "columns": [],
            "missing_columns": [], "contexts": [], "evidence_families": [],
            "table_rows": None, "stored_any_rows": None, "missing_unknown_rows": None,
            "stored_fraction": None, "false_cells": None, "pair_records": None,
            "assayed_pair_denominator": None, "unattributed_pair_records": None, "unmeasured_records": None,
            "unmapped_records": None, "below_detection_records": None,
            "failed_qc_records": None, "unknown_missing_records": None,
            "observed_records": None, "explicit_negative_records": None,
            "raw_available": None, "processed_available": None, "access_status": "not_checked",
            "url": "", "note": "", "rejection_reason": "", "question": ""}


def _projection(dataset, column, catalog):
    code = dataset.organism
    # Some legacy sources are registered in the parasite's *exploration* scope
    # but explicitly store their outputs in a host protein table. Keep both
    # addresses rather than projecting host proteins onto parasite genes.
    for host, filename in organisms.HOST_TABLES.items():
        if dataset.path == "starplast/data/" + filename:
            code = host
            break
    if column.startswith((slots.EDGE_PREFIX, "bridge:")):
        return code, "pair"
    if code in organisms.HOST_TABLES:
        return code, "protein"
    frame = pd.DataFrame(columns=[column])
    matched = {s.unit for s in catalog if s.organism == code and slots.declared_columns(frame, s)}
    # The existing catalog's host_gene rows describe host references within the
    # parasite UI scope. A missing host address cannot be guessed from a symbol.
    if "host_gene" in matched:
        return code, "unresolved_host"
    if len(matched) > 1:
        raise ValueError(f"Conflicting slot units for {dataset.key}:{column}: {sorted(matched)}")
    return code, next(iter(matched), "gene")


def build_inventory(tables=None, graphs=None, bridges=None, sources=None, catalog=None,
                    raw_paths=None, refusals=(), observations=(), processed_paths=None) -> pd.DataFrame:
    """Reconcile registry declarations with explicitly qualified installed tables.

    ``tables`` keys are (organism, storage unit), with host references keyed as
    protein rather than gene. ``graphs`` and ``bridges`` are keyed by organism.
    Missing raw paths indicate no located local source, not remote inaccessibility;
    no network requests are made. Optional typed Observation records distinguish
    known missingness causes; their counts are records, not independent genes.
    Refusals are dictionaries with source_id, organism, unit, question and reason.
    ``processed_paths`` separately records located processed files; their presence
    does not establish raw-assay availability or an entity mapping.
    """
    tables, graphs, bridges = tables or {}, graphs or {}, bridges or {}
    sources = datasets.registry() if sources is None else list(sources)
    catalog = slots.all_slots() if catalog is None else tuple(catalog)
    raw_paths = raw_paths or {}
    for (code, unit), table in tables.items():
        if code not in organisms.SPACES and code not in organisms.HOST_TABLES:
            raise ValueError(f"Unknown table organism {code!r}")
        identifier = "host_id" if unit == "protein" else "gene_id" if unit == "gene" else "metabolite"
        if identifier in table:
            if table[identifier].isna().any() or table[identifier].duplicated().any():
                raise ValueError("Table denominator requires unique, nonmissing entity identifiers")
            if unit == "gene" and code in organisms.SPACES:
                if not table[identifier].map(organisms.get(code).matches).all():
                    raise ValueError("Table accessions do not match the explicit organism")
    if len({d.key for d in sources}) != len(sources):
        raise ValueError("Dataset keys must be unique")
    records = []
    for dataset in sources:
        partitions = {}
        for column in dataset.columns:
            partitions.setdefault(_projection(dataset, column, catalog), []).append(column)
        if not partitions:
            partitions[(dataset.organism, "source_records")] = []
        for (code, unit), columns in sorted(partitions.items()):
            row = _base(dataset.key, code, unit)
            row.update(registry_organism=dataset.organism, columns=columns, url=dataset.url or "",
                       note=dataset.note, raw_available=bool(raw_paths.get(dataset.key)))
            if processed_paths is not None:
                row['processed_available'] = bool(processed_paths.get(dataset.key))
            declared = pd.DataFrame(columns=columns)
            matching = [s for s in catalog if s.organism == dataset.organism
                        and (slots.declared_columns(declared, s)
                             or set(s.patterns) & set(columns))]
            row["contexts"] = sorted({" / ".join(s.context_path) or s.context for s in matching})
            row["evidence_families"] = sorted({" / ".join(s.evidence_path) for s in matching
                                               if s.evidence_path}) or [dataset.level + " / " + dataset.kind]
            if unit == "pair":
                total, unattributed, missing = 0, 0, []
                graph = graphs.get(code)
                for column in columns:
                    if column.startswith(slots.EDGE_PREFIX):
                        name = column[len(slots.EDGE_PREFIX):]
                        endpoints = (name + "__a", name + "__b")
                        if graph is None or any(key not in graph.files for key in endpoints):
                            missing.append(column)
                        else:
                            a, b = (graph[key] for key in endpoints)
                            if a.ndim != 1 or b.ndim != 1 or len(a) != len(b):
                                raise ValueError(f"Incomplete paired edge records for {name}")
                            total += len(a)
                    else:
                        table = bridges.get(code)
                        if table is None or "bridge" not in table:
                            missing.append(column)
                        else:
                            # Different sources may share the same bridge label;
                            # only source-tagged rows can be attributed to a paper.
                            source_column = "source"
                            if source_column not in table:
                                missing.append(column)
                                unattributed += int((table.bridge == column.split(":", 1)[1]).sum())
                            else:
                                selected = (table.bridge == column.split(":", 1)[1]) & (table[source_column] == dataset.key)
                                total += int(selected.sum())
                row.update(pair_records=total if len(missing) < len(columns) else None,
                           unattributed_pair_records=unattributed or None,
                           missing_columns=missing,
                           status="unattributed" if unattributed else "unavailable" if len(missing) == len(columns) else "partial" if missing else "installed")
            elif unit == "source_records":
                row["status"] = "source_only" if row["raw_available"] or row['processed_available'] else "unavailable"
            else:
                table = tables.get((code, unit))
                if table is None:
                    row["missing_columns"] = columns
                else:
                    present = [c for c in columns if c in table]
                    missing = [c for c in columns if c not in table]
                    mask = pd.Series(False, index=table.index)
                    false_cells = 0
                    for column in present:
                        observed = _present(table[column])
                        mask |= observed
                        if pd.api.types.is_bool_dtype(table[column]):
                            false_cells += int((observed & table[column].eq(False)).sum())
                    count = int(mask.sum()) if present else None
                    row.update(table_rows=len(table), stored_any_rows=count,
                               missing_unknown_rows=len(table) - count if count is not None else None,
                               stored_fraction=count / len(table) if count is not None and len(table) else None,
                               false_cells=false_cells if present else None, missing_columns=missing,
                               status="unavailable" if not present else "partial" if missing else "installed")
            records.append(row)
    for refusal in refusals:
        required = {"source_id", "organism", "unit", "question", "reason"}
        if set(refusal) != required or not all(isinstance(v, str) and v.strip() for v in refusal.values()):
            raise ValueError("A refusal needs a source, organism, unit, question and reason")
        row = _base(refusal["source_id"], refusal["organism"], refusal["unit"])
        row.update(status="rejected", rejection_reason=refusal["reason"], question=refusal["question"])
        records.append(row)
    counts = {}
    for observation in observations:
        key = observation.source_id, observation.organism, observation.entity_type
        counts.setdefault(key, Counter())[observation.missing_state] += 1
        if (observation.missing_state == "observed" and observation.evidence_status == "measured"
                and observation.value is False):
            counts[key]["explicit_negative"] += 1
    mapping = {"not_assayed": "unmeasured_records", "ambiguous_mapping": "unmapped_records",
               "below_detection": "below_detection_records", "failed_qc": "failed_qc_records",
               "unknown": "unknown_missing_records", "observed": "observed_records"}
    mapping["explicit_negative"] = "explicit_negative_records"
    for row in records:
        key = row["source_id"], row["organism"], row["unit"]
        if key in counts and row["status"] != "rejected":
            for state, field in mapping.items():
                row[field] = counts[key][state]
    addressed = {(r["source_id"], r["organism"], r["unit"]) for r in records if r["status"] != "rejected"}
    if set(counts) - addressed:
        raise ValueError("Observations contain sources, organisms or units outside this inventory")
    return pd.DataFrame(records)
