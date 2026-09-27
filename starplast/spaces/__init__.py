"""Shared validation for building one organism's data space without silent data loss."""
from pathlib import Path

import numpy as np
import pandas as pd


def validate_nodes(space, nodes: pd.DataFrame, previous: pd.DataFrame | None = None) -> None:
    """Require canonical unique IDs; refuse lost old IDs, columns or measured cells."""
    if not nodes.columns.is_unique or "gene_id" not in nodes or nodes.empty:
        raise ValueError("a space needs rows, unique columns and gene_id")
    ids = nodes["gene_id"]
    if ids.isna().any() or not ids.is_unique or not ids.map(space.matches).all():
        raise ValueError(f"{space.code}: gene IDs must be unique, nonmissing and canonical")
    if previous is None or previous.empty:
        return
    # Validate the baseline too; an ambiguous index must not hide a lost measurement.
    validate_nodes(space, previous)
    missing_columns = set(previous) - set(nodes)
    missing_ids = set(previous.gene_id) - set(ids)
    if missing_columns or missing_ids:
        raise ValueError(f"space build loses {len(missing_ids)} IDs or columns {sorted(missing_columns)}")
    old = previous.set_index("gene_id")
    aligned = nodes.set_index("gene_id").reindex(old.index)
    lost = [c for c in old if (old[c].notna() & aligned[c].isna()).any()]
    if lost:
        raise ValueError(f"space build loses measured values: {lost}")


def validate_graph(nodes: pd.DataFrame, graph_path) -> None:
    """Require finite display coordinates, exact gene order and in-range edge indices."""
    with np.load(graph_path, allow_pickle=False) as graph:
        if not {"xyz", "gene_ids"} <= set(graph.files):
            raise ValueError("graph needs xyz and gene_ids")
        if not np.array_equal(graph["gene_ids"].astype(str), nodes.gene_id.astype(str)):
            raise ValueError("graph gene order differs from nodes")
        xyz = graph["xyz"]
        if xyz.shape != (len(nodes), 3) or not np.isfinite(xyz).all():
            raise ValueError("graph needs three finite coordinates per gene")
        for key in graph.files:
            if key.endswith(("__a", "__b")):
                indices = graph[key]
                if (indices.ndim != 1 or not np.issubdtype(indices.dtype, np.integer)
                        or (indices < 0).any() or (indices >= len(nodes)).any()):
                    raise ValueError(f"invalid graph indices: {key}")
                prefix = key.rsplit("__", 1)[0]
                if not {prefix + "__a", prefix + "__b", prefix + "__w"} <= set(graph.files):
                    raise ValueError(f"incomplete edge layer: {prefix}")
                for suffix in ("__a", "__b", "__w", "__r"):
                    if prefix + suffix in graph.files:
                        values = graph[prefix + suffix]
                        if values.shape != indices.shape or not np.isfinite(values).all():
                            raise ValueError(f"invalid edge array: {prefix + suffix}")


def build_pack(space, directory, version: str, column_licenses: dict, destination,
               previous=None) -> str:
    """Validate a prepared space and create its distributable pack; return the archive SHA256.

    Callers supply the per-column license decisions. Unreviewed (V) and local-only (X) sources
    cannot enter a distributed pack. Raw-source acquisition remains the organism builder's job.
    """
    from .. import packs
    directory = Path(directory)
    nodes = pd.read_parquet(directory / space.nodes)
    old = pd.read_parquet(previous) if previous is not None else None
    validate_nodes(space, nodes, old)
    if set(column_licenses) != set(nodes.columns):
        raise ValueError("every node column needs exactly one license declaration")
    files = {space.nodes: directory / space.nodes}
    if space.graph:
        validate_graph(nodes, directory / space.graph)
        files[space.graph] = directory / space.graph
    return packs.build(space, version, files, column_licenses, destination)
