#!/usr/bin/env python3
"""Rebuild bundled display coordinates while preserving every existing edge array."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def rebuild(directory, prefix=""):
    """Atomically replace one layout, recording the node snapshot and actual execution."""
    from starplast.build_graph import embed
    directory = Path(directory)
    node_path = directory / f"{prefix}nodes.parquet"
    graph_path = directory / f"{prefix}graph.npz"
    nodes = pd.read_parquet(node_path)
    with np.load(graph_path, allow_pickle=False) as previous:
        arrays = {key: previous[key] for key in previous.files}
    if len(arrays["xyz"]) != len(nodes):
        raise ValueError("graph and node table differ in length")
    if "gene_ids" in arrays and not np.array_equal(arrays["gene_ids"], nodes.gene_id.astype(str)):
        raise ValueError("graph gene order differs from node table; rebuild edges first")
    coords, record = embed(nodes, return_metadata=True)
    if record["executed_method"] != "umap":
        raise RuntimeError(f"refusing to ship a fallback layout: {record['fallback']}")
    record["node_file_sha256"] = hashlib.sha256(node_path.read_bytes()).hexdigest()
    arrays.update(xyz=np.asarray(coords), gene_ids=nodes.gene_id.to_numpy(dtype=str),
                  layout_metadata=np.array(json.dumps(record)))
    temporary = graph_path.with_suffix(".tmp.npz")
    np.savez_compressed(temporary, **arrays)
    os.replace(temporary, graph_path)
    print(f"{graph_path.name}: {len(nodes):,} genes, {len(record['features'])} features, "
          f"{record['executed_method']} via {record['backend']}", flush=True)


def main():
    """Regenerate both organism layouts through the public recipe builder."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path(__file__).resolve().parents[1]/"starplast/data")
    args = parser.parse_args()
    for prefix in ("", "pf_"):
        rebuild(args.data, prefix)


if __name__ == "__main__":
    main()
