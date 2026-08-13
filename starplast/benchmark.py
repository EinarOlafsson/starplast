#!/usr/bin/env python3
"""Building the same map twice, on the CPU and on the GPU, and showing what differs.

Two questions a switch labelled "GPU acceleration" raises and does not answer: how much faster, and
is it the same map. The second matters more here. cuml's UMAP is not umap-learn -- different
implementation, different random stream -- so turning the switch on does not speed a map up, it
produces **a different map of the same data**. That is a fine thing to work with and a terrible
thing to discover halfway through a comparison, so this runs both and puts them side by side.

The number that says how different: **k-nearest-neighbour overlap**. Two embeddings can look alike
and place different genes together, and they can look unlike -- rotated, mirrored, differently
scaled -- while preserving every neighbourhood. Coordinates cannot be compared directly; who is
next to whom can, and that is what the map is read for.
"""
from __future__ import annotations

import os
import time

import numpy as np


def knn_overlap(A, B, k: int = 15) -> float:
    """Fraction of each point's k nearest neighbours shared between two embeddings.

    1.0 is the same neighbourhoods throughout; 0.0 is no gene keeping a single neighbour. Rotation,
    reflection and scale do not touch it, which is what makes it the right comparison for two
    embeddings that were never going to have the same coordinates.
    """
    from scipy.spatial.distance import squareform, pdist
    A, B = np.asarray(A, dtype=float), np.asarray(B, dtype=float)
    n = len(A)
    if n < 3 or len(B) != n:
        return float("nan")
    k = int(min(max(k, 1), n - 1))
    da, db = squareform(pdist(A)), squareform(pdist(B))
    np.fill_diagonal(da, np.inf)
    np.fill_diagonal(db, np.inf)
    na = np.argsort(da, axis=1)[:, :k]
    nb = np.argsort(db, axis=1)[:, :k]
    shared = [len(set(na[i]).intersection(nb[i])) for i in range(n)]
    return float(np.mean(shared) / k)


def _embed_with(nodes, spec, want_gpu: bool, log):
    """One embedding, with the GPU forced on or off for the duration, and the seconds it took."""
    from . import embedding, gpu
    before = os.environ.get(gpu.ENV_GPU)
    os.environ[gpu.ENV_GPU] = "1" if want_gpu else "0"
    try:
        t0 = time.perf_counter()
        Y, names, rows = embedding.embed(nodes, spec, log=log)
        seconds = time.perf_counter() - t0
        return {"coords": np.asarray(Y), "seconds": seconds, "backend": gpu.backend()["umap"],
                "genes": int(np.sum(rows)), "features": len(names)}
    finally:
        if before is None:
            os.environ.pop(gpu.ENV_GPU, None)
        else:
            os.environ[gpu.ENV_GPU] = before


def compare(nodes, spec=None, sample: int = 2000, seed: int = 42, log=print) -> dict:
    """Build the same embedding both ways and report the times and the difference.

    Bounded by `sample` because this is a demonstration, not a run: the point is to see the two maps
    beside each other and the two clocks beside each other, and waiting four minutes for the full
    proteome to answer "are these different" is not seeing it sooner.
    """
    from . import gpu
    from .embedding import EmbeddingSpec
    spec = spec or EmbeddingSpec()
    if sample and sample < len(nodes):
        idx = np.sort(np.random.default_rng(seed).choice(len(nodes), sample, replace=False))
        nodes = nodes.iloc[idx]
    out = {"sample": int(len(nodes)), "gpu_available": bool(gpu.available()["cuml"])}

    log(f"comparing on {len(nodes):,} genes")
    cpu = _embed_with(nodes, spec, False, log)
    out["cpu"] = cpu
    log(f"  CPU: {cpu['backend']} in {cpu['seconds']:.1f}s")

    if not out["gpu_available"]:
        # Said plainly rather than shown as a missing bar: torch and cupy do array work and neither
        # implements UMAP, so with only those installed there is no second map to compare.
        out["note"] = ("no cuml on this machine, so there is no GPU UMAP to compare against -- "
                       "pip install starplast-gpu")
        return out
    gpu_run = _embed_with(nodes, spec, True, log)
    out["gpu"] = gpu_run
    log(f"  GPU: {gpu_run['backend']} in {gpu_run['seconds']:.1f}s")
    out["speedup"] = cpu["seconds"] / max(gpu_run["seconds"], 1e-9)
    out["knn_overlap"] = knn_overlap(cpu["coords"], gpu_run["coords"])
    out["note"] = (f"{out['speedup']:.1f}x faster, and {out['knn_overlap']:.0%} of each gene's 15 "
                   f"nearest neighbours are the same in both. The maps are not the same map: they "
                   f"are two maps of the same data.")
    log("  " + out["note"])
    return out
