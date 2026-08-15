# 32 — Finish the GPU path: kmeans, DBSCAN, t-SNE

**Status: complete (v0.31.0, 2026-08-14).** cuML resolvers and GPU-first dispatch now cover k-means,
DBSCAN and 2-D t-SNE; 3-D t-SNE and failures fall back loudly. A search pins its backend for its
whole lifetime and records the snapshot in its manifest.

Measured on the RTX 3090 for 100 identical clusterings of 8,140 three-dimensional points: k-means
fell from 37.7670 s on scikit-learn to 13.2551 s on cuML (2.85x; 73% peak utilization), and DBSCAN
from 3.7874 s to 0.3644 s (10.40x; 29% peak utilization). Both paths therefore remain shipped.

## The evidence

A ten-task search was run on the RTX 3090 through `starplast-discover`. Measured while it ran:

    nvidia-smi: 16% GPU utilisation, 2.7 GB of 24 GB used

Sixteen per cent, on a job whose whole purpose is to run a hundred embeddings and clusterings. The
card builds a UMAP in about three seconds and then waits while scikit-learn clusters 8,140 points on
one core.

This is not a limitation, it is a gap in the wiring, and it matters for this workload specifically:
**every winning configuration so far chose kmeans or DBSCAN.**

| run | winning algorithm | clusters | runs on |
|---|---|---|---|
| `bigA_00_guilt_compartment_best` | kmeans | 60 | CPU |
| `bigA_01_guilt_lopit_unified` | kmeans | 60 | CPU |
| `bigA_02_auprc_compartment_best` | DBSCAN | 488 | CPU |
| `bigB_00_guilt_fit_invitro_hff` | DBSCAN | 573 | CPU |

t-SNE is CPU-only in every configuration, and it was chosen by two of the five group-A winners.

## The cause

In `starplast/clustering.py`, `cluster()` returns from the CPU branches **before** the GPU block is
ever consulted. The order today is:

    if algorithm == "dbscan":        -> sklearn DBSCAN, returns
    if algorithm == "kmeans":        -> sklearn KMeans, returns
    if algorithm == "agglomerative": -> sklearn Agglomerative, returns
    on_gpu = gpu.hdbscan_class()     <- only reached for hdbscan

The three CPU branches were added in 0.27.0 above an existing GPU block that only ever handled
HDBSCAN. Nothing is wrong with the branches; they are simply in front of the door.

`starplast/embedding.py` has the same shape for `method == "tsne"`: it goes straight to
`sklearn.manifold.TSNE` with no GPU attempt at all.

## What to do

**1. `starplast/gpu.py`** — add resolvers beside `umap_class()` and `hdbscan_class()`:

    kmeans_class()   -> cuml.cluster.KMeans or None
    dbscan_class()   -> cuml.cluster.DBSCAN or None
    tsne_class()     -> cuml.manifold.TSNE or None

Same contract as the existing two: return None when cuML is absent, never raise, and let
`backend()` report what resolved. Check what this cuML actually exposes before writing the names —
the environment has `cuml-cu12 26.8.0`, and the API has moved between releases.

**2. `starplast/clustering.py`** — try the GPU first for every algorithm, fall back on the CPU, and
say which ran. The existing HDBSCAN block already has the right shape, including the loud warning
on fallback ("a silent fallback is how 'the GPU switch does nothing' becomes impossible to
diagnose") and the `len(X) >= 1000` floor, which exists because cuML below a thousand points is
slower than sklearn once the copy is counted. Keep both. Do not restructure the function into
something clever; a small dispatch table keyed by algorithm is enough.

**3. `starplast/embedding.py`** — the same for t-SNE. Note that cuML's t-SNE has historically
supported only `n_components=2`; if that still holds, fall back to the CPU for 3 components rather
than silently returning a 2-D map, and log the reason.

**4. Say it out loud, as the UMAP path already does.** `embed()` logs which library ran and warns
that a cuML UMAP is a different map from the CPU one rather than the same map faster. The same
sentence is true of cuML kmeans and DBSCAN — different implementations, not accelerated copies —
and a search that mixes them across configurations is comparing libraries as well as
hyperparameters. Either pin the backend for the length of a run and record it in the manifest
(`searches.from_climb` already writes a manifest, add the backend to it), or accept the mixing and
document it. Pinning is the better answer.

## Acceptance

* `nvidia-smi` shows utilisation well above 16% during a `starplast-discover` sweep that picks
  kmeans or DBSCAN. Measure it, do not assume it.
* A timed comparison in the commit message: one identical configuration, CPU vs GPU, wall clock for
  a hundred clusterings of 8,140 points. If a GPU path turns out to be **slower**, remove it rather
  than ship it behind a switch — that decision was already made once for the parts of the GPU work
  that measured slower, and it is the right one.
* `gpu.backend()` reports every resolved class, and the manifest of a saved search records which
  backend produced it.
* 100% coverage on the changed modules, no `pragma`. The GPU branches are testable the way the
  HDBSCAN one already is: monkeypatch the resolver to return a fake class, and separately to raise,
  and assert the fallback logs and returns CPU labels.

## Do not

Do not re-run the ten-task search to "see if it is faster" before the timing test exists. The
searches already saved under `~/.cache/starplast/searches/bigA_*`, `bigB_*` were produced by the
CPU clustering path and are the baseline; overwriting them loses the comparison.
