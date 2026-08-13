#!/usr/bin/env python3
"""Hyperparameter walks, saved embeddings, and importing your own table.

Three jobs that belong together because they are what a user does before trusting a map:

* **`walk_umap_iter`** — search `n_neighbors` x `min_dist` on a seeded random subsample, so a sweep costs
  seconds rather than the minutes a full 8,140-gene UMAP takes, and **yield each configuration as it is
  computed** rather than a table at the end. `walk_umap` is that collected and ranked, for callers who
  want the answer rather than the process. The subsample is drawn with a fixed seed so two sweeps are
  comparable; the seed is recorded in every row.
* **`EmbeddingStore`** — save a coordinate set together with the full `EmbeddingSpec` that produced it.
  An embedding without its recipe cannot be compared with another or reported in a methods section.
* **`import_table`** — read a user's CSV and align it to gene ids, with regex repair when the identifier
  column is formatted differently.

**On scoring a UMAP.** There is no ground truth for "a good embedding", so the walk reports several
quantities that measure different things and lets the user decide, rather than manufacturing one number:

    trustworthiness      does the map preserve the original neighborhoods (0-1, higher better)
    continuity_proxy     correlation between original and embedded pairwise distances
    n_clusters_hdbscan   how many clusters fall out, and how much is noise

A setting that maximises trustworthiness often produces one blob; a setting that produces many clusters
often distorts. Both numbers are shown because that trade-off is the actual decision.
"""
from __future__ import annotations

import json
import os
import re
from dataclasses import asdict, dataclass

import numpy as np
import pandas as pd

from .embedding import RENAMED_BLOCKS, EmbeddingSpec, build_matrix, normalize

DEFAULT_SEED = 42
#: The min_cluster_size the walk's own cluster check uses. Named because the interface rebuilds a
#: walk configuration and clusters it again to show what the row's `n_clusters_hdbscan` counted --
#: with a different value there, the map on screen would not be the map the number describes.
WALK_MIN_CLUSTER_SIZE = 15


# --------------------------------------------------------------------------- UMAP walk
def _subsample(n_rows: int, size: int, seed: int) -> np.ndarray:
    if size >= n_rows:
        return np.arange(n_rows)
    return np.random.default_rng(seed).choice(n_rows, size=size, replace=False)


def _quality(X: np.ndarray, Y: np.ndarray, n_neighbors: int) -> dict:
    from sklearn.manifold import trustworthiness
    out = {}
    k = int(min(max(n_neighbors, 5), (len(X) - 1) // 2))
    try:
        out["trustworthiness"] = float(trustworthiness(X, Y, n_neighbors=k))
    except Exception:
        out["trustworthiness"] = np.nan
    # Distance correlation on a bounded random subset -- the full matrix is quadratic and unnecessary.
    # The distances go to the GPU when there is one and the matrix is big enough to pay for the
    # copy: this runs once per configuration, so a 288-point sweep runs it 288 times.
    try:
        from scipy.stats import spearmanr
        from . import gpu
        idx = _subsample(len(X), min(600, len(X)), DEFAULT_SEED)
        a, b = _condensed(X[idx]), _condensed(Y[idx])
        out["continuity_proxy"] = float(spearmanr(a, b).statistic)
    except Exception:
        out["continuity_proxy"] = np.nan
    return out


def _condensed(A):
    """The upper triangle of A's distance matrix, from the GPU where that is worth it."""
    from . import gpu
    if gpu.worth_it(np.asarray(A)):
        D = gpu.pairwise_distances(A)
        iu = np.triu_indices(len(D), k=1)
        return D[iu]
    from scipy.spatial.distance import pdist
    return pdist(A)


@dataclass
class WalkStep:
    """One configuration of a walk, complete: its scores and the map that produced them.

    The walk used to hand back a table when the whole sweep had finished, so a 288-configuration
    sweep -- half an hour -- showed nothing at all until it ended, and the embeddings themselves were
    computed, scored and thrown away. A step carries the coordinates as well as the numbers, which is
    what lets configuration 12 be looked at while 13 is still computing.

    `genes` is a boolean mask over the node table rather than a count, because a walk embeds a
    subsample: without it there is no way to say which gene each coordinate belongs to, and a map
    whose points cannot be named is a picture rather than a map.
    """
    index: int                   # 1-based, so it reads as "7 of 20" without arithmetic
    total: int                   # configurations this walk expects to run, at most
    row: dict                    # the scores, exactly as they go into the table
    coords: np.ndarray           # (n_sampled, n_components), at the same scale as any other map
    genes: np.ndarray            # boolean mask over the node table: which genes have a position
    spec: EmbeddingSpec          # the full recipe, carrying THIS configuration's hyperparameters
    name: str = ""               # key it was stored under, or "" if the walk was given no store

    @property
    def label(self) -> str:
        """The configuration in one line, for a thumbnail caption or a status message."""
        return f"n_neighbors={self.spec.n_neighbors}  min_dist={self.spec.min_dist:g}"


def walk_umap_iter(nodes: pd.DataFrame, spec: EmbeddingSpec,
                   n_neighbors_values=(5, 15, 25, 50, 100),
                   min_dist_values=(0.0, 0.1, 0.25, 0.5),
                   sample_size: int = 2000, seed: int = DEFAULT_SEED,
                   cluster_check: bool = True, store=None, should_stop=None, log=print):
    """Sweep UMAP hyperparameters, yielding each configuration as it finishes.

    This is the walk; `walk_umap` is this collected into a table. Emitting per configuration is what
    makes a sweep watchable -- rows and thumbnails appear one at a time rather than all at once at
    the end -- and it costs nothing, because the embedding was being built anyway.

    Given a `store`, every configuration is saved with its full recipe as it is computed, so a walk
    that is stopped half way still leaves behind everything it had finished.

    `should_stop` is asked before each configuration and simply ends the sweep: the caller keeps
    every step already yielded. A stop is a decision that enough has been seen rather than an error,
    so nothing is raised and nothing is discarded.
    """
    X, names, rows = build_matrix(nodes, spec, log=lambda *a: None)
    # Sorted, so the coordinates line up with the node table. `rng.choice` returns its picks in
    # random order, which made row i of the embedding an arbitrary gene: any attempt to say which
    # gene a point is -- a mask, a gene_id list, clicking a point -- silently named the wrong one.
    # Sorting picks the same genes, in the order the table has them. search.py already does this.
    take = np.sort(_subsample(X.shape[0], sample_size, seed))
    Xs = X[take]
    # Node-table positions of the sampled genes, so a coordinate can be traced back to a gene.
    where = np.arange(len(nodes))[rows][take]
    genes = np.zeros(len(nodes), dtype=bool)
    genes[where] = True
    gene_ids = (nodes.gene_id.to_numpy()[where] if "gene_id" in nodes.columns else None)
    total = len(n_neighbors_values) * len(min_dist_values)
    log(f"walk_umap: {len(Xs):,} of {X.shape[0]:,} genes (seed {seed}), "
        f"{total} settings, {X.shape[1]} features")

    try:
        import umap
    except ImportError:
        log("umap-learn not installed")
        return

    i = 0
    for nn in n_neighbors_values:
        if nn >= len(Xs):
            continue
        for md in min_dist_values:
            if should_stop is not None and should_stop():
                log(f"  stopped after {i} of {total} configurations -- each one is already saved")
                return
            Y = np.asarray(umap.UMAP(n_components=spec.n_components, n_neighbors=nn, min_dist=md,
                                     metric=spec.metric,
                                     random_state=spec.random_state).fit_transform(Xs))
            row = {"n_neighbors": nn, "min_dist": md, "seed": seed,
                   "sample_size": len(Xs), **_quality(Xs, Y, nn)}
            if cluster_check:
                from .clustering import cluster, NOISE
                lab = cluster(Y, algorithm="hdbscan", min_cluster_size=WALK_MIN_CLUSTER_SIZE)
                row["n_clusters_hdbscan"] = int(len(set(lab[lab != NOISE])))
                row["noise_frac"] = float((lab == NOISE).mean())
            i += 1
            # Scored on the raw output and displayed from the normalized one. Normalizing is a
            # uniform move-and-scale, so it cannot change trustworthiness or the clustering, but
            # computing the numbers first keeps that guarantee obvious rather than argued.
            used = EmbeddingSpec(**{**asdict(spec), "n_neighbors": int(nn), "min_dist": float(md)})
            name = ""
            if store is not None:
                name = f"walk_nn{nn}_md{md:g}_seed{seed}"
                store.save(name, Y, used, gene_ids=gene_ids, features=names,
                           extra={"walk": True, "scores": row})
            log(f"  n_neighbors={nn:4d} min_dist={md:<5} "
                f"trust={row['trustworthiness']:.3f} "
                f"clusters={row.get('n_clusters_hdbscan', '-')}")
            yield WalkStep(index=i, total=total, row=row, coords=normalize(Y), genes=genes,
                           spec=used, name=name)


def walk_umap(nodes: pd.DataFrame, spec: EmbeddingSpec,
              n_neighbors_values=(5, 15, 25, 50, 100),
              min_dist_values=(0.0, 0.1, 0.25, 0.5),
              sample_size: int = 2000, seed: int = DEFAULT_SEED,
              cluster_check: bool = True, store=None, on_step=None,
              should_stop=None, log=print) -> pd.DataFrame:
    """Sweep UMAP hyperparameters on a seeded subsample of the genes, and rank the settings.

    `on_step` is called with each `WalkStep` the moment it is computed, which is how the interface
    fills its table and its gallery a configuration at a time. The ranked table is still returned,
    because ranking needs the whole sweep and a caller that only wants the answer should not have to
    accumulate it.
    """
    out = []
    for step in walk_umap_iter(nodes, spec, n_neighbors_values=n_neighbors_values,
                               min_dist_values=min_dist_values, sample_size=sample_size,
                               seed=seed, cluster_check=cluster_check, store=store,
                               should_stop=should_stop, log=log):
        if on_step is not None:
            on_step(step)
        out.append(step.row)
    if not out:
        # Same shape as a populated result. Built from [], the frame has no columns at all and sorting
        # raises KeyError -- so a grid where every setting was skipped crashed instead of reporting
        # that nothing was runnable.
        return pd.DataFrame(columns=["n_neighbors", "min_dist", "seed", "sample_size",
                                     "trustworthiness", "continuity_proxy"])
    return pd.DataFrame(out).sort_values("trustworthiness", ascending=False)


# --------------------------------------------------------------------------- saved embeddings
class EmbeddingStore:
    """Named embeddings on disk: coordinates plus the recipe that produced them."""

    def __init__(self, root: str):
        self.root = root
        os.makedirs(root, exist_ok=True)

    def _paths(self, name: str):
        safe = re.sub(r"[^A-Za-z0-9._-]", "_", name)
        return (os.path.join(self.root, f"{safe}.npz"),
                os.path.join(self.root, f"{safe}.json"))

    def save(self, name: str, coords: np.ndarray, spec: EmbeddingSpec,
             gene_ids=None, features=None, extra: dict | None = None) -> str:
        """Store coordinates with the full recipe that produced them."""
        npz, meta = self._paths(name)
        arrays = {"xyz": np.asarray(coords, dtype=np.float32)}
        if gene_ids is not None:
            arrays["gene_id"] = np.asarray(gene_ids, dtype=object)
        np.savez_compressed(npz, **arrays)
        json.dump({"name": name, "spec": asdict(spec), "n_genes": int(len(coords)),
                   "features": list(features or []), **(extra or {})},
                  open(meta, "w"), indent=1)
        return npz

    def load(self, name: str):
        """Reload a stored embedding, its gene ids and its recipe."""
        npz, meta = self._paths(name)
        if not os.path.exists(npz):
            raise FileNotFoundError(name)
        z = np.load(npz, allow_pickle=True)
        m = json.load(open(meta)) if os.path.exists(meta) else {}
        spec = EmbeddingSpec.from_dict(m.get("spec", {})) if m.get("spec") else None
        # The gene ids are saved into the npz and were not returned, so a reopened embedding came back
        # as coordinates with no way to say which gene each row is. That makes it useless for anything
        # gene-specific -- including turning a recovered structure into named predictions, which is the
        # reason for saving embeddings at all.
        if "gene_id" in z.files:
            m["gene_ids"] = [str(g) for g in z["gene_id"]]
        return z["xyz"], spec, m

    def list(self) -> pd.DataFrame:
        """Every stored embedding, newest first."""
        rows = []
        for f in sorted(os.listdir(self.root)):
            if not f.endswith(".json"):
                continue
            m = json.load(open(os.path.join(self.root, f)))
            s = m.get("spec", {})
            rows.append({"name": m.get("name", f[:-5]), "n_genes": m.get("n_genes"),
                         # Through the rename table, so the list says what a block is called now
                         # rather than what it was called on the day each file was written -- one
                         # block under two spellings reads as two blocks.
                         "blocks": ",".join(RENAMED_BLOCKS.get(b, b) for b in s.get("blocks", [])),
                         "na_policy": s.get("na_policy"), "scaling": s.get("scaling"),
                         "n_neighbors": s.get("n_neighbors"), "min_dist": s.get("min_dist")})
        return pd.DataFrame(rows)


# --------------------------------------------------------------------------- import
GENE_RX = re.compile(r"TGME49_\d{5,6}", re.I)
# What import_table EXTRACTS. Wider than GENE_RX on purpose: a table keyed on TGGT1_ or TGVEG_ has to
# reach `resolve` to be mapped forward, and extracting with the ME49-only pattern discarded those rows
# before the identity layer ever saw them -- silently, and in exactly the case the docstring promises
# is handled. TGGT1_ accessions are more common than TGME49_ in published supplements.
IMPORT_RX = re.compile(r"TG(?:ME49|GT1|VEG)_\d{5,6}", re.I)
# Deliberately loose: catches TgME49.208830, TGME49-208830, tgme49 208830, bare 6-digit accessions.
# Used only to *suggest* a column, never to parse one -- a malformed column is precisely the case where
# the user needs a suggestion, so refusing to guess is unhelpful.
LOOSE_RX = re.compile(r"(tgme49|tggt1|tgveg)\s*[._\- ]?\s*\d{5,6}", re.I)


def suggest_gene_column(df: pd.DataFrame):
    """Return (column, fraction_matching, needs_repair) for the likeliest identifier column."""
    best = (None, 0.0, False)
    for c in df.columns:
        s = df[c].astype(str)
        exact = float(s.str.contains(GENE_RX, regex=True, na=False).mean())
        loose = float(s.map(lambda x: bool(LOOSE_RX.search(str(x)))).mean())
        score = max(exact, loose)
        if score > best[1]:
            best = (c, score, exact < loose)
    return best if best[1] > 0.1 else (None, 0.0, False)


def suggest_repair(values) -> str | None:
    """A regex/replacement pair that would turn a malformed identifier column into canonical form."""
    s = pd.Series(values).astype(str)
    for pat, rep in ((r"(?i)tgme49\s*[.\- ]\s*(\d{5,6})", r"TGME49_\1"),
                     (r"(?i)\btgme49(\d{5,6})\b", r"TGME49_\1"),
                     (r"(?i).*?(tgme49_\d{5,6}).*", r"\1")):
        fixed = s.str.replace(pat, rep, regex=True)
        if fixed.str.contains(GENE_RX, regex=True, na=False).mean() > 0.5:
            return pat, rep
    return None


def import_table(path_or_df, gene_column: str | None = None, pattern: str | None = None,
                 replacement: str = r"\g<0>", resolve=None, prefix: str = "imported_",
                 log=print) -> pd.DataFrame:
    """Read a user table and return it keyed by canonical gene id.

    `pattern` is applied to the identifier column before matching, so a column formatted
    `TgME49.208830` or `gene|TGME49_208830|v2` can be repaired without editing the file:

        import_table(df, "id", pattern=r"TgME49\\.(\\d+)", replacement=r"TGME49_\\1")

    `resolve` should be `identity.GeneIndex`-backed, so previous and strain accessions map to current
    ones -- the failure that cost a published screen every one of its rows.
    """
    df = path_or_df if isinstance(path_or_df, pd.DataFrame) else (
        pd.read_csv(path_or_df, sep=None, engine="python"))
    col = gene_column or suggest_gene_column(df)[0]
    if col is None:
        raise ValueError("no column looks like a gene identifier; pass gene_column=")

    raw = df[col].astype(str)
    if pattern:
        raw = raw.str.replace(pattern, replacement, regex=True)
    ids = raw.str.extract(f"({IMPORT_RX.pattern})", expand=False, flags=re.I)
    ids = ids.str.upper().str.replace("TGME49_", "TGME49_", regex=False)
    if resolve is not None:
        ids = ids.map(lambda x: resolve(x) if isinstance(x, str) else x)

    matched = ids.notna()
    log(f"import: {matched.sum():,} of {len(df):,} rows resolved to a gene id "
        f"(column {col!r}" + (f", pattern {pattern!r}" if pattern else "") + ")")
    if not matched.any():
        log("  nothing matched -- check the column, or supply a regex to reformat it")

    out = df.loc[matched].copy()
    out.insert(0, "gene_id", ids[matched].values)
    num = [c for c in out.columns
           if c != "gene_id" and pd.api.types.is_numeric_dtype(out[c])]
    out = out.groupby("gene_id")[num].mean()
    out.columns = [f"{prefix}{c}" for c in out.columns]
    return out
