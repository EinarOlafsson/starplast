#!/usr/bin/env python3
"""Searching for structures that recover a held-out label.

This is the point of the whole apparatus. Everything else — identity, missing-value policy, block
weighting, clustering walks — exists so that this question can be asked without fooling ourselves:

    Is there a combination of datasets and hyperparameters whose structure recovers a label that was
    never given to it?

If a map built from expression and fitness alone puts 97% of apicoplast proteins in one cluster, then
compartment is *predictable from that structure*, and the unlabelled genes in that cluster acquire a
testable prediction. The same question for cell-cycle or stage tells you when a hypothetical gene is
likely active. That is inference, not description — and it is only valid because the target was held out.

Three rules make the search trustworthy, and each of them was a bug before it was a rule:

1. **The target and everything that restates it are excluded from the embedding.** Naming the target
   column is not enough: `compartment` fed one embedding and its twin `lopit_map` was duly reported as
   the top discovery at V = 0.96. `excluded_for` resolves the full set by measured association.
2. **Precision and recall are scored separately and both are kept.** A cluster that is 100% apicoplast
   but holds 5% of apicoplast proteins is useless for inference; one holding 97% at 90% precision is the
   goal. A single blended number hides exactly that difference, so the objective is per-label F1 with
   both parts stored beside it.
3. **Every run records what it would take to reproduce it** — the full `EmbeddingSpec`, the seed, the
   clustering parameters, the excluded columns, and the resulting scores. A hit nobody can rebuild is
   not a result.

The search is deliberately exhaustive rather than clever. This space has not been mapped, so the useful
first move is to cover it and look, not to optimise into one corner of it.
"""
from __future__ import annotations

import itertools
import json
import os
import time
from dataclasses import asdict

import numpy as np
import pandas as pd

from .clustering import NOISE, _association_with_inputs, cluster
from .embedding import EmbeddingSpec, build_matrix, columns_for

DEFAULT_SEED = 42

# Targets worth asking about, and the columns that carry them.
TARGETS = {
    "compartment": "compartment",           # hyperLOPIT, 3,827 labelled
    "compartment_best": "compartment_best",  # + orthoLOPIT transfer
    "lopit_unified": "lopit_unified",       # 12 cross-species categories
    "attention_depth": "attention_depth",
}


def excluded_for(nodes: pd.DataFrame, target: str, threshold=0.8) -> set:
    """The target plus every column that substantially restates it.

    The threshold is deliberately stricter than the 0.95 used to flag derived columns after the fact.
    Here we are choosing what an embedding may *see*, and a 0.85-associated column leaks nearly as much
    as an identical one -- `lopit_mcmc` is a second inference over the same experiment at 0.74.
    """
    assoc = _association_with_inputs(nodes, {target})
    out = {target} | {c for c, v in assoc.items() if v >= threshold}
    return out


def _spec_without(spec: EmbeddingSpec, nodes: pd.DataFrame, banned: set) -> EmbeddingSpec:
    """Drop any block or categorical that would feed a banned column into the embedding."""
    keep_blocks = []
    for b in spec.blocks:
        cols = columns_for(nodes, EmbeddingSpec(blocks=(b,))).get(b, [])
        if cols and not (set(cols) & banned):
            keep_blocks.append(b)
    keep_cat = tuple(c for c in spec.categorical if c not in banned)
    return EmbeddingSpec(**{**asdict(spec), "blocks": tuple(keep_blocks), "categorical": keep_cat})


def score_recovery(labels: np.ndarray, truth: pd.Series, min_label=15) -> tuple:
    """How well does the clustering recover `truth`? Per-label precision, recall and F1.

    For each label the best single cluster is taken -- the question is whether the structure isolates
    the label somewhere, not whether the clustering happens to use the same number of groups.
    """
    ok = (labels != NOISE) & truth.notna().to_numpy()
    if ok.sum() < 50:
        return {}, pd.DataFrame()
    lab, tru = labels[ok], truth[ok].astype(str).to_numpy()
    rows = []
    for t in pd.unique(tru):
        n_t = int((tru == t).sum())
        if n_t < min_label:
            continue
        best = None
        for cl in np.unique(lab):
            inter = int(((lab == cl) & (tru == t)).sum())
            if not inter:
                continue
            prec = inter / int((lab == cl).sum())
            rec = inter / n_t
            f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
            if best is None or f1 > best["f1"]:
                best = {"label": t, "cluster": int(cl), "n_label": n_t,
                        "n_in_cluster": inter, "precision": prec, "recall": rec, "f1": f1}
        if best:
            rows.append(best)
    per = pd.DataFrame(rows)
    if per.empty:
        return {}, per
    # Weighted by label size, so a structure that nails one tiny class does not outrank one that
    # organises the whole proteome.
    w = per.n_label / per.n_label.sum()
    summary = {"mean_f1": float((per.f1 * w).sum()),
               "best_f1": float(per.f1.max()),
               "best_label": str(per.loc[per.f1.idxmax(), "label"]),
               "n_labels_scored": int(len(per)),
               "n_labels_recovered": int((per.f1 >= 0.5).sum())}
    return summary, per


def _flushing(log):
    """A long search redirected to a file shows nothing until it ends unless output flushes."""
    def out(m):
        log(m)
        try:
            import sys; sys.stdout.flush()
        except Exception:
            pass
    return out


def search(nodes: pd.DataFrame, target: str = "compartment",
           block_sets=None, na_policies=("median",), scalings=("rank",),
           n_neighbors_values=(15, 50), min_dist_values=(0.0, 0.25),
           min_cluster_sizes=(25, 60), sample_size: int | None = None,
           seed: int = DEFAULT_SEED, store=None, save_above: float | None = None,
           log=print) -> tuple:
    """Walk combinations of datasets and hyperparameters, scoring recovery of a held-out target."""
    log = _flushing(log)
    truth_col = TARGETS.get(target, target)
    if truth_col not in nodes.columns:
        raise ValueError(f"target {truth_col!r} not in the table")
    banned = excluded_for(nodes, truth_col)
    log(f"target {truth_col!r}: {int(nodes[truth_col].notna().sum()):,} labelled genes")
    log(f"  excluded from every embedding ({len(banned)}): {', '.join(sorted(banned))}")

    if block_sets is None:
        base = ["expression_summary", "expression_raw", "fitness_screens",
                "published_screens", "protein_features", "interactions"]
        base = [b for b in base if columns_for(nodes, EmbeddingSpec(blocks=(b,))).get(b)]
        block_sets = [tuple(c) for r in (1, 2, 3) for c in itertools.combinations(base, r)]
        log(f"  {len(block_sets)} dataset combinations from {len(base)} blocks")

    rows, per_label, runs = [], [], 0
    total = (len(block_sets) * len(na_policies) * len(scalings)
             * len(n_neighbors_values) * len(min_dist_values) * len(min_cluster_sizes))
    log(f"  {total} runs")
    t0 = time.time()
    sub = None
    if sample_size and sample_size < len(nodes):
        sub = np.random.default_rng(seed).choice(len(nodes), sample_size, replace=False)
        sub.sort()

    for blocks, pol, sc in itertools.product(block_sets, na_policies, scalings):
        spec0 = EmbeddingSpec(blocks=blocks, na_policy=pol, scaling=sc, random_state=seed)
        spec0 = _spec_without(spec0, nodes, banned)
        if not spec0.blocks:
            continue
        try:
            X, names, keep = build_matrix(nodes, spec0, log=lambda *a: None)
        except ValueError:
            continue
        idx = np.arange(len(nodes))[keep]
        if sub is not None:
            sel = np.isin(idx, sub)
            X, idx = X[sel], idx[sel]
        if len(X) < 200:
            continue
        try:
            import umap
        except ImportError:                                    # pragma: no cover
            log("umap-learn not installed"); return pd.DataFrame(), pd.DataFrame()

        for nn, md in itertools.product(n_neighbors_values, min_dist_values):
            if nn >= len(X):
                continue
            Y = np.asarray(umap.UMAP(n_components=3, n_neighbors=nn, min_dist=md,
                                     metric="euclidean", random_state=seed).fit_transform(X))
            for mcs in min_cluster_sizes:
                lab = cluster(Y, algorithm="hdbscan", min_cluster_size=mcs)
                summary, per = score_recovery(lab, nodes[truth_col].iloc[idx])
                runs += 1
                if not summary:
                    continue
                row = {"target": truth_col, "blocks": "+".join(spec0.blocks),
                       "na_policy": pol, "scaling": sc, "n_neighbors": nn, "min_dist": md,
                       "min_cluster_size": mcs, "seed": seed, "n_genes": int(len(X)),
                       "n_features": X.shape[1],
                       "n_clusters": int(len(set(lab[lab != NOISE]))),
                       "noise_frac": float((lab == NOISE).mean()), **summary}
                rows.append(row)
                per = per.assign(**{k: row[k] for k in
                                    ("blocks", "na_policy", "scaling", "n_neighbors",
                                     "min_dist", "min_cluster_size", "target")})
                per_label.append(per)
                # Default to saving the best runs relative to what this search actually found,
                # not an absolute bar. A fixed 0.6 threshold saved nothing at all on the first
                # real search, whose best was 0.592 -- the gate silently discarded the answer.
                keep_at = save_above if save_above is not None else 0.0
                if store is not None and summary["best_f1"] >= keep_at:
                    name = (f"{truth_col}_{'+'.join(spec0.blocks)}_{pol}_{sc}"
                            f"_nn{nn}_md{md}_mcs{mcs}")
                    store.save(name, Y, spec0, gene_ids=nodes.gene_id.iloc[idx], features=names,
                               extra={"target": truth_col, "excluded": sorted(banned),
                                      "clustering": {"algorithm": "hdbscan",
                                                     "min_cluster_size": mcs},
                                      "scores": summary})
        if runs and runs % 40 == 0:
            log(f"    {runs}/{total} runs, {time.time() - t0:.0f}s")

    R = pd.DataFrame(rows).sort_values("mean_f1", ascending=False) if rows else pd.DataFrame()
    if store is not None and not R.empty and save_above is None:
        log(f"  saved embeddings for every run; the top result is the first row of the table")
    P = pd.concat(per_label, ignore_index=True) if per_label else pd.DataFrame()
    log(f"  done: {runs} runs in {time.time() - t0:.0f}s")
    if not R.empty:
        b = R.iloc[0]
        log(f"  best: mean F1 {b.mean_f1:.3f} (best label {b.best_label} at F1 {b.best_f1:.3f}) "
            f"from {b.blocks} nn={b.n_neighbors} md={b.min_dist} mcs={b.min_cluster_size}")
    return R, P


def predictions(nodes: pd.DataFrame, labels: np.ndarray, truth: pd.Series, gene_index,
                min_precision=0.8, min_cluster_labelled=10) -> pd.DataFrame:
    """Turn a recovered structure into predictions for the unlabelled genes in each cluster.

    Only clusters that are already purely one label are used, and the precision achieved on the
    *labelled* members is carried through as the prediction's confidence -- it is the only honest
    estimate of how often the prediction will be right.
    """
    out = []
    tru = truth.to_numpy()
    for cl in np.unique(labels[labels != NOISE]):
        m = labels == cl
        known = m & pd.notna(tru)
        if known.sum() < min_cluster_labelled:
            continue
        vals, counts = np.unique(tru[known].astype(str), return_counts=True)
        top, n_top = vals[counts.argmax()], counts.max()
        prec = n_top / known.sum()
        if prec < min_precision:
            continue
        for gi in np.asarray(gene_index)[m & pd.isna(tru)]:
            out.append({"gene_id": gi, "predicted": top, "cluster": int(cl),
                        "cluster_precision": float(prec),
                        "n_labelled_in_cluster": int(known.sum())})
    return pd.DataFrame(out)
