#!/usr/bin/env python3
"""Leakage-closed GPU UMAP/clustering walk with nested label holdouts.

The fixed proteome is treated transductively: all genes contribute their non-target features to each
unsupervised embedding, while labels are used only to learn cluster calls. Hyperparameters are selected
on inner validation labels and evaluated once on untouched outer labels. The full selection procedure is
also repeated after label permutation, controlling the optimism from searching many structures.
"""
from __future__ import annotations

import argparse
import json
import os
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np
import pandas as pd

from starplast import gpu, paths, slots
from starplast.clustering import cluster
from starplast.embedding import EmbeddingSpec, embed
from starplast.holdout_cv import (choose_final_structure, inference_table, nested_structure_cv,
                                  pooled_category_scores)
from starplast.search import excluded_for


def feature_sets(nodes: pd.DataFrame, banned: set) -> dict[str, tuple[str, ...]]:
    """Independent slot collections offered to the walk, with no banned column in any slot."""
    eligible = []
    coverage = {}
    for slot in slots.all_slots("Tg"):
        columns = slots.source_columns(nodes, slot)
        if slot.role != "feature" or not columns or set(columns) & banned:
            continue
        eligible.append(slot)
        coverage[slot.key] = float(nodes[list(columns)].notna().any(axis=1).mean())
    candidates = {
        "all_independent": tuple(slot.key for slot in eligible),
        "measured_only": tuple(slot.key for slot in eligible
                               if slot.evidence_path and slot.evidence_path[0] in
                               ("molecular measurements", "perturbational measurements")),
        "broad_coverage": tuple(slot.key for slot in eligible if coverage[slot.key] >= .20),
        "non_transcription": tuple(slot.key for slot in eligible
                                   if not slot.biology_path or slot.biology_path[0] != "gene expression"),
    }
    out = {}
    seen = set()
    for name, block_set in candidates.items():
        if block_set and block_set not in seen:
            out[name] = block_set
            seen.add(block_set)
    return out


def walk_structures(nodes, target, scope, *, quick=False, log=print):
    """Build GPU embeddings and clusterings, returning labels plus an auditable config table."""
    banned = excluded_for(nodes, target, scope=scope)
    offered = feature_sets(nodes, banned)
    neighbors = (25, 60) if quick else (15, 30, 60, 120)
    min_dists = (.0, .2)
    hdbscan_grid = ((30, None, "eom"), (75, 10, "eom")) if quick else tuple(
        (mcs, ms, method) for mcs in (20, 50, 100, 200)
        for ms in (None, 10) for method in ("eom", "leaf"))
    kmeans_grid = (20, 40) if quick else (20, 40, 80)
    structures, rows, coords_by_embedding = {}, [], {}
    total_embeddings = len(offered) * len(neighbors) * len(min_dists)
    embedding_number = 0
    started = time.time()
    for set_name, blocks in offered.items():
        for n_neighbors in neighbors:
            for min_dist in min_dists:
                embedding_number += 1
                embedding_id = f"e{embedding_number:03d}"
                spec = EmbeddingSpec(
                    name=f"{target}__{scope}__{set_name}", blocks=blocks,
                    na_policy="median", scaling="rank", method="umap", n_components=3,
                    n_neighbors=n_neighbors, min_dist=min_dist, random_state=42)
                log(f"[{target}/{scope}] embedding {embedding_number}/{total_embeddings}: "
                    f"{set_name}, nn={n_neighbors}, min_dist={min_dist}")
                coords, feature_names, kept = embed(nodes, spec, log=lambda message: log(f"  {message}"))
                if not np.asarray(kept).all():
                    raise RuntimeError("the validation walk requires one coordinate per gene")
                coords_by_embedding[embedding_id] = np.asarray(coords, dtype=np.float32)
                common = {"embedding_id": embedding_id, "feature_set": set_name,
                          "n_slots": len(blocks), "n_features": len(feature_names),
                          "n_neighbors": n_neighbors, "min_dist": min_dist,
                          "embedding_seconds": time.time() - started}
                for mcs, min_samples, method in hdbscan_grid:
                    structure_id = (f"{embedding_id}__hdb_mcs{mcs}_ms{min_samples or 0}_{method}")
                    labels = cluster(coords, algorithm="hdbscan", min_cluster_size=mcs,
                                     min_samples=min_samples, cluster_selection_method=method)
                    structures[structure_id] = np.asarray(labels, dtype=np.int32)
                    real = labels[labels >= 0]
                    rows.append({"structure_id": structure_id, **common, "algorithm": "hdbscan",
                                 "min_cluster_size": mcs, "min_samples": min_samples,
                                 "selection_method": method,
                                 "n_clusters": len(np.unique(real)),
                                 "noise_fraction": float((labels < 0).mean())})
                for k in kmeans_grid:
                    structure_id = f"{embedding_id}__kmeans_k{k}"
                    labels = cluster(coords, algorithm="kmeans", n_clusters=k, random_state=42)
                    structures[structure_id] = np.asarray(labels, dtype=np.int32)
                    rows.append({"structure_id": structure_id, **common, "algorithm": "kmeans",
                                 "n_clusters": k, "noise_fraction": 0.0})
                try:
                    import cupy
                    cupy.get_default_memory_pool().free_all_blocks()
                except Exception:
                    pass
    return structures, pd.DataFrame(rows), coords_by_embedding, sorted(banned), offered


def permutation_test(structures, truth, splits, observed, *, permutations=25, seed=907):
    """Repeat nested model selection after shuffling labels; embeddings are legitimately reused."""
    rng = np.random.default_rng(seed)
    values = truth.astype("object").to_numpy()
    known = np.unique(np.concatenate([np.r_[train, validation, test]
                                      for train, validation, test in splits]))
    null = []
    for permutation in range(permutations):
        shuffled = values.copy()
        shuffled[known] = rng.permutation(shuffled[known])
        outer, _, _, _ = nested_structure_cv(
            structures, pd.Series(shuffled, index=truth.index, name=truth.name),
            folds=len(splits), seed=42, split_labels=splits)
        null.append(float(outer.macro_f1.mean()))
    null = np.asarray(null)
    p = float((1 + (null >= observed).sum()) / (1 + len(null)))
    return null, p


def run_experiment(nodes, target, scope, output, *, quick=False, permutations=25, log=print):
    """Run, save and summarize one target/scope experiment."""
    exp = output / f"{target}__{scope}"
    exp.mkdir(parents=True, exist_ok=True)
    structures, configs, coordinates, banned, offered = walk_structures(
        nodes, target, scope, quick=quick, log=log)
    truth = nodes[target]
    outer, per, validation, splits = nested_structure_cv(structures, truth, folds=5, seed=42)
    category = pooled_category_scores(per)
    final_id = choose_final_structure(validation)
    observed = float(outer.macro_f1.mean())
    null, permutation_p = permutation_test(
        structures, truth, splits, observed, permutations=permutations)
    final_embedding = final_id.split("__", 1)[0]
    predictions = inference_table(structures[final_id], truth, nodes.gene_id, category)

    configs.to_csv(exp / "structures.csv", index=False)
    outer.to_csv(exp / "outer_folds.csv", index=False)
    per.to_csv(exp / "outer_per_category.csv", index=False)
    category.to_csv(exp / "category_summary.csv", index=False)
    validation.to_csv(exp / "inner_validation_scores.csv", index=False)
    predictions.to_csv(exp / "inferred_genes.csv", index=False)
    pd.DataFrame({"permuted_nested_macro_f1": null}).to_csv(exp / "permutation_null.csv",
                                                            index=False)
    np.savez_compressed(exp / "best_structure.npz", coordinates=coordinates[final_embedding],
                        clusters=structures[final_id], gene_id=nodes.gene_id.astype(str).to_numpy())
    manifest = {
        "target": target, "scope": scope, "validation": "5-fold nested stratified CV",
        "embedding_mode": "transductive unsupervised fixed-proteome embedding",
        "selection_labels": "inner validation only", "reported_labels": "locked outer folds",
        "cluster_call_min_support": 5, "cluster_call_min_probability": .55,
        "permutations": permutations, "permutation_p": permutation_p,
        "nested_macro_f1_mean": observed,
        "nested_macro_f1_sd": float(outer.macro_f1.std()),
        "nested_coverage_mean": float(outer.coverage.mean()),
        "final_structure": final_id, "n_structures": len(structures),
        "n_predictions": len(predictions), "excluded_columns": banned,
        "feature_sets": {name: list(blocks) for name, blocks in offered.items()},
        "backend": gpu.backend(), "seed": 42,
    }
    (exp / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default="results/headless_gpu_cv_2026_08_14")
    parser.add_argument("--quick", action="store_true", help="small smoke-test grid")
    parser.add_argument("--permutations", type=int, default=25)
    args = parser.parse_args()
    os.environ["STARPLAST_GPU"] = "1"
    backend = gpu.backend()
    if not gpu.available()["cuml"] or "cuml" not in backend["umap"]:
        raise SystemExit("cuML GPU UMAP/HDBSCAN is required; refusing a silent CPU fallback")
    nodes = pd.read_parquet(paths.cache_file("nodes.parquet"))
    output = Path(args.output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    experiments = (("compartment", "biology"), ("compartment", "target_family"),
                   ("cellcycle_phase", "biology"), ("cellcycle_phase", "target_family"))
    manifests = []
    for target, scope in experiments:
        manifests.append(run_experiment(nodes, target, scope, output, quick=args.quick,
                                        permutations=args.permutations, log=lambda x: print(x, flush=True)))
    summary = pd.DataFrame([{key: value for key, value in manifest.items()
                             if not isinstance(value, (list, dict))} for manifest in manifests])
    summary.to_csv(output / "summary.csv", index=False)
    (output / "manifest.json").write_text(json.dumps({"backend": backend,
                                                       "experiments": manifests}, indent=2) + "\n")
    print(summary.to_string(index=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
