#!/usr/bin/env python3
"""One question, several methods, one control judging all of them.

Instruction 47's deliverable. Every method answers the SAME recipe -- same leakage closure, same
seed, same holdout, same validation control -- so the numbers below differ because the methods do,
not because each was written up on its own terms.

**Read the F1 column, not the gene count.** Both are honest and only one is comparable. A
classifier's partition is its own predicted classes, so a label is concentrated within its group by
construction and the naming step is far more permissive there than it is over a clustering's 60-odd
groups. Mean F1 is measured the same way for every method -- out of fold for the supervised ones --
and is the column that says which method actually recovers the label.

    python scripts/compare_methods.py [--limit N] [--methods umap+hdbscan,logistic]
"""
from __future__ import annotations

import argparse
import os
import sys
from dataclasses import replace

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(_ROOT, "instructions", "done", "47_method_comparison.csv")
NODES = os.path.join(_ROOT, "starplast", "data", "nodes.parquet")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--methods", default="umap+hdbscan,logistic")
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()

    from starplast.questions import as_recipe, load, shipped
    from starplast.recipes import run

    nodes = pd.read_parquet(NODES)
    catalogue = shipped(nodes, load())
    if args.limit:
        catalogue = catalogue[:args.limit]
    wanted = [m.strip() for m in args.methods.split(",") if m.strip()]
    rows = []
    for i, q in enumerate(catalogue, 1):
        base = as_recipe(q)
        for method in wanted:
            result = run(nodes, replace(base, method=method),
                         tune=(method == "umap+hdbscan"), log=lambda *a: None)
            corroborated = 0
            if len(result.inference) and "control_agrees" in result.inference:
                corroborated = int(result.inference["control_agrees"].fillna(False).sum())
            rows.append({
                "question": base.question, "axis": base.axis, "method": method,
                "ran": result.ok, "stopped_because": result.stopped_because,
                "labels_scored": result.summary.get("n_labels_scored", 0),
                "mean_f1": result.summary.get("mean_f1", float("nan")),
                "best_f1": result.summary.get("best_f1", float("nan")),
                "best_label": result.summary.get("best_label", ""),
                "genes_named": len(result.inference), "corroborated": corroborated,
                "groups": result.quality.get("clusters", 0),
            })
            pd.DataFrame(rows).to_csv(args.out, index=False)
            last = rows[-1]
            print(f"[{i}/{len(catalogue)}] {method:14} mean F1 {last['mean_f1']:.3f} "
                  f"genes {last['genes_named']:5} corroborated {last['corroborated']:5}", flush=True)

    table = pd.DataFrame(rows)
    print(f"\n{len(table)} runs -> {args.out}\n")
    if len(table):
        summary = (table[table.ran].groupby("method")
                   .agg(questions=("question", "nunique"), mean_f1=("mean_f1", "mean"),
                        genes=("genes_named", "sum"), corroborated=("corroborated", "sum"))
                   .reset_index())
        summary["corroboration_rate"] = summary.corroborated / summary.genes.replace(0, pd.NA)
        print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
