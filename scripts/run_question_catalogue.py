#!/usr/bin/env python3
"""Run every shipped question end to end, and write what each one actually answered.

Instruction 45's acceptance is not "the twenty exist" but "the twenty run end to end and return
genes", so this runs them and records the outcome per question: the map that was built, what
recovered, how many genes were named, and whether the independent control corroborated any of them.

It is slow on purpose -- each question tunes its own UMAP and its own clustering, because a recipe
that reused another question's map would be answering with a map chosen for a different question.
Expect minutes per question. Results are written after EACH question rather than at the end, so a run
that is interrupted keeps what it has.

    python scripts/run_question_catalogue.py [--limit N] [--no-tune]
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(_ROOT, "instructions", "done", "45_question_results.csv")
NODES = os.path.join(_ROOT, "starplast", "data", "nodes.parquet")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--no-tune", action="store_true")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--pdf-dir", default="", help="write a paper-ready PDF per question here")
    ap.add_argument("--alternatives", type=int, default=0,
                    help="extra maps to build and show in the PDF, so the winner has a comparison")
    args = ap.parse_args()

    from starplast.questions import as_recipe, load, shipped
    from starplast.recipes import run

    nodes = pd.read_parquet(NODES)
    catalogue = shipped(nodes, load())
    if args.limit:
        catalogue = catalogue[:args.limit]
    print(f"{len(catalogue)} questions to run\n")

    rows = []
    for i, q in enumerate(catalogue, 1):
        recipe = as_recipe(q)
        t0 = time.time()
        print(f"[{i}/{len(catalogue)}] {recipe.question}")
        result = run(nodes, recipe, tune=not args.no_tune,
                     alternatives=args.alternatives if args.pdf_dir else 0, log=lambda *a: None)
        if args.pdf_dir:
            from starplast.report import recipe_pdf
            safe = "".join(c if c.isalnum() else "_" for c in recipe.question)[:60].strip("_")
            recipe_pdf(result, os.path.join(args.pdf_dir, f"{i:02d}_{safe}.pdf"))
        agreed = 0
        if len(result.inference) and "control_agrees" in result.inference:
            agreed = int(result.inference["control_agrees"].fillna(False).sum())
        rows.append({
            "question": recipe.question, "axis": recipe.axis,
            "expect_refusal": recipe.expect_refusal,
            "ran": result.ok, "stopped_because": result.stopped_because,
            "clusters": result.quality.get("clusters", 0),
            "settings": json.dumps(result.settings),
            "labels_scored": result.summary.get("n_labels_scored", 0),
            "mean_f1": result.summary.get("mean_f1", float("nan")),
            "best_f1": result.summary.get("best_f1", float("nan")),
            "best_label": result.summary.get("best_label", ""),
            "genes_named": len(result.inference),
            "genes_corroborated": agreed,
            "note": result.inference_note,
            "seconds": round(time.time() - t0, 1),
        })
        pd.DataFrame(rows).to_csv(args.out, index=False)
        last = rows[-1]
        print(f"    {'ran' if last['ran'] else 'STOPPED: ' + last['stopped_because']}"
              f" | {last['genes_named']} genes, {last['genes_corroborated']} corroborated"
              f" | {last['seconds']}s\n")

    table = pd.DataFrame(rows)
    answered = int((table.genes_named > 0).sum())
    print(f"\n{len(table)} questions -> {args.out}")
    print(f"  ran to completion : {int(table.ran.sum())}")
    print(f"  returned genes    : {answered}")
    print(f"  with corroboration: {int((table.genes_corroborated > 0).sum())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
