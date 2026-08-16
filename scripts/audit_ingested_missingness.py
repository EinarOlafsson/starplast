#!/usr/bin/env python3
"""Measure whether missing coverage separates each newly exposed assay block.

The output is a reproducible audit table, not a pass/fail threshold: missingness often tracks real
abundance biology, but its map-scale centroid gap has to travel with any result from a partial assay.
Run after rebuilding the cache::

    python scripts/audit_ingested_missingness.py --out results/data_ingestion_2026_08_14.csv
"""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import embedding, paths  # noqa: E402


SLOTS = (
    "Tg_transcription_in_vivo_brain_acute",
    "Tg_transcription_in_vivo_brain_chronic",
    "Tg_transcription_purified_bradyzoite_in_vivo",
    "Tg_transcription_under_stress_conversion",
    "Tg_transcription_under_tf_or_chromatin_perturbation",
    "Tg_protein_abundance_tachyzoite",
    "Tg_protein_abundance_other_life_stages",
    "Tg_phosphorylation_quantitative",
    "Tg_host_transcriptional_effect_per_effector",
)


def audit(nodes: pd.DataFrame, method: str = "umap", log=print) -> pd.DataFrame:
    """Return one missingness row per source column for every newly exposed slot."""
    rows = []
    for key in SLOTS:
        spec = embedding.EmbeddingSpec(name=key, blocks=(key,), method=method,
                                       n_components=3, random_state=42)
        coords, _names, kept = embedding.embed(nodes, spec, log=log)
        source = embedding.columns_for(nodes, spec).get(key, [])
        subset = nodes.iloc[kept]
        report = embedding.missingness_leak(coords, subset, source)
        checked = set(report.column) if not report.empty else set()
        skipped = []
        for column in source:
            if column not in checked:
                missing = int(subset[column].isna().sum())
                skipped.append({"column": column, "missing_frac": missing / len(subset),
                                "centroid_gap": float("nan"),
                                "status": "too few measured or missing genes (<30)"})
        if not report.empty:
            report["status"] = "measured"
        report = pd.concat([report, pd.DataFrame(skipped)], ignore_index=True)
        if not report.empty:
            report.insert(0, "slot", key)
            rows.append(report)
    columns = ["slot", "column", "missing_frac", "centroid_gap", "status"]
    return pd.concat(rows, ignore_index=True)[columns] if rows else pd.DataFrame(columns=columns)


def main(argv=None) -> int:
    """Run the audit against the shipped cache and save a CSV."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="results/data_ingestion_2026_08_14.csv")
    parser.add_argument("--method", choices=embedding.METHODS, default="umap")
    args = parser.parse_args(argv)
    nodes = pd.read_parquet(paths.cache_file("nodes.parquet"))
    result = audit(nodes, method=args.method)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    result.to_csv(args.out, index=False)
    print(f"{len(result)} missingness checks -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
