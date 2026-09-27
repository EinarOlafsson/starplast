#!/usr/bin/env python3
"""Split the legacy host cache without losing identifiers or measurements.

The source is kept intact. Run with --write to create the separate Hs/Mm caches and a provenance
record; otherwise only validate and report. Identity-only rows are resolved against the local
species-specific UniProt mapping files and the existing human interaction bridge sources.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from starplast import datasets, host, organisms, paths  # noqa: E402


def identity_memberships(frame, data_dir, dataset_root):
    """Resolve identity-only rows from local sources, recording their exact file hashes."""
    columns = [c for c in frame if c not in ("host_id", "host_name")]
    wanted = set(frame.loc[~frame[columns].notna().any(axis=1), "host_id"])
    memberships = {code: set() for code in host.HOST_TABLES}
    records = []
    sources = [(organisms.HUMAN, Path(data_dir) / "host_bridges.parquet"),
               (organisms.HUMAN, Path(data_dir) / "pf_host_bridges.parquet")]
    for species, code in (("human", organisms.HUMAN), ("mouse", organisms.MOUSE)):
        sources.append((code, Path(dataset_root) / host.UNIPROT_ROOT
                        / host.UNIPROT_IDMAP[species][0]))
    for code, path in sources:
        if not path.exists():
            continue
        if path.suffix == ".parquet":
            found = wanted & set(pd.read_parquet(path, columns=["host_id"]).host_id)
        else:
            found = set()
            with gzip.open(path, "rt") as fh:
                for line in fh:
                    accession = line.split("\t", 1)[0].split("-")[0]
                    if accession in wanted:
                        found.add(accession)
        memberships[code].update(found)
        records.append({"source": path.name, "organism": code,
                        "sha256": datasets.digest(str(path)),
                        "identity_only_rows": sorted(found)})
    return memberships, records


def migrate(source, data_dir, dataset_root, write=False):
    """Validate a lossless split and optionally write it; return its provenance report."""
    source, data_dir = Path(source), Path(data_dir)
    original = pd.read_parquet(source)
    membership, records = identity_memberships(original, data_dir, dataset_root)
    tables = host.split_legacy_table(original, membership)
    before = original.set_index("host_id")
    all_ids = set()
    for code, table in tables.items():
        if all_ids & set(table.host_id):
            raise ValueError("the split assigned a host identifier to two organisms")
        all_ids.update(table.host_id)
        after = table.set_index("host_id")
        pd.testing.assert_frame_equal(before.loc[after.index, after.columns], after)
    if all_ids != set(original.host_id):
        raise ValueError("the split lost host identifiers")
    for column in original:
        if column not in ("host_id", "host_name"):
            count = sum(int(t[column].notna().sum()) for t in tables.values() if column in t)
            if count != int(original[column].notna().sum()):
                raise ValueError(f"the split lost measurements in {column}")
    report = {"source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
              "rows_before": len(original), "identity_sources": records,
              "tables": {code: {"file": host.HOST_TABLES[code], "rows": len(t),
                                "coverage": {c: int(t[c].notna().sum()) for c in t}}
                         for code, t in tables.items()}}
    if write:
        # Check every existing output first: a rerun must never replace later acquisitions.
        for code, table in tables.items():
            output = data_dir / host.HOST_TABLES[code]
            if output.exists():
                pd.testing.assert_frame_equal(pd.read_parquet(output), table)
        for code, table in tables.items():
            output = data_dir / host.HOST_TABLES[code]
            temporary = output.with_suffix(".tmp.parquet")
            table.to_parquet(temporary, index=False)
            pd.testing.assert_frame_equal(table, pd.read_parquet(temporary))
            os.replace(temporary, output)
        (data_dir / "host_species_migration.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    """Validate the old host cache, optionally writing the two species tables."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "starplast/data")
    parser.add_argument("--source", type=Path)
    parser.add_argument("--dataset-root", default=paths.dataset_root())
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    report = migrate(args.source or args.data / "host_proteins.parquet", args.data,
                     args.dataset_root, write=args.write)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
