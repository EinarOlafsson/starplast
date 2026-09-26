#!/usr/bin/env python3
"""Merge the derived deposit tables into the shipped caches.

`starplast.deposits` derives each deposit into `starplast/data/deposit_<key>.tsv`; this puts those
columns into `nodes.parquet` (Toxoplasma), `pf_nodes.parquet` (Plasmodium) and
`host_proteins.parquet` (host), the same columns a full `build_graph` would produce. Idempotent: a
column is replaced from its table each time, never appended twice.

It refuses to write when a merge would LOSE something -- a column that had values and now has
fewer, or a host protein that disappears -- because a merge that silently drops data exits 0 and
looks like one that worked.

Run:  python scripts/derive_deposits.py     # raw deposits -> deposit tables (+ the notebook)
      python scripts/add_deposits.py         # deposit tables -> shipped caches
"""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from starplast import deposits, host, identity  # noqa: E402


def toxo_resolver(data: str, ids):
    """GT1 and old accessions to current ME49, through the identity layer."""
    ix = identity.build_index(ids, os.path.join(data, "toxodb_identity.tsv"),
                              log=lambda *a: None)
    identity.add_strain_accessions(
        ix, {"GT1": os.path.join(data, "toxodb_strain_gt1.tsv"),
             "VEG": os.path.join(data, "toxodb_strain_veg.tsv")}, log=lambda *a: None)

    def resolve(acc):
        hit = ix.lookup.get(identity.norm(acc))
        return hit[0] if hit else None
    return resolve


def merge_parasite(path: str, organism: str, base: str, resolve=None, log=print,
                   dry_run: bool = False) -> int:
    nodes = pd.read_parquet(path)
    ids = nodes["gene_id"].astype(str)
    cols = deposits.parasite_columns(base, organism, ids, resolve=resolve, log=log)
    if cols.empty or not len(cols.columns):
        log(f"{organism}: no deposit columns")
        return 0
    lost = [c for c in cols.columns if c in nodes.columns
            and cols[c].notna().sum() < nodes[c].notna().sum()]
    if lost:
        log(f"REFUSED ({organism}): would lose values in {lost}")
        return 2
    empty = [c for c in cols.columns if cols[c].notna().sum() == 0]
    if empty:
        log(f"REFUSED ({organism}): empty columns {empty} -- is the identity layer resolving?")
        return 2
    before = nodes.shape[1]
    for c in cols.columns:
        nodes[c] = cols[c].to_numpy()
    for c in cols.columns:
        n = int(nodes[c].notna().sum())
        log(f"  {c:<40} {n:>6,} genes ({100 * n / len(nodes):.0f}%)")
    log(f"{organism}: {before} columns -> {nodes.shape[1]}")
    if not dry_run:
        nodes.to_parquet(path, index=False)
    return 0


def merge_host(path: str, base: str, log=print, dry_run: bool = False) -> int:
    new = deposits.host_columns(base)
    if new.empty:
        log("host: no deposit columns")
        return 0
    existing = pd.read_parquet(path) if os.path.exists(path) else pd.DataFrame()
    merged = host.merge_tissue(existing, new)
    if len(existing):
        gone = set(existing["host_id"]) - set(merged["host_id"])
        named_before = int(existing["host_name"].notna().sum())
        named_after = int(merged.set_index("host_id").loc[list(existing["host_id"]),
                                                          "host_name"].notna().sum())
        if gone or named_after < named_before:
            log(f"REFUSED (host): {len(gone)} proteins lost, names {named_before} -> "
                f"{named_after}")
            return 2
    for c in new.columns:
        if c not in ("host_id", "host_name"):
            log(f"  {c:<40} {int(merged[c].notna().sum()):>6,} host proteins")
    log(f"host: {len(existing):,} -> {len(merged):,} proteins, {merged.shape[1]} columns")
    if not dry_run:
        merged.to_parquet(path, index=False)
    return 0


def main(argv=None, log=print) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--base", default=ROOT)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args(argv)
    data = os.path.join(args.base, "starplast", "data")
    tg = os.path.join(data, "nodes.parquet")
    status = merge_parasite(tg, "Tg", args.base,
                            resolve=toxo_resolver(data, pd.read_parquet(tg, columns=["gene_id"])
                                                  ["gene_id"].astype(str)),
                            log=log, dry_run=args.dry_run)
    pf = os.path.join(data, "pf_nodes.parquet")
    if os.path.exists(pf) and any(d.organism == "Pf" for d in deposits.DEPOSITS):
        status = max(status, merge_parasite(pf, "Pf", args.base, log=log, dry_run=args.dry_run))
    status = max(status, merge_host(os.path.join(data, "host_proteins.parquet"), args.base,
                                    log=log, dry_run=args.dry_run))
    return status


if __name__ == "__main__":
    raise SystemExit(main())
