#!/usr/bin/env python3
"""Add the verified mass-spectrometry columns to the cached node table.

`build_graph` is the proper home for this and it needs the upstream `toxonet` interim tables, which
are not on every machine that has the cache. So this does the one merge on its own: read the node
table, ask `proteomics.load_all` for its columns, write them back.

Idempotent by construction -- a column is recomputed from the deposit each time rather than appended
-- so running it twice is running it once. It prints the coverage of every column it writes, because
a merge that silently wrote a column of NaN would look exactly like a merge that worked.
"""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import proteomics  # noqa: E402


def main(argv=None) -> int:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--nodes", default=os.path.join(root, "starplast", "data", "nodes.parquet"))
    p.add_argument("--base", default=root)
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args(argv)

    nodes = pd.read_parquet(args.nodes)
    before = nodes.shape[1]
    columns = proteomics.load_all(args.base, nodes.gene_id.astype(str), log=print)
    if columns.empty:
        print("no verified deposits found; nothing to add")
        return 1
    for column in columns.columns:
        nodes[column] = columns[column].to_numpy()
    got = {c: int(nodes[c].notna().sum()) for c in columns.columns}
    print(f"{before} columns -> {nodes.shape[1]}")
    for column, n in sorted(got.items()):
        print(f"  {column:<28} {n:>6,} genes measured ({100 * n / len(nodes):.0f}%)")
    if min(got.values()) == 0:
        # A column of nothing is worse than a missing column: the slot flips to filled and the map
        # gains a feature that is entirely absent, which reads as a measurement nobody made.
        print("REFUSED: a column came back empty", file=sys.stderr)
        return 2
    if args.dry_run:
        print("dry run, not written")
        return 0
    nodes.to_parquet(args.nodes, index=False)
    print(f"written to {args.nodes}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
