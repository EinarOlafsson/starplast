#!/usr/bin/env python3
"""Add newly verified columns to the cached node table.

`build_graph` is the proper home for these and it needs the upstream `toxonet` interim tables, which
are not on every machine that has the cache. So this does the merge on its own: read the node table,
ask each source for its columns, write them back. Every source here is also in `build_graph`'s own
loader list, so a full rebuild produces the same table -- this is a shortcut, not a second pipeline.

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

from starplast import (chromatin, codons, expression, identity, iedb,  # noqa: E402
                       palmitome, proteomics, toxodb_evidence, variation)


def main(argv=None) -> int:
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--nodes", default=os.path.join(root, "starplast", "data", "nodes.parquet"))
    p.add_argument("--base", default=root)
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args(argv)

    nodes = pd.read_parquet(args.nodes)
    before = nodes.shape[1]
    ids = nodes.gene_id.astype(str)

    # Accessions in a deposit are whatever was current when it was submitted, so they go through the
    # identity layer rather than being matched as strings -- the same routing build_graph uses, and
    # the reason a table of TGGT1_ ids joins anything at all.
    data = os.path.join(args.base, "starplast", "data")
    ix = identity.build_index(ids, os.path.join(data, "toxodb_identity.tsv"),
                              log=lambda *a: None)
    identity.add_strain_accessions(
        ix, {"GT1": os.path.join(data, "toxodb_strain_gt1.tsv"),
             "VEG": os.path.join(data, "toxodb_strain_veg.tsv")}, log=lambda *a: None)

    def resolve(acc):
        hit = ix.lookup.get(identity.norm(acc))
        return hit[0] if hit else None

    columns = proteomics.load_all(args.base, ids, log=print, resolve=resolve)

    for table in (expression.gse245775_differentiation_ribosome_profiling(
                      args.base, resolve=resolve, log=print),
                  expression.gse223620_bfd2_rip(args.base, resolve=resolve, log=print),
                  variation.strain_snps(args.base, resolve=resolve, log=print),
                  codons.codon_usage(args.base, resolve=resolve, log=print),
                  chromatin.chromatin_signals(args.base, resolve=resolve, log=print),
                  palmitome.palmitome(args.base, resolve=resolve, log=print),
                  toxodb_evidence.evidence(args.base, resolve=resolve, log=print),
                  toxodb_evidence.enzyme_classification(args.base, resolve=resolve,
                                                        log=print),
                  iedb.bcell_epitopes(args.base, resolve=resolve, log=print)):
        if table.empty:
            continue
        aligned = table.reindex(pd.Index(ids))
        for column in table.columns:
            columns[column] = aligned[column].to_numpy()

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
