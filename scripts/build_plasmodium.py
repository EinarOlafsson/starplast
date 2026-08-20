#!/usr/bin/env python3
"""Rebuild the Plasmodium cache, and refuse to write a build that lost something.

`plasmodium.build_all` exists so the table is reproducible rather than the product of whatever was
typed at a prompt -- and until now the INVOCATION was exactly that, typed fresh each time, with the
comparison against the previous cache left to whoever remembered to make it. The comparison is the
part that matters: a build that exits 0 is not evidence that a build was correct. Two silent losses
have already been caught only by making it, and both were in columns that still existed:

* four genes of a screen, lost when a strain accession stopped resolving;
* 1,306 values of `best_model_agreement`, wiped while its column stayed in the table.

So this script does the whole dangerous step in one place: build, diff against the shipped cache on
columns GAINED, columns LOST and per-column coverage in BOTH directions, and write only if nothing
went backwards. `--allow-loss WHY` is the override, and it requires a sentence, because a loss that
someone decided to accept and a loss nobody noticed must not look the same in the history.

    python scripts/build_plasmodium.py                 # build, diff, write if nothing was lost
    python scripts/build_plasmodium.py --dry-run       # build and diff, write nothing
    python scripts/build_plasmodium.py --allow-loss "the deposit was withdrawn upstream"
"""
from __future__ import annotations

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import host, paths, pf_graph, plasmodium  # noqa: E402

#: A coverage change smaller than this in either direction is not worth a line of output. It is ONE
#: gene: the 4-gene regression that motivated this check would have printed.
NOTABLE = 1


def diff(old: pd.DataFrame, new: pd.DataFrame, key: str = "gene_id") -> dict:
    """What changed between two builds, in both directions.

    `key` is the column that identifies a row, because this is now run over the host table too and
    its rows are proteins. Writing a host table with the gene key would compare two empty sets and
    report every dropped protein as nothing at all.
    """
    old_cols, new_cols = set(old.columns), set(new.columns)
    shared = sorted(old_cols & new_cols)
    coverage = {c: (int(old[c].notna().sum()), int(new[c].notna().sum())) for c in shared}
    return {"gained": sorted(new_cols - old_cols),
            "lost": sorted(old_cols - new_cols),
            "rows": (len(old), len(new)),
            "key": key,
            "genes_lost": sorted(set(old.get(key, [])) - set(new.get(key, []))),
            "genes_gained": sorted(set(new.get(key, [])) - set(old.get(key, []))),
            "fell": {c: v for c, v in coverage.items() if v[1] - v[0] <= -NOTABLE},
            "rose": {c: v for c, v in coverage.items() if v[1] - v[0] >= NOTABLE}}


def report(d: dict, log=print) -> None:
    """Print the diff. Coverage going UP is printed too -- it is how two identity bugs surfaced."""
    log(f"rows: {d['rows'][0]:,} -> {d['rows'][1]:,}")
    log(f"columns: +{len(d['gained'])} / -{len(d['lost'])}")
    for c in d["gained"]:
        log(f"  + {c}")
    for c in d["lost"]:
        log(f"  - {c}   <-- LOST")
    noun = "genes" if d.get("key", "gene_id") == "gene_id" else "rows"
    for c, (before, after) in sorted(d["fell"].items()):
        log(f"  ! {c}: {before:,} -> {after:,} {noun}   <-- coverage FELL")
    for c, (before, after) in sorted(d["rose"].items()):
        log(f"  ^ {c}: {before:,} -> {after:,} {noun}")
    if d["genes_lost"]:
        log(f"  {noun} dropped: {len(d['genes_lost'])} ({', '.join(d['genes_lost'][:5])} ...)")
    if d["genes_gained"]:
        log(f"  {noun} added: {len(d['genes_gained'])}")


def regressions(d: dict) -> list:
    """Everything that went backwards, as sentences. Empty means the build may be written."""
    out = [f"column lost: {c}" for c in d["lost"]]
    out += [f"coverage fell in {c}: {before:,} -> {after:,}"
            for c, (before, after) in sorted(d["fell"].items())]
    if d["genes_lost"]:
        noun = "genes" if d.get("key", "gene_id") == "gene_id" else "rows"
        out.append(f"{len(d['genes_lost'])} {noun} dropped from the table")
    return out


def main(argv=None, log=print) -> int:
    """Build, diff, and write unless something went backwards."""
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dataset-root", default=None, help="defaults to the resolved dataset root")
    ap.add_argument("--dry-run", action="store_true", help="build and diff, write nothing")
    ap.add_argument("--allow-loss", default=None, metavar="WHY",
                    help="write even though something was lost, and say why")
    args = ap.parse_args(argv)

    root = args.dataset_root or paths.dataset_root()
    table = paths.cache_file(plasmodium.TABLE)
    graph = paths.cache_file(pf_graph.GRAPH)
    # The host table first, because a host slot is graded against it. It is not a Plasmodium file
    # -- a mouse macrophage surfaceome answers a Toxoplasma slot -- but this is the one build that
    # already assembles it, and splitting it out would give two scripts one output. Merged rather
    # than overwritten: `build_graph` writes the Toxoplasma pulldown's columns into the same file,
    # and the owners share a key rather than a column.
    tissue = host.tissue_references(root, log=log)
    if not tissue.empty:
        host_path = paths.cache_file("host_proteins.parquet")
        existing = pd.read_parquet(host_path) if os.path.exists(host_path) else pd.DataFrame()
        merged = host.merge_tissue(existing, tissue)
        # Diffed like the gene table, and for the same reason: this merge silently dropped 311 host
        # symbols the first time it ran, and nothing said so because nothing was looking.
        if len(existing):
            host_diff = diff(existing, merged, key="host_id")
            log("host table:")
            report(host_diff, log=log)
            problems = regressions(host_diff)
            if problems and not args.allow_loss:
                log("REFUSING to write the host table. " + "; ".join(problems))
                return 2
        if not args.dry_run:
            merged.to_parquet(host_path, index=False)
            log(f"wrote {host_path}  ({len(merged):,} host proteins, "
                f"{len(merged.columns)} columns)")

    new = plasmodium.build_all(root, log=log)
    if new.empty:
        log("build produced nothing -- is the PlasmoDB gene report under reference/plasmodb?")
        return 1

    if os.path.exists(table):
        d = diff(pd.read_parquet(table), new)
        report(d, log=log)
        problems = regressions(d)
        if problems and not args.allow_loss:
            log("REFUSING to write. " + "; ".join(problems))
            log("Re-run with --allow-loss \"the reason\" if this is deliberate.")
            return 2
        if problems:
            log(f"writing anyway, on the stated ground: {args.allow_loss}")
    else:
        log(f"no previous cache at {table} -- nothing to diff against")

    if args.dry_run:
        log("--dry-run: not written")
        return 0
    new.to_parquet(table, index=False)
    log(f"wrote {table}  ({len(new):,} genes, {len(new.columns)} columns)")
    pf_graph.save(new, graph, log=log, dataset_root=root)
    log(f"wrote {graph}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
