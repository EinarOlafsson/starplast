#!/usr/bin/env python3
"""Discover and archive pan-Apicomplexan source files. Instruction 38.

Discovery is resolved rather than remembered: the site-to-clade table is the only thing written
down, and the release, the organisms and the files are read from each site's own index at run time.

    python scripts/build_archive.py --discover-only        # enumerate, write the catalogue, fetch nothing
    python scripts/build_archive.py --limit 20             # fetch the first 20 sources
    python scripts/build_archive.py --organisms 3          # a small walk, for checking the pipeline

Nothing is overwritten: a file already on disk is recorded as `present` and left alone. Failures are
recorded with their reason rather than omitted, so the next session does not rediscover the same
dead ends.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from dataclasses import asdict

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ARCHIVE = "/mnt/firecuda2/Claude/toxoplasma_projects/datasets/pan_apicomplexan"
CATALOG = os.path.join(_ROOT, "instructions", "done", "38_archive_catalog.json")
MANIFEST = os.path.join(_ROOT, "instructions", "done", "38_archive_manifest.csv")


def main() -> int:
    from starplast import archive
    ap = argparse.ArgumentParser()
    ap.add_argument("--discover-only", action="store_true")
    ap.add_argument("--limit", type=int, default=0, help="fetch at most this many sources")
    ap.add_argument("--organisms", type=int, default=0, help="organisms per site, 0 for all")
    ap.add_argument("--root", default=ARCHIVE)
    ap.add_argument("--catalog", default=CATALOG)
    ap.add_argument("--manifest", default=MANIFEST)
    args = ap.parse_args()

    sources = archive.discover_veupathdb(limit_organisms=args.organisms)
    os.makedirs(os.path.dirname(args.catalog), exist_ok=True)
    with open(args.catalog, "w") as fh:
        json.dump([asdict(s) for s in sources], fh, indent=1)
    print(f"\n{len(sources)} sources -> {args.catalog}")
    by_clade = {}
    for s in sources:
        by_clade[s.scope] = by_clade.get(s.scope, 0) + 1
    for clade, n in sorted(by_clade.items()):
        print(f"  {n:5}  {clade}")
    if args.discover_only:
        return 0
    if args.limit:
        sources = sources[:args.limit]
    print(f"\nfetching {len(sources)} into {args.root}")
    archive.build(sources, args.root, args.manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
