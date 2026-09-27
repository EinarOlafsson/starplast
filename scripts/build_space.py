#!/usr/bin/env python3
"""Validate a prepared organism space and package it without losing installed measurements."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


def main(argv=None):
    """Build a versioned data pack from tables and an explicit per-column license mapping."""
    from starplast import organisms, packs
    from starplast.spaces import build_pack
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--organism", required=True, choices=organisms.codes())
    parser.add_argument("--data", required=True, type=Path)
    parser.add_argument("--version", required=True)
    parser.add_argument("--licenses", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--previous", type=Path, help="Earlier node table; defaults to the installed table")
    args = parser.parse_args(argv)
    space = organisms.get(args.organism)
    previous = args.previous
    if previous is None:
        installed = packs.installed(space.code)
        try:
            candidate = installed / space.nodes if installed else Path(organisms.nodes_path(space.code))
        except FileNotFoundError:
            candidate = None
        if candidate is not None and candidate.is_file():
            previous = candidate
    checksum = build_pack(space, args.data, args.version, json.loads(args.licenses.read_text()),
                          args.output, previous)
    print(f"{checksum}  {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
