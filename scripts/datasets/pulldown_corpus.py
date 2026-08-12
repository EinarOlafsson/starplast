#!/usr/bin/env python3
"""Pulldown corpus

55 IP-MS / co-IP studies with a tagged Toxoplasma protein

    level / kind : post_translation / IPMS
    provides     : (edges or build inputs only)
    coverage     : 29 studies with data, 140 files
    local path   : datasets/post_translation/IPMS/

Quirks that cost time once:
    Downloaded and indexed; NOT yet parsed into edges.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pulldown_corpus.py
"""
from _common import run

KEY = "pulldown_corpus"

if __name__ == "__main__":
    run(KEY)
