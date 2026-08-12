#!/usr/bin/env python3
"""IP-MS of tagged baits

Replicated pulldown vs untagged control

    level / kind : post_translation / IPMS
    provides     : n_ipms_partners
    coverage     : 64 pairs / 48 genes
    accession    : PXD043808, PXD065585

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/ipms_baits.py
"""
from _common import run

KEY = "ipms_baits"

if __name__ == "__main__":
    run(KEY)
