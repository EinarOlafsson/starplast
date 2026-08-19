#!/usr/bin/env python3
"""Plasmodium host interaction degree (COMPUTED)

How many host proteins a gene was crosslinked to, where it was looked at

    level / kind : reference / crosslink_MS
    provides     : n_host_targets
    coverage     : 117 genes seen, 10 with a host partner
    local path   : starplast/data/pf_nodes.parquet
    derived from : bridge:host

Quirks that cost time once:
    THREE states, and the middle one is why this is not a fillna(0). The Toxoplasma column of
    the same name is 0 everywhere without a curated host target, which is right there because
    its source is a curated table covering the literature. This source is ONE experiment: a gene
    it never detected has not been shown to lack host partners. So the 117 genes seen in the
    crosslink data carry a count -- zero included, because being crosslinked only to parasite
    proteins is a real observation -- and the other 5,603 stay missing. Classification is three-
    way for a reason a test found: a PF3D7 accession the node table does not carry is a PARASITE
    protein with no row, and reading it as host inflated this count. The shipped numbers were
    unaffected, because every accession in this file is in the table, but the fix is what stops
    the next file from being wrong.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_host_degree.py
"""
from _common import run

KEY = "pf_host_degree"

if __name__ == "__main__":
    run(KEY)
