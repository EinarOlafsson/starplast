#!/usr/bin/env python3
"""IEDB epitopes mapped to genes (via ToxoDB)

How many IEDB epitopes ToxoDB maps to this gene

    level / kind : reference / immunity
    provides     : iedb_epitope_count
    coverage     : 221 genes
    accession    : ToxoDB / IEDB
    url          : https://toxodb.org/toxo/service/record-types/transcript/searches/GenesWithEpitopes/reports/attributesTabular
    local path   : starplast/data/toxodb_epitopes.tsv

Quirks that cost time once:
    ToxoDB's own join of IEDB against the ME49 proteome, at all three confidence levels. NOT
    split by epitope type -- the integration does not expose that -- so the count is T-cell and
    B-cell epitopes alike, which is why the column is named iedb_ and not t_cell_. Verified by
    what comes out on top: SRS29B (SAG1) with 45, then GRA6, GRA7, GRA2 and ROP18. Those are the
    canonical Toxoplasma serology antigens, in the order a serologist would put them.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/toxodb_epitopes.py
"""
from _common import run

KEY = "toxodb_epitopes"

if __name__ == "__main__":
    run(KEY)
