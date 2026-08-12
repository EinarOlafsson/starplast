#!/usr/bin/env python3
"""OrthoMCL orthogroups

Orthogroup assignment and cross-species bridge

    level / kind : reference / orthology
    provides     : orthogroup
    coverage     : 16,793 groups
    accession    : OrthoMCL release 6.21
    url          : https://orthomcl.org/common/downloads/release-6.21/groups_OrthoMCL-6.21.txt.gz
    local path   : datasets/MASTER_parasite_wide_by_orthogroup.csv

Quirks that cost time once:
    The RELEASE is load-bearing, not decoration: all 16,793 groups reproduce exactly from 6.21,
    while 6.20 differs in 180 cells and Current_Release (v7) renumbers every group into an OG7_
    namespace that matches nothing here. The shipped CSV is derived from this file -- filter to
    tgon/pfal/cpar/tbrt, then pivot wide. Note the T. brucei taxon code is tbrt, though the CSV
    column is tbru.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/orthomcl.py
"""
from _common import run

KEY = "orthomcl"

if __name__ == "__main__":
    run(KEY)
