#!/usr/bin/env python3
"""OrthoMCL orthogroups

Orthogroup assignment and cross-species bridge

    level / kind : reference / orthology
    provides     : orthogroup
    coverage     : 16,793 groups
    accession    : OrthoMCL
    local path   : datasets/MASTER_parasite_wide_by_orthogroup.csv

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
