#!/usr/bin/env python3
"""In-vivo fitness of hyperLOPIT-unassigned proteins

Two targeted libraries tested during mouse infection

    level / kind : DNA / CRISPR_screen
    provides     : fit_hyperlopit_unassigned_invivo_lib1, fit_hyperlopit_unassigned_invivo_lib2
    coverage     : measured at build time
    citation     : Tachibana Y et al., CRISPR screens identify genes essential for in vivo virulence among proteins of hyperLOPIT-unassigned localization. mBio 2024
    PMID         : 39082802
    accession    : GSE253884;GSE253885
    url          : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253884
    local path   : datasets/toxoplasma_acquisition_2026_08_14/

Quirks that cost time once:
    Only the newly measured in-vivo fitness values enter; copied comparator columns in the
    summary workbook are not duplicated.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gse253884_5.py
"""
from _common import run

KEY = "gse253884_5"

if __name__ == "__main__":
    run(KEY)
