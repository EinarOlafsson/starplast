#!/usr/bin/env python3
"""Differentiation reporter CRISPR screen (COMPUTED)

Guide enrichment in reporter-positive parasites against the bulk population

    level / kind : DNA / CRISPR_screen
    provides     : diff_reporter_log2_mNG_over_bulk
    coverage     : 235 genes
    accession    : GSE132237
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132237/suppl/GSE132237_RAW.tar
    local path   : datasets/quarantine/2026_08_16_pride/Tg/essentiality_in_a_second_background/GSE132237_RAW.tar

Quirks that cost time once:
    COMPUTED here: the deposit publishes guide counts, not the ratio. Guides are summed per
    gene, scaled to a common library size, and log2((mNG+ + 1)/(bulk + 1)) is averaged over the
    two reporter lines at 10 days. It is the comparison the authors' design names, and the
    column says ratio rather than phenotype so nobody mistakes it for a number they reported. A
    targeted screen against nucleic-acid binding proteins, so 240 genes is its full extent and
    not a coverage failure. Which member is which sample comes from the series matrix, never
    from the file name.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/differentiation_screen.py
"""
from _common import run

KEY = "differentiation_screen"

if __name__ == "__main__":
    run(KEY)
