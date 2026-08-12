#!/usr/bin/env python3
"""Oocyst sporulation series

Unsporulated / sporulating / sporulated, 2 replicates (6 columns)

    level / kind : transcription / RNAseq
    provides     : expr_sporulated
    coverage     : 7,974 (98.0%)
    accession    : GSE206344
    url          : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206344

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.expression.load_all()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gse206344.py
"""
from _common import run

KEY = "gse206344"

if __name__ == "__main__":
    run(KEY, sheet='1', normalised_by='expression.load_all()')
