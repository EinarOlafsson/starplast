#!/usr/bin/env python3
"""Pru tachyzoite / 72-hour bradyzoite stage array

Matched tachyzoite and alkaline-induced bradyzoite expression

    level / kind : transcription / microarray
    provides     : rna22258_tachyzoite, rna22258_bradyzoite
    coverage     : 7,253 genes
    accession    : GSE22258
    url          : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE22258
    local path   : datasets/stagetranscriptome_GSE22258_series_matrix.txt.gz

Quirks that cost time once:
    Already keyed by TGME49 accessions. Kept separate from the newer RNA-seq stage series; it is
    not averaged as though microarray intensity were FPKM.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gse22258.py
"""
from _common import run

KEY = "gse22258"

if __name__ == "__main__":
    run(KEY)
