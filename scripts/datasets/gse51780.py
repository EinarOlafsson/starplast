#!/usr/bin/env python3
"""Feline merozoite transcriptome

Merozoite expression with matched tachyzoite comparators

    level / kind : transcription / microarray
    provides     : rna51780_tachy_r1, rna51780_tachy_r2, rna51780_mero_r3, rna51780_mero_r4, rna51780_mero_r5, rna51780_mero_r6
    coverage     : measured at build time
    citation     : Behnke MS et al., Toxoplasma gondii merozoite gene expression analysis with comparison to the life cycle. BMC Genomics 2014;15:350
    PMID         : 24885521
    accession    : GSE51780
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE51nnn/GSE51780/matrix/GSE51780_series_matrix.txt.gz
    local path   : datasets/toxoplasma_acquisition_2026_08_14/GSE51780_series_matrix.txt.gz

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gse51780.py
"""
from _common import run

KEY = "gse51780"

if __name__ == "__main__":
    run(KEY)
