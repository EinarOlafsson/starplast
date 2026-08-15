#!/usr/bin/env python3
"""Stage transcriptome

Tachyzoite, day 3/5/7, in vivo tissue cyst (12 columns)

    level / kind : transcription / RNAseq
    provides     : expr_tachy, expr_cyst, expr_max, rna108740_Tachyzoites_T2_FPKM, rna108740_Tachyzoites_T4_FPKM, rna108740_Tissue_cysts_A_FPKM, rna108740_Tissue_cysts_B_FPKM, rna108740_Day3_CS4_FPKM, rna108740_Day3_CS5_FPKM, rna108740_Day3_CS6_FPKM, rna108740_Day5_CS2_FPKM, rna108740_Day5_K4_14_FPKM, rna108740_Day7_CS7_FPKM, rna108740_Day7_CS8_FPKM, rna108740_Day7_CS9_FPKM
    coverage     : 7,739 (95.1%)
    accession    : GSE108740
    url          : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108740

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.expression.load_all()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gse108740.py
"""
from _common import run

KEY = "gse108740"

if __name__ == "__main__":
    run(KEY, normalized_by='expression.load_all()')
