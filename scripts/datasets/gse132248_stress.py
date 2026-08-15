#!/usr/bin/env python3
"""Alkaline-stress differentiation transcriptome

Unstressed tachyzoites and alkaline-stressed bradyzoites

    level / kind : transcription / RNAseq
    provides     : stress_5-1_unstress_red_24, stress_5-2_unstress_red_24, stress_5-24_stress_green_48, stress_5-27_stress_green_48, stress_6-1_unstress_red_24, stress_6-4_unstress_red_24, stress_6-16_stress_green_48, stress_6-18_stress_green_48, stress_7-1_unstress_red_24, stress_7-18_stress_green_48
    coverage     : 7,880 (96.8%)
    citation     : Waldman BS et al., Identification of a Master Regulator of Differentiation in Toxoplasma. Cell 2020;180:359-372.e16
    PMID         : 31955846
    accession    : GSE132248
    url          : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE132248
    local path   : toxo_stage_atlas/data/transcriptomics/GSE132248_STAR_counts_matrix.tsv

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.expression.load_all()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gse132248_stress.py
"""
from _common import run

KEY = "gse132248_stress"

if __name__ == "__main__":
    run(KEY, normalized_by='expression.load_all()')
