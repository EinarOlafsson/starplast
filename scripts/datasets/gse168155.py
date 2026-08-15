#!/usr/bin/env python3
"""CPSF4 RNA-processing perturbation transcriptome

RNA response at 7, 24 and 48 hours after CPSF4 depletion

    level / kind : transcription / RNAseq
    provides     : cpsf4rna168155_ut_1, cpsf4rna168155_ut_2, cpsf4rna168155_iaa_7h_1, cpsf4rna168155_iaa_7h_2, cpsf4rna168155_iaa_24h_1, cpsf4rna168155_iaa_24h_2, cpsf4rna168155_iaa_48h_1, cpsf4rna168155_iaa_48h_2
    coverage     : measured at build time
    citation     : Farhat DC et al., A plant-like mechanism coupling m6A reading to polyadenylation safeguards transcriptome integrity. eLife 2021;10:e68312
    PMID         : 34263725
    accession    : GSE168155
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE168nnn/GSE168155/suppl/GSE168155_Matrix_table_processed_data.xlsx
    local path   : datasets/toxoplasma_acquisition_2026_08_14/GSE168155_Matrix_table_processed_data.xlsx

Quirks that cost time once:
    This is a perturbation-response transcriptome, not a direct gene-wise m6A map.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gse168155.py
"""
from _common import run

KEY = "gse168155"

if __name__ == "__main__":
    run(KEY)
