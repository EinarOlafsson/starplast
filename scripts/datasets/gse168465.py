#!/usr/bin/env python3
"""Primary brain-cell parasite differentiation time course

Parasite base mean and log2 fold-change at days 1, 2, 4, 7 and 14

    level / kind : transcription / RNAseq
    provides     : brain168465_1d_base_mean, brain168465_1d_lfc, brain168465_2d_base_mean, brain168465_2d_lfc, brain168465_4d_base_mean, brain168465_4d_lfc, brain168465_7d_base_mean, brain168465_7d_lfc, brain168465_14d_base_mean, brain168465_14d_lfc
    coverage     : measured at build time
    citation     : Mouveaux T et al., Primary brain cell infection by Toxoplasma gondii reveals spontaneous bradyzoite differentiation and modification of neuron biology
    PMID         : 34610266
    accession    : GSE168465
    url          : https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168465
    local path   : datasets/stagetranscriptome_GSE168465_DESeq2-Toxo-all-time-points.xlsx

Quirks that cost time once:
    Dual host-parasite RNA-seq; only the workbook explicitly containing Toxoplasma gene results
    enters this map. p-values remain evidence metadata, not features.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gse168465.py
"""
from _common import run

KEY = "gse168465"

if __name__ == "__main__":
    run(KEY)
