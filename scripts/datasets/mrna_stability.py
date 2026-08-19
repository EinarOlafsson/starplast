#!/usr/bin/env python3
"""mRNA stability after actinomycin D

Proportion of transcript remaining after five hours of transcription block

    level / kind : transcription / RNAseq
    provides     : mrna_remaining_5h_actinomycin
    coverage     : 412 genes
    PMID         : 39899594
    accession    : PLoS Pathogens 1012857 Table S12
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC11801735/supplementaryFiles
    local path   : starplast/data/mrna_stability.tsv

Quirks that cost time once:
    A direct measurement: block transcription, wait, see what is left. Untreated parasites at
    five hours, so the column is stability and not the iron response the paper is about. BIASED
    BY CONSTRUCTION and the bias is worth stating -- the table is the 426 transcripts that fell
    below 75% remaining, so it describes the unstable tail and a gene absent from it is stable
    OR was not measured, which the column cannot distinguish. Consistent with that: ribosomal-
    protein transcripts, which are classically stable, are under-represented among the
    responders at odds 0.37.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/mrna_stability.py
"""
from _common import run

KEY = "mrna_stability"

if __name__ == "__main__":
    run(KEY)
