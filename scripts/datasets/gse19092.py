#!/usr/bin/env python3
"""Synchronized tachyzoite cell-cycle transcriptome

Two replicates across blocked, asynchronous and hourly release states

    level / kind : transcription / microarray
    provides     : cellcycle19092_async_r1, cellcycle19092_async_r2, cellcycle19092_blocked_r1, cellcycle19092_blocked_r2, cellcycle19092_1h_r1, cellcycle19092_1h_r2, cellcycle19092_2h_r1, cellcycle19092_2h_r2, cellcycle19092_3h_r1, cellcycle19092_3h_r2, cellcycle19092_4h_r1, cellcycle19092_4h_r2, cellcycle19092_5h_r1, cellcycle19092_5h_r2, cellcycle19092_6h_r1, cellcycle19092_6h_r2, cellcycle19092_7h_r1, cellcycle19092_7h_r2, cellcycle19092_8h_r1, cellcycle19092_8h_r2, cellcycle19092_9h_r1, cellcycle19092_9h_r2, cellcycle19092_10h_r1, cellcycle19092_10h_r2, cellcycle19092_11h_r1, cellcycle19092_11h_r2, cellcycle19092_12h_r1, cellcycle19092_12h_r2
    coverage     : measured at build time
    citation     : Behnke MS et al., Coordinated progression through two subtranscriptomes underlies the tachyzoite cycle of Toxoplasma gondii. PLoS ONE 2010;5:e12354
    PMID         : 20865045
    accession    : GSE19092
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE19nnn/GSE19092/matrix/GSE19092_series_matrix.txt.gz
    local path   : datasets/toxoplasma_acquisition_2026_08_14/GSE19092_series_matrix.txt.gz

Quirks that cost time once:
    Legacy GPL7186 probes are mapped through the platform ToxoDB field and the project's
    previous-ID resolver.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gse19092.py
"""
from _common import run

KEY = "gse19092"

if __name__ == "__main__":
    run(KEY)
