#!/usr/bin/env python3
"""Single-parasite transcriptional atlas (cell cycle)

Measured cell-cycle phase per gene, and pseudotime cluster

    level / kind : transcription / scRNAseq
    provides     : cellcycle_phase, cellcycle_pseudotime
    coverage     : 873 genes phased, 7,499 clustered
    citation     : Xue Y et al. eLife 2020;9:e54129
    PMID         : 32065584
    url          : https://cdn.elifesciences.org/articles/54129/elife-54129-supp3-v2.csv
    local path   : datasets/transcription/scRNAseq/32065584/cellcycle_phase_RH.csv

Quirks that cost time once:
    Tab-separated despite the .csv extension. RH files use TGGT1_ accessions and the Pru files
    in the same supplement use TGME49_; the prefix is the only thing that distinguishes them.
    The only MEASURED discrete cell-cycle label in the project.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.cellcycle.add_all()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/xue_singlecell.py
"""
from _common import run

KEY = "xue_singlecell"

if __name__ == "__main__":
    run(KEY, sep='\t', normalised_by='cellcycle.add_all()')
