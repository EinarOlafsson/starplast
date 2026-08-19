#!/usr/bin/env python3
"""m6A methylome (MeRIP peaks)

How many m6A peaks the authors called on this gene in tachyzoites

    level / kind : transcription / MeRIP
    provides     : n_m6a_peaks
    coverage     : 837 genes (10%)
    PMID         : 34324585
    accession    : PLoS Pathogens 1009335 Table S3A
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8354455/supplementaryFiles
    local path   : starplast/data/m6a_peaks.tsv

Quirks that cost time once:
    GEO carries no MeRIP for Toxoplasma -- its deposit for this study is depletion RNA-seq,
    which says which transcripts DEPEND on m6A and not which CARRY it. The peaks are in the
    paper. This slot had already been written up as unservable in instruction 41 when the
    supplement turned up, and that entry is now struck through rather than deleted. 866 of 8,922
    genes carry a peak, which is the right order for m6A. Verified against the paper's own
    second dataset: marked genes are enriched among those responding to METTL3 depletion, odds
    1.35, p = 2e-03 -- modest because removing a writer has broad indirect effects, but the
    direction a writer's own substrates have to take.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/m6a_peaks.py
"""
from _common import run

KEY = "m6a_peaks"

if __name__ == "__main__":
    run(KEY)
