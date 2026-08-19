#!/usr/bin/env python3
"""HDAC3 occupancy (CUT&TAG, COMPUTED)

Mean HDAC3 CUT&TAG coverage over the promoter, relative to the genome mean

    level / kind : DNA / CUTandTAG
    provides     : cuttag_hdac3_promoter_ut
    coverage     : 8,140 genes (100%)
    accession    : GSE277553
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277553/suppl/GSE277553_RAW.tar
    local path   : datasets/quarantine/2026_08_16_unverified/Tg/chromatin_accessibility/

Quirks that cost time once:
    COMPUTED the same way as the ATAC column, from the untreated arm; the deposit's other arm is
    an AP2XII-5 knockout. Two replicates and not three: GSM8524430_UT_2.bw begins with eight
    0xFF bytes and is not a bigWig, identically whether taken from the series tar or fetched
    from GEO as a sample file, so the corruption is in the deposit. It is skipped with a message
    rather than silently, because a replicate dropped without saying so makes the mean smaller
    than the note beside the column claims.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gse277553_cuttag.py
"""
from _common import run

KEY = "gse277553_cuttag"

if __name__ == "__main__":
    run(KEY)
