#!/usr/bin/env python3
"""Promoter accessibility (ATAC-seq, COMPUTED)

Mean ATAC coverage over the promoter, relative to the genome mean

    level / kind : DNA / ATACseq
    provides     : atac_promoter_ut
    coverage     : 7,988 genes (98%)
    accession    : GSE313048
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE313nnn/GSE313048/suppl/GSE313048_ATACseq_GCN5b-KD_UT.bw
    local path   : datasets/quarantine/2026_08_16_unverified/Tg/acetylation/

Quirks that cost time once:
    COMPUTED here: GEO serves this as a bigWig with no peak calls and no per-gene table, so the
    summary is the mean coverage 1 kb either side of the transcription start, over the genome
    mean. No peak calling and no thresholds -- those would be modelling choices invented here
    rather than taken from the authors. Only the UNTREATED arm is read; the deposit is a GCN5b
    knockdown and the depleted arm answers a different question. Verified against expression:
    rho +0.50, and the top decile of expressed genes carries 2.0 log2 more promoter signal than
    the bottom decile. Correlates with fitness at rho +0.02, so it is not merely tracking
    essentiality.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gse313048_atac.py
"""
from _common import run

KEY = "gse313048_atac"

if __name__ == "__main__":
    run(KEY)
