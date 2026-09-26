#!/usr/bin/env python3
"""Genome-wide mRNA decay after actinomycin D

Wild-type mRNA remaining after 4 h of transcription block, relative to the median

    level / kind : transcription / RNAseq
    provides     : mrna_log2_remaining_4h_actinomycin
    coverage     : 5,944 (73%)
    citation     : Giuliano CJ et al., Convergent evolution of metabolic regulation governs redox adaptation in Toxoplasma. Cell 2026
    PMID         : 42580337
    accession    : GSE329845
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE329nnn/GSE329845/suppl/GSE329845_Processed_Data.csv.gz
    local path   : datasets/transcription/RNAseq/GSE329845/GSE329845_Processed_Data.csv.gz

Quirks that cost time once:
    Raw counts, one 4 h time point, three replicates of vehicle and actinomycin D in wild type
    (the TgPRO knockout arms are not used). Median-of-ratios size factors, a moderated linear
    model with replicate blocking (starplast/deposits.py). No spike-in, so a global loss of RNA
    is normalized away: 0 means as stable as the typical transcript, not fully stable, and this
    is a relative decay, never a half-life. Replicates agree at rho 0.96; ribosomal-protein
    mRNAs are stable (+1.4 against -0.15, p = 6e-14). Genome-wide where mrna_stability is the
    412-gene unstable tail of another study; the two do not correlate (rho -0.03 on 286 shared
    genes), which that tail's selection explains, so they are kept as separate columns in
    different units.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/mrna_decay_gse329845.py
"""
from _common import run

KEY = "mrna_decay_gse329845"

if __name__ == "__main__":
    run(KEY)
