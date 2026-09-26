#!/usr/bin/env python3
"""Human fibroblast response to Toxoplasma infection

Infected against uninfected HFF, per host gene: moderated log2 change and FDR

    level / kind : reference / transcription
    provides     : hff_tg_infection_log2fc, hff_tg_infection_padj
    coverage     : 10,631 human genes
    accession    : GSE335016
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE335nnn/GSE335016/suppl/GSE335016_245740-HOMO-gene_fpkm.txt.gz
    local path   : datasets/host/fibroblast_infection/GSE335016/GSE335016_245740-HOMO-gene_fpkm.txt.gz

Quirks that cost time once:
    No paper is attached to the series (Europe PMC and GEO searched 2026-09-26), so the citation
    is the accession. FPKM only; the contrast is computed here (starplast/deposits.py):
    moderated t on log2(FPKM + 1), three replicates per arm, genes at a mean FPKM of 1, wild-
    type infection against uninfected. The TGGT1_245740 knockout arm is not shipped -- it
    answers a question about one parasite gene. Textbook response: CXCL8 +5.4, IL6 +2.7, CXCL10
    +2.3, ISG15 +2.1; GAPDH and ACTB flat. 315 genes at padj < 0.05.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_hff_tg_infection.py
"""
from _common import run

KEY = "host_hff_tg_infection"

if __name__ == "__main__":
    run(KEY)
