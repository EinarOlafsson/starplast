#!/usr/bin/env python3
"""Transcription in a drug-tolerant latent state

Change in each gene in latent parasites, and the paper's 200-gene latency classifier

    level / kind : transcription / scRNAseq
    provides     : latency_log2fc, latency_classifier_member
    coverage     : 4,887 genes (85%)
    citation     : Aryal A et al., Environmental stress promotes entry into a pre-existing latent state in Plasmodium falciparum. bioRxiv 2026, doi:10.64898/2026.09.13.751295 (preprint)
    accession    : bioRxiv 10.64898/2026.09.13.751295 Supplementary Tables 7 and 24
    url          : https://www.biorxiv.org/content/10.64898/2026.09.13.751295v1.supplementary-material
    local path   : datasets/transcription/scRNAseq/pf_latency_2026/ST07_media-4.xlsx

Quirks that cost time once:
    Dormancy is how this parasite survives artemisinin, and the slot for it was empty. Two
    things the deposit offers are deliberately NOT shipped. Its FDR: cells were treated as
    replicates, so 87% of genes reach padj 0.05 and the column would call almost everything
    significant. And absence as evidence: ribosomal protein genes and most var and rifin genes
    are missing from the table entirely, so a gene with no value was not shown to be unchanged.
    No GEO or SRA accession is given in the preprint, which is why the accession field names the
    supplement.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_latency_transcriptome.py
"""
from _common import run

KEY = "pf_latency_transcriptome"

if __name__ == "__main__":
    run(KEY)
