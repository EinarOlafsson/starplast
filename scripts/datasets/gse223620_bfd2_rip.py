#!/usr/bin/env python3
"""BFD2-bound transcriptome (RIP-seq, COMPUTED)

Enrichment of each transcript in the BFD2 immunoprecipitation

    level / kind : transcription / RIPseq
    provides     : bfd2_rip_log2_ip_over_input
    coverage     : 7,463 genes (92%)
    accession    : GSE223620
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223620/suppl/GSE223620_ProcessedDataFile_BFD2.RIPseq.xls.gz
    local path   : datasets/quarantine/2026_08_16_unverified/Tg/rna_binding_protein_targets/

Quirks that cost time once:
    COMPUTED: the deposit publishes read counts for the IP and the input, not the ratio. Both
    are scaled to a common library size and the log2 taken, at a floor of 20 reads across the
    pair -- deliberately low, because the point of a RIP is the enriched tail and a stricter
    floor would drop the genes the slot asks about. Verified against the published mechanism:
    BFD2 binds and stabilises the BFD1 transcript, and BFD1 comes out at +3.52, rank 32 of
    8,090, the top 0.4%. BFD2's own transcript is unremarkable at +0.31, which is what says the
    enrichment is not an artefact of the tagged locus. One column and not a set: the slot holds
    one protein's targets, and a second RIP would sit beside this rather than be averaged into
    it.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gse223620_bfd2_rip.py
"""
from _common import run

KEY = "gse223620_bfd2_rip"

if __name__ == "__main__":
    run(KEY)
