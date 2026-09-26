#!/usr/bin/env python3
"""5' UTR architecture from reannotated transcripts

Length, upstream AUGs and ORFs, and start-context strength of each 5' UTR

    level / kind : transcription / RNAseq
    provides     : utr5_length, utr5_n_uaugs, utr5_n_uorfs, utr5_n_oorfs, utr5_n_inframe_ext, utr5_kozak_score
    coverage     : 5,992 (73.6%)
    citation     : Peters ML et al., 5' untranslated regions tune Toxoplasma translation. bioRxiv 2025, doi:10.1101/2025.07.14.664749 (preprint)
    accession    : GSE302108 / bioRxiv 10.1101/2025.07.14.664749 Supplementary Data 4
    url          : https://www.biorxiv.org/content/10.1101/2025.07.14.664749v1.supplementary-material
    local path   : datasets/translation/riboseq/GSE302108/SupplementaryData4_TE_info_and_UTRfeatures.xlsx

Quirks that cost time once:
    Sequence properties of the untranslated region, from a reannotation built on long reads and
    ribosome footprints -- not a measurement of the gene's behaviour, which is why they are
    their own slot rather than part of translation. They predict translation the way they
    should: upstream AUGs against efficiency at rho -0.47, Kozak strength with it at +0.22. The
    study's reporter assay is NOT here: its 30,235 scored sequences are variants of twelve
    endogenous UTRs, so it describes sequences, not genes.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gse302108_utr5.py
"""
from _common import run

KEY = "gse302108_utr5"

if __name__ == "__main__":
    run(KEY)
