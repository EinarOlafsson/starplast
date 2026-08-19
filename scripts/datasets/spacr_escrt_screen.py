#!/usr/bin/env python3
"""Host ESCRT recruitment screen (UNPUBLISHED)

Per-gene effect on host TSG101 recruitment to the vacuole, by two models

    level / kind : DNA / CRISPR_screen
    provides     : escrt_recruitment_xgboost, escrt_recruitment_maxvit
    coverage     : 13 and 8 genes
    citation     : Olafsson EB et al., A pooled image-based CRISPR screen identifies EAF1 as a T. gondii modulator of ESCRT subversion. bioRxiv 2026 (under submission)
    accession    : spaCR screen, bioRxiv 10.64898/2026.07.08.737057
    url          : https://doi.org/10.64898/2026.07.08.737057
    local path   : datasets/quarantine/2026_08_16_unverified/Tg/escrt_screen/

Quirks that cost time once:
    A pooled image-based screen of secretory proteins, deconvolved to gene effects by
    regression. The phenotype is host ESCRT recruitment and NOT invasion or egress, so it fills
    a slot of its own rather than the invasion slot whose high-content-imaging context it
    matches -- a slot names a question, not a method. Both deconvolution models are shipped
    because their agreement is the verification, and the agreement is close: TGGT1_244480 is
    rank 1 in each, TGGT1_409250 and GRA14 fill out the top three of each, and EAF1 -- which the
    PRIDE deposit PXD080696 identifies as TGGT1_225160 -- is rank 7 in BOTH. Two independent
    deconvolutions landing the same gene at the same rank is a stronger statement than any
    single ordering. An earlier version of this note called TGGT1_244480 EAF1 and said EAF1 was
    rank 1; both were wrong, and the deposit's own abstract is what settled it. The MYR1 host
    bridge added the same day is verified by ESCRT machinery topping it, so two unrelated
    datasets in this map now point at the same biology.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/spacr_escrt_screen.py
"""
from _common import run

KEY = "spacr_escrt_screen"

if __name__ == "__main__":
    run(KEY)
