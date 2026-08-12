#!/usr/bin/env python3
"""In vitro CRISPR fitness (HFF)

Competitive growth in fibroblasts

    level / kind : DNA / CRISPR_screen
    provides     : fit_invitro_hff
    coverage     : 7,325 (90.0%)
    citation     : A Genome-wide CRISPR Screen in Toxoplasma Identifies Essential Apicomplexan Genes (Sidik et al. 2016)
    PMID         : 27594426
    url          : https://ars.els-cdn.com/content/image/1-s2.0-S0092867416310704-mmc3.xlsx

Quirks that cost time once:
    Competitive growth, NOT essentiality. Predicted from protein features at R2 = 0.453, while
    the other screens are predicted at -0.105 to +0.102.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.screens.crispr_screens()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/crispr_invitro.py
"""
from _common import run

KEY = "crispr_invitro"

if __name__ == "__main__":
    run(KEY, normalised_by='screens.crispr_screens()')
