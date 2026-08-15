#!/usr/bin/env python3
"""GRA12 strains and mouse subspecies

Median L2FC in vitro and in vivo, DISCO score; two screens

    level / kind : DNA / CRISPR_screen
    provides     : crispr_gra12s1_l2fc_invitro, crispr_gra12s1_l2fc_invivo, crispr_gra12s1_disco, crispr_gra12s2_l2fc_invitro, crispr_gra12s2_l2fc_invivo, crispr_gra12s2_disco
    coverage     : 236 / 232
    citation     : GRA12 is a common virulence factor across Toxoplasma gondii strains and mouse subspecies
    PMID         : 40240328
    url          : https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-025-58876-2/MediaObjects/41467_2025_58876_MOESM5_ESM.xlsx
    local path   : datasets/DNA/CRISPR_screen/40240328/

Quirks that cost time once:
    The two screens are NOT replicates: in-vivo L2FC correlate at r = 0.41.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.screens.crispr_screens()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gra12.py
"""
from _common import run

KEY = "gra12"

if __name__ == "__main__":
    run(KEY, normalized_by='screens.crispr_screens()')
