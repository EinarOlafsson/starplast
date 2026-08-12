#!/usr/bin/env python3
"""Young 2019 in vivo screen

In vivo fitness

    level / kind : DNA / CRISPR_screen
    provides     : fit_invivo_young2019
    coverage     : 115
    citation     : Young J et al., A CRISPR platform for targeted in vivo screens identifies Toxoplasma gondii virulence factors in mice. Nat Commun 2019;10:3963
    PMID         : 31481656
    url          : https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-11855-w/MediaObjects/41467_2019_11855_MOESM6_ESM.xlsx

Quirks that cost time once:
    Same paper as invivo_platform, a different supplementary table (MOESM6 vs MOESM5). TARGETED,
    not genome-wide: the libraries are 200, 800 and 3200 gRNAs, which is why this covers 115
    genes rather than the proteome.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.screens.crispr_screens()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/crispr_young2019.py
"""
from _common import run

KEY = "crispr_young2019"

if __name__ == "__main__":
    run(KEY, normalised_by='screens.crispr_screens()')
