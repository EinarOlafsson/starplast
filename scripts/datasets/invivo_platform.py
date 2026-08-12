#!/usr/bin/env python3
"""In vivo CRISPR platform

Mean log fold-change across replicates

    level / kind : DNA / CRISPR_screen
    provides     : crispr_invivo_platform_lfc
    coverage     : 168
    citation     : A CRISPR platform for targeted in vivo screens identifies Toxoplasma gondii virulence factors in mice
    PMID         : 31481656
    url          : https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-019-11855-w/MediaObjects/41467_2019_11855_MOESM5_ESM.xlsx
    local path   : datasets/DNA/CRISPR_screen/31481656/

Quirks that cost time once:
    Cites pre-2012 TGME49_0xxxxx accessions for every gene; must go through identity.py.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.screens.crispr_screens()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/invivo_platform.py
"""
from _common import run

KEY = "invivo_platform"

if __name__ == "__main__":
    run(KEY, normalised_by='screens.crispr_screens()')
