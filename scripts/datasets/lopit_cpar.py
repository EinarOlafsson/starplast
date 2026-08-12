#!/usr/bin/env python3
"""C. parvum hyperLOPIT

Donor labels for orthoLOPIT transfer

    level / kind : post_translation / LOPIT
    provides     : (edges or build inputs only)
    coverage     : 1,107 usable
    citation     : Guerin et al. 2023
    url          : https://ars.els-cdn.com/content/image/1-s2.0-S1931312823001051-mmc4.xlsx
    local path   : datasets/lopit_cryptosporidium_parvum_MEASURED_Guerin2023.csv

Quirks that cost time once:
    MASTER_parasite_wide_by_orthogroup.csv has an empty cpar_lopit_native column; this data is
    joined from source instead.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.localisation.lopit_labels()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/lopit_cpar.py
"""
from _common import run

KEY = "lopit_cpar"

if __name__ == "__main__":
    run(KEY, normalised_by='localisation.lopit_labels()')
