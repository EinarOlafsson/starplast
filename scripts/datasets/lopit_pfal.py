#!/usr/bin/env python3
"""P. falciparum LOPIT

Donor labels for orthoLOPIT transfer

    level / kind : post_translation / LOPIT
    provides     : (edges or build inputs only)
    coverage     : 1,646 usable
    local path   : datasets/lopit_plasmodium_falciparum_3D7.csv

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.localisation.lopit_labels()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/lopit_pfal.py
"""
from _common import run

KEY = "lopit_pfal"

if __name__ == "__main__":
    run(KEY, normalised_by='localisation.lopit_labels()')
