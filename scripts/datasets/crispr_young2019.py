#!/usr/bin/env python3
"""Young 2019 in vivo screen

In vivo fitness

    level / kind : DNA / CRISPR_screen
    provides     : fit_invivo_young2019
    coverage     : 115

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
