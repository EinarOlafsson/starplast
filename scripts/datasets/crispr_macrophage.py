#!/usr/bin/env python3
"""Macrophage CRISPR screens

Naive BMDM and IFN-gamma survival

    level / kind : DNA / CRISPR_screen
    provides     : fit_naive_bmdm, fit_ifng
    coverage     : 7,402 (90.9%)
    citation     : Forward genetics screens using macrophages to identify Toxoplasma gondii genes important for resistance to IFN-gamma (2015) -- CONFIRM this is the source
    PMID         : 25867017

Quirks that cost time once:
    Sign convention is INVERTED relative to the other screens.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.screens.crispr_screens()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/crispr_macrophage.py
"""
from _common import run

KEY = "crispr_macrophage"

if __name__ == "__main__":
    run(KEY, normalised_by='screens.crispr_screens()')
