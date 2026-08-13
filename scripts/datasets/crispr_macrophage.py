#!/usr/bin/env python3
"""Macrophage CRISPR screens

Naive BMDM and IFN-gamma survival

    level / kind : DNA / CRISPR_screen
    provides     : fit_naive_bmdm, fit_ifng
    coverage     : 7,402 (90.9%)
    citation     : Wang Y et al., Genome-wide screens identify Toxoplasma gondii determinants of parasite fitness in IFN-gamma-activated murine macrophages. Nat Commun 2020;11:5258
    PMID         : 33067458
    url          : https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-020-18991-8/MediaObjects/41467_2020_18991_MOESM5_ESM.xlsx

Quirks that cost time once:
    Sign convention is INVERTED relative to the other screens. THE PREVIOUS CITATION WAS WRONG
    and its own 'CONFIRM this is the source' warning was justified: PMID 25867017 is a 2015 JoVE
    video protocol using CHEMICAL mutagenesis, verified against PubMed -- not a CRISPR screen,
    and it predates the first one. The entry had copied that protocol's title verbatim. A video
    protocol cannot be the source of 7,402 per-gene fitness scores; Wang 2020 is genome-wide in
    IFN-gamma-activated macrophages, which is exactly what these two columns are.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.screens.crispr_screens()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/crispr_macrophage.py
"""
from _common import run

KEY = "crispr_macrophage"

if __name__ == "__main__":
    run(KEY, normalized_by='screens.crispr_screens()')
