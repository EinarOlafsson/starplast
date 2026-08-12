#!/usr/bin/env python3
"""Phosphosite counts

Count of phosphosites per protein, no positions

    level / kind : post_translation / phosphoproteomics
    provides     : n_phosphosites
    coverage     : 1,175 (14.4%)
    citation     : Treeck M et al. 2011 -- CONFIRM against the file on disk
    url          : https://ars.els-cdn.com/content/image/1-s2.0-S1931312811002885-mmc2.xls

Quirks that cost time once:
    The URL downloads a real phosphoproteomics table, but that it is the source of THIS column
    is inference rather than verification; confirm before citing. Missing for 85.6% of genes;
    effectively an indicator of having been in a phosphoproteomics experiment.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.screens.proteomics()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/phosphosites.py
"""
from _common import run

KEY = "phosphosites"

if __name__ == "__main__":
    run(KEY, normalised_by='screens.proteomics()')
