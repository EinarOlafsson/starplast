#!/usr/bin/env python3
"""Pru proteome and IP abundance

Median log2 iBAQ across replicates

    level / kind : translation / proteomics
    provides     : protein_ibaq_log2
    coverage     : 748 (9.2%)
    accession    : PXD043808, PXD065585
    url          : https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD065585
    local path   : toxo_stage_atlas/data/proteomics/

Quirks that cost time once:
    Immunoprecipitation experiments of 424 and 594 proteins. Enrichment, NOT a deep proteome; do
    not report as proteome-wide.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.screens.proteomics()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/proteome_pru.py
"""
from _common import run

KEY = "proteome_pru"

if __name__ == "__main__":
    run(KEY, normalised_by='screens.proteomics()')
