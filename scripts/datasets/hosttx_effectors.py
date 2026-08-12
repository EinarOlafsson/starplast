#!/usr/bin/env python3
"""Host-transcription effector screen

Hotelling T2 statistic per effector, adjusted p

    level / kind : DNA / CRISPR_screen
    provides     : hosttx_T2, hosttx_padj
    coverage     : 252
    citation     : High-throughput identification of Toxoplasma gondii effector proteins that target host cell transcription
    PMID         : 37827122
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12033024/supplementaryFiles
    local path   : datasets/DNA/CRISPR_screen/37827122/

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.screens.crispr_screens()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/hosttx_effectors.py
"""
from _common import run

KEY = "hosttx_effectors"

if __name__ == "__main__":
    run(KEY, normalised_by='screens.crispr_screens()')
