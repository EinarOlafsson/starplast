#!/usr/bin/env python3
"""S-nitrosylation (iodoTMT)

S-nitrosylation sites reported per gene

    level / kind : post_translation / proteomics
    provides     : n_nitrosylation_sites
    coverage     : 660 genes measured
    accession    : PXD046083
    url          : https://www.ebi.ac.uk/pride/archive/projects/PXD046083
    local path   : datasets/quarantine/2026_08_16_pride/Tg/S_nitrosylation/

Quirks that cost time once:
    Counted only from the iodoTMT tables. Counting the whole txt folder put 90% of the proteome
    in this slot, which is what a modification measured on nearly every gene should always look
    like: a bug.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pride_nitrosylation.py
"""
from _common import run

KEY = "pride_nitrosylation"

if __name__ == "__main__":
    run(KEY)
