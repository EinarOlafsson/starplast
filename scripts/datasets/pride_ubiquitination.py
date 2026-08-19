#!/usr/bin/env python3
"""Ubiquitination / SUMOylation (GlyGly)

GlyGly sites reported per gene

    level / kind : post_translation / proteomics
    provides     : n_ubiquitination_sites
    coverage     : 128 genes measured
    accession    : PXD042937
    url          : https://www.ebi.ac.uk/pride/archive/projects/PXD042937
    local path   : datasets/quarantine/2026_08_16_pride/Tg/ubiquitination_SUMOylation/

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pride_ubiquitination.py
"""
from _common import run

KEY = "pride_ubiquitination"

if __name__ == "__main__":
    run(KEY)
