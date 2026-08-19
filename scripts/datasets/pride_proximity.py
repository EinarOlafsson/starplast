#!/usr/bin/env python3
"""Proximity labelling

Proximity partners reported per gene

    level / kind : post_translation / proteomics
    provides     : n_proximity_partners
    coverage     : 1,734 genes measured
    accession    : PXD059579
    url          : https://www.ebi.ac.uk/pride/archive/projects/PXD059579
    local path   : datasets/quarantine/2026_08_16_pride/Tg/interaction_proximity_labelling/

Quirks that cost time once:
    mzIdentML rather than MaxQuant, and shipped as a lone .gz -- which is one compressed file
    and not an archive, a distinction that read as an empty deposit until it was handled.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pride_proximity.py
"""
from _common import run

KEY = "pride_proximity"

if __name__ == "__main__":
    run(KEY)
