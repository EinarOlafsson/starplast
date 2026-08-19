#!/usr/bin/env python3
"""Lysine lactylome

Lactylation sites reported per gene

    level / kind : post_translation / proteomics
    provides     : n_lactylation_sites
    coverage     : 515 genes measured
    accession    : PXD031526
    url          : https://www.ebi.ac.uk/pride/archive/projects/PXD031526
    local path   : datasets/quarantine/2026_08_16_pride/Tg/lactylation/

Quirks that cost time once:
    The second lactylation deposit tried. PXD022700 ships a RAR that bsdtar cannot open, and the
    0 genes that produced was a fact about the reader, not about a study that names 537
    proteins. This one reads, and its `La (K)Sites` table is MaxQuant's lactylation search.
    Keyed entirely on TGGT1_ accessions -- 515 of its 524 genes reach the map through the
    identity layer and would reach none without it.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pride_lactylation.py
"""
from _common import run

KEY = "pride_lactylation"

if __name__ == "__main__":
    run(KEY)
