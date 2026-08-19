#!/usr/bin/env python3
"""Plasmodium to human contacts (crosslinking MS)

Parasite protein to erythrocyte protein, measured as a crosslink

    level / kind : reference / crosslink_MS
    provides     : bridge:host
    coverage     : 10 pairs, 10 parasite genes, 7 human proteins
    PMID         : 41966402
    accession    : Cell Rep mmc1 sheet D
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13200099/supplementaryFiles
    local path   : starplast/data/pf_host_bridges.parquet

Quirks that cost time once:
    The half of the crosslink file the edge layer throws away, and it is a BRIDGE rather than an
    edge for the reason instruction 39 gives: the pair's two ends live in different tables and a
    human protein has no index in this one. Second bridge in the project after the Toxoplasma
    host IP-MS one, and the first thing it needed was a species-aware bridge lookup -- both arms
    key their bridge `host`, because both cross to a human protein, so the name cannot say whose
    contacts these are and only the parasite end can. Validated on an interaction that is in the
    textbooks: MESA (PF3D7_0500800) crosslinks to erythrocyte ankyrin, and the rest of the human
    side is stomatin, calpain, actin and spectrin beta -- the membrane skeleton, which is what
    an exported parasite protein should be touching.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_host_bridge_xlms.py
"""
from _common import run

KEY = "pf_host_bridge_xlms"

if __name__ == "__main__":
    run(KEY)
