#!/usr/bin/env python3
"""Plasmodium crosslinking MS contacts

Protein pairs joined by a measured crosslink

    level / kind : reference / crosslink_MS
    provides     : edge:xlms
    coverage     : 79 parasite-parasite pairs
    PMID         : 41966402
    accession    : Cell Rep mmc1 sheet D
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13200099/supplementaryFiles
    local path   : datasets/reference/plasmodb/crosslink/41966402/mmc1.xlsx

Quirks that cost time once:
    Two filters, neither optional. The experiment crosslinked PARASITE INSIDE ERYTHROCYTE, so a
    third of the 106 protein pairs have a human protein at one or both ends -- spectrin, band 3,
    protein 4.2. Those are real contacts and they are a HOST BRIDGE rather than a parasite-
    parasite edge, so they are dropped from this layer instead of being indexed against a table
    with no row for them. And a pair whose ends resolve to one gene is a homomeric crosslink:
    evidence the protein self-associates, not an edge, and drawing it would put a zero-length
    line in the graph. Validated on complexes that have to be there: EXP2, PTEX150 and HSP101 --
    three subunits of the PTEX translocon -- crosslink to one another, prohibitin 1 to
    prohibitin 2, and RAP1 to RAP2. A contact map that missed those would not be measuring
    contacts. Which end is host and which is parasite goes through the ACCESSION INDEX and not
    through a pattern: the first version matched PF3D7_ in the mapping field and shipped 73
    edges, but the source writes some rows with a UniProt symbol instead -- `sp|Q6ZMA7|Pfs16` is
    PF3D7_0406200, a parasite gene -- so six real contacts were dropped as host-at-one-end and a
    parasite protein was on its way into a host bridge. Resolving recovers all 79, and a test
    fails if the count ever drops back.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_crosslink_ms.py
"""
from _common import run

KEY = "pf_crosslink_ms"

if __name__ == "__main__":
    run(KEY)
