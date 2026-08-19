#!/usr/bin/env python3
"""Crosslinking MS interactome

How many proteins this one crosslinks to

    level / kind : post_translation / XLMS
    provides     : n_crosslink_partners
    coverage     : 494 proteins
    PMID         : 40874616
    accession    : mBio 02159-25 supplementary file s0004
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12505969/supplementaryFiles
    local path   : starplast/data/crosslink_partners.tsv

Quirks that cost time once:
    395 high-confidence protein pairs, counted per protein. Verified against the two largest
    obligate complexes any cell has: 21 of 32 proteasome subunits are in the interactome (odds
    30.8, p = 1e-18) and 57 of 158 ribosomal proteins (odds 9.7, p = 4e-30). Crosslinking finds
    stable abundant complexes, and if it did not find those two it would not be finding
    complexes. Absent is absent: a protein with no partner here may be in no complex or may
    simply not have crosslinked.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/crosslink_interactome.py
"""
from _common import run

KEY = "crosslink_interactome"

if __name__ == "__main__":
    run(KEY)
