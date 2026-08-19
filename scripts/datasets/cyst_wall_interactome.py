#!/usr/bin/env python3
"""Cyst wall interactome

Strongest bait signal and how many baits saw the protein

    level / kind : post_translation / IPMS
    provides     : cyst_wall_max_spectral, cyst_wall_n_baits
    coverage     : 56 proteins
    PMID         : 32019789
    accession    : PMC7002340 Data Set S1
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7002340/supplementaryFiles
    local path   : datasets/quarantine/2026_08_16_unverified/Tg/cyst_wall/

Quirks that cost time once:
    Two numbers because the bait count is the more honest one: a protein found by thirteen
    independent pulldowns is in the cyst wall in a way a single strong hit is not. Verified by
    what comes out on top -- MAG1 and MAG2, the canonical cyst matrix proteins, with MCP3, MCP4
    and SRS44 beside them. 57 of the table's 265 rows name a Toxoplasma accession; the rest are
    human, because the pulldowns were done on infected cultures and the table lists everything
    identified.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/cyst_wall_interactome.py
"""
from _common import run

KEY = "cyst_wall_interactome"

if __name__ == "__main__":
    run(KEY)
