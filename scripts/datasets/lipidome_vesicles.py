#!/usr/bin/env python3
"""Membrane lipid composition of parasite vesicles

Lipid species abundance, and its proportion against the host cell

    level / kind : reference / lipidomics
    provides     : lipid_ev_level_log2, lipid_ev_vs_host_clr
    coverage     : 194 lipid species
    PMID         : 41716462
    accession    : Front Cell Infect Microbiol 1745625 Tables 1-3
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12913473/supplementaryFiles
    local path   : starplast/data/metabolites.parquet

Quirks that cost time once:
    The SECOND study in the metabolite table, and the name-matching cost the entry above warns
    about is not paid here: lipid shorthand and polar compound names are different naming
    systems and none of the 194 species collides with the 1,102 compounds. Rows are appended,
    not joined. These are read as PARASITE lipids because the vesicles came from post-egress
    tachyzoites in host-cell-free medium, and because the composition does not move when the
    host does: across four host backgrounds the host cells differ in 1,018-1,362 species and the
    vesicles the same parasite released in them differ in 0-4. A lipidome of an infected culture
    would not have supported this slot at all -- most of that lipid is host. The `vs_host`
    column is COMPUTED here, sample-centred so it compares proportion rather than amount; the
    archive's own EV-minus-cell column is not used because its transform could not be reproduced
    to better than 3 log units.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/lipidome_vesicles.py
"""
from _common import run

KEY = "lipidome_vesicles"

if __name__ == "__main__":
    run(KEY)
