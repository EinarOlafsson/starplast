#!/usr/bin/env python3
"""PlasmoDB gene identity

Symbols, previous IDs, product descriptions for the Plasmodium arm

    level / kind : reference / identity
    provides     : gene_id, product
    coverage     : 5,791 P. falciparum 3D7 transcripts
    accession    : PlasmoDB 3D7
    url          : https://plasmodb.org/plasmo/service/record-types/transcript/searches/GenesByTaxon/reports/attributesTabular
    local path   : starplast/data/plasmodb_identity.tsv

Quirks that cost time once:
    The Plasmodium arm's identity layer, fetched 2026-08-18 through the same WDK report as the
    ToxoDB tables and kept in its own file -- two identifier spaces in one index is the merge
    this project refuses everywhere else. It exists because the first Pf source keyed on
    anything but current accessions joined ZERO rows: the 2014 ribosome-profiling deposit
    reports the pre-2012 chromosome-based ids (`PFE0630c`, `PF13_0222`), and a string join found
    none of its 3,629 genes. 9,106 previous ids resolve; 66 are claimed by two current genes
    each -- a gene model SPLIT, seen from the other side -- and those are withdrawn rather than
    assigned to whichever came first, the same rule the Toxoplasma layer applies to 153 strings.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/plasmodb_identity.py
"""
from _common import run

KEY = "plasmodb_identity"

if __name__ == "__main__":
    run(KEY)
