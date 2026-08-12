#!/usr/bin/env python3
"""ToxoDB gene identity

Symbols, previous IDs, product descriptions

    level / kind : reference / identity
    provides     : gene_id, product
    coverage     : 8,843 ME49 genes
    accession    : ToxoDB ME49
    url          : https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByTaxon/reports/attributesTabular
    local path   : starplast/data/toxodb_identity.tsv

Quirks that cost time once:
    Retrieved 2026-08-11 via the REST API; strain tables for GT1 and VEG alongside.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/toxodb_identity.py
"""
from _common import run

KEY = "toxodb_identity"

if __name__ == "__main__":
    run(KEY)
