#!/usr/bin/env python3
"""Plasmodium complexes from crosslinking MS

Which crosslink-derived complex a gene belongs to, and whether it reaches the host

    level / kind : reference / crosslink_MS
    provides     : complex_id, complex_size, complex_spans_host
    coverage     : 128 genes in 42 complexes
    PMID         : 41966402
    accession    : Cell Rep mmc5 Clusters
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13200099/supplementaryFiles
    local path   : datasets/reference/plasmodb/complexes/41966402/mmc5.xlsx

Quirks that cost time once:
    `complex_spans_host` is the informative column and is kept rather than dropped. Seven of the
    47 clusters contain human proteins as well as parasite ones, which is not contamination --
    the experiment crosslinked parasite inside erythrocyte, so a complex reaching into the host
    is a finding. But a parasite gene in one of those has partners this table cannot name, and a
    reader taking `complex_size` at face value would over-count its parasite neighbours. Only
    parasite members get a row; the host members belong to a bridge table. Same study as the
    crosslink edge layer and a different question: that one is which pairs touch, this one is
    which assembly a protein sits in.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_complexes.py
"""
from _common import run

KEY = "pf_complexes"

if __name__ == "__main__":
    run(KEY)
