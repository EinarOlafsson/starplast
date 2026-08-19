#!/usr/bin/env python3
"""N-myristoylated proteome

The authors' confidence that this protein is myristoylated, 3 high to 1 low

    level / kind : post_translation / proteomics
    provides     : myristoylation_confidence
    coverage     : 65 substrates
    PMID         : 32618271
    accession    : eLife 57861 supplementary file 4
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7373427/supplementaryFiles
    local path   : starplast/data/myristoylome.tsv

Quirks that cost time once:
    A category the catalog did not have. N-myristoylation is co-translational and irreversible,
    has its own enzyme in NMT and its own drug programme, and its substrates are published.
    Verified against chemistry rather than against annotation: myristoylation happens on an
    N-terminal glycine, and all 65 of 65 substrates have glycine at position 2 against 5.8% of
    every other gene (Fisher p = 4e-79). No other column in the map can be checked that cleanly.
    Taken from the paper rather than from PXD019677, its PRIDE deposit, which ships MaxQuant
    archives of 250-340 MB apiece; the answer is a 65-row table in supplementary file 4.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/myristoylome.py
"""
from _common import run

KEY = "myristoylome"

if __name__ == "__main__":
    run(KEY)
