#!/usr/bin/env python3
"""Plasmodium N-myristoylome (NMT-inhibitor sensitive)

Proteins whose click-chemistry capture drops when N-myristoyltransferase is blocked

    level / kind : post_translation / myristoylation
    provides     : is_myristoylated
    coverage     : 16 substrates of 609 assayed
    PMID         : 34695132
    accession    : PLoS Biol 3001408 S11
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8544853/supplementaryFiles
    local path   : datasets/post_translation/myristoylome/34695132/pbio.3001408.s011.xlsx

Quirks that cost time once:
    The evidence for a substrate is not being pulled down -- background comes down too -- but
    coming down LESS when the transferase is inhibited, so the loader requires significance AND
    a negative difference. A positive difference under a blocked transferase would be a protein
    that came down MORE without it, which is not what a substrate does. THREE states, not two:
    16 substrates, 593 assayed and not substrates, and 5,111 genes never in the pulldown, which
    stay missing -- collapsing the last two would claim the whole proteome had been tested for
    myristoylation by one experiment that saw 609 proteins. Sparse because the biology is:
    Plasmodium has roughly thirty predicted NMT substrates. The list validates itself -- GAP45,
    ARO, CDPK1, Rab-5B, ARF1 and ISP3 are the canonical N-myristoylated families in
    apicomplexans and all are present.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_myristoylome.py
"""
from _common import run

KEY = "pf_myristoylome"

if __name__ == "__main__":
    run(KEY)
