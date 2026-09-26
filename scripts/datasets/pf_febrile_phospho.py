#!/usr/bin/env python3
"""Phosphorylation under febrile heat stress

How far each protein's phosphorylation moves at 39 degrees, and at how many sites

    level / kind : post_translation / proteomics
    provides     : febrile_phospho_max_abs_log2fc, febrile_phospho_n_sites_up, febrile_phospho_n_sites_down
    coverage     : 1,874 proteins (33%)
    citation     : Jones D et al., Physiological febrile heat stress increases cytoadhesion through increased protein trafficking. Elife 2026;14:RP107860
    PMID         : 42126964
    accession    : PXD073843 / eLife Supplementary File 1
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13171106/supplementaryFiles
    local path   : datasets/post_translation/phosphosites/pf_febrile_2026/elife-107860-supp1-v1.xlsx

Quirks that cost time once:
    A fever is a condition this parasite actually meets, which is what makes this a different
    question from the phosphoproteome already in the map. THIRTY-NINE degrees, not forty: the
    methods give the temperature the cultures were held at, and the preprint's abstract does
    not. The paper's 143 up and 53 down reproduce exactly, but those count ROWS, and a row is a
    site at one phosphorylation multiplicity; collapsed to unique sites, each taking its
    strongest change, it is 128 and 51, which is what these columns count. Specific rather than
    global: 60% of proteins with a rising site are exported to the host cell against 7% of all
    quantified proteins. A site whose peptide could belong to more than one gene is dropped
    rather than assigned to the first.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_febrile_phospho.py
"""
from _common import run

KEY = "pf_febrile_phospho"

if __name__ == "__main__":
    run(KEY)
