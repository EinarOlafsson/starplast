#!/usr/bin/env python3
"""StarPath crosslink MS

Measured physical proximity; residue-level crosslinks and Chai-1 complexes

    level / kind : post_translation / XLMS
    provides     : n_xlink_partners, best_model_agreement
    coverage     : 2,842 pairs / 1,630 genes
    citation     : Mapping a Toxoplasma gondii interactome by crosslinking mass spectrometry and machine learning (2025)
    PMID         : 40874616
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12505969/supplementaryFiles
    local path   : starpath_crosslinks.json, starpath_dump/cifs/

Quirks that cost time once:
    RH88 accessions do NOT map to ME49 by suffix; use the alias column. 60% of predicted
    complexes place no crosslink within reach; only 162 pairs are trustworthy.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/starpath_xlms.py
"""
from _common import run

KEY = "starpath_xlms"

if __name__ == "__main__":
    run(KEY)
