#!/usr/bin/env python3
"""Plasmodium lysine lactylome (resolved from NF54)

Lactylated lysines per gene, reported against NF54 and resolved to 3D7

    level / kind : post_translation / lactylation
    provides     : n_lactylsites, has_lactyl
    coverage     : 144 genes
    PMID         : 41417877
    accession    : PLoS Genet 1011991 S1
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12742760/supplementaryFiles
    local path   : datasets/post_translation/lactylome/41417877/pgen.1011991.s014.xlsx

Quirks that cost time once:
    Reported against the NF54 annotation, not 3D7, so joining on the accession string would have
    dropped all 186 genes without a word -- the exact failure the Toxoplasma identity layer
    exists to prevent, met here for the first time on this arm. `plasmodium.strain_map` resolves
    it through orthogroups holding exactly one gene on each side, which drops the paralogous
    surface families rather than guessing which member a measurement belongs to. Orthology is a
    claim about ancestry, so it is CHECKED against one about identity: 96.8% of the 4,310 pairs
    have exactly the same protein length and 99.2% are within 5%, which is what it should look
    like when one line was cloned from the other, and a test fails if a future release breaks
    it. 144 of 186 genes resolve; the rest are in multi-gene groups. Site counts use a 0.75
    localisation cut and the flag does not, the same split as acetylation.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_lactylome.py
"""
from _common import run

KEY = "pf_lactylome"

if __name__ == "__main__":
    run(KEY)
