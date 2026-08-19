#!/usr/bin/env python3
"""Plasmodium lysine acetylome

Acetylated lysines per gene, and whether the gene was seen acetylated at all

    level / kind : post_translation / acetylation
    provides     : n_acetylsites, has_acetyl
    coverage     : 1,145 genes, 2,163 localised sites
    PMID         : 26813983
    accession    : Sci Rep 19722 S2
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4728587/supplementaryFiles
    local path   : datasets/post_translation/acetylome/26813983/srep19722-s2.xls

Quirks that cost time once:
    Two columns built to two standards, because identifying an acetylated PEPTIDE and localising
    the acetyl group to a particular lysine are different claims. The FLAG uses every
    identification; the COUNT uses only sites with an Ascore of 0.75 or better, since a site
    count is meaningless if you do not know which lysine. The list is titled Final and is not
    pre-filtered on localisation -- Ascores run down to 0 -- so taking its length as a site
    count would have been wrong by about a quarter. Self-validating: fourteen histones appear,
    and the most heavily acetylated proteins are the PHD finger proteins, the MYST
    acetyltransferase and the coactivator ADA2, which is to say the acetylation machinery
    itself.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_acetylome.py
"""
from _common import run

KEY = "pf_acetylome"

if __name__ == "__main__":
    run(KEY)
