#!/usr/bin/env python3
"""Plasmodium falciparum abstract corpus (COMPUTED layer)

Who is named in the malaria literature, how deeply, and which genes appear together

    level / kind : reference / literature
    provides     : n_publications, n_papers_focal, n_papers_substantive, attention_depth, edge:comention
    coverage     : 43,482 abstracts; 732 genes named in the first 10,000
    accession    : PubMed
    url          : https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi
    local path   : .claude/skills/plasmodium-scientist/corpus/pubmed_plasmodium.jsonl

Quirks that cost time once:
    Fetched by `scripts/fetch_pubmed_corpus.py`, which slices the query by YEAR because NCBI
    stops at ten thousand twice over: esearch will not page past it and efetch answers 400 for a
    retstart beyond it. Both limits are silent -- the first version of that script fetched 9,989
    abstracts of 43,482 and reported success. The scan itself is the Toxoplasma arm's,
    unchanged: `identity` for who is named, `corpus` for what a document is, `literature` for
    the counting and the attention correction, so the two arms' attention numbers mean the same
    thing. What is organism-specific is the accession shapes -- this literature cites `PF3D7_`,
    `PFA0110w`, `PF13_0222` and `MAL1P4.01` in the same paragraph -- and the `Pf` symbol prefix,
    both of which are now arguments to `identity.build_index` rather than constants in it. Hard-
    coded to Toxoplasma they registered 9 of 9,106 previous accessions and the corpus read as
    one that never mentions a gene. Abstracts only: there is no `incidental` tier, since that
    means a mention in a body or a caption, and no full-text corpus is loaded for this arm.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_literature.py
"""
from _common import run

KEY = "pf_literature"

if __name__ == "__main__":
    run(KEY)
