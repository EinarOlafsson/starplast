#!/usr/bin/env python3
"""PubMed Central open-access full texts

Sectioned JATS XML

    level / kind : reference / literature
    provides     : n_fulltext
    coverage     : 6,667 articles
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML
    local path   : /mnt/wd4tb/skill_corpora/toxoplasma-scientist/

Quirks that cost time once:
    A biased subset: only what publishers deposited open access. ASSEMBLED HERE: for each PMID
    in the pubmed corpus that has a PMCID, fetch the JATS from the URL above. Machine-local by
    size, which is why the built cache is what ships.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pmc_oa.py
"""
from _common import run

KEY = "pmc_oa"

if __name__ == "__main__":
    run(KEY)
