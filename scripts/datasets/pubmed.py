#!/usr/bin/env python3
"""PubMed abstracts

Titles and abstracts for co-mention and attention

    level / kind : reference / literature
    provides     : n_publications, attention_depth
    coverage     : 33,924 records
    url          : https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=Toxoplasma&retmax=100000
    local path   : .claude/skills/toxoplasma-scientist/corpus/pubmed_toxoplasma.jsonl

Quirks that cost time once:
    ASSEMBLED HERE from an E-utilities query rather than downloaded as a file: the esearch above
    returns the PMID set, and efetch retrieves each record. The count grows over time, so a
    rebuild will not reproduce 33,924 exactly -- record the date. This corpus was built 2026.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pubmed.py
"""
from _common import run

KEY = "pubmed"

if __name__ == "__main__":
    run(KEY)
