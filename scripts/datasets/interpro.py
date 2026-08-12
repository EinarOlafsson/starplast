#!/usr/bin/env python3
"""InterPro domains

Domain identity and count

    level / kind : reference / domains
    provides     : n_interpro, interpro_id, interpro_desc, pfam_id
    coverage     : 8,140
    url          : https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByTaxon/reports/attributesTabular?organism=%5B%22Toxoplasma%20gondii%20ME49%22%5D&reportConfig=%7B%22attributes%22%3A%5B%22gene_source_id%22%2C%22interpro_id%22%2C%22interpro_description%22%2C%22pfam_id%22%2C%22pfam_description%22%5D%2C%22includeHeader%22%3Atrue%2C%22attachmentType%22%3A%22plain%22%7D
    local path   : datasets/interpro_tgon.csv

Quirks that cost time once:
    Domain annotation partly records study effort, not conserved architecture.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/interpro.py
"""
from _common import run

KEY = "interpro"

if __name__ == "__main__":
    run(KEY)
