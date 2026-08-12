#!/usr/bin/env python3
"""InterPro domains

Domain identity and count

    level / kind : reference / domains
    provides     : n_interpro, interpro_id, interpro_desc, pfam_id
    coverage     : 8,140
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
