#!/usr/bin/env python3
"""AlphaFold DB

Per-gene mean pLDDT; coordinates fetched on demand

    level / kind : reference / structure
    provides     : mean_plddt
    coverage     : 6,480 (79.6%)
    citation     : Jumper et al. 2021, AlphaFold Protein Structure Database
    accession    : AlphaFold DB

Quirks that cost time once:
    Missing for the largest proteins, which here are disproportionately secreted effectors.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/alphafold.py
"""
from _common import run

KEY = "alphafold"

if __name__ == "__main__":
    run(KEY)
