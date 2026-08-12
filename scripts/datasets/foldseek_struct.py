#!/usr/bin/env python3
"""Foldseek structural similarity

TM-align over Toxoplasma AlphaFold models, TM >= 0.7

    level / kind : post_translation / structure
    provides     : n_struct_similar
    coverage     : 11,684 pairs / 2,338 genes

Quirks that cost time once:
    Needs no orthology, so it reaches lineage-specific effectors homology edges cannot.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/foldseek_struct.py
"""
from _common import run

KEY = "foldseek_struct"

if __name__ == "__main__":
    run(KEY)
