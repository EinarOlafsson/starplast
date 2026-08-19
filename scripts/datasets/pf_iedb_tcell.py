#!/usr/bin/env python3
"""Plasmodium T-cell epitopes (IEDB)

Distinct T-cell epitope sequences per gene

    level / kind : reference / immunity
    provides     : n_tcell_epitopes
    coverage     : 44 antigens, 1,542 distinct epitopes
    accession    : IEDB tcell_search
    url          : https://query-api.iedb.org/tcell_search?parent_source_antigen_source_org_name=ilike.*Plasmodium%20falciparum*
    local path   : datasets/reference/plasmodb/iedb_pf_tcell_epitopes.tsv

Quirks that cost time once:
    Kept apart from the antibody half on purpose, and the numbers show why they are not
    interchangeable: 434 antigens carry an antibody epitope and only 44 carry a T-cell one.
    Pooling them, or filling either slot with the other, would answer one question with the
    other's number. Same UniProt keying and same distinct-sequence counting as the antibody
    table, and the loader reads whichever halves are present so one fetch failing does not cost
    the other column.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_iedb_tcell.py
"""
from _common import run

KEY = "pf_iedb_tcell"

if __name__ == "__main__":
    run(KEY)
