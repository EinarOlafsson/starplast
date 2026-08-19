#!/usr/bin/env python3
"""Plasmodium antibody epitopes (IEDB)

Distinct antibody epitope sequences per gene

    level / kind : reference / immunity
    provides     : n_bcell_epitopes
    coverage     : 434 antigens, 7,366 distinct epitopes
    accession    : IEDB bcell_search
    url          : https://query-api.iedb.org/bcell_search?parent_source_antigen_source_org_name=ilike.*Plasmodium%20falciparum*
    local path   : datasets/reference/plasmodb/iedb_pf_bcell_epitopes.tsv

Quirks that cost time once:
    DISTINCT sequences, not assay records. MSP1 alone carries 1,739 epitopes out of 14,610
    records, so counting records would rank antigens by how many groups have studied them rather
    than by how much of the protein antibodies recognise. Reached through UniProt rather than
    through product descriptions -- the Toxoplasma arm has to match descriptions because IEDB's
    Toxoplasma antigen names are verbatim ToxoDB text, while the falciparum names carry the
    accession, which is a better key. The 209 PlasmoDB accessions naming more than one gene are
    dropped: an epitope belongs to a protein, and attaching it to whichever paralogue sorted
    first would be inventing the answer. 434 of 444 antigens resolve. Absent is absent and not
    zero, because IEDB records what somebody tested. Validated on the history of the field: MSP1
    is the top antigen and CSP, the RTS,S vaccine antigen, is present.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_iedb_bcell.py
"""
from _common import run

KEY = "pf_iedb_bcell"

if __name__ == "__main__":
    run(KEY)
