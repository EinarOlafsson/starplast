#!/usr/bin/env python3
"""Antibody epitopes (IEDB)

Distinct antibody epitope sequences per gene

    level / kind : reference / immunity
    provides     : n_bcell_epitopes
    coverage     : 34 genes, 222 distinct epitopes
    accession    : IEDB bcell_search
    url          : https://query-api.iedb.org/bcell_search?parent_source_antigen_source_org_name=ilike.*Toxoplasma*
    local path   : starplast/data/iedb_bcell_epitopes.tsv

Quirks that cost time once:
    From IEDB directly, NOT through ToxoDB, because ToxoDB's epitope integration is not split by
    type and antigenicity is a question about antibodies. 189 of the 221 genes in the ToxoDB
    epitope column have no antibody record at all, so the two are genuinely different
    measurements filling different slots. Distinct epitope SEQUENCES and not assay records: a
    protein studied by twenty groups accumulates twenty records for one peptide, and counting
    records would rank antigens by fashion. Verified by what tops it -- SRS29B (SAG1) 47, GRA6
    29, GRA1 26, GRA4 20, GRA7 16, MIC3 11, which is the panel commercial Toxoplasma
    serodiagnostic kits are built from. Antigens are mapped by their product description,
    because IEDB names them verbatim from ToxoDB; the trailing-symbol route resolves ten fewer
    and loses SRS29B, the most studied antigen in the organism.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/iedb_bcell.py
"""
from _common import run

KEY = "iedb_bcell"

if __name__ == "__main__":
    run(KEY)
