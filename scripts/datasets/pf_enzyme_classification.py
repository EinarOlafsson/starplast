#!/usr/bin/env python3
"""Plasmodium enzyme classification (PlasmoDB)

EC numbers reported directly and inferred from orthologs, kept apart

    level / kind : reference / annotation
    provides     : ec_number, has_ec, ec_number_orthology
    coverage     : Declared coverage: 1,220 direct annotations; orthology-derived field separate
    accession    : PlasmoDB ec_numbers and ec_numbers_derived
    url          : https://plasmodb.org/plasmo/service/record-types/transcript/searches/GenesByTaxon/reports/attributesTabular
    local path   : datasets/reference/plasmodb/plasmodb_pf3d7_ec.tsv

Quirks that cost time once:
    PlasmoDB reports direct EC annotations separately from EC numbers derived from orthologs.
    Gene-level curation/prediction provenance and the assignment release remain unresolved for
    the installed table. `has_ec` reports a direct assignment; experimental enzyme activity has
    not been independently verified. The slot selects the direct field by default. Orthology-
    derived annotations require an explicit selection and retain their transfer lineage.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_enzyme_classification.py
"""
from _common import run

KEY = "pf_enzyme_classification"

if __name__ == "__main__":
    run(KEY)
