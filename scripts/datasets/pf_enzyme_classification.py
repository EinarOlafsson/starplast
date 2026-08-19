#!/usr/bin/env python3
"""Plasmodium enzyme classification (PlasmoDB)

EC number per gene, curated and orthology-derived kept apart

    level / kind : reference / annotation
    provides     : ec_number, has_ec, ec_number_orthology
    coverage     : 1,220 curated, 335 more from orthology
    accession    : PlasmoDB ec_numbers and ec_numbers_derived
    url          : https://plasmodb.org/plasmo/service/record-types/transcript/searches/GenesByTaxon/reports/attributesTabular
    local path   : datasets/reference/plasmodb/plasmodb_pf3d7_ec.tsv

Quirks that cost time once:
    PlasmoDB serves two EC fields and they are different KINDS of evidence -- one curated for
    this organism, one inferred from the gene's OrthoMCL group -- so they are separate columns
    and `has_ec` counts only the curated one. Merged they would be 1,584 genes with no way to
    tell which 335 were never annotated here at all, which is inference standing where
    annotation should. The slot lists the curated column first and its policy is `one`, so the
    leading candidate wins and the derived field is there to be chosen deliberately rather than
    by default.

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
