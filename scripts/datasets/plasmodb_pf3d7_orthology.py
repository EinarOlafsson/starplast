#!/usr/bin/env python3
"""Plasmodium falciparum 3D7 orthology and paralogy

OrthoMCL group, ortholog and paralog counts from the same PlasmoDB report

    level / kind : reference / orthology
    provides     : ortholog_number, orthogroup, paralog_number, has_paralog
    coverage     : 5,720 P. falciparum genes
    accession    : PlasmoDB GenesByTaxon attributesTabular
    url          : https://plasmodb.org/plasmo/service/record-types/transcript/searches/GenesByTaxon/reports/attributesTabular
    local path   : starplast/data/pf_nodes.parquet

Quirks that cost time once:
    Split from the attribute report on 2026-09-25 (see plasmodb_pf3d7_attributes): orthology is
    computed from sequence across species, not measured on this gene, and grouping it with the
    piggyBac screen closed one over the other for nothing.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/plasmodb_pf3d7_orthology.py
"""
from _common import run

KEY = "plasmodb_pf3d7_orthology"

if __name__ == "__main__":
    run(KEY)
