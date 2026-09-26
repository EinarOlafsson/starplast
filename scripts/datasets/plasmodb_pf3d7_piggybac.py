#!/usr/bin/env python3
"""Plasmodium falciparum piggyBac saturation mutagenesis

Mutagenesis index and fitness score from a genome-saturating transposon screen

    level / kind : DNA / insertion screen
    provides     : piggybac_mis, piggybac_mfs
    coverage     : 5,720 P. falciparum genes
    citation     : Zhang M et al. Uncovering the essential genes of the human malaria parasite Plasmodium falciparum by saturation mutagenesis. Science 2018
    PMID         : 29724925
    accession    : PlasmoDB GenesByTaxon attributesTabular
    url          : https://plasmodb.org/plasmo/service/record-types/transcript/searches/GenesByTaxon/reports/attributesTabular
    local path   : starplast/data/pf_nodes.parquet

Quirks that cost time once:
    PMID resolved by E-utilities on 2026-09-25 (esearch 'Plasmodium falciparum saturation
    mutagenesis piggyBac essential genes 2018' returns only 29724925; esummary confirms the
    title and Zhang M as first author). Served through the PlasmoDB attribute report and split
    from it the same day. The direction is the trap: a LOW mutagenesis index means the gene
    resists disruption and is therefore essential, and inverting it would swap the essential and
    dispensable genomes without crashing, so the test checks it against biology -- ribosomal
    proteins come out at median MIS 0.15 and the var, rifin and stevor families at 0.94.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/plasmodb_pf3d7_piggybac.py
"""
from _common import run

KEY = "plasmodb_pf3d7_piggybac"

if __name__ == "__main__":
    run(KEY)
