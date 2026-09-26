#!/usr/bin/env python3
"""Plasmodium falciparum strain variation

SNP counts across sequenced strains, from the same PlasmoDB report

    level / kind : reference / variation
    provides     : snp_total_all_strains, snp_nonsynonymous, snp_synonymous, snp_noncoding, snp_stop_codon, snp_nonsyn_syn_ratio
    coverage     : 5,720 P. falciparum genes
    accession    : PlasmoDB GenesByTaxon attributesTabular
    url          : https://plasmodb.org/plasmo/service/record-types/transcript/searches/GenesByTaxon/reports/attributesTabular
    local path   : starplast/data/pf_nodes.parquet

Quirks that cost time once:
    Split from the attribute report on 2026-09-25 (see plasmodb_pf3d7_attributes): population
    sequencing is its own experiment, and its columns are the only ones of the report that are
    associated with each other.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/plasmodb_pf3d7_snps.py
"""
from _common import run

KEY = "plasmodb_pf3d7_snps"

if __name__ == "__main__":
    run(KEY)
