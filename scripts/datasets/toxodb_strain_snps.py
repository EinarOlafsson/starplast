#!/usr/bin/env python3
"""Strain variation (ToxoDB HTS SNPs)

SNPs per gene across every sequenced strain, split by effect

    level / kind : reference / variation
    provides     : snp_total_all_strains, snp_nonsynonymous, snp_synonymous, snp_noncoding, snp_stop_codon
    coverage     : 8,140 genes (100%)
    accession    : ToxoDB ME49
    url          : https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByTaxon/reports/attributesTabular
    local path   : starplast/data/toxodb_strain_snps.tsv

Quirks that cost time once:
    Retrieved 2026-08-16 through the same REST report as the identity table, asking for the five
    gene_hts_*_snps attributes. Verified against known biology rather than against a metadata
    field: nonsynonymous SNPs per kb come out at 112 for the SRS surface antigens, 55 for the
    ROP5/ROP18/GRA15 virulence loci, 30 across all genes and 2.9 for ribosomal proteins. That
    ordering -- what immunity sees, then the strain-typing markers, then the conserved core --
    is the check. Zero is a measurement here, not a gap: 690 genes carry no SNP in any sequenced
    strain.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/toxodb_strain_snps.py
"""
from _common import run

KEY = "toxodb_strain_snps"

if __name__ == "__main__":
    run(KEY)
