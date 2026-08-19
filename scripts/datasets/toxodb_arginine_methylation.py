#!/usr/bin/env python3
"""Monomethylarginine proteome (via ToxoDB)

Monomethylarginine sites reported per gene

    level / kind : post_translation / proteomics
    provides     : n_arginine_methylation_sites
    coverage     : 368 genes
    accession    : ToxoDB Yakubu monomethylarginine
    url          : https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByPTM/reports/attributesTabular
    local path   : starplast/data/toxodb_arginine_methylation.tsv

Quirks that cost time once:
    Yakubu et al.'s RH proteomics, through ToxoDB's PTM search. Verified by substrate class
    rather than by a metadata field: RNA-binding, RRM and helicase proteins are enriched
    2.9-fold among the methylated (Fisher p = 0.002) and transporters and membrane proteins are
    depleted at odds 0.43. That is the PRMT substrate profile -- RG and RGG motifs sit in RNA-
    binding proteins. The slot it fills did not exist before: arginine methylation has its own
    writers, its own substrate class and its own deposit, and its absence from the catalog was a
    gap rather than a lack of data.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/toxodb_arginine_methylation.py
"""
from _common import run

KEY = "toxodb_arginine_methylation"

if __name__ == "__main__":
    run(KEY)
