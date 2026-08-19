#!/usr/bin/env python3
"""Secreted-fraction partition

How a secreted protein splits between the soluble and vesicular fractions

    level / kind : post_translation / proteomics
    provides     : secretome_soluble_over_vesicle_log2
    coverage     : 165 proteins
    PMID         : 40874616
    accession    : ToxoDB Ramirez-Flores vesicles
    url          : https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByProteomicstgonGT1_quantitativeMassSpec_Ramirez_Flores_Vesicles_RSRC/reports/attributesTabular
    local path   : starplast/data/secretome_partition.tsv

Quirks that cost time once:
    REFUSED earlier in the same session and then accepted, which is worth stating rather than
    reversing quietly. It was first asked which proteins are enriched in secreted VESICLES, and
    the answer was incoherent -- dense granule proteins depleted, ribosomal proteins the most
    enriched thing in it. The two fractions it compares are BOTH secreted material, so the
    question it can answer is how a secreted protein partitions between them, and asked that way
    it behaves: micronemes, which dominate classical excretory-secretory antigen preparations,
    at +3.72, and the GPI-anchored surface antigens at -0.68. Its limit is that there is no
    negative list -- 171 proteins were seen in secreted material and nothing says what was
    looked for and missed, so absence is not evidence.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/secretome_partition.py
"""
from _common import run

KEY = "secretome_partition"

if __name__ == "__main__":
    run(KEY)
