#!/usr/bin/env python3
"""Antimalarial target engagement (thermal profiling)

How many of 25 antimalarials measurably engage each protein, and how many tested it

    level / kind : post_translation / proteomics
    provides     : engaged_n_compounds_tested, engaged_n_compounds_hit
    coverage     : 3,126 proteins (55%)
    citation     : Pazicky S et al., Thermal proteome profiling identifies new drug targets in Plasmodium falciparum parasites. bioRxiv 2026, doi:10.64898/2026.01.30.702724 (preprint)
    accession    : PXD048737-PXD048772 / bioRxiv 10.64898/2026.01.30.702724 Table S3
    url          : https://www.biorxiv.org/content/10.64898/2026.01.30.702724v1.supplementary-material
    local path   : datasets/post_translation/thermal/pf_tpp_2026/SupplTable3_media-6.csv

Quirks that cost time once:
    COUNTS, not the continuous response, and the reason is the reason: the continuous response
    agrees between replicates at r 0.0-0.55, so a per-compound number would be mostly noise,
    while the hit call reproduces the paper (99 stabilised hits; cladosporine-KRS,
    MMV665915-ACS10 and KAF156-prohibitin all in the top thirteen; PfATP4 with cipargamin; DHODH
    with the DSM compounds). The denominator ships with the count, because a protein hit twice
    out of five assays is not the evidence of one hit twice out of 25. NOT a melting
    temperature: this design holds temperature fixed and varies dose, so it cannot give one --
    the melting points are the separate MAP-X meltome, a different experiment from the same
    laboratory.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_target_engagement.py
"""
from _common import run

KEY = "pf_target_engagement"

if __name__ == "__main__":
    run(KEY)
