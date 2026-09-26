#!/usr/bin/env python3
"""Chromatin-state proximity proteomes

How enriched each protein is near heterochromatin, active marks and the centromere

    level / kind : post_translation / BioID
    provides     : chromprox_hp1_log2fc, chromprox_hp1_hit, chromprox_h3k27ac_log2fc, chromprox_h3k27ac_hit, chromprox_h3k4me3_log2fc, chromprox_h3k4me3_hit, chromprox_centromere_log2fc, chromprox_centromere_hit
    coverage     : 2,020 proteins (35%)
    citation     : Ramon-Zamorano G et al., Protein landscape of the chromatin states in the malaria parasite Plasmodium falciparum. bioRxiv 2025, doi:10.1101/2025.09.23.678001 (preprint)
    accession    : bioRxiv 10.1101/2025.09.23.678001 Tables S1, S5, S7
    url          : https://www.biorxiv.org/content/10.1101/2025.09.23.678001v1.supplementary-material
    local path   : datasets/post_translation/BioID/pf_chromatin_proxiome_2026/TableS1_HP1_proximity_proteome_media-3.xlsx

Quirks that cost time once:
    Two empty slots at once -- chromatin state, and proximity labelling -- for an organism whose
    gene regulation is largely chromatin. Four of the seven baits ship. The thresholds are the
    paper's own and differ per bait, because the baits differ in signal; one threshold for all
    of them would be tidier and wrong. The two euchromatin counts reproduce exactly (H3K27ac 99,
    H3K4me3 48); the HP1 column gives 91, which is NOT the paper's 61 -- that is its combined
    heterochromatin group over two HP1 baits and this is one of them, so they are not the same
    quantity. Controls land correctly: HP1, HDA2, GDV1 and AP2-HC in heterochromatin, CENH3 at
    the centromere, BDP1 and GCN5 in active chromatin. SIR2A is in none, because no bait
    detected it -- absence of evidence. A protein group spanning two genes is dropped, not
    assigned to the first.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_chromatin_proxiome.py
"""
from _common import run

KEY = "pf_chromatin_proxiome"

if __name__ == "__main__":
    run(KEY)
