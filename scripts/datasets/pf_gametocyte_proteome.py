#!/usr/bin/env python3
"""Mature gametocyte proteome and translatome

What a stage V gametocyte contains, and which proteins it is still making

    level / kind : translation / proteomics
    provides     : gametocyte_proteome_log2, gametocyte_newly_made
    coverage     : 2,544 proteins (44%)
    citation     : Alves E et al., The translatome of quiescent Plasmodium falciparum gametocytes reveals parasite pyridoxal kinase as a target. bioRxiv 2026, doi:10.64898/2026.03.24.713170 (preprint)
    accession    : PXD075878 / bioRxiv 10.64898/2026.03.24.713170 Extended Data Tables 1, 3
    url          : https://www.biorxiv.org/content/10.64898/2026.03.24.713170v1.supplementary-material
    local path   : datasets/translation/proteomics/pf_gametocyte_2026/Extended Datasets.xlsx

Quirks that cost time once:
    The transmission stage, and the first protein abundance in the map for any stage other than
    the blood stage. Newly made means found with the click-chemistry label and in no control --
    the deposit's group 4, which is the paper's 705 proteins; taking the whole sheet would call
    1,179 newly made, because the other groups also appear without the label or when synthesis
    is blocked. Shipped as MEMBERSHIP rather than an enrichment on purpose: there is one pooled
    sample per condition and most labelled proteins are absent from the controls, so a ratio
    would be a ratio to nothing. A protein the label found but the abundance run did not
    quantify keeps its flag and has no abundance.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_gametocyte_proteome.py
"""
from _common import run

KEY = "pf_gametocyte_proteome"

if __name__ == "__main__":
    run(KEY)
