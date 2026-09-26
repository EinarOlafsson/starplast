#!/usr/bin/env python3
"""Host genes required for rhoptry discharge

Whether knocking out a human gene stops Toxoplasma discharging its rhoptries

    level / kind : reference / CRISPR_screen
    provides     : rhoptry_discharge_score, rhoptry_discharge_beta, rhoptry_discharge_fdr
    coverage     : 18,739 human genes
    citation     : Valleau D et al., Clustering of host N-glycans licenses Toxoplasma rhoptry discharge. bioRxiv 2025, doi:10.1101/2025.10.16.682961 (preprint)
    accession    : bioRxiv 10.1101/2025.10.16.682961 v2 Table S1
    url          : https://www.biorxiv.org/content/10.1101/2025.10.16.682961v2.supplementary-material
    local path   : datasets/host/k562_rhoptry_screen/v2_media-1.xlsx

Quirks that cost time once:
    The first host-gene screen in the map, and it opens the host side of invasion: positive
    means the knockout blocks discharge, so the gene is required for it. Recovers its own
    biology -- SLC35A2, the transporter the paper is about, ranks 1 of 20,010, and a twenty-gene
    N-glycan panel has median rank 135 with 85% in the top tenth, while B4GALT1, FUT8 and
    ST6GAL1, which the paper argues are not involved, are not hits. The Wald FDR is shipped
    rather than the permutation one, which is quantised into a few values and would read as
    ties. Keyed by reviewed UniProt through gene symbol; 155 microRNA loci and other
    unresolvable symbols are dropped. A preprint.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_k562_rhoptry_screen.py
"""
from _common import run

KEY = "host_k562_rhoptry_screen"

if __name__ == "__main__":
    run(KEY)
