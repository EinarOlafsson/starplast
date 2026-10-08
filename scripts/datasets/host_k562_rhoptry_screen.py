#!/usr/bin/env python3
"""Host genes required for rhoptry discharge

Whether knocking out a human gene stops Toxoplasma discharging its rhoptries

    level / kind : reference / CRISPR_screen
    provides     : rhoptry_discharge_score, rhoptry_discharge_beta, rhoptry_discharge_fdr
    coverage     : 18,700 unambiguous protein projections; all 20,010 gene scores retained; 39 ambiguous symbols withheld
    citation     : Valleau D et al., Clustering of host N-glycans by the microneme MIC1/4/6 complex licenses Toxoplasma rhoptry discharge. EMBO J 2026, doi:10.1038/s44318-026-00911-z
    PMID         : 42791346
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
    ties. Keyed by reviewed UniProt through gene symbol. The 2026 journal Dataset EV1 matches
    all 20,010 source genes and all three score fields exactly; the original v2 source file and
    all gene scores are retained. First-publication lineage is bioRxiv
    doi:10.1101/2025.10.16.682961 (2025), not a new experiment in 2026. Mapping correction:
    18,700 unambiguous protein projections; 39 ambiguous symbols are withheld and 1,271 unmapped
    symbols remain gene evidence. All 20,010 original gene rows, including mapping alternatives,
    are preserved in deposit_host_k562_rhoptry_gene_evidence.parquet. The source-verified ledger
    host_projection_withdrawals.json records the 117 corrected protein cells and separate
    source/legacy decimal encodings. Withheld does not mean absent or a negative experiment. See
    results/host_symbol_mapping_2026_10_08/.

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
