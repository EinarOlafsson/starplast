#!/usr/bin/env python3
"""AP2XII-1/AP2XI-2 perturbation total proteome

Replicate abundance and log2 fold-change during pre-sexual conversion

    level / kind : translation / proteomics
    provides     : proteome_log2_normalized_and_imputated_abundances_UT_R1, proteome_log2_normalized_and_imputated_abundances_UT_R2, proteome_log2_normalized_and_imputated_abundances_UT_R3, proteome_log2_normalized_and_imputated_abundances_T_24h_R1, proteome_log2_normalized_and_imputated_abundances_T_24h_R2, proteome_log2_normalized_and_imputated_abundances_T_24h_R3, proteome_log2_normalized_and_imputated_abundances_T_32h_R1, proteome_log2_normalized_and_imputated_abundances_T_32h_R2, proteome_log2_normalized_and_imputated_abundances_T_32h_R3, proteome_log2_normalized_and_imputated_abundances_T_48h_R1, proteome_log2_normalized_and_imputated_abundances_T_48h_R2, proteome_log2_normalized_and_imputated_abundances_T_48h_R3, proteome_lfc_UT_Vs_T_24h_log2_fold_change, proteome_lfc_UT_Vs_T_32h_log2_fold_change_1, proteome_lfc_UT_Vs_T_48h_log2_fold_change_2
    coverage     : 3,005 (36.9%)
    citation     : Antunes AV et al., In vitro production of cat-restricted Toxoplasma pre-sexual stages. Nature 2024;625:366-376
    PMID         : 38093015
    accession    : PXD039400, PXD042658
    url          : https://www.ebi.ac.uk/pride/archive/projects/PXD039400
    local path   : toxo_stage_atlas/data/proteomics/PXD039400_PXD042658_supp_SupplTable3_total_proteome.xlsx

Quirks that cost time once:
    An actual total proteome, not IP enrichment. The abundance and fold-change columns are
    distinct quantification types and expression.total_proteome normalizes them separately
    before they enter the cache.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.expression.load_all()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/proteome_total.py
"""
from _common import run

KEY = "proteome_total"

if __name__ == "__main__":
    run(KEY, normalized_by='expression.load_all()')
