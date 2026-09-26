#!/usr/bin/env python3
"""Serum-restriction CRISPR screens (10% vs 1% FBS)

Fitness in lipid-rich and lipid-limited medium, and the lipid-dependence difference

    level / kind : DNA / CRISPR_screen
    provides     : fit_lipid_rich_p8, fit_lipid_limited_p8, fit_lipid_rich_p4p5, fit_lipid_limited_p4p5, fit_serum_differential_p8, fit_serum_differential_p4p5
    coverage     : 7,395 (90.8%)
    citation     : Bitew MA et al., A genome-wide CRISPR screen identifies GRA38 as a key regulator of lipid homeostasis during Toxoplasma gondii adaptation to lipid-rich conditions. Nat Commun 2025;16:11177
    PMID         : 41407671
    accession    : Nat Commun 2025 Supplementary Data 2 (MOESM4)
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12711892/supplementaryFiles
    local path   : datasets/DNA/CRISPR_screen/41407671/41467_2025_66137_MOESM4_ESM.xlsx

Quirks that cost time once:
    Two independent genome-wide screens in HFF (Exp1 at passage 8; Exp2 at passages 4, 5, 8),
    each in 10% and 1% serum. Replicate screens of one passage are averaged. Sign as published
    and verified on the data: negative fitness = depleted = needed (ribosomal proteins median
    -7.8 against -4.9, AUC 0.79-0.84); a NEGATIVE differential means needed in lipid-RICH medium
    -- GRA38 (TGGT1_312420), the paper's gene, is -6.66 in 10% and +0.33 in 1%. Fitness in
    either serum IS fibroblast fitness again (rho 0.79-0.85 with fit_invitro_hff) and is grouped
    with it for leakage; the differential is orthogonal to it (rho 0.06) and is the new axis.
    The workbook's lipidomics sheets are not gene-level and are not used.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/crispr_serum_restriction.py
"""
from _common import run

KEY = "crispr_serum_restriction"

if __name__ == "__main__":
    run(KEY)
