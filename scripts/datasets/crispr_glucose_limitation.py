#!/usr/bin/env python3
"""CRISPR screen without glucose or glutamine

Fitness on glucose alone, on glutamine alone, and which carbon source a gene needs

    level / kind : DNA / CRISPR_screen
    provides     : fit_complete_medium_2025, fit_no_glutamine, fit_no_glucose, fit_glucose_dependence, fit_glucose_dependence_fdr
    coverage     : 7,393 (90.8%)
    citation     : Uboldi AD et al., Differentiation of Toxoplasma into latent forms is linked to central carbon metabolism and requires a GID/CTLH-type E3 ligase. bioRxiv 2025, doi:10.1101/2025.07.27.667068 (preprint)
    accession    : bioRxiv 10.1101/2025.07.27.667068 Table S1
    url          : https://www.biorxiv.org/content/10.1101/2025.07.27.667068v1.supplementary-material
    local path   : datasets/DNA/CRISPR_screen/glucose_limitation_2025/TableS1.xlsx

Quirks that cost time once:
    An eighth genome-wide knockout screen, and the first that changes the carbon source. Arms
    are measured against the post-selection library, so they are not directly comparable with
    the input-referenced screens; the dependence column is the authors' own contrast. Verified
    against the paper's biology: GDH1 ranks 1 of 8,155 and PEPCK 2 among genes needed without
    glucose, with BFD1 and all four GID subunits on the other side. The differential is NOISY --
    replicates agree at rho 0.10-0.17 and only 175 genes reach FDR 0.05 -- so the FDR ships
    beside it and the column is a screen, not a measurement of one gene. A preprint: cite the
    journal version once it exists.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/crispr_glucose_limitation.py
"""
from _common import run

KEY = "crispr_glucose_limitation"

if __name__ == "__main__":
    run(KEY)
