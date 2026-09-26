#!/usr/bin/env python3
"""Bradyzoite subtypes in the mouse brain (single-cell)

Average expression in each of five subtypes of in vivo bradyzoite

    level / kind : transcription / scRNAseq
    provides     : bzsub_A_expr, bzsub_B_expr, bzsub_C_expr, bzsub_D_expr, bzsub_E_expr
    coverage     : 7,739 (95%)
    citation     : Ulu A et al., Bradyzoite subtypes rule the crossroads of Toxoplasma development. Nat Commun 2026;17:1783
    PMID         : 41580398
    accession    : GSE311669 / Supplementary Data 1
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12917143/supplementaryFiles
    local path   : datasets/transcription/scRNAseq/41580398/41467_2026_68489_MOESM3_ESM.xlsx

Quirks that cost time once:
    Five columns because a tissue cyst is not one transcriptional state, which is the paper's
    finding. All five are bradyzoites -- BAG1, LDH2, ENO1, SRS9 and CST1 above the 92nd
    percentile in every group, SAG1 below the 41st -- and Group B carries SRS22A at 3.27 against
    0.34-0.66 elsewhere, the subtype signature. Mean expression agrees with the shipped in vivo
    bradyzoite column at rho 0.75.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/brady_subtypes.py
"""
from _common import run

KEY = "brady_subtypes"

if __name__ == "__main__":
    run(KEY)
