#!/usr/bin/env python3
"""CRISPR screen at low and high parasite density

Fitness at low and high infection density, and which genes high density needs

    level / kind : DNA / CRISPR_screen
    provides     : fit_density_low, fit_density_high, fit_density_dependence, fit_density_dependence_log10padj, fit_density_dim
    coverage     : 7,461 (91.7%)
    citation     : Giuliano CJ, Kalluraya CA, Kloehn J, Sloan MA, Bunkofske ME, Hunter CA, Soldati-Favre D, Harding CR, Lourido S. Convergent evolution of metabolic regulation governs redox adaptation in Toxoplasma. Cell 2026 Aug, doi:10.1016/j.cell.2026.07.029
    PMID         : 42580337
    accession    : Cell 2026 Supplementary Table 1 (mmc2)
    url          : https://ars.els-cdn.com/content/image/1-s2.0-S0092867426008275-mmc2.xlsx
    local path   : datasets/DNA/CRISPR_screen/42580337/mmc2.xlsx

Quirks that cost time once:
    A genome-wide knockout library selected for four passages at MOI 1, then split for four more
    at MOI 0.3 (low density) and MOI 3 (high). The arms are the authors' passage-8 gene scores
    (mean gRNA log2 fold change to input; negative = needed) and are fibroblast fitness again: r
    = 0.995 between them, as the paper states, and rho 0.70 with fit_invitro_hff, so they are
    held out with it. The dependence column is the authors' contrast on barcoded gRNA-UMI
    clones, log2(high / low): NEGATIVE = needed at high density. It is orthogonal to bulk
    fitness (rho -0.08) and is the new axis; -log10 of its Bonferroni-adjusted t-test p ships
    beside it, and fit_density_dim marks the paper's 31 density-inhibited mutants (NMNAT, NAD
    synthetase, NAD kinase, nicotinamidase, the glucose transporter, TgPRO) among the ~6,150
    genes the contrast scored. Verified: r = 0.9951, 31 hits within the 32 genes at adj. p <
    0.05, the 12 high-confidence hits exactly from the stated rule. The 266-gene targeted
    follow-up screen (with nicotinamide) is not shipped. GT1 accessions, resolved to ME49
    through the identity layer. The article is not open access; PMC blocks scripted download,
    the publisher CDN serves the table.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/crispr_parasite_density.py
"""
from _common import run

KEY = "crispr_parasite_density"

if __name__ == "__main__":
    run(KEY)
