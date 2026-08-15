#!/usr/bin/env python3
"""Oocyst-versus-tachyzoite phosphoproteome

Measured-site counts and strongest up/down phosphosite ratios

    level / kind : post_translation / phosphoproteomics
    provides     : phospho_up_sites, phospho_up_ratio, phospho_down_sites, phospho_down_ratio, phospho_sites_measured
    coverage     : 1,603 (19.7%)
    citation     : Wang Z-X et al., Comparative Phosphoproteomic Analysis of Sporulated Oocysts and Tachyzoites of Toxoplasma gondii Reveals Stage-Specific Patterns. Molecules 2022;27:1109
    PMID         : 35164288
    accession    : PXD017032
    url          : https://www.ebi.ac.uk/pride/archive/projects/PXD017032
    local path   : toxo_stage_atlas/data/proteomics/PXD017032_supp_TableS1_upregulated_phosphosites_oocyst_vs_tachy.xlsx

Quirks that cost time once:
    A per-gene cache cannot retain residue positions. It carries how many sites were measured,
    how many moved each way and the median ratios; the source workbooks remain the residue-level
    record.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.expression.load_all()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/phospho_quantitative.py
"""
from _common import run

KEY = "phospho_quantitative"

if __name__ == "__main__":
    run(KEY, normalized_by='expression.load_all()')
