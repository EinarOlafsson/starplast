#!/usr/bin/env python3
"""Oocyst developmental-stage iTRAQ proteome

iTRAQ abundance ratios across oocyst developmental stages

    level / kind : translation / proteomics
    provides     : oocyst_itraq_115_113, oocyst_itraq_116_113, oocyst_itraq_115_114, oocyst_itraq_116_114, oocyst_itraq_115_113_1, oocyst_itraq_116_113_1, oocyst_itraq_115_114_1, oocyst_itraq_116_114_1
    coverage     : 2,079 (25.5%)
    citation     : Possenti A et al., Proteomic Differences between Developmental Stages of Toxoplasma gondii Revealed by iTRAQ-Based Quantitative Proteomics. Front Microbiol 2017;8:1732
    PMID         : 28626452
    accession    : PXD003765
    url          : https://www.ebi.ac.uk/pride/archive/projects/PXD003765
    local path   : toxo_stage_atlas/data/proteomics/PXD003765_supp_mmc2_iTRAQ_ratios_2095proteins.xls

Quirks that cost time once:
    The legacy XLS needs conversion before pandas can read it. Ratios are retained as ratios
    rather than logged or centered, because moving their reference changes the measurement.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.expression.load_all()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/oocyst_itraq.py
"""
from _common import run

KEY = "oocyst_itraq"

if __name__ == "__main__":
    run(KEY, normalized_by='expression.load_all()')
