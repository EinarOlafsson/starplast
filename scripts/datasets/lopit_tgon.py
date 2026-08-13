#!/usr/bin/env python3
"""T. gondii hyperLOPIT

Subcellular compartment, MAP and MCMC, with posteriors

    level / kind : post_translation / LOPIT
    provides     : compartment, lopit_map, lopit_mcmc, lopit_prob_map, lopit_prob_mcmc
    coverage     : 3,827 (47.0%)
    citation     : A Comprehensive Subcellular Atlas of the Toxoplasma Proteome via hyperLOPIT (Barylyuk et al. 2020)
    PMID         : 33053376
    url          : https://ars.els-cdn.com/content/image/1-s2.0-S193131282030514X-mmc5.xls
    local path   : datasets/lopit_toxoplasma_gondii_ME49.csv

Quirks that cost time once:
    MAP and MCMC disagree for 980 of 3,827 (26%). Assignment tracks abundance, so the unassigned
    half is biased toward low-abundance proteins.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.localization.lopit_labels()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/lopit_tgon.py
"""
from _common import run

KEY = "lopit_tgon"

if __name__ == "__main__":
    run(KEY, normalized_by='localization.lopit_labels()')
