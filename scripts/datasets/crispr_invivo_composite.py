#!/usr/bin/env python3
"""In vivo CRISPR composite scores

Differential fitness x significance in six mouse tissues after acute infection

    level / kind : DNA / CRISPR_screen
    provides     : fit_invivo_PE, fit_invivo_lung, fit_invivo_liver, fit_invivo_spleen, fit_invivo_heart, fit_invivo_brain
    coverage     : 7,395 (90.8%)
    citation     : Giuliano CJ et al., CRISPR-based functional profiling of the Toxoplasma gondii genome during acute murine infection. Nat Microbiol 2024;9:2323-2343
    PMID         : 38977907
    accession    : Nat Microbiol 2024 Supplementary Data 5, sheet Genome-Wide Differential
    url          : https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-024-01754-2/MediaObjects/41564_2024_1754_MOESM3_ESM.zip
    local path   : datasets/DNA/CRISPR_screen/38977907/Supplementary Data 5 (Genome-wide mouse screen scores).xlsx

Quirks that cost time once:
    CITATION CORRECTED 2026-09-26. The four original tissues were fetched from ToxoDB's
    tgonGt1CrisprFunc* tracks and cited to the 2019 platform paper (PMID 31481656), whose own
    supplements are 36-620-gene libraries. They are this study's genome-wide composite scores:
    Spearman 0.9991-0.9993 and a maximum absolute difference of 5e-8 against the supplement,
    each tissue matching itself best among all 32 numeric columns. Heart and brain come from the
    same sheet (scripts/derive_deposits.py; starplast/deposits.py) and are the only in vivo
    heart and brain CRISPR fitness there is; both are bottleneck-limited (medians -7.9 and -7.1
    raw), so most genes carry little information there. The ToxoDB attributes are gone from
    release 71, so the supplement is now the durable source. GT1 accessions; the 68 loci GT1
    splits into A/B and ME49 does not are averaged onto the one ME49 gene.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.screens.crispr_screens()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/crispr_invivo_composite.py
"""
from _common import run

KEY = "crispr_invivo_composite"

if __name__ == "__main__":
    run(KEY, normalized_by='screens.crispr_screens()')
