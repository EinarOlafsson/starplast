#!/usr/bin/env python3
"""In vivo CRISPR composite scores

Peritoneum, lung, liver, spleen composite scores

    level / kind : DNA / CRISPR_screen
    provides     : fit_invivo_PE, fit_invivo_lung, fit_invivo_liver, fit_invivo_spleen
    coverage     : 7,395 (90.8%)
    PMID         : 31481656
    accession    : ToxoDB tgonGt1CrisprFunc*
    url          : https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByTaxon/reports/attributesTabular?organism=%5B%22Toxoplasma%20gondii%20GT1%22%5D&reportConfig=%7B%22attributes%22%3A%5B%22primary_key%22%2C%22tgonGt1CrisprMeanPhenotype%22%2C%22tgonGt1CrisprFuncPE%22%2C%22tgonGt1CrisprFuncLung%22%2C%22tgonGt1CrisprFuncLiver%22%2C%22tgonGt1CrisprFuncSpleen%22%5D%2C%22includeHeader%22%3Atrue%2C%22attachmentType%22%3A%22plain%22%7D

Quirks that cost time once:
    Corresponds to the in vivo CRISPR platform paper; confirm before citing. The URL pulls the
    tgonGt1CrisprFunc* tracks straight from ToxoDB, keyed on GT1.

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
