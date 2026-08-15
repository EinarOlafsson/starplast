#!/usr/bin/env python3
"""MORC depletion and BFD1 perturbation transcriptome

MORC knockdown, BFD1 knockout and BFD1 stabilization series

    level / kind : transcription / RNAseq
    provides     : morc_MORC_UT1, morc_MORC_UT2, morc_MORC_UT3, morc_MORC_IAA1, morc_MORC_IAA2, morc_MORC_IAA3, morc_MORC-KD-BFD1-KO_UT1, morc_MORC-KD-BFD1-KO_UT2, morc_MORC-KD-BFD1-KO_UT3, morc_MORC-KD-BFD1-KO_IAA1, morc_MORC-KD-BFD1-KO_IAA2, morc_MORC-KD-BFD1-KO_IAA3, morc_DD-BFD1-Ty_UT1, morc_DD-BFD1-Ty_UT2, morc_DD-BFD1-Ty_UT3, morc_DD-BFD1-Ty_Shield1, morc_DD-BFD1-Ty_Shield2, morc_DD-BFD1-Ty_Shield3
    coverage     : 7,841 (96.3%)
    accession    : PXD058095
    url          : https://www.ebi.ac.uk/pride/archive/projects/PXD058095
    local path   : toxo_stage_atlas/data/proteomics/PXD058095_supp_DatasetEV1_MORC_RNAseq_counts_TPM.xlsx

Quirks that cost time once:
    The accession and file location sit in a proteomics collection, but the shipped workbook is
    explicitly RNA-seq counts/TPM and is normalized as transcription. The originating
    publication still needs confirmation before citation.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.expression.load_all()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/morc_depletion.py
"""
from _common import run

KEY = "morc_depletion"

if __name__ == "__main__":
    run(KEY, normalized_by='expression.load_all()')
