#!/usr/bin/env python3
"""In vivo brain-stage transcriptome

Tachyzoites, acute/chronic whole brain, and purified bradyzoites

    level / kind : transcription / RNAseq
    provides     : invivo_TZ_1, invivo_TZ_2, invivo_WholeBrain_Acute_1, invivo_WholeBrain_Acute_2, invivo_WholeBrain_Acute_3, invivo_WholeBrain_Chronic_1, invivo_WholeBrain_Chronic_2, invivo_WholeBrain_Chronic_3, invivo_BZ_28DPI_1, invivo_BZ_28DPI_3, invivo_BZ_90DPI_1, invivo_BZ_90DPI_2, invivo_BZ_120DPI_2, invivo_BZ_120DPI_3
    coverage     : 7,663 (94.1%)
    citation     : Garfoot AL et al., Proteomic and transcriptomic analyses of early and late-chronic Toxoplasma gondii infection shows novel and stage specific transcripts. BMC Genomics 2019;20:859
    PMID         : 31726967
    url          : https://static-content.springer.com/esm/art%3A10.1038%2Fs12864-019-6213-0/MediaObjects/12864_2019_6213_MOESM4_ESM.csv
    local path   : datasets/translation/proteomics/31726967/12864_2019_6213_MOESM4_ESM.csv

Quirks that cost time once:
    These columns are FPKM-derived transcript abundance, despite the mixed
    transcriptome/proteome paper and the legacy proteomics directory. They belong to
    transcription slots, never fitness.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.expression.load_all()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/invivo_brain_transcriptome.py
"""
from _common import run

KEY = "invivo_brain_transcriptome"

if __name__ == "__main__":
    run(KEY, normalized_by='expression.load_all()')
