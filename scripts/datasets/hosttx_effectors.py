#!/usr/bin/env python3
"""Host-transcription effector screen

Hotelling T2 plus full per-effector host-response signature

    level / kind : DNA / CRISPR_screen
    provides     : hosttx_T2, hosttx_padj, hosttx_signature_pc01, hosttx_signature_pc02, hosttx_signature_pc03, hosttx_signature_pc04, hosttx_signature_pc05, hosttx_signature_pc06, hosttx_signature_pc07, hosttx_signature_pc08, hosttx_signature_pc09, hosttx_signature_pc10, hosttx_signature_pc11, hosttx_signature_pc12, hosttx_signature_pc13, hosttx_signature_pc14, hosttx_signature_pc15, hosttx_signature_pc16, hosttx_signature_pc17, hosttx_signature_pc18, hosttx_signature_pc19, hosttx_signature_pc20, hosttx_signature_norm, hosttx_signature_n_de
    coverage     : 252 screened / 22 full signatures
    citation     : High-throughput identification of Toxoplasma gondii effector proteins that target host cell transcription
    PMID         : 37827122
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12033024/supplementaryFiles
    local path   : datasets/DNA/CRISPR_screen/37827122/

Quirks that cost time once:
    The 737,726-row host differential-expression table is represented by 20 PCA coordinates, its
    L2 norm and substantial-DE count; PCA is a dimensional summary, not a host-gene measurement.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.screens.crispr_screens() + screens.host_transcription_signatures()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/hosttx_effectors.py
"""
from _common import run

KEY = "hosttx_effectors"

if __name__ == "__main__":
    run(KEY, normalized_by='screens.crispr_screens() + screens.host_transcription_signatures()')
