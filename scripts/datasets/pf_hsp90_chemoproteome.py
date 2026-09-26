#!/usr/bin/env python3
"""Blood-stage proteome and Hsp90 dependence

Protein abundance in a DMSO control, what two Hsp90 inhibitors do to it, and the paper's chaperone-dependent call

    level / kind : translation / proteomics
    provides     : proteome_blood_log2, hsp90_inhibition_ga_log2fc, hsp90_inhibition_xl_log2fc, hsp90_dependent
    coverage     : 3,049 proteins (53%)
    citation     : Ibrasheva N et al., Chemoproteomic profiling of Plasmodium falciparum Hsp90 inhibition reveals functional link to DNA replication pathways. bioRxiv 2026, doi:10.64898/2026.08.28.747854 (preprint)
    accession    : PXD079493 / bioRxiv 10.64898/2026.08.28.747854 Tables S1-S2
    url          : https://www.biorxiv.org/content/10.64898/2026.08.28.747854v1.supplementary-material
    local path   : datasets/translation/proteomics/pf_hsp90_2026/TableS1_media-2.xlsx

Quirks that cost time once:
    Two empty slots and one new question. The DMSO arm is the first protein abundance in the map
    for the asexual blood stage -- the stage this parasite spends its life in, confirmed from
    the PRIDE record rather than assumed -- and the two inhibitors are the first abundance under
    stress. Chaperone dependence is the new question: which proteins need Hsp90 to stay folded,
    measured with two chemically unrelated inhibitors so the answer is the chaperone's and not
    one compound's -- and that conjunction is doing real work, because the two continuous
    responses agree only at rho 0.17. The paper's 131 hits reproduce EXACTLY from its own rule
    (down by at least 0.5 log2 at p < 0.05 under both), all 131 of them; 124 survive into the
    table, because seven protein groups name more than one gene and are dropped rather than
    assigned to the first. UniProt accessions are mapped to genes through the deposit's own
    Spectronaut report, not an external lookup. The PRIDE description says 133 hits where the
    table lists 131; the reproducible number ships and the discrepancy is recorded. A preprint.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_hsp90_chemoproteome.py
"""
from _common import run

KEY = "pf_hsp90_chemoproteome"

if __name__ == "__main__":
    run(KEY)
