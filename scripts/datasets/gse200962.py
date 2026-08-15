#!/usr/bin/env python3
"""Bradyzoite restriction-checkpoint transcriptome

Cyclin perturbations in tachyzoite and bradyzoite conditions

    level / kind : transcription / RNAseq
    provides     : restriction200962_7th_gt_2_kk31_kim_s1_l001, restriction200962_7th_gt_20_kk27_kim_s14_l001, restriction200962_7th_gt_22_kk29_kim_s15_l001, restriction200962_7th_gt_4_kk32_kim_s2_l001, restriction200962_kk17_p2wo_ph8_2_s5, restriction200962_kk18_p2wo_ph8_2_s1, restriction200962_kk19_p2w_ph8_2_s22, restriction200962_kk20_p2w_ph8_2_s16, restriction200962_kk21_p5wo_ph8_2_s13, restriction200962_kk22_p5wo_ph8_2_s15, restriction200962_kk23_p5w_ph8_2_s12, restriction200962_kk24_p5w_ph8_2_s11, restriction200962_kk5_p2wo_ph7_4_s6, restriction200962_kk6_p2wo_ph7_4_s19, restriction200962_kk7_p2w_ph7_4_s14, restriction200962_kk8_p2w_ph7_4_s23
    coverage     : measured at build time
    accession    : GSE200962
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/GSE200962/suppl/GSE200962_gene_count_matrix_geo.csv.gz
    local path   : datasets/toxoplasma_acquisition_2026_08_14/GSE200962_gene_count_matrix_geo.csv.gz

Quirks that cost time once:
    No publication is linked from GEO; raw sample identifiers are retained verbatim in column
    names rather than assigned conditions by guesswork.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/gse200962.py
"""
from _common import run

KEY = "gse200962"

if __name__ == "__main__":
    run(KEY)
