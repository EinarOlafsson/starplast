#!/usr/bin/env python3
"""Plasmodium falciparum life-stage and polysomal RNA

Transcript abundance across seven life stages, and what is on ribosomes

    level / kind : transcription / RNAseq
    provides     : expr_ring, expr_early_trophozoite, expr_late_trophozoite, expr_schizont, expr_gametocyte_ii, expr_gametocyte_v, expr_ookinete, expr_asexual_blood, expr_oocyst, expr_sporozoite, polysomal_ring, polysomal_trophozoite, polysomal_schizont, steady_state_ring, steady_state_trophozoite, steady_state_schizont, protein_stage_share_ring, protein_stage_share_trophozoite, protein_stage_share_schizont, antisense_asexual_blood, antisense_oocyst, antisense_sporozoite
    coverage     : 5,720 P. falciparum genes
    accession    : PlasmoDB: Su seven stages, Bunnik polysomal IDC, Gomez-Diaz mosquito stages
    url          : https://plasmodb.org/plasmo/service/record-types/transcript/searches/GenesByTaxon/reports/attributesTabular
    local path   : starplast/data/pf_nodes.parquet

Quirks that cost time once:
    Three studies in one report, split across slots so that no column answers two questions: the
    seven-stage study answers the individual stages, and the polysomal study is an IDC time
    course at 0h, 18h and 36h whose steady-state arm answers cell-cycle phase while its
    polysomal arm answers translation. Polysome-associated RNA is what is ON ribosomes rather
    than a transcript level, and keeping its steady-state partner is what makes that distinction
    measurable instead of assumed. Stage labels are checked against marker genes rather than
    trusted, because the columns are matched by substring out of one wide report and a
    mislabelling would silently shift a stage: CSP peaks in sporozoite, MSP1 in schizont, Pfs16
    in gametocyte II. Pfs25 is the informative one -- its TRANSCRIPT peaks in gametocyte V
    rather than in the ookinete where the protein acts, which is the textbook translational-
    repression stockpile and would look like an off-by-one error to anyone who checked the
    protein instead. The antisense columns sit beside their sense partners rather than being
    reduced to a ratio, because the denominator is what makes the ratio interpretable. Adding
    them exposed a silent matcher fault worth recording: `sense - asexual blood stages` is a
    SUBSTRING of `antisense - asexual blood stages`, so plain containment made each sense entry
    match two headers, fail its one-match test and vanish -- fetching antisense DELETED sense
    and the slot count fell by two while a slot was being added. The matcher now requires the
    label to start the header or follow a non-alphanumeric character, and a test asserts no
    declared sample goes missing. The TMT proteome in the same report is COMPOSITIONAL and is
    named for it: PlasmoDB serves the channels row-normalised, so a gene's three values sum to a
    constant and the columns anti-correlate by construction. They say which stage a protein sits
    in, not how much there is, so `protein abundance · asexual blood stage` is left EMPTY rather
    than filled with a share. The tell would have been well hidden: ring protein correlates
    -0.25 with ring mRNA, which reads as a biological puzzle and is only the normalisation
    showing through.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/plasmodb_pf3d7_expression.py
"""
from _common import run

KEY = "plasmodb_pf3d7_expression"

if __name__ == "__main__":
    run(KEY)
