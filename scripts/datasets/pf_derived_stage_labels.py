#!/usr/bin/env python3
"""Plasmodium peak expression and stage label (DERIVED)

Maximum expression across stages, and which stage a gene belongs to

    level / kind : transcription / RNAseq
    provides     : expr_max, stage_enriched_derived, stage_margin_derived
    coverage     : 5,720 genes for the maximum, 310 labelled
    local path   : starplast/data/pf_nodes.parquet
    derived from : expr_ring, expr_early_trophozoite, expr_late_trophozoite, expr_schizont, expr_gametocyte_ii, expr_gametocyte_v, expr_ookinete, expr_oocyst, expr_sporozoite

Quirks that cost time once:
    COMPUTED from the stage columns and declaring it, so leakage closure excludes them together.
    The stage call reuses `cellcycle.stage_enrichment` rather than reimplementing it, which is
    deliberate: if the two arms' stage labels are ever compared, a difference should mean the
    biology differs and not that one z-scored and the other did not. Only 310 of 5,720 genes are
    labelled, because a gene is left unlabelled unless one stage leads the next by half a z-unit
    -- a label that is really a coin toss looks like a measurement in every table it reaches.
    Read the class counts with the same caveat the Toxoplasma arm carries: ookinete takes 181 of
    the 310 not because it uses more genes but because ring, trophozoite and schizont are highly
    correlated with one another and rarely win by a margin, while the mosquito stages are
    separable. The margin rule is working; the interpretation is what needs care.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_derived_stage_labels.py
"""
from _common import run

KEY = "pf_derived_stage_labels"

if __name__ == "__main__":
    run(KEY)
