#!/usr/bin/env python3
"""Life-cycle stage enrichment (DERIVED)

Which stage a gene's own expression is highest in

    level / kind : transcription / RNAseq
    provides     : stage_enriched_derived, stage_margin_derived
    coverage     : 1,911 of 8,140 genes called
    derived from : expr_tachy, expr_cyst, expr_sporulated

Quirks that cost time once:
    DERIVED, not measured: computed here from expr_tachy / expr_cyst / expr_sporulated by
    z-scoring each and taking the argmax where it leads by 0.5 z. It is a restatement of those
    columns, so holding it out against an embedding built on them is circular by construction.
    Left unlabelled where no stage leads clearly.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.cellcycle.add_all()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/stage_enriched.py
"""
from _common import run

KEY = "stage_enriched"

if __name__ == "__main__":
    run(KEY, normalised_by='cellcycle.add_all()')
