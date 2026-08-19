#!/usr/bin/env python3
"""Plasmodium export prediction (ExportPred)

Predicted export to the erythrocyte, as an ordinal confidence tier

    level / kind : reference / annotation
    provides     : export_pred_tier, is_exported
    coverage     : 440 genes called at some threshold, 191 at the default
    accession    : PlasmoDB GenesByExportPrediction
    url          : https://plasmodb.org/plasmo/service/record-types/transcript/searches/GenesByExportPrediction/reports/attributesTabular
    local path   : datasets/reference/plasmodb/exportpred/

Quirks that cost time once:
    PREDICTED, not measured, which is why it answers `export / PEXEL trafficking` and not
    `exposure to host cytosol` -- that slot wants a measured exportome and a sequence model
    filling it would be a model answering for an experiment. PlasmoDB serves ExportPred as a
    search with a score threshold rather than as a per-gene attribute, so the score is recovered
    by asking at several thresholds and keeping the highest a gene survives; the scale
    saturates, since asking for 20 returns nothing, so 10 is the algorithm's own default and the
    top tier rather than an arbitrary cut. It is a TIER and not a boolean because the default
    loses real biology: MESA and PfEMP3 are exported by any textbook and both fall below 10,
    while KAHRP and the FIKK kinases sit above it. Absence is a real negative here and not a gap
    -- a sequence model was evaluated on every protein, so its silence is a prediction of not-
    exported, which is the opposite of the screen columns.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/plasmodb_pf3d7_exportpred.py
"""
from _common import run

KEY = "plasmodb_pf3d7_exportpred"

if __name__ == "__main__":
    run(KEY)
