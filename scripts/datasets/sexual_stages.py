#!/usr/bin/env python3
"""Sexual development in the cat (single-cell atlas)

Enrichment at 8 days post-infection, when gametogony happens

    level / kind : transcription / scRNAseq
    provides     : sexual_stage_8dpi_log2fc
    coverage     : 4,463 genes (55%)
    PMID         : 41929010
    accession    : PMC13042011 supplementary media-2
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13042011/supplementaryFiles
    local path   : starplast/data/sexual_stages.tsv

Quirks that cost time once:
    The only stage in the map that is both inside a cat and sexual; the enteroepithelial column
    holds the asexual stages that precede it. Verified by what rises: oocyst wall protein sits
    at the 97th percentile, and the oocyst wall is built at the end of the sexual cycle, while
    ribosomal housekeeping genes sit at the 23rd. Accessions arrive as `DEAD/DEAHboxhelicase-
    TGME49-220860` -- product description glued to the accession with hyphens for underscores --
    so they are extracted and normalised rather than matched.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/sexual_stages.py
"""
from _common import run

KEY = "sexual_stages"

if __name__ == "__main__":
    run(KEY)
