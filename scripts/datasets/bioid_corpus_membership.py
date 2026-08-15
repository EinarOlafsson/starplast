#!/usr/bin/env python3
"""BioID/TurboID supplement membership corpus

Number of downloaded proximity-labeling studies whose supplement names each gene

    level / kind : post_translation / proximity_labelling
    provides     : n_bioid_studies
    coverage     : measured at build time
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/supplementaryFiles
    local path   : datasets/post_translation/BioID/

Quirks that cost time once:
    Membership is not enrichment and is never converted to an interaction edge.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.interaction_studies.parse_studies()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/bioid_corpus_membership.py
"""
from _common import run

KEY = "bioid_corpus_membership"

if __name__ == "__main__":
    run(KEY, normalized_by='interaction_studies.parse_studies()')
