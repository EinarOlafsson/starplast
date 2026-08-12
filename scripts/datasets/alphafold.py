#!/usr/bin/env python3
"""AlphaFold DB

Per-gene mean pLDDT; coordinates fetched on demand

    level / kind : reference / structure
    provides     : mean_plddt
    coverage     : 6,480 (79.6%)
    citation     : Varadi et al. 2024 NAR (database); Jumper et al. 2021 Nature (method)
    accession    : UP000001529 (taxid 508771), AlphaFold DB
    url          : https://alphafold.ebi.ac.uk/api/prediction/{acc}

Quirks that cost time once:
    Missing for the largest proteins, which here are disproportionately secreted effectors. NO
    ANONYMOUS BULK DOWNLOAD: the per-proteome tar exists at gs://public-datasets-deepmind-
    alphafold-v4/proteomes/proteome-tax_id-508771-0_v4.tar but plain HTTPS returns 403, so it
    needs `gcloud storage cp` and a Google account. The per-accession API above is the
    credential-free route and is what structures.py uses. The database paper is Varadi et al.,
    not Jumper et al. -- Jumper is the method, and there is no paper by Jumper titled after the
    database.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/alphafold.py
"""
from _common import run

KEY = "alphafold"

if __name__ == "__main__":
    run(KEY)
