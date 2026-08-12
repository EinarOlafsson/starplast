#!/usr/bin/env python3
"""P. falciparum LOPIT

Donor labels for orthoLOPIT transfer

    level / kind : post_translation / LOPIT
    provides     : (edges or build inputs only)
    coverage     : 1,646 usable
    citation     : Chisholm SA et al., The spatial proteome of the Plasmodium falciparum schizont. Nat Commun 2026;17:6192 -- CONFIRM against the file on disk
    PMID         : 42218142
    url          : https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-026-73664-2/MediaObjects/41467_2026_73664_MOESM3_ESM.xlsx
    local path   : datasets/lopit_plasmodium_falciparum_3D7.csv

Quirks that cost time once:
    The URL downloads and is the right kind of data, but this file predates the registry entry,
    so that this paper is the source of THIS csv is inference, not verification. Confirm against
    the file before citing.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.localisation.lopit_labels()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/lopit_pfal.py
"""
from _common import run

KEY = "lopit_pfal"

if __name__ == "__main__":
    run(KEY, normalised_by='localisation.lopit_labels()')
