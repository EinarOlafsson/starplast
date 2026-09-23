#!/usr/bin/env python3
"""Local AlphaFold 3 protein models

sequence-verified AF3 confidence, coverage and confident-region geometry

    level / kind : post_translation / computed structure
    provides     : af3_sequence_coverage, af3_mean_plddt, af3_plddt_q25, af3_confident_sequence_fraction, af3_very_confident_sequence_fraction, af3_low_confidence_modelled_fraction, af3_confident_rg_angstrom, af3_confident_contacts_per_residue, af3_ptm
    coverage     : 1,210 T. gondii genes with exactly matched local AF3 sequences
    citation     : Einar Olafsson, unpublished local AF3 model collection (2026)
    url          : https://github.com/EinarOlafsson/starplast/blob/main/docs/structures.md
    local path   : starplast/data/af3_features.parquet

Quirks that cost time once:
    Computed predictions, not experimental measurements. Fragments retain residue ranges;
    confidence summaries combine overlapping fragments per residue. Geometry only from complete
    models. Coordinates remain on the local shared drive; see source manifest.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.structure_catalog.attach_features()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/local_af3.py
"""
from _common import run

KEY = "local_af3"

if __name__ == "__main__":
    run(KEY, normalized_by='structure_catalog.attach_features()')
