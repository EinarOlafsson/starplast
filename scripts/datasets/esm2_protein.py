#!/usr/bin/env python3
"""Frozen ESM-2 protein representations

320 sequence representation coordinates in a separate gene-keyed feature table

    level / kind : reference / computed sequence features
    provides     : (edges or build inputs only)
    coverage     : 8,064 T. gondii genes with valid bundled CDS translations
    citation     : Lin et al., Science (2023); frozen ESM-2 8M model, local sequence encoding
    url          : https://huggingface.co/facebook/esm2_t6_8M_UR50D
    local path   : starplast/data/esm_features.parquet

Quirks that cost time once:
    Computed representation, not an experimentally measured trait. Model revision and sequence
    hashes are in esm_manifest.parquet. Long proteins use overlapping windows with residue-
    complete pooling. Features are loaded for prediction rather than adding 320 uninterpretable
    colour controls to the display map.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.sequence_features.encode()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/esm2_protein.py
"""
from _common import run

KEY = "esm2_protein"

if __name__ == "__main__":
    run(KEY, normalized_by='sequence_features.encode()')
