#!/usr/bin/env python3
"""Organelle-surface proximity proteomes

Enrichment near the cytosolic face of the apicoplast, mitochondrion and ER

    level / kind : post_translation / BioID
    provides     : surface_apicoplast_log2fc, surface_apicoplast_stringent, surface_mitochondrion_log2fc, surface_mitochondrion_stringent, surface_er_log2fc, surface_er_stringent
    coverage     : 742 proteins
    citation     : Parker KV, Huet D. A proximity biotinylation approach for the identification of membrane contact site proteins in Toxoplasma gondii. bioRxiv 2026, doi:10.64898/2026.08.05.743015 (preprint)
    accession    : bioRxiv 10.64898/2026.08.05.743015 Table S1
    url          : https://www.biorxiv.org/content/10.64898/2026.08.05.743015v1.supplementary-material
    local path   : datasets/post_translation/BioID/organelle_surface_turboid_2026/media-1.xlsx

Quirks that cost time once:
    Baits anchored in each outer membrane with their tail in the cytosol, so this is the OUTSIDE
    of an organelle -- a different question from hyperLOPIT, which says which organelle a
    protein is in. Checked against it anyway: stringent mitochondrial hits are mitochondrial at
    odds 13.8 and ER hits ER at 7.6, but the apicoplast bait shows no enrichment for apicoplast
    proteins (odds 1.2, p = 0.8), so that arm is recorded as proximity and not as a location. A
    protein a bait never detected keeps a MISSING flag, not a zero: that bait did not test it.
    One accession in the deposit is a backtick, repaired from its product text only because that
    text names exactly one gene. A preprint.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/organelle_surface_turboid.py
"""
from _common import run

KEY = "organelle_surface_turboid"

if __name__ == "__main__":
    run(KEY)
