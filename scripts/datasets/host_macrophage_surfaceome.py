#!/usr/bin/env python3
"""Mouse bone-marrow macrophage cell-surface repertoire

Which host proteins are EXPOSED on the surface of the macrophage a tachyzoite invades

    level / kind : reference / proteomics
    provides     : bmdm_surface_detected, bmdm_surface_intensity
    coverage     : 1,296 mouse surface proteins, 150 of them on primary BMDM
    PMID         : 25894527
    accession    : PLoS ONE 0121314 S1 File
    url          : https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0121314.s002&type=supplementary
    local path   : datasets/host/cspa/25894527/S1_File.xlsx

Quirks that cost time once:
    Cell-surface capture, so this is a REPERTOIRE and not a proteome: a protein here is on the
    outside of an intact cell, which is the surface the parasite meets, rather than merely
    present somewhere in it. The row space is every protein the atlas saw on the surface of ANY
    mouse cell type, which makes a `False` a real negative -- the same capture ran on
    macrophages and did not find it -- and that is what instruction 39 means by a repertoire
    slot being FILLED rather than averaged. The shipped flag is a NULLABLE boolean: the host
    table also holds human red cell rows where the question was never asked, and a plain bool
    would turn those into `False`. Of the 41 human and 31 mouse cell types in the atlas, this is
    the only one either parasite lives in -- there is no erythrocyte, no hepatocyte and no
    primary fibroblast in it, so it fills one slot and not six. The lab's own copy of the file
    is a git-lfs pointer served as the spreadsheet (132 bytes that open as nothing); the
    journal's supplement is the same bytes and is what this fetches. The deposit's two matrices
    disagree for twelve proteins, nine with a macrophage intensity and no detection mark and
    three the other way; either sheet counts as the authors having measured it there, and the
    count of disagreements is logged rather than smoothed. Self-validating: the strongest
    signals are Emr1 (F4/80), Siglec1 (CD169), Itgb2, Cd47 and H2-K1.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_macrophage_surfaceome.py
"""
from _common import run

KEY = "host_macrophage_surfaceome"

if __name__ == "__main__":
    run(KEY)
