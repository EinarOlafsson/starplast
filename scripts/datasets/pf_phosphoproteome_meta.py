#!/usr/bin/env python3
"""Plasmodium phosphosites (re-analysis of all public data)

Distinct phosphorylated residues per gene, pooled across every public study

    level / kind : post_translation / phosphoproteomics
    provides     : n_phosphosites, has_phospho
    coverage     : 16,318 sites over 2,503 genes
    accession    : PXD046874
    url          : https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/11/PXD046874/
    local path   : datasets/reference/plasmodb/phosphosites/

Quirks that cost time once:
    A re-analysis of every public Plasmodium phosphoproteomics dataset through one pipeline,
    which is what makes a per-gene count meaningful: the same serine found by three groups is
    one site rather than three. Counted as distinct (gene, position) pairs and NOT as rows -- a
    site-centric table still carries one row per peptidoform and per source run, so summing rows
    would count how often a protein was looked at instead of how many sites it has, and the
    files hold millions of rows for 16,318 sites. Accessions are stripped of their transcript
    and product suffix (`PF3D7_1346300.1-p1`), or a gene with two products counts its sites
    twice. Where a source has a merged table the merged one is used and its per-run siblings are
    skipped. Missingness mirrors the Toxoplasma arm: the COUNT stays missing where nothing was
    detected, because how many sites a protein has is genuinely unknown if mass spectrometry
    never saw it, while the FLAG is False, because whether it was ever observed phosphorylated
    is a question about the evidence and the answer is no. Validated on orderings rather than
    totals: 67% of kinases carry a site against 44% of genes at large, and site count rises with
    protein length at rho +0.46.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_phosphoproteome_meta.py
"""
from _common import run

KEY = "pf_phosphoproteome_meta"

if __name__ == "__main__":
    run(KEY)
