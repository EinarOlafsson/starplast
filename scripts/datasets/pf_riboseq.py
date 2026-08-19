#!/usr/bin/env python3
"""Plasmodium ribosome profiling across the asexual cycle

Ribosome-footprint and mRNA density per gene at five points of the blood-stage cycle

    level / kind : translation / riboseq
    provides     : riboseq_rpf_ring, riboseq_rpf_early_trophozoite, riboseq_rpf_late_trophozoite, riboseq_rpf_schizont, riboseq_rpf_merozoite, riboseq_mrna_ring, riboseq_mrna_early_trophozoite, riboseq_mrna_late_trophozoite, riboseq_mrna_schizont, riboseq_mrna_merozoite
    coverage     : 3,501 genes (61%), 2,182 at the ring and 1,174 at the merozoite
    PMID         : 25493618
    accession    : GSE58402
    url          : https://ftp.ncbi.nlm.nih.gov/geo/series/GSE58nnn/GSE58402/suppl/
    local path   : datasets/translation/riboseq/25493618/

Quirks that cost time once:
    The first MEASURED translation on this arm: the slot was previously answerable only through
    polysome-associated RNA, which is what is on ribosomes rather than how much ribosome is on
    it. Both arms of the experiment ship as conditions and the ratio between them does NOT,
    which is the second time this project has computed a Plasmodium translation efficiency and
    refused it -- and this time on a different instrument, a different strain and a different
    decade. Ribosomal proteins carry far more footprint than the rest at every stage (median
    log1p 5.6 to 7.9 against 3.5 to 4.1, p <= 2e-18), which is the check that says the
    measurement behaves; but the ratio puts them BELOW the rest at the schizont, and correlates
    negatively with codon adaptation at every stage (rho -0.01 to -0.12), where the textbook
    expectation is positive. Two checks disagreeing is the contradictory case, whose answer is
    to ship the conditions. What the replication adds is where to look: `codon_cai_ribosomal`
    does not separate the very ribosomal proteins it is built from in this genome (0.710 against
    0.717, p = 0.37) while it does in Toxoplasma (0.771 against 0.714, p = 4e-22), so the
    quantity that fails to behave is the codon index, not the footprints. Stage labels are the
    deposit's own and are CHECKED against the independent PlasmoDB stage series: ring, early
    trophozoite, late trophozoite and schizont each correlate highest with their own stage (rho
    0.65 to 0.77). The merozoite arm has no counterpart there and lands on the ring, which is
    the neighbouring point of the cycle rather than a contradiction. Strain W2, not 3D7, so the
    surface-antigen families are the place to distrust it. Keyed on pre-2012 accessions and
    resolved through `plasmodb_identity`; the deposit's `-a`/`-b` split entries are dropped
    rather than summed, since RPKM is already length-normalised.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_riboseq.py
"""
from _common import run

KEY = "pf_riboseq"

if __name__ == "__main__":
    run(KEY)
