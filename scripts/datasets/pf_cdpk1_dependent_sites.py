#!/usr/bin/env python3
"""Plasmodium CDPK1-dependent phosphosites

Phosphosites per gene that are lost when PfCDPK1 is knocked down

    level / kind : post_translation / phosphoproteomics
    provides     : cdpk1_dependent_sites
    coverage     : 62 genes, 73 sites
    PMID         : 28680058
    accession    : Nat Commun 00053 Supplementary Data 2a
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5498596/supplementaryFiles
    local path   : datasets/post_translation/kinase_substrate/28680058/41467_2017_53_MOESM3_ESM.xls

Quirks that cost time once:
    Keyed by SEQUENCE, because nothing else in the file can be resolved: the sites are numbered
    in a 2017 annotation (`3885720(S422)`) that neither current accessions nor PlasmoDB's
    previous-id list carries. The 15-residue window around each site is an identifier when it
    occurs in exactly one protein -- 73 of 79 match one gene, 6 match none, none matches two --
    and the translation comes from `codons.translate` over the CDS table this project already
    ships, so the genetic code is not written down twice. The mapping is then CHECKED against a
    field it did not use: 69 of 73 rows agree with the current product description, and the four
    that do not are re-annotations rather than wrong genes (a `conserved membrane protein` now
    named basal complex protein bleb, a `formin 2` now an Eps15-like protein, and two that my
    word matcher split on a digit). The set then reproduces the paper's own conclusion from the
    other side: 14.5% of the 62 genes are motor, IMC or invasion machinery against 1.6% of the
    proteome, and GAP45, myosin A, actin I and IMC1c/1g are all in it. Named for DEPENDENCE and
    not for substrate -- a site lost under knockdown may be phosphorylated by this kinase or by
    something downstream of it, and the file cannot tell them apart. The paper's PfPKA-R result
    is not in this sheet and is not claimed here.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_cdpk1_dependent_sites.py
"""
from _common import run

KEY = "pf_cdpk1_dependent_sites"

if __name__ == "__main__":
    run(KEY)
