#!/usr/bin/env python3
"""Enzyme classification (ToxoDB)

EC number per gene, and whether it has one

    level / kind : reference / annotation
    provides     : ec_number, has_ec
    coverage     : 1,313 enzymes of 8,140 genes
    accession    : ToxoDB ME49
    url          : https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByTaxon/reports/attributesTabular
    local path   : starplast/data/toxodb_ec_numbers.tsv

Quirks that cost time once:
    An annotation rather than a measurement, with the same standing as the InterPro domains that
    fill `domain content`. It matters because it is the ONLY gene-indexed metabolic datum there
    is -- every other metabolism question in the catalog is about metabolites, and a metabolite
    is not a gene. Verified against conservation: enzymes have a Plasmodium ortholog 59.6% of
    the time against 30.4% for other genes (odds 3.39, p = 8e-88), are lineage-specific a third
    as often, and are more costly to lose in vitro. `has_ec` is 0 and not missing where ToxoDB
    reports no EC: the whole proteome was asked, so no assignment is an answer about the gene.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/toxodb_ec_numbers.py
"""
from _common import run

KEY = "toxodb_ec_numbers"

if __name__ == "__main__":
    run(KEY)
