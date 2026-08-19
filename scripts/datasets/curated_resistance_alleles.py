#!/usr/bin/env python3
"""Validated resistance-conferring mutations (CURATED)

Mutations shown to CAUSE drug resistance by putting them back into a clean background

    level / kind : reference / literature
    provides     : resistance_allele_count, resistance_compound_count, resistance_substitutions, resistance_compounds
    coverage     : 1 gene, 3 substitutions, 3 compounds
    PMID         : 24533298
    accession    : Int J Parasitol Drugs Drug Resist, PMIDs 24533298 and 25941623
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3862444/fullTextXML

Quirks that cost time once:
    The ONLY curated source in the map: Toxoplasma has no resistome to parse, so this is read
    out of papers one allele at a time and every row carries its paper, substitution, compound
    and how causality was shown. The bar is that the mutation was put BACK -- introduced into a
    clean background and shown to produce the resistance -- and that bar is why this holds one
    gene rather than ten. Nine genes were curated from the artemisinin and auranofin in-vitro-
    evolution studies and then DROPPED: the artemisinin table is titled 'Mutations found in
    candidate genes', and the auranofin paper names SOD2 as its likeliest locus and then reports
    that SOD2 L201P was not sufficient to confer resistance when introduced into wild-type
    parasites. Three well-known alleles are deliberately absent and the module says why: DHFR-TS
    pyrimethamine alleles (primary text pre-PMC and unreachable), DHODH N302S (primary not open
    access), and cytochrome b atovaquone alleles (mitochondrially encoded, so no row exists in a
    table of nuclear genes -- the eighteen nuclear cytochrome b hits are b-c1 subunits and would
    be the wrong gene). Absence here is ignorance rather than a negative result: nobody selected
    resistance in most genes, so only curated genes carry a value.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/curated_resistance_alleles.py
"""
from _common import run

KEY = "curated_resistance_alleles"

if __name__ == "__main__":
    run(KEY)
