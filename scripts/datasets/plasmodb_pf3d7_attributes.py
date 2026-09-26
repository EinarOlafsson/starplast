#!/usr/bin/env python3
"""Plasmodium falciparum 3D7 gene attributes

The second species: gene structure and the protein's sequence properties

    level / kind : reference / annotation
    provides     : length, molecular_weight, isoelectric_point, transcript_length, exon_count, n_tm, is_tm, has_signal_peptide
    coverage     : 5,720 P. falciparum genes
    accession    : PlasmoDB GenesByTaxon attributesTabular
    url          : https://plasmodb.org/plasmo/service/record-types/transcript/searches/GenesByTaxon/reports/attributesTabular
    local path   : starplast/data/pf_nodes.parquet

Quirks that cost time once:
    The FIRST parasite table beyond Toxoplasma, and it is its own table -- nothing is merged,
    because PF3D7 and TGME49 identifiers do not align and neither does the data behind them, so
    a merged frame would encode which species was convenient to work on as though it were
    biology. Column names deliberately match the Toxoplasma table where the quantity is the
    same, so a slot pattern reads the same on both arms, but the Pf slots name their own
    patterns rather than inheriting Toxoplasma's -- inheriting them would have claimed
    falciparum slots with gondii numbers. The report is served by the TRANSCRIPT record type, so
    5,791 rows describe 5,720 genes and the rows are collapsed on the longest transcript; taking
    the row count at face value double-weights 71 genes. One download carries five unrelated
    sources, which were one registry entry until 2026-09-25: its columns' median association
    with each other was 0.13, with 70% of pairs under 0.2 -- the loosest 'experiment' in either
    table -- so holding out the piggyBac screen also removed codon usage, orthology and SNPs.
    They are five entries now, and the closure's 'shared experiment' rule means what it says.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/plasmodb_pf3d7_attributes.py
"""
from _common import run

KEY = "plasmodb_pf3d7_attributes"

if __name__ == "__main__":
    run(KEY)
