#!/usr/bin/env python3
"""Metabolome and isotope labelling under iron deprivation

Steady-state metabolite levels and the fraction labelled from glucose or glutamine

    level / kind : reference / metabolomics
    provides     : metabolite_level_log2fc_iron_depleted, metabolite_level_padj, labelled_fraction_glucose, labelled_fraction_glutamine
    coverage     : 1,102 metabolites
    PMID         : 41925342
    accession    : mBio 03788-25 Tables S3 and S5
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13170339/supplementaryFiles
    local path   : starplast/data/metabolites.parquet

Quirks that cost time once:
    Rows are COMPOUNDS, not genes -- the first table in the project that is not the node table,
    and the shape instruction 39 describes for host tables. The flux column is one minus the
    unlabelled isotopologue, which is the only labelling readout comparable between molecules of
    different carbon number. The archive's two labelling sheets differ in whether they carry a
    title row, and assuming they did not silently dropped the glucose arm: the sheet read fine
    and had no column called `Metabolite`. Joining a second study means matching compound NAMES,
    which is lossy; that cost is unpaid with one study and is the first thing to fix when a
    second arrives.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/metabolome_iron.py
"""
from _common import run

KEY = "metabolome_iron"

if __name__ == "__main__":
    run(KEY)
