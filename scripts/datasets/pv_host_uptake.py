#!/usr/bin/env python3
"""Host proteins at the vacuole

How enriched a host protein is at the parasitophorous vacuole

    level / kind : post_translation / proteomics
    provides     : pv_enrichment_log2
    coverage     : 12 host proteins
    PMID         : 34898650
    accession    : PLoS Pathogens 1010138 supplementary table
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC8700025/supplementaryFiles
    local path   : starplast/data/host_proteins.parquet

Quirks that cost time once:
    A property OF a host protein rather than a bridge, because the bait is the compartment and
    not a named parasite gene -- a bridge needs a parasite gene at one end. Averaged over three
    infection contexts: tachyzoite-infected fibroblast, bradyzoite-infected fibroblast and
    neuron. The sheet lists parasite and host proteins together, which is how the authors show
    the experiment worked -- the dense granule proteins top it -- and only the host rows are
    kept. Top of those: PDCD6/ALG-2 at +5.90, VPS37C at +4.47, then CHMP4B, PEF1 and VPS28. That
    is the FOURTH independent dataset in this map to put ALG-2 at the host-parasite interface,
    after the MYR1, EAF1 and GRA35 pulldowns.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pv_host_uptake.py
"""
from _common import run

KEY = "pv_host_uptake"

if __name__ == "__main__":
    run(KEY)
