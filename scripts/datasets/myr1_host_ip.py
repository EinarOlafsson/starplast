#!/usr/bin/env python3
"""MYR1 host interactome (bridge)

Host proteins co-immunoprecipitating with the parasite protein MYR1

    level / kind : post_translation / IPMS
    provides     : bridge:host
    coverage     : 219 host proteins, 1 parasite gene
    PMID         : 32075880
    accession    : PXD016383
    url          : https://www.ebi.ac.uk/pride/archive/projects/PXD016383
    local path   : starplast/data/host_bridges.parquet

Quirks that cost time once:
    The first BRIDGE in the project: pairs whose two ends are in different tables, which
    `graph.npz` cannot hold because its edges are index pairs into the parasite table. Verified
    twice over. The parasite side: MYR1 is rank 1 of 325 in its own pulldown at +34.4, with MYR3
    at 47 and GRA44, GRA7 and GRA9 in the top ten. The host side against independent literature:
    PDCD6, which is ALG-2, is rank 1 of 219, and a 2026 paper reports Toxoplasma GRA8 engaging
    host ALG-2 at the vacuole; its partner ALIX is rank 12 and VPS28 rank 62, so ESCRT as a
    class sits at p = 0.006 against the rest of the host proteins. Two filters do the work and
    both were got wrong first: a group is a contaminant group if ANY entry in it is one --
    MaxQuant prefixes only on the leading entry, so keratin hides mid-group and is otherwise the
    four most enriched host proteins -- and two unique peptides are required in BOTH bait
    replicates, which takes 674 host groups to 219.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/myr1_host_ip.py
"""
from _common import run

KEY = "myr1_host_ip"

if __name__ == "__main__":
    run(KEY)
