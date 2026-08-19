#!/usr/bin/env python3
"""Plasmodium extracellular vesicle proteome

Parasite proteins found in extracellular vesicles, and how many preparations found them

    level / kind : post_translation / proteomics
    provides     : ev_studies
    coverage     : 184 proteins, 53 of them in both preparations
    PMID         : 28944300
    accession    : Wellcome Open Res 11910 S2 (PRIDE PXD006925)
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5583745/supplementaryFiles
    local path   : datasets/post_translation/secretome/28944300/5c097a1c-efe5-4ed8-b97b-f9ba656268a6.xlsx

Quirks that cost time once:
    The columns say `extracellular vesicle` and not `secreted`, because the second word would
    assert a route this measurement does not establish. The sheet read is the paper's own
    compilation: the union of two independent EV preparations with a membership column each, so
    'how many studies saw this' is a fact in the file rather than a join. Its other columns are
    seroreactivity and antibody-array results from unrelated studies -- claims about immunity,
    not about vesicles -- and are deliberately not read. The deposit itself is raw-only (24 RAW
    files, no RESULT), which is the usual shape here: per-protein numbers come from the
    supplement. Checked in the direction a secretome should go -- 20.1% carry a signal peptide
    against 10.2% of the proteome (p = 8e-05), exported proteins run 6.0% against 3.3% (p =
    0.06, same direction and not significant at 184 genes) -- and RESA, KAHRP, MSP1 and Ag332
    are all present. The confound that must travel with the column is abundance: EV genes have a
    median blood-stage expression of 71.3 against 12.2 for the rest (p = 8e-31), so this is what
    mass spectrometry found in a vesicle preparation and not a list of what the parasite
    exports. Same caveat as hyperLOPIT assignment on the other arm, and for the same reason.
    Absence is unknown and stays missing. A companion boolean completed with False would have
    read as 5,720 genes tested and 5,536 negative, and graded the slot A at 100% for an
    experiment that identified 184 proteins -- so this one ships a single column whose presence
    is the evidence.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_secretome.py
"""
from _common import run

KEY = "pf_secretome"

if __name__ == "__main__":
    run(KEY)
