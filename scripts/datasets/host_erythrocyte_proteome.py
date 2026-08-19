#!/usr/bin/env python3
"""Human red blood cell proteome, by fraction

Which human proteins are present in the cell the blood stage lives in

    level / kind : reference / proteomics
    provides     : rbc_membrane_psms, rbc_cytoplasm_psms
    coverage     : 5,264 human proteins: 4,777 membrane, 2,350 cytoplasmic
    PMID         : 41654503
    accession    : Sci Data 06792 Supplementary Table S1
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12992542/supplementaryFiles
    local path   : datasets/host/erythrocyte/41654503/41597_2026_6792_MOESM2_ESM.xlsx

Quirks that cost time once:
    The first HOST TISSUE reference in this project, and the thing instruction 39's host slots
    have been waiting for: a proteome of the cell the parasite lives in answers a question no
    pulldown can, since a pulldown says what a bait touched and this says what was there to
    touch. The two fractions are kept apart because they are different measurements -- a protein
    in the membrane extract is at the surface the merozoite invades through, one in the
    cytoplasm is in the haemoglobin around it. Self-validating as a fractionation should be:
    spectrin beta heads the membrane list and haemoglobin alpha the cytoplasmic one. Rows are
    HUMAN proteins keyed by UniProt accession and live in the host table, never in a parasite
    one; a row can name several genes (`HBA1; HBA2`) and the string is kept as given rather than
    one of them chosen.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/host_erythrocyte_proteome.py
"""
from _common import run

KEY = "host_erythrocyte_proteome"

if __name__ == "__main__":
    run(KEY)
