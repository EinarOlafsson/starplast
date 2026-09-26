#!/usr/bin/env python3
"""Male and female fertility, transferred from P. berghei

Whether a knockout loses male or female fertility, measured in the rodent parasite

    level / kind : DNA / CRISPR_screen
    provides     : fertility_female, fertility_male
    coverage     : 1,121 orthologs (19.6%)
    citation     : Sayers C et al., Systematic screens for fertility genes essential for malaria parasite transmission reveal conserved aspects of sex in a divergent eukaryote. Cell Syst 2024;15:1075-1091.e6
    PMID         : 39541984
    accession    : Cell Systems Table S1 sheet A
    url          : https://www.cell.com/cms/10.1016/j.cels.2024.10.008/attachment/mmc2.xlsx/mmc2.xlsx
    local path   : datasets/DNA/transmission_screen/39541984/mmc2.xlsx

Quirks that cost time once:
    Barcoded P. berghei knockouts crossed and scored through the mosquito, so the sexes are
    measured SEPARATELY on the same mutant -- which is the point, and which the controls
    confirm: HAP2, P48/45, P230, CDPK4 and MAPK2 fail in males only, P47, NEK4 and DMC1 in
    females only, and the redundant P25/P28 pair is correctly not called. The union of reduced
    mutants is 348, the paper's own number. The table names the falciparum ortholog itself, so
    no orthology call is made here; a gene named by two mutants, or a mutant covering two genes,
    is dropped rather than averaged. A TRANSFER, grouped with orthology for leakage: an
    ortholog's phenotype is not an independent measurement of this gene.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_pb_fertility_transfer.py
"""
from _common import run

KEY = "pf_pb_fertility_transfer"

if __name__ == "__main__":
    run(KEY)
