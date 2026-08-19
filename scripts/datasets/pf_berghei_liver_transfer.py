#!/usr/bin/env python3
"""P. berghei liver-stage fitness, transferred to falciparum

How a berghei knockout fares through the liver, carried onto its falciparum ortholog

    level / kind : DNA / CRISPR_screen
    provides     : pb_transferred_liver_log2fc, pb_transferred_liver_reduced
    coverage     : 754 falciparum genes; 180 reduced
    PMID         : 31730853
    accession    : Cell 2019 Table S2 (PlasmoGEM liver stage)
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6904910/supplementaryFiles
    local path   : datasets/reference/plasmodb/pb_transfer/31730853/mmc2.xlsx

Quirks that cost time once:
    Two borrowings, both deliberate. The falciparum id comes from the BLOOD-STAGE screen's table
    -- the same consortium pairing the same mutants -- rather than from an orthology this
    project derived. And the value is the authors' own blood-stage-CORRECTED figure for the
    salivary-gland-to-blood transition, because that transition ends in blood and the
    uncorrected column would call every blood-essential gene liver-essential. The sheet has two
    header rows and repeats `Log2-FC / SD / Power` per transition, so columns are read by
    position and the loader refuses the file if the transition is not where it expects it. 507
    genes are dropped for `no power`: too few barcodes to say anything, which is not a
    measurement of no effect. Validated on the genes the field would name -- LISP1, the
    UIS/ETRAMP early transcribed membrane proteins and perforin-like protein 1 all come out
    reduced, which is the textbook set for liver development and hepatocyte egress. The same
    file's two MOSQUITO transitions are NOT shipped: the markers available to check them (P25,
    P28, SOAP, chitinase) are the redundant ones, so nothing in the data confirms the direction,
    and a transmission slot filled on an unchecked axis is what this campaign refuses.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_berghei_liver_transfer.py
"""
from _common import run

KEY = "pf_berghei_liver_transfer"

if __name__ == "__main__":
    run(KEY)
