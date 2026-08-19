#!/usr/bin/env python3
"""Plasmodium co-immunoprecipitation interactome (EPIC)

Which parasite proteins came down with each tagged bait, against its own control

    level / kind : post_translation / IP-MS
    provides     : n_ip_ms_partners
    coverage     : 98 edges over 65 partners and 3 baits
    PMID         : 28691708
    accession    : Nat Commun 16044 Supplementary Tables 1-5
    url          : https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5508133/supplementaryFiles
    local path   : datasets/reference/plasmodb/ip_ms/28691708/ncomms16044-s1.pdf

Quirks that cost time once:
    Published as PAGES OF A PDF rather than as a table, and its captions are wrong: the table
    headed `parasite interacting proteins` holds human ones and the table headed `human` holds
    parasite ones. So no caption is read for what a table contains -- the bait comes from the
    column header, which names the pulldown and its control, and the organism from the
    identifier space. The check that licenses the assignment is that a BAIT MUST TOP ITS OWN
    TABLE: PV1 leads the PV1 pulldown at 97 and 120 spectra, PV2 leads PV2's, EXP3 leads EXP3's,
    and a page whose first row is not its bait is a continuation rather than a new table. Read
    two ways before it was trusted -- pypdf and poppler independently give 38, 44, 13, 11 and 37
    rows on the five table pages. The count layout is the trap: eight numbers are TWO
    experiments of four (bait, bait, control, control), so splitting them down the middle
    compares experiment 1 with experiment 2 and reports an enriched partner with 145 spectra in
    the untagged line. Caught because that is impossible; read correctly, no pair has more
    spectra in its control than in its bait and the median control is 0. Two things are
    deliberately not carried: the 38 rows of the PfEMP1B pulldown, whose bait is a var-gene
    transgene with no accession in the table, and one PV2 row (PIESP2) that carries seven counts
    instead of eight, since the missing number could be either arm. Degree is missing outside
    the experiment -- four pulldowns are not a survey, and a zero would say `nothing binds this`
    about a protein nobody tested.

Fetches the source, reads it, resolves its accessions to current ToxoDB ME49, and reports the
coverage that resolution achieves. The shipped columns are assembled by
`starplast.build_graph.load_nodes()`; this script is the per-dataset view of the same source, for
inspecting or re-fetching one dataset without running the whole build.

Run:  python scripts/datasets/pf_ip_ms.py
"""
from _common import run

KEY = "pf_ip_ms"

if __name__ == "__main__":
    run(KEY)
