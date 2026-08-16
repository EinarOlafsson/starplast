# 41 — Fill as many slots as possible: Toxoplasma, Plasmodium, and host

**Status: open. Requested 2026-08-15.** This is the acquisition campaign. 39 defines the tables and
40 the window to check them in; this one puts data in.

## Where it stands

239 slots. Measured from the shipped node table and `slots.declared_columns`, not from the CSV:

| arm | filled | candidate only | empty |
|---|---|---|---|
| *T. gondii* | 59 | 42 | 11 |
| *Plasmodium* | 0 | 25 | 102 |
| host | — | — | not yet defined (see 39) |

**25 of the 25 Plasmodium candidates carry no data accession.** All 193 citations across both tables
resolve and every claimed title matches the real one, so nothing is fabricated — but a citation
nobody can download is not a filled slot, and the atlas deliberately refuses to let that read as
progress.

## The rules, before any of the work

1. **Never type a PMID or an accession.** Resolve it, from E-utilities or the repository's own API,
   and record what the query was. This rule has held for all 193 citations; it is the reason the
   tables can be trusted at all.
2. **A slot is filled when a column exists in the node table**, not when a paper is cited. `filled`
   comes from `slots.declared_columns`; nothing else may set it.
3. **Register every dataset** in `datasets.REGISTRY` with its columns, its accession, and
   `derived_from` where it is a computation over other columns. The leakage closure reads that field
   and cannot protect what is not declared.
4. **Do not merge species.** One table per parasite species, orthology as a bridge slot — see 39,
   including the rule that `target_family` closure spans species.
5. **Host genes never become rows in a parasite table.**

## Fix the candidate search first

The Plasmodium candidate list contains 16 citations that are not about a malaria parasite —
red-cell physiology, endothelial transporters, *Anopheles* immunity, essential-oil antimicrobial
screens. That is a fault in the query specified in the last handoff, not in the resolver: it used
`malaria[Title/Abstract]` as a standalone term, which sweeps in vector biology, host physiology and
natural-product pharmacology.

Replace it with a filter that requires all three:

* a **parasite species term** (`Plasmodium falciparum`, `P. berghei`, `P. vivax`, `P. knowlesi`,
  `P. yoelii`) — not the disease name;
* an **assay term** from `ASSAY_TERMS` for that slot's axis;
* an **accession pattern** in the abstract or the linked data availability
  (`GSE\d+`, `PRJ[EDN][ABN]\d+`, `PXD\d+`, `E-MTAB-\d+`, `SRP\d+`, `MTBLS\d+`), or a hit in the
  repository search rather than in PubMed at all.

Then re-run it over every Plasmodium slot and replace the current candidates wholesale. Keep the
16 off-target ones nowhere — a wrong candidate is worse than an empty slot, because an empty slot
is honest.

## Priorities

**Toxoplasma — highest value first, because the tables are already there.**

1. `relation` — 3 filled, **9 empty**, the weakest axis: crosslink MS, IP-MS parasite–parasite,
   proximity labelling, structural similarity, shared-complex membership. 39 makes pair-indexed data
   first-class, so this axis is about to carry much more weight than its coverage suggests.
2. `regulation` — **0 filled** across six slots that all have candidates: ChIP-seq per factor,
   histone marks, chromatin accessibility, m6A, splicing, RNA half-life. Candidates exist; nothing
   has been ingested.
3. The 20 citation-only slots that already carry an accession — these are the cheapest wins in the
   whole table and need only fetching and normalising.
4. The 11 empty slots: find candidates or record why none exists.

**Plasmodium — acquisition, per 39's tiers.** *P. falciparum* (3D7) and *P. berghei* (ANKA) first;
between them they carry nearly all the genome-scale data, and they are the pair where transfer is
most useful. Start with the axes where genome-scale data certainly exists: transcription (IDC time
course, single-cell atlas, gametocyte, liver stage, sporozoite), fitness (piggyBac saturation
mutagenesis for *falciparum*, PlasmoGEM barcoded knockouts for *berghei*), protein abundance,
phosphoproteome, chromatin. Resolve each through the repository, not from memory.

**Host — per 39.** Human and mouse first, *Anopheles* named by species. Tissue references for the
context nodes 39 lists: erythrocyte, hepatocyte, dermis, brain, fibroblast, monocyte, midgut,
salivary gland. Sources named there (HPA, GTEx, Tabula Muris, VectorBase, PRIDE, Ensembl/UniProt).

## Per-slot procedure

    resolve  -> candidate with an accession, recorded with its query
    fetch    -> raw file into the dataset archive, checksummed, release pinned
    map      -> to the table's gene identifiers; VEuPathDB names must be RESOLVED against the
                release index, never derived -- see the fetcher note in the skills corpus, where
                deriving-and-hoping silently skipped organisms and still printed DONE
    register -> datasets.REGISTRY entry: columns, accession, derived_from, note
    verify   -> coverage recomputed from the node table; the slot flips to filled only here
    regenerate -> scripts/generate_slot_table.py, and check the atlas

## Prune before scoping

The Plasmodium `fitness` axis has 33 empty slots because the stage expansion crossed 16 generic
fitness questions with the full stage list. Some of those combinations nobody has measured or
plausibly will — fitness in the ookinete, for one. Prune by hand before the acquisition target is
scoped against the count, or the campaign will be measured against a denominator that is partly
fictional.

## Acceptance

* Every newly filled slot has a `datasets.REGISTRY` entry and a column that
  `slots.declared_columns` returns.
* The column partition still holds per table: none claimed twice, none claimed by nobody.
* Every candidate in both tables resolves, its title matches, and it names the right organism —
  assert it in a test that can be re-run offline against a cached response.
* No slot is marked filled from a citation.
* Coverage back to 100% on the code, no `pragma`; it is currently 99%.
* The atlas regenerates and its counts change in the direction the work went.


## What the 2026-08-16 pass established, and the trap it found twice

**Filled: 65 of 105 Toxoplasma feature slots, up from 60.** Five came from PRIDE deposits verified
by reading the files: acetylation (3,921 genes), kinase substrates (2,575), proximity labelling
(1,734), S-nitrosylation (660), ubiquitination (128). Eight more were never empty -- they are
`unit=pair` slots filled by graph edges, and every census had been asking only about node columns.

**Of 102 datasets proposed from GEO, 8 survived adversarial verification.** That is the number to
plan around: a GEO query that returns a plausible hit for every slot is not finding data, it is
finding the nearest sequencing study, and roughly nine in ten of those do not survive being read.

### The trap, in its clearest form

`GSE313582` was proposed for SIX Toxoplasma slots: histone marks, interaction with host proteins,
strain variation, transcription in naive macrophage, transcription in IFN-gamma macrophage, and m6A.
Its sample columns are `RH_GCN5_KD_UT_1 ... RH_GCN5_KD_IAA` -- an auxin-induced GCN5 knockdown
RNA-seq in RH tachyzoites. It is none of those six. It is a chromatin-perturbation transcriptome,
which is a slot Toxoplasma already has filled, so ingesting it would have added nothing and claimed
six things.

The same shape appeared in the Plasmodium arm, where `GSE270631` passed a first reviewer for
"maximum observed across stages" and was refuted on reading the deposit: a PfSET10 conditional
knockout, two separately DESeq2-normalised matrices whose column totals differ by 48%, and
blood-stage samples only.

**So: no GEO candidate is ingested without reading its sample column names.** A title and an
organism are not enough -- both of these had the right organism and a plausible title.

### Where the remaining 40 actually stand

* **19 GEO-servable** (regulation 6, transcription 5, fitness 5, translation 3): candidates exist
  and are downloaded, but the six above were all one wrong dataset. Verify sample columns first.
* **3 metabolism** need MetaboLights, **2 immunity** need IEDB, **3 chemistry** are mostly
  supplementary tables in papers. None of these three repositories is wired.
* **4 PTM/abundance** (glycosylation, palmitoylation, crosslink MS, secretome) returned honestly
  zero from PRIDE for Toxoplasma. They may have no deposit at all, and an empty slot is the correct
  answer if so.


### Every Toxoplasma GEO candidate, read by its sample titles (2026-08-16)

Reading `!Sample_title` from each downloaded series matrix settles all nineteen at once, and it is
the check that should have been run before any of them were proposed. Two datasets account for ten
of the assignments:

| dataset | what its samples actually are | slots it was proposed for |
|---|---|---|
| GSE313582 | `Illumina-GCN5b-KD-UT/IAA` -- GCN5b knockdown RNA-seq | m6A, RNA half-life, TF binding, histone marks, host interaction, strain variation, IFN-gamma macrophage |
| GSE313048 | `ATACseq-GCN5b-KD-*` -- ATAC-seq of the same knockdown | glycosylation, palmitoylation, protein turnover |
| GSE287334 | `Pru-MORC-KD-*` -- a MORC knockdown | T-cell epitope content |
| GSE300509 | growth on different host cell types | seroreactivity |
| GSE175919 | `DMSO` vs `MC1742` -- transcriptional response to an HDAC inhibitor | drug sensitivity per gene, resistance-conferring mutation |

None of those fits the slot it was proposed for. An ATAC-seq is not a palmitoylome, and a knockdown
transcriptome is not a measure of strain variation.

**But three are real data in the wrong slot, and should be REASSIGNED rather than discarded:**

* `GSE132237` -- "CRISPR-Cas9 screens to identify regulators of differentiation", samples
  `L1 input library, L1 p4, L1 p5, L1 p6`. That is a genuine pooled screen with passages, and it was
  proposed for `essentiality in a second background`, which is where it belongs. Verify and ingest.
* `GSE245775` -- `parental_tachyzoite_RIBOseq_rep1...`. Ribosome profiling, proposed for
  `stage-conversion phenotype`, which it is not. It belongs in `translation · bradyzoite` or
  `translation · per cell-cycle phase`, both of which are empty.
* `GSE253885` -- "In vivo CRISPR screens for hyperLOPIT-unassigned proteins", samples
  `unassigned_2_lib, _P4, _WT_01`. A real in vivo screen, proposed for `secretome / excreted`, which
  it is not. It belongs in a fitness slot.

So of nineteen GEO-servable Toxoplasma slots, the honest yield is **three datasets, none of them in
the slot it was proposed for**. That is the shape of the problem: the resolver finds real
Toxoplasma data and assigns it by keyword proximity, and the keyword is nearly always wrong.

**The fix for the next pass** is to propose from the sample titles rather than from the study
abstract: `esummary` gives them cheaply, and "does this series contain samples of the kind this slot
asks about" is answerable from `!Sample_title` alone in the great majority of these cases.

### The registry catches a second kind of trap: the same study entering twice

`PXD017032` was ingested as `n_kinase_substrate_sites` and had to be removed. It is Wang 2022's
**sporulated-oocyst vs tachyzoite phosphoproteome**, already in the registry as `phospho_quantitative`
from the authors' own supplement. Counted a second time from the raw deposit it produced a column
correlating at **rho = 0.729** with `phospho_sites_measured` over 1,592 shared genes, and it sat in a
slot whose own citations name CDPK1 and CDPK7 substrate deposits — neither of which it is.

Two separate faults, and the second is the one worth the note: the first (wrong slot) is the keyword
trap already documented above, appearing this time in work done here rather than in the resolver. The
second is new — **a dataset already in the map arriving again through a different door**, as counts
from the raw deposit rather than as the published table. An embedding would then weight that one
experiment twice.

`datasets.REGISTRY` is what caught it, and only because the accession was written down. So:

> **Check a new accession against `{d.accession for d in datasets.REGISTRY}` before ingesting it.**
> The registry's invariant tests do the rest — level must be one of the five on-disk names, every
> entry needs a download URL or a declared derivation, and every entry needs a generated fetch script.

Filling the five PRIDE deposits and the differentiation screen into the registry is what surfaced
this; the five had been merged into the node table without registry entries, and the acceptance rule
requiring one exists precisely so that cannot happen quietly.

### The PRIDE round, and what it cost to verify six deposits (2026-08-16)

Six deposits examined for six empty Toxoplasma slots. **One was ingested.** The other five each
failed in a different way, and the ways are worth more than the yield:

| deposit | proposed for | verdict |
|---|---|---|
| PXD031526 | lactylation | **INGESTED**, 515 genes. `La (K)Sites.txt` inside a readable RAR |
| PXD008574 | exposure to host cytosol | refused — label swaps do not reproduce each other |
| PXD044588 | secretome / excreted | refused — the mzid cannot distinguish secreted from detected |
| PXD056853 | glycosylation | refused — the deposit is ***Dictyostelium discoideum*** |
| PXD028969 | secretome / excreted | refused — the deposit is ***Cryptosporidium*** and human |
| PXD033642 | protein turnover | real Toxoplasma CETSA, but only identifications are published |

**Two of the candidates in the Toxoplasma table are not Toxoplasma deposits at all.** The instruction
above records 16 off-target citations in the *Plasmodium* arm and treats it as a fault in that query;
it is not confined to that arm. Every PRIDE candidate in both tables needs its `organisms` field
checked against the slot's organism, and that check costs one API call each.

**PXD008574 is the instructive refusal.** It is a genuine N-terminomics study of WT versus ASP5
knockout, SILAC with a proper label swap, and every metadata field says so. Orienting both swaps to
WT/KO and correlating them gives r = +0.09 (n = 70), +0.44 (n = 342) and +0.007 (n = 53). A label
swap that does not reproduce itself is measuring noise, and a column built from it would rank genes
by noise while looking like a measurement. **No metadata field could have caught this — only
computing the number and checking it against itself.**

Note also that it is N-terminomics, so even had it verified it would not fill `exposure to host
cytosol`: ASP5 cleaves in the Golgi, before export. It would have needed a new
`proteolytic processing / N-terminome` slot on the PTM axis, which remains undefined and should be.

### A silent-drop bug the refuse-if-empty guard caught

`proteomics.column_for` joined deposit accessions to the node index as strings. The lactylation
deposit is keyed entirely on `TGGT1_` and reported **0 of its 524 genes** — indistinguishable from a
study with no coverage. It now routes through the identity layer, as `build_graph` always has.

The guard that surfaced it is `scripts/add_verified_columns.py` refusing to write when any column
comes back empty. That refusal has now paid for itself twice; keep it.
