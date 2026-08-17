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

### The ToxoDB round: eight datasets, five kept, and three kinds of failure (2026-08-16)

ToxoDB does not only host datasets — it runs them through its own pipelines and exposes the per-gene
result through searches whose report can be asked for the value. That is the richest single source
found so far, and it is the only route to several of these: the epitope mapping is ToxoDB's join of
IEDB against the ME49 proteome; the Foe palmitome's own paper is not open access and PMC blocks
automated download of its supplements.

| dataset | slot | verdict |
|---|---|---|
| Foe palmitome | palmitoylation | **kept**, 470 genes |
| IEDB epitopes | T-cell epitope content | **kept**, 221 genes |
| Ramakrishnan enteroepithelial | transcription · in vivo enteric | **kept**, 7,739 genes |
| Hakimi/Ali H4 acetylation | chromatin state · histone marks | **kept**, 7,515 genes |
| Saeij 29 strains (ME49 arm) | transcription · in naive macrophage | **kept**, 8,140 genes |
| Stuart/Ralph nanopore | splicing / isoform use | **kept**, 798 genes |
| Einstein H3K4me1 | chromatin state · histone marks | refused — backwards |
| Ramirez-Flores vesicles | secretome / excreted | refused — measures something else |
| Gregory sense/antisense | noncoding and antisense | refused — does not reproduce |

**The three failures are three different failures, and none is visible in metadata.**

1. **Backwards.** H3K4me1's "marked" genes have LESS accessible promoters than unmarked ones
   (−0.51 against +0.39, p = 3e-50). H4 acetylation from the same site, same assay type and same
   query shape gives +0.43 with ATAC, which is what says the refusal is about the data.
2. **Measures something else.** In the vesicle proteome the dense granule proteins sit at −0.85 and
   the microneme proteins at −3.72 — depleted from the vesicle fraction — while ribosomal proteins
   are the most enriched thing in it. Probably a correct measurement of vesicle partitioning; still
   not an answer to what the parasite secretes.
3. **Does not reproduce.** The sense/antisense coupling was independent of expression, which looked
   like a good sign, and the ME49 and GT1 time courses of the same analysis share 12 of their top
   200 genes where chance gives 28.

**Every one of these was caught by comparing the candidate column against something else** — another
column already in the map, or the same measurement computed a second way. That is the check to
budget for; it is not expensive and nothing else substitutes for it.

### A fact about the API that is worth not rediscovering

`fold_change_avg` is **comp/ref**, not ref/comp, and nothing in the API says which. Established
against a case with one possible answer: tachyzoites as reference against tissue cysts puts BAG1 at
+21.6 and LDH2 at +16.1, both bradyzoite-specific, and SAG1 at −14.3. Check it against a known
comparison whenever a new fold-change search is added — getting it backwards raises nothing, it
silently inverts a column.

### Paper supplements beat repositories when what you want is a study's analysis (2026-08-16)

Four slots filled from open-access supplementary tables through EuropePMC, after repository sweeps
had reported all four as unserved:

| slot | source | verification |
|---|---|---|
| fitness · oxidative stress | PMC8216390 Data Sheet 1 | catalase at −6.15, the extreme of the screen |
| cyst wall composition | PMC7002340 Data Set S1 | MAG1 and MAG2 top the interactome |
| target engagement / thermal shift | PMC9436416 supp. file 3 | CAM1 and CAM2 at the 98th percentile |
| **N-myristoylation** (new slot) | PMC7373427 supp. file 4 | **65 of 65 substrates start with glycine** |

**Two of these had been examined and refused earlier the same day, from their PRIDE deposits.**
`PXD033642` is the thermal-profiling study and publishes only identifications; the ED scores are in
the paper. `PXD019677` is the myristoylation study and ships MaxQuant archives of 250–340 MB apiece;
the answer is a 65-row table in supplementary file 4. A third, glycosylation, was reported here as
"returned honestly zero from PRIDE" — the deposit exists and is indexed under *glycoprotein*, not
*glycosylation*.

> **When the wanted quantity is a study's ANALYSIS — a score, a fit, a curated list — the paper's
> supplement is the primary source and the repository is the fallback.** Repositories hold the
> evidence a study was built from, not the conclusions it reached. Search PubMed, resolve the PMCID,
> and pull `europepmc/webservices/rest/{pmcid}/supplementaryFiles`, which returns a zip.

### A third species turned up in the Toxoplasma candidate table

`PMID 38747635`, proposed for `exposure to host cytosol`, is a ***Theileria annulata*** paper — its
supplement is keyed on `TaC12_001700` and `TA19380`. With PXD056853 (*Dictyostelium*) and PXD028969
(*Cryptosporidium*), that is three wrong-organism candidates found by opening the file. The organism
field must be checked for every candidate in both tables, from the record and not from the title.

Two other named candidates for `phosphorylation · kinase-substrate` fail differently and are worth
recording so nobody re-fetches them: `PXD019677` is the myristoylation study and not CDPK1
substrates, and `PXD019655` is CDPK7 but ships only Proteome Discoverer `.pdResult` and `.msf` files
of 3–8 GB, which nothing outside that program reads.

## Where the Toxoplasma arm actually stands, slot by slot (2026-08-16, end of the campaign)

**97 of 105 gene slots filled**, plus 8 of 10 pair slots and 0 of 3 metabolite slots. The eight
empty gene slots are listed below with what was searched, because "empty" without that is
indistinguishable from "nobody looked".

| slot | searched | finding |
|---|---|---|
| transcription · in IFN-gamma macrophage | GEO with `"Toxoplasma gondii"[Organism]` + interferon/macrophage | **0 hits.** The one promising paper (PMID 34928716) is *Cryptosporidium* |
| translation · per cell-cycle phase | GEO, PubMed | no cell-cycle-resolved ribosome profiling exists |
| protein turnover | PRIDE, PubMed for pulse-SILAC / pulse-chase / half-life | none for Toxoplasma. Thermal stability is now its own slot and is filled |
| drug sensitivity | GEO, PRIDE, PMC supplements | no genome-wide chemogenomic screen published |
| RNA modification · m6A / 5mC | GEO | no MeRIP-seq. GSE178355 is depletion RNA-seq of METTL3/WTAP/YTH1 — which transcripts *depend on* m6A, not which *carry* it |
| resistance-conferring mutation | PubMed, PMC supplements | only single-gene selections (ROP1, PRELID); nothing genome-wide and per-gene |
| invasion and egress phenotype | GEO organism-filtered, PubMed | 73 GEO hits, all transcriptomes; no per-gene imaging screen with a table |
| fitness · in vivo gut | GEO, PubMed | the in vivo screens are mouse peritoneum/lung/liver/spleen. No feline enteric screen has been run |

### Two slots were removed or moved rather than left empty

* `drug sensitivity per gene` (chemistry) **merged into** `drug sensitivity` (fitness). They asked
  one question — does disrupting this gene change survival under a compound — on two axes, and
  whichever got data first would leave the other permanently and misleadingly empty.
* `metabolite levels`, `metabolic flux`, `lipid composition` are now `unit=metabolite`; they were
  never answerable from a table whose rows are genes.

### What worked, in the order the yield came

1. **Paper supplements through EuropePMC**, once it was clear that repositories hold the evidence a
   study was built from and not the conclusions it reached. Seven fills.
2. **ToxoDB's own search reports** — it runs datasets through its pipelines and will return the
   per-gene value. Six fills.
3. **Re-reading data already on disk.** Four fills, including a screen sitting in a directory named
   for the slot it filled.
4. **Re-asking a refused source a better question.** Two fills: the antisense LEVEL after the
   antisense CHANGE was refused, and the secretome PARTITION after "enriched in vesicles" was.

### And the counts that say what verification is worth

Of everything examined, **six sources were refused after being read**: H3K4me1 (backwards against
ATAC), the sense/antisense change (did not reproduce), PXD008574 (label swaps disagreed), the vesicle
proteome under its first framing, PXD017032 (already in the map under another name), and the
GSE302108 MPRA (keyed on synthetic sequences, not genes). **Four candidates were the wrong organism**
— *Dictyostelium*, *Cryptosporidium*, *Theileria*, and *Cryptosporidium* again — each found by
opening the file rather than by reading the title.

### The seven that remain, after a second and third search pass

The audit above was written after one pass and was wrong within the hour: m6A fell to a paper
supplement. So the remaining seven were searched again, and this records what the later passes
added rather than repeating the earlier table.

| slot | second/third pass | still absent because |
|---|---|---|
| transcription · in IFN-gamma macrophage | PubMed all-fields, GEO organism-filtered | every hit is HOST response. The parasite's own transcriptome in an activated macrophage has not been sequenced |
| translation · per cell-cycle phase | PubMed for polysome/ribosome + cell cycle | the 2024 single-cell atlas (PMC11358496) gives RNA peak time and ATAC per phase, and no translation |
| protein turnover | PubMed for proteasome/degradation proteomes | the 2026 iron paper (PMC13170339) is translational remodelling, not turnover rates |
| drug sensitivity | GEO, PRIDE, PMC | still no genome-wide chemogenomic screen |
| resistance-conferring mutation | PubMed for in vitro evolution + sequencing | single-gene selections only (ROP1, PRELID, auranofin) |
| invasion and egress phenotype | GEO organism-filtered, PubMed for phenotypic/secondary screens | the splitCas9 phenotypic screen (PMID 35538310) has no PMC record and no accessible table |
| fitness · in vivo gut | GEO, PubMed | the in vivo screens are mouse; PMC13230962 is a real 309-gene in vivo screen but mouse and reference-strain, so it is a second dataset for a filled slot rather than this one |

**Two of these are one dataset away and worth watching**: a splitCas9 phenotypic screen exists for
invasion/egress but is not deposited anywhere reachable, and the in vivo screens would answer the gut
slot if anyone ran one in a cat.

**Two are unlikely to move**: no genome-wide chemogenomic screen and no pulse-SILAC proteome has been
published for this organism at all.

### What the campaign changed about the catalog itself

Five slots were re-specified rather than filled, and each was a defect the emptiness had been hiding:

* three metabolism slots asked gene-indexed questions about metabolites (`unit=metabolite` now);
* `predicted complex membership` was written around AlphaFold-Multimer, a method with no Toxoplasma
  dataset, when the question has a measured answer;
* `protein turnover / stability` bundled two different properties behind a slash;
* `drug sensitivity per gene` and `drug sensitivity` were one question on two axes;
* and two new PTM categories — arginine methylation, N-myristoylation — plus thermal stability were
  missing from the tree entirely while their data sat published.

A slot that stays empty is worth re-reading as a question about the SLOT, not only about the data.

## Final state of the Toxoplasma arm

| unit | filled | of | what the remainder needs |
|---|---|---|---|
| gene | 98 | 105 | seven slots, searched three times each (table above) |
| pair | 9 | 10 | `interaction · with host proteins` needs the host tables — instruction 39 |
| metabolite | 2 | 3 | `lipid composition` needs a parasite lipidome |
| **total** | **109** | **118** | |

### `lipid composition`, searched and recorded

Three candidate routes, none of which is a Toxoplasma membrane lipidome:

* **the iron metabolomics already in the table** carries 61 lipid species, and filling this slot from
  them would be the levels column claimed by two slots. Refused as leakage, not for lack of data.
* **PMC12913473** (2026, extracellular vesicles) measures four HOST cell lines — fibroblast, IPEC,
  myotube, Vero — and EVs from infected cultures. Host membrane, not parasite.
* **PMC10033509** (2023, ester- and ether-linked phosphatidylethanolamine) publishes figure source
  data — plaque sizes, vacuole counts, parasite lengths — and not a table keyed by lipid species.

The 2007 lipidomic analysis of Toxoplasma (PMID 17988103) predates supplementary data tables.

### The metabolite table exists now, and the host tables can follow the same shape

Filling two metabolism slots required the first non-gene table in the project, and the plumbing it
needed generalises: `slots.resolve` takes the unit the caller holds, `slots.is_filled` looks in the
table matching a slot's unit, and the slot table grades against that table's own denominator.
`UNIT_TABLES` is where a unit declares it has rows.

**That is most of what instruction 39 needs for the host tables.** What remains there is the bridge —
a pair slot whose two ends are in different tables — which is exactly what
`interaction · with host proteins` is waiting for, and which no amount of downloading will supply.
