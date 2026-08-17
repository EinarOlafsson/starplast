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

### The host bridge fell after being written off twice — and what that cost the audit above

`interaction · with host proteins` is filled. It was declared blocked twice in the same session:
once as "structural, not data" and once as "the host side does not verify". Both were wrong, and the
second was wrong for an avoidable reason — I read the top of the host list, saw tubulin and filamin
and HSP90, and stopped before asking what the FIRST entry was. It was PDCD6, which is ALG-2, which a
2026 paper independently reports Toxoplasma engaging at the vacuole. ALIX follows at rank 12 and
VPS28 at 62; ESCRT as a class sits at p = 0.006.

So the audit above should be read as what was found, not as what exists. Two of its entries have
already been overturned — m6A and this one.

### The last two, checked again after that

* **`invasion and egress phenotype`.** A pooled image-based CRISPR screen for exactly this exists:
  `PPR1275982`, 2026, *A pooled image-based CRISPR screen identifies EAF1 as a T. gondii modulator*.
  It is a preprint with `isOpenAccess=N` and `hasSuppl=N` — no reachable table. Together with the
  splitCas9 screen (PMID 35538310, no PMC record) that is two screens that would fill this slot and
  neither is deposited anywhere a reader can get at.
* **`transcription · in IFN-gamma macrophage`.** `GSE229505` is a dual perturb-seq in IFN-gamma
  stimulated cells and looks like the answer until you read it: it publishes raw single-cell
  matrices, and its readout is the HOST transcriptome per parasite knockout, which this map already
  carries as `hosttx_signature_`. The parasite's own transcriptome in an activated macrophage is
  still unsequenced.

**Eight slots remain.** Four need an experiment nobody has run in this organism, two need a
deposited screen that exists only as an inaccessible preprint, one needs a parasite lipidome, and one
needs cell-cycle-resolved ribosome profiling. Given the record above, that list should be re-tested
rather than trusted — the two overturned entries were both overturned within hours.

## The ESCRT-recruitment corpus, downloaded 2026-08-16

Eight sources for host-ESCRT engagement at the Toxoplasma vacuole, in
`datasets/quarantine/2026_08_16_escrt/`. None is ingested yet — this is acquisition.

| source | what it holds | usable per-gene table? |
|---|---|---|
| **PXD080696** | EAF1 and GRA35 affinity purification, DIA. `EAF1_II_vs_WT_III.csv` and `combined_results.csv` | **yes** — differential abundance vs wild type |
| **PMC9426488** | GRA64 IP in tachyzoite and bradyzoite, plus TurboID in HFF and neuron | **yes** — four summary sheets |
| **PMC8700025** | Toxoplasma exploiting host ESCRT; results summary with log2 fold change vs control | **yes** |
| **PMC11559087** | Ulp1 TurboID with an enrichment-fold column, 162 filtered proteins | **yes** |
| **PMC11377541** | DCS1, DCS2, PP2A-B2, PP2A-C2 co-IPs — ESCRT-adjacent abscission machinery | **yes**, four co-IP datasets |
| PXD051495 | the same abscission study's raw MGF and Mascot `.dat` | no processed table |
| PXD024491 | mitochondrial SPOT shedding, LFQ and TMT designs in `.7z` | designs only, data compressed |
| PMC12453208, PMC12669045 | non-canonical ESCRT activation; VIP1 at the PVM | figure source data only |

**Five carry a per-gene table and are ready to read.** The obvious next step is a
`host ESCRT engagement` layer built from them together rather than one at a time: five independent
baits converging on the same host machinery is a far stronger claim than any one IP, and it is the
kind of agreement that the MYR1 bridge and the spaCR screen already show separately.

Note the convergence already in the map: the MYR1 host bridge is verified by ESCRT topping it
(PDCD6/ALG-2 rank 1, ALIX 12, VPS28 62, p = 0.006), and the spaCR screen's top hits are EAF1 and
GRA14. Three unrelated datasets, one biology.

## Fourth pass, 2026-08-16: the metabolite axis closes

**112 of 119.** `lipid composition` is filled and the metabolite unit is complete at 3 of 3. Seven
gene slots remain.

### What filled it, and why three earlier passes missed it

The slot was never short of Toxoplasma lipidomes — there are several, and passes one through three
found them. It was short of a lipidome *of the parasite*. A lipidome taken from an infected culture
is mostly host cell, and a lipid measured there is not thereby a parasite lipid, so each candidate
was correctly refused on the same ground and the ground was recorded as though it applied to the
whole question. It applied to those experiments.

The fill (PMID 41716462) works because the vesicles were collected from tachyzoites **after egress,
in host-cell-free medium**, and because the study grew the parasite in four host backgrounds. The
host cells differ from each other in 1,018–1,362 lipid species; the vesicles the same parasite
released in them differ in 0–4. Composition set by the host would move when the host is replaced.
That contrast, not the paper's wording, is what licenses the column, so `tests/test_lipids.py`
recomputes both numbers from the shipped archive: if it stops holding, the column comes out.

Two refusals inside the accepted source, both worth keeping:

- The archive's own EV-minus-cell column is **not** used. It correlates 0.92 with every
  reconstruction attempted and matches none to better than 3 log units, so its transform is unknown.
  A number whose definition cannot be reproduced cannot be verified. The column shipped instead is
  computed here, sample-centred, and the module states the rule.
- The raw difference would have been wrong even if reproducible: a vesicle holds far less total
  material than a cell, so every species reads as depleted and the column measures sample size. The
  word in the slot is *composition*, which is a question about proportion.

Sign was then checked against biology the study does not itself argue: cholesteryl ester,
sphingomyelin, HexCer and cardiolipin come out depleted, PI enriched — what a parasite that
scavenges cholesterol and carries a GPI-rich surface should show. An inverted convention would have
failed that.

### The pattern this adds

A fifth working pattern, and the one with the best yield per hour so far: **when a slot has been
refused three times on the same sentence, suspect the sentence.** "A lipidome of an infected culture
is mostly host" is true and was the right refusal each time; it is not the same statement as "the
parasite's lipid composition cannot be measured". Re-asking what would make the measurement valid —
*host-cell-free material, and a control that varies the host* — named the experiment, and the
experiment existed. This generalises the fourth pattern (re-asking a refused source a better
question) from a source to a slot.

### Four candidates opened and refused this pass

| candidate | for | why not |
|---|---|---|
| PMID 35976251, eLife 80336 | protein turnover | The review draft lists it as "temporal/thermal proteome profiling", and the temporal half was read as turnover. Every one of its seven supplements is CETSA or phospho: `mineCETSA_curve_fits`, AUC, melting curves. There are no half-lives in it. The two thermal columns already shipped are all it holds. |
| PMC11510713, BONCAT-iTRAQ | protein turnover | Measures **newly synthesised** protein, which is synthesis and not degradation, and measures it under pyrimethamine rather than at steady state. 220 proteins. Wrong quantity, not merely thin. |
| PMC9167752, kinome HiT screen | invasion and egress phenotype | Title promises a regulator of invasion and egress and delivers it — by deep follow-up on one gene, SPARK. The per-gene screen columns are `Microscopy phenotype` (Cell Division I/II, Doublets, Singlets, Accumulated IMC1) and a lytic clearance call. Division phenotypes and monolayer clearance, neither of which is this slot, and both of which belong to slots already filled. |
| ToxoDB "Lipidome and palmitoylome" | lipid composition | Named for this slot and is not it: a chemical-proteomic study of the lipid**ated** proteome, gene-indexed, already in the map as the palmitome. |

A ToxoDB pass this round also enumerated all 180 datasets and all 234 gene searches against the
seven remaining slots. `GenesByPhenotypeEvidence` looked like curated phenotype annotation and takes
no parameters — it is dataset filtering, not a vocabulary.

### Where the seven stand

Searched again this pass, nothing found: parasite transcriptome in IFN-γ-activated macrophages
(the naive-macrophage arm exists and is in; the activated one is not published as parasite-side
data), cell-cycle-resolved ribosome profiling, proteome-wide turnover, per-gene drug sensitivity,
and an enteric fitness screen. `resistance-conferring mutation` is the one of the seven that is
**curatable rather than blocked** — Toxoplasma in vitro evolution papers do exist (ROP1 P207S with
ROP8 and TGGT1_237700 for KG8; TgMAPK1 L162Q/I171N for BKIs; artemisinin resistance by serial
passage) but each contributes a handful of genes and no systematic resistome exists for this
organism the way Cowell 2018 does for *Plasmodium*. That is a curation job of roughly eight papers
read end to end, and it yields a deliberately sparse column where sparseness is the biology. It is
the obvious next move on this axis and has not been started.

### `resistance-conferring mutation`: worked, and refused on what the papers actually say

This was named above as the one of the seven that is curatable rather than blocked. It was then
curated, and the curation is what refused it — so the slot stays empty for a reason that is now
specific rather than "not found".

The corpus is real and was read from primary documents, not from summaries. From the artemisinin
in-vitro-evolution study (PMID 31806760), Table 1 and the supplement give nine genes with coding
changes and allele frequencies tracked across 8 µM, 16 µM and 100 µM. All nine accessions resolve in
the current annotation, and each one's product matches the label the paper gives it — a per-row check
that would catch a mistyped accession, and it caught something else first: the paper's body text
calls Ark1 `TGME49_239240` while its table says `TGME49_239420`. The table is right. `239420` resolves
to "protein kinase"; `239240` does not exist in the current annotation at all.

Two findings then stopped the ingest:

1. **The one case where causality was tested, it failed.** The auranofin study (PMID 33816332) ranks
   its variants and names `TGGT1_316330` (SOD2) and `TGGT1_294640` (RNR) as "most likely resistance
   conferring". Read further and the same paper says it "did not reveal a consensus resistance
   locus", and that SOD2 L201P "was **not sufficient to confer resistance** when introduced into
   wild-type parasites". A ranked candidate table with the top candidate experimentally excluded is
   not a resistance column. Its table also carries synonymous changes (L606L, P66P), which are
   background variants and not phenotype at all.
2. **The artemisinin table is titled "Mutations found in candidate genes".** Candidate is the
   authors' word. The evidence for selection is genuinely strong — two independent lines converge on
   the same two genes, and on the same residue in Ark1 (Cys274Arg in F4, Cys274Phe in B2) — but
   convergence under selection is evidence of selection, not a demonstration that the mutation
   confers the resistance.

So the honest column here would be *mutation rising to fixation under drug pressure*, which is not
what the slot says, and shipping it under this slot's name would put a causal claim on the map that
its sources decline to make. The slot's name is not the problem to fix either: "resistance-conferring"
is the right question, and the map should be able to say that Toxoplasma has almost no validated
answer to it.

**What would fill it**, and does exist in principle: mutations validated by introducing them into a
clean background and recovering the resistance. TgMAPK1 L162Q and I171N for bumped kinase inhibitors
are the clearest published case, and the classical DHFR-TS and cytochrome b alleles are others.
That is a small set — probably under fifteen genes — and sparseness would be the biology rather than
a coverage failure. It needs the primary text of each, and the two routes tried here both failed:
`PMC8092512` and `PMC12172689` have no full text in EuropePMC, and the publisher PDF returns 403.
The remaining route is the accepted manuscripts or the authors.

Note for whoever picks this up: an LLM summary of these papers gets this wrong in both directions.
Asked for the mutations, it returned the auranofin candidates without the sentence excluding them,
and it attached the artemisinin quote for Ark1 to the DegP2 row. Both papers were then read as PDF
and XML directly, which is the only reason the refusal is trustworthy.

## Fifth pass, 2026-08-16: the egress screen, and the Plasmodium arm opens

**Toxoplasma 113 of 119. Plasmodium 7 of 103. Combined 120 of 222.**

### `invasion and egress phenotype`, filled on its egress half

Two earlier passes read this as blocked, and both times for the same reason: the Toxoplasma papers
whose titles promise invasion and egress deliver it by characterising *one gene*. The kinome HiT
screen was refused last pass on exactly that. What was missing was the search that asks for the
screen rather than the finding, and it exists — PMID 35538310, an arrayed splitCas9 screen that
disrupted 319 genes one per well and scored the images by eye.

The category codes had to be **earned**. The workbook gives each gene codes like `E3` and `F1 A2`
and ships no legend; the legend is in a figure that is an image. So the mapping rests on two things
rather than on the initials:

* The paper names its categories in three independent places — supplementary discussion, figure
  legend, abstract — and names the same four each time: replication, apicoplast, F-actin, egress.
* The check that does not use initials at all: it names exactly two genes as the egress mutants it
  went on to characterise, CGP (`TGGT1_240380`) and SLF (`TGGT1_208420`), and **both carry an `E`**.
  The test asserts it, so the mapping is a claim the suite can lose.

The **subscript is deliberately not read**. `E3` and `E4` differ in something no accessible text
defines, and a severity column invented out of a digit is a number with no measurement behind it.

Missingness carries half the meaning here and is encoded on purpose: 319 genes were looked at and 99
had a phenotype, so a screened gene with no egress call is `False` and the other 7,800 stay missing.
Collapsing those would tell the map that nearly every Toxoplasma gene has been checked for an egress
defect and passed.

Invasion is still not covered, and the `detail` says so: the screen's own figure legend calls it a
screen for actin dynamics, apicoplast segregation and egress, and invasion is a property its hits
were shown to have *afterwards* rather than a category anything was scored into.

### The Plasmodium arm, 0 -> 7

Instruction 39's per-species design was settled and unbuilt, so all 103 Pf slots read empty for want
of anywhere to look rather than for want of data. One PlasmoDB report now fills seven at grade A:
sequence basics, domain content, conservation breadth, paralogy, strain variation, membrane topology
and asexual-blood-stage fitness. That single call does the work of five Toxoplasma acquisitions,
because PlasmoDB curates into gene attributes what for Toxoplasma had to be found a paper at a time.

The design point that cost the most thought: the Pf columns are named the SAME as the Toxoplasma
ones wherever the quantity is the same, so the two arms can be read side by side — which broke the
existing leakage guard, since it was resting on the pattern strings being disjoint. The guard moved
into the resolution layer instead. A table now reports its own species from what its accessions look
like, derived rather than declared, and both `declared_columns` and `resolve` refuse a slot whose
organism disagrees.

### Where the six remaining Toxoplasma slots stand

Searched again this pass, and this pass exhausted the catalogues rather than the phrasings:
**every one of ToxoDB's 180 datasets, all 234 gene searches, and all 56 RNA-seq datasets** were
enumerated against these six. The only macrophage transcriptome ToxoDB serves is the naive 29-strain
panel already in the map; there is no IFN-γ-activated arm. GEO's nine Toxoplasma ribosome-profiling
series were opened one by one and none is cell-cycle resolved. PRIDE returns nothing for turnover.
The CRISPR phenotype searches carry the in-vitro and in-vivo fitness arms already ingested and no
drug arm.

`resistance-conferring mutation` remains the one that is blocked on **access** rather than on
existence: it needs mutations validated by reintroduction, and the primary text for the two
candidate papers is absent from EuropePMC and 403 from the publisher.

## Sixth pass, 2026-08-16: `resistance-conferring mutation`, filled by raising the bar

**Toxoplasma 114 of 119.** Five gene slots remain.

Last pass curated this slot and refused it, on the grounds that the Toxoplasma in-vitro-evolution
studies produce *candidate* loci and that the one candidate ever tested for causality failed. That
refusal was right about those studies and wrong about the slot. The question "which genes carry a
mutation that confers resistance" has an answer in this organism; it is just not in the evolution
papers. It is in the papers that **put the mutation back**.

`TgMAPK1` / `TgMAPKL-1` (`TGME49_312570`) carries L162Q and I171N for 1NM-PP1, with cross-resistance
to 3BrB-PP1 and 3MB-PP1 (PMID 24533298), and the gatekeeper S191Y (PMID 25941623) — where one
residue swapped two ways in one background gives a sensitive clone and a resistant one, which is as
clean as causal evidence gets. Both papers are open access and were read directly.

So the slot is filled by **one gene**, and one gene is the honest number under that bar. Nine genes
were curated from the artemisinin and auranofin studies and then dropped.

### Three alleles deliberately not carried, recorded as data rather than omitted

`resistance.DOCUMENTED_ELSEWHERE` names them so the next pass does not rediscover the same dead ends,
and a test asserts none has been added without a source:

* **DHFR-TS** (`TGME49_249180`) pyrimethamine alleles — established well enough that the mutant is
  the field's standard selectable marker, but the primary text is pre-PMC and unreachable, and a
  review's paraphrase is not the measurement.
* **DHODH** (`TGME49_210790`) N302S — primary paper not open access.
* **Cytochrome b** M129L, I254L for atovaquone — excluded for a *structural* reason rather than an
  access one. Toxoplasma cytochrome b is mitochondrially encoded and has no row in a table of 8,140
  nuclear genes; the eighteen nuclear "cytochrome b" hits are b-c1 subunits and b5-domain proteins,
  and putting an atovaquone allele on one of those would be a plain error.

### What makes a hand-typed table safe

This is the only source in the map with no parser to fail and no archive to disagree with, so the
usual safety net is gone. What replaces it is a per-row check that costs nothing: every allele
records the gene's product **as the current annotation gives it**, and a test asserts they match. A
transposed digit either names nothing or names a different protein, and either way the product stops
matching — the same check that caught the artemisinin paper writing its own kinase accession two
ways. A second test asserts no row's validation text reads as *candidate*, *associated with* or
*correlated*, so the bar cannot erode by someone adding a weaker row later.

### The five that remain

`transcription · in IFN-gamma macrophage`, `translation · per cell-cycle phase`, `protein turnover`,
`drug sensitivity`, `fitness · in vivo gut`. All five were re-searched this pass at catalogue level
rather than by phrasing — every ToxoDB dataset, search and RNA-seq series, every Toxoplasma
ribosome-profiling series in GEO, and PRIDE. Each needs an experiment that has not been done in this
organism.

## Seventh pass, 2026-08-16: the campaign closes at 114, and says so in the table

**Toxoplasma 114 of 119.** The five that remain were re-searched once more, this time by chasing
four specific candidates rather than re-running phrasings. All four were opened and all four missed:

| candidate | for | why not |
|---|---|---|
| ToxoDB "4 mouse cell types" | IFN-γ macrophage | The four are neurons, skeletal muscle, astrocytes and fibroblasts. No macrophage, no activation. |
| GSE230866 | IFN-γ macrophage | Profiles IFN-γ-activated cells infected with three parasite strains — and is `taxon: Homo sapiens`. Only the host side was sequenced. |
| PMID 41407671, GRA38 | drug sensitivity | A genome-wide CRISPR screen under a perturbation, and the perturbation is 1% vs 10% serum. Lipid limitation is a nutrient, not a compound. |
| PMID 39082802 | fitness in vivo gut | In vivo CRISPR screens of 600 hyperLOPIT-unassigned proteins, scored for systemic virulence. Not enteric. |

### What changed instead: the table now distinguishes two kinds of empty

A dash said a slot had no data and said nothing about whether anyone had looked, so a question nobody
has searched read exactly like one seven passes had exhausted. Those want opposite next actions, and
losing the difference has a concrete cost that was paid twice during this campaign: the same searches
were re-run because the previous refusal had not recorded where it looked.

Every empty Toxoplasma slot now carries three fields in the generated table — `blocked_by`,
`searched`, `would_fill_it` — and the verdict is one of two words. **`missing`** means the
measurement has not been made in this organism. **`unreachable`** means it has been made and the data
cannot be got at. All five currently read `missing`; `resistance-conferring mutation` would have read
`unreachable` before it was filled, which is exactly the distinction that made it worth another
attempt while the other five were not.

Three tests hold this up: every empty slot must carry a verdict, every verdict must name where
someone looked and what would fill the slot, and a slot that gets filled must lose its verdict so a
stale one cannot outlive the gap it described. A new empty slot fails until someone writes down what
they searched, which is the cheapest possible moment to write it.

### The five, and the experiments they are waiting for

* `transcription · in IFN-gamma macrophage` — dual RNA-seq inside IFN-γ-activated macrophages with
  the parasite reads kept.
* `translation · per cell-cycle phase` — ribosome profiling of synchronised or FUCCI-sorted
  tachyzoites. The FUCCI probes now exist (PMID 40590555), so this one has become newly feasible.
* `protein turnover` — pulse-SILAC or a cycloheximide chase with proteome-wide degradation rates.
* `drug sensitivity` — a genome-wide CRISPR screen under compound pressure.
* `fitness · in vivo gut` — a pooled screen through the enteroepithelial stages.

Four of the five are ordinary experiments that simply have not been run in Toxoplasma. That is the
campaign's real finding and it is now visible in the artifact rather than only in this file.

## Eighth pass: the Plasmodium transcriptome, 7 -> 17

**Combined 131 of 222.** Toxoplasma stays at 114 of 119; the five that remain need experiments
nobody has run, and this pass did not pretend otherwise. The work went where slots could still be
filled from data that exists.

One PlasmoDB report carrying three studies fills ten more slots. The split between them was the
design decision, and it is the leakage rule rather than convenience that dictated it:

| study | answers |
|---|---|
| Su seven stages | ring, trophozoite, schizont, gametocyte, ookinete |
| Gomez-Diaz mosquito stages | asexual blood stage, oocyst, sporozoite |
| Bunnik polysomal IDC (0h, 18h, 36h) | **steady-state arm** → transcription per cell-cycle phase; **polysomal arm** → translation |

Handing the stage slots and the cell-cycle slot the same columns would have been one measurement
claimed twice. Two studies, two slots, no shared column.

### Stage labels are checked against biology, not trusted

The columns come out of PlasmoDB as one wide report and are matched by substring, so a mislabelling
would silently shift a stage and nothing downstream would notice. The test asserts marker genes peak
where a century of malaria biology puts them: CSP in sporozoite, MSP1 in schizont, Pfs16 in
gametocyte II.

**Pfs25 is the informative one.** Its transcript peaks in **gametocyte V**, not in the ookinete where
the protein does its work — the textbook translational-repression stockpile. Anyone validating
against the protein literature would read that as an off-by-one error and "fix" a correct column.
The test pins the transcript behaviour and the docstring says why.

### A test that was passing by accident

`test_every_registered_dataset_contributes_at_least_one_column` checks each dataset against the
tables it could land in, and it did not know about the species tables. The first Plasmodium dataset
passed it anyway — because the Pf column names are deliberately mirrored from the Toxoplasma ones,
and `length` exists in both. The second dataset, whose stage names are Plasmodium-specific, is what
exposed it. Mirroring the names buys comparability between the arms and costs exactly this: checks
that key on a column name can no longer tell the arms apart, and have to be told which table to look
in. That is now the third place this session where that cost has come due, after `declared_columns`
and `resolve`.

## Ninth pass: two catalogue-complete negatives, and a slot refused after it was nearly filled

**Toxoplasma 114 of 119, Plasmodium 17 of 103, combined 131 of 222.**

### The two Toxoplasma slots that were still "maybe" are now definite

Both were re-tested by enumerating the whole catalogue rather than sampling keywords, which is the
difference between "I did not find it" and "it is not there":

* **`protein turnover`** — all **201** Toxoplasma deposits in PRIDE listed and scanned. Four mention
  stability or a drug; none measures degradation. (The four are the DegP2 thermal-stability set,
  two zaprinast phosphoproteomes, and an HDAC-inhibitor antigen study.)
* **`drug sensitivity`** — **100** Toxoplasma screen papers enumerated from EuropePMC and scanned.
  The only CRISPR drug-resistance screen among them is *Leishmania*.

Both verdicts in the slot table now say what was enumerated, so the next pass can see the difference
between a search and a sweep.

### `Pf_protein abundance` was built, validated, and then refused

PlasmoDB serves a TMT proteome of ring, trophozoite and schizont as three gene attributes. They were
fetched, merged, and the slot filled — and then the validation caught it.

Protein correlated with its own transcript at **−0.25** in ring. Plasmodium's proteome is known to
lag its transcriptome, so the first guess was that the lag would put the peak off-diagonal; it did
not, and the structure made no sense as biology either way. The actual answer was in the shape of the
data: every gene's three values **sum to 12.07 ± 0.20**, and the three columns anti-correlate with
one another (−0.46, −0.74, −0.14). PlasmoDB serves that study **row-normalised**. The numbers are a
protein's distribution across the cycle, not its abundance.

So the columns are named `protein_stage_share_*` and `protein abundance · asexual blood stage` stays
empty. What makes this worth recording is how well the error was hidden: named `protein_ring` it
would have been claimed by an abundance slot automatically, and the only symptom was a −0.25
correlation that reads as a biological puzzle rather than as an artifact. Two tests now pin it — one
asserts the columns really are compositional (so a future PlasmoDB release serving true abundances
fails rather than keeping a wrong label), and one asserts no slot claims a share.

That is the same refusal as the vesicle EV-minus-cell column and the auranofin candidates: the
question is never whether a number is available, it is whether it is the quantity the slot names.

## Tenth pass: the Plasmodium graph, and the third species guard

**Toxoplasma 114 of 119, Plasmodium 20 of 103, combined 134 of 222.**

A second graph file, because an edge is a pair of indices into a table and a *falciparum* gene has no
index in the Toxoplasma one. Three layers — `orthogroup`, `domain`, `coexpression` — built with the
Toxoplasma constructions copied rather than re-invented, so that a difference between the arms means
the biology differs and not that the edges were drawn by different rules.

### The duplicate the Toxoplasma arm never showed

Two proteins can share more than one InterPro domain, and emitting the pair once per shared domain
draws the same edge repeatedly — which reads as repeated evidence. The naive construction produced
**10,764 duplicate emissions among 34,887** here, nearly all inside the var, rifin and stevor
families that share whole multi-domain architectures.

The Toxoplasma arm builds this layer the same way and emits **no duplicates at all** — checked, not
assumed — which is why the fault was invisible for the life of the project. It would acquire it the
moment its annotation gained a pair sharing two domains. The count is now the weight, which is both
correct and strictly more informative than the 1.0 it replaced: a pair sharing eleven domains says
so.

### The species guard, third instance

`orthogroup`, `domain` and `coexpression` are named identically in both graphs — deliberately. So a
Plasmodium pair slot handed the Toxoplasma graph finds every layer it asked for and **reads as
filled by another organism's edges**. Three did, the day the Plasmodium graph was built, and it
surfaced only because an unrelated assertion (`slots.coverage("Pf", …)["filled"] == 0`) still
encoded the old world where Plasmodium had nothing.

A graph carries no accessions to identify itself, so it is identified by the table it arrives with:
`nodes` and `graph` describe one species or the caller has mixed two caches. Unknown stays permissive
so callers holding only a graph still work.

That is now **three places** where mirroring the column and layer names between arms has cost a
name-keyed check — `declared_columns`, `resolve`, and now `is_filled` — plus the dataset-contribution
test that was passing by accident. The mirroring is still right: it is what lets the arms be compared
at all. But the rule it implies should be written down before a third species arrives: **nothing may
identify an organism by the name of a column or a layer. Only the table's own accessions do that.**
