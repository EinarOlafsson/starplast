# 41 — Fill as many slots as possible: Toxoplasma, Plasmodium, and host

**Status: open. Requested 2026-08-15.**

## The work queue, as it stands 2026-08-18: seventeen slots are past discovery

Checked rather than assumed, and it changes what the next session should do first. Of the 54 empty
Plasmodium gene slots, **seventeen already carry resolved candidate studies** from the earlier
candidate passes -- accession, title, year and journal, merged onto the slot by
`_merge_candidate_references`. For these the remaining work is fetch, verify, key to Pf accessions
and load; discovery is done.

| slot | candidates | first candidate |
|---|---|---|
| `Pf_resistance-conferring mutation` | 8 | 24533298: Identification of mutations in TgMAPK1 of Toxoplasma gondii conferring resista |
| `Pf_invasion and egress phenotype` | 7 | 35538310: A splitCas9 phenotypic screen in Toxoplasma gondii identifies proteins involve |
| `Pf_essentiality in a second background` | 7 | eight passages in the differentiation reporter strain | 41686849: LAMP-coupled CRISPR-Ca |
| `Pf_metabolite levels` | 6 | 38928131: Mechanisms Underlying the Effects of Chloroquine on Red Blood Cells Metabolism |
| `Pf_metabolic flux` | 6 | 35064153: Metabolic adjustments of blood-stage Plasmodium falciparum in response to subl |
| `Pf_lipid composition` | 6 | 41964222: Deoxy-Piezo1 hyperactivity elevates pump-leak fluxes and lactate production in |
| `Pf_target engagement / thermal shift` | 6 | 42451739: Essential Oil of Symplocos chinensis (Lour.) Druce: Chemical Composition, Anti |
| `Pf_RNA-binding protein targets` | 6 | 33207342: The Route of Infection Influences the Contribution of Key Immunity Genes to An |
| `Pf_transcription · liver stage` | 6 | 38657074: Autonomous circadian rhythms in the human hepatocyte regulate hepatic drug met |
| `Pf_transcription · mosquito stages` | 6 | 41510253: A divergent Plasmodium NEK4 acts as a key regulator driving the early events o |
| `Pf_transcription · dormancy / recrudescence` | 6 | 42374402: Mapping the intellectual landscape of malaria drug repurposing: a systematic a |
| `Pf_fitness · liver stage` | 6 | 39714137: Autophagy protein Atg7 is essential for maintaining malaria parasite cellular  |
| `Pf_fitness · transmission` | 6 | 42457405: Engineered promoter system enables high-efficiency transgenic CRISPR editing i |
| `Pf_antigenic variation family expression` | 6 | 34389510: Plasmodium falciparum SET2 domain is allosterically regulated by its PHD-like  |
| `Pf_host receptor binding` | 6 | 42291314: MAHRP2 is required for tether formation and cytoadherence in Plasmodium falcip |
| `Pf_field variation and resistance markers` | 6 | 42557356: Genomic surveillance reveals co-occurrence of Plasmodium falciparum drug resis |
| `Pf_host ESCRT recruitment` | 1 | pooled image-based CRISPR screen of secretory proteins |

The other 37 empty Pf slots have no resolved candidate and still need the search.

**The dangerous step is the load, not the download.** `build_graph` exits 0 while losing data, so
every fill must diff the rebuilt cache against the previous one on columns gained, columns LOST, and
per-column coverage in BOTH directions -- presence alone has already missed a 4-gene regression and a
1,306-value wipe in columns that still existed, and coverage going UP is how two identity-layer bugs
surfaced.
 This is the acquisition campaign. 39 defines the tables and
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

## Eleventh pass: export prediction, and a threshold that would have lost the textbook

**Toxoplasma 114 of 119, Plasmodium 21 of 103, combined 135 of 222.**

PlasmoDB's gene attributes carry no PEXEL, export or PTM annotation — checked across all 3,070 — so
the remaining Plasmodium PTM and localisation slots need deposits one at a time, the same acquisition
work Toxoplasma took. One exception was sitting in the search list rather than the attribute list:
**ExportPred**, a hidden Markov model of the signal sequence and PEXEL motif.

### Why it is a tier and not a boolean

PlasmoDB serves it as a search with a score threshold, so the score was recovered by asking at
several thresholds and keeping the highest each gene survives. The scale saturates — asking for 20
returns nothing — so 10 is the algorithm's own default and the top tier rather than an arbitrary cut.

Taking that default as a boolean would have been wrong in a way that a composition check nearly
missed. At score ≥ 10 the 191 genes are exactly the textbook exportome: *Plasmodium* exported
proteins, rifin, PfEMP1, stevor, the FIKK kinases, KAHRP. It looks right. But **MESA and PfEMP3 —
exported by any textbook — both fall below 10** and only appear at the permissive threshold. A
boolean at the default would have called two of the best-known exported proteins in the organism not
exported, and the composition check would still have passed.

Two tests hold the distinction: one asserts the known exportome is found at the default, the other
asserts MESA and PfEMP3 are present *below* it.

### The other judgement, written into the assembly

Absence from ExportPred is a **real negative**, not a gap: it is a sequence model evaluated on every
protein, so its silence is a prediction of "not exported". That is the opposite of the screen columns,
where absence inside the screened set means normal and outside it means nobody looked. The rule now
lives in `plasmodium.build_all` rather than in whatever was typed at a prompt, which is also what
made the whole table reproducible for the first time.

`exposure to host cytosol` is left empty on purpose. That slot wants a measured exportome, and a
prediction filling it would be a model answering for an experiment.

## Twelfth pass: the pooled phosphoproteome

**Toxoplasma 114 of 119, Plasmodium 22 of 103, combined 136 of 222.**

PlasmoDB serves no PTM annotation as gene attributes — checked across all 3,070 and all 325 searches
— so the Plasmodium PTM slots need deposits one at a time. A sweep of all **251 *P. falciparum*
PRIDE deposits** found 26 carrying a PTM, and one of them is worth more than the rest:
**PXD046874**, a re-analysis of every public Plasmodium phosphoproteomics dataset through a single
pipeline. That is what makes a per-gene count meaningful — the same serine found by three groups is
one site rather than three.

**16,318 distinct sites over 2,503 genes.**

### The counting trap

The deposit's tables are *site-centric* and still carry **one row per peptidoform and per source
run**. Summing rows would count how often a protein was looked at rather than how many sites it has,
and the difference is not subtle: the files hold millions of rows for 16,318 sites. Sites are
counted as distinct `(gene, position)` pairs, accessions are stripped of their transcript-and-product
suffix (`PF3D7_1346300.1-p1`, or a gene with two products counts twice), and where a source has a
merged table its per-run siblings are skipped.

### Validated on orderings, not totals

A count like this cannot be checked against a published number, because pooling changes it. It can be
checked against things that must be true whatever the total: **67% of kinases carry a site against
44% of genes at large** (activation loops and autophosphorylation), and **site count rises with
protein length at rho +0.46**. Both survive the numbers being re-derived.

### Missingness, mirrored from the Toxoplasma arm

The **count** stays missing where nothing was detected — how many sites a protein has is genuinely
unknown if mass spectrometry never saw it. The **flag** is False, because whether it was ever
observed phosphorylated is a question about the evidence, and the answer is no. That is the
Toxoplasma convention exactly, checked against it rather than reinvented, so the arms read the same.

## Thirteenth pass: the palmitome, and a sheet that lies about itself

**Toxoplasma 114 of 119, Plasmodium 23 of 103, combined 137 of 222.**

The 26 Plasmodium PTM deposits found in the PRIDE sweep mostly carry search-engine output — `.msf`,
MaxQuant internals, `.mzid` — rather than per-gene tables, which is the situation instruction 41
already has a rule for: **paper supplements beat repositories when the wanted quantity is the study's
analysis**. Applied here it went straight to a usable table.

### The sheet named for the wrong thing

PMID 36250062's workbook has a sheet called **`nrPalmitoylatedProteins`** with 3,105 rows. Its name
says palmitoylated. Its contents are the **union of palmitoyl-ABLE — a motif prediction over 2,902
proteins — and the 503 actually observed.**

Reading it by name would have called **54% of the proteome palmitoylated**, against published
palmitomes of 400 to 500. What gave it away was not the name but the size, and then the first rows:
PfEMP1 and rifin, the families that dominate any cysteine-presence prediction. The loader reads the
observed column of the `Palmitome` sheet and the docstring says why, with a test asserting the
predicted column cannot leak in.

This is the fourth time this campaign that a source's own label was the thing to distrust — after the
compositional TMT proteome served as "abundance", the archive delta whose transform could not be
reproduced, and the ExportPred default that drops MESA and PfEMP3.

### Validated on substrates and mechanism, not on a total

GAP45 and CDPK1 — the canonical Plasmodium palmitoylation substrates — are both present, and membrane
proteins are enriched **2.1-fold** among the palmitoylated (44% against 27%, p = 7e-15), which is
what a membrane-anchoring modification has to do. ARO is a known miss; no palmitome is complete, and
absence here means not observed.

## Fourteenth pass: a contrast deliberately not computed

**Toxoplasma 114 of 119, Plasmodium 24 of 103, combined 138 of 222.**

The Sir2 knockout microarray is the chromatin perturbation this arm has: wild type against *sir2a*
and *sir2b* knockouts, at ring, trophozoite and schizont.

The obvious column is knockout minus wild type, and it is **not** computed. The independent check on
it came out ambiguous. Sir2a silences subtelomeric *var* genes, so *var* should rise in the *sir2a*
knockout — and it does in ring (+0.135, p = 4e-12) and schizont (+0.130, p = 2e-38), but **falls in
trophozoite** (−0.240, p = 3e-22), with all effects small against a spread of 0.7.

That is consistent with the canonical result being subset-specific and with *var* probes
cross-hybridising across sixty paralogues. It is not a clean confirmation. A derived column carrying
an unexplained sign flip would state more confidence than there is, so the nine **conditions** ship
instead — those are unambiguous: PlasmoDB names them, the values are log intensities, and the array
medians align within 0.1. A test asserts that alignment, because it is the precondition that makes
differencing them meaningful at all, and it is the thing that would break first if PlasmoDB
renormalised.

This is a different kind of restraint from the four refusals before it. Those were cases where a
number was not the quantity the slot named. Here the number probably *is* the right quantity and the
evidence for it is merely weaker than a column implies — so the answer is not to refuse the data but
to ship it one step further back, and say why.

## Fifteenth pass: `Pf_chromatin accessibility` built and refused

**No change to the totals: Toxoplasma 114 of 119, Plasmodium 24 of 103, combined 138 of 222.**

A Plasmodium ATAC-seq study (PMID via PMC13032594) ships a peak table already annotated to genes,
with the paper's own `Promoter` calls — no coordinate pipeline needed, 3,798 genes with a
promoter-proximal peak. It was built and then refused on three checks that all point the same way.

| check | Toxoplasma benchmark | this data |
|---|---|---|
| rho(promoter ATAC, mRNA) | **+0.495** (`atac_promoter_ut`) | **+0.067** |
| expressed vs silent promoters | — | 161.7 vs 155.7, no separation |
| heterochromatic families | should be LESS accessible | **MORE accessible**, 178.1 vs 157.5, p = 3e-17 |

The third is the one that settles it. *var*, *rifin* and *stevor* sit in heterochromatin, and their
promoters come out **more** accessible than the genome at large. That is backwards, and it is the
signature of read pile-up in subtelomeric multigene families whose members are near-identical —
mappability, not chromatin. Normalising by peak width or taking the maximum peak instead of the mean
changes the correlation by less than 0.005, so it is not a summarisation choice.

This is the same refusal as H3K4me1 in the Toxoplasma arm, and deliberately decided the same way:
there, marked genes had *less* accessible promoters than unmarked ones, and the H4-acetylation data
from the same site and assay behaved correctly — which proved the method was fine and the data was
not. Here the Toxoplasma ATAC column plays that role: the same measurement, in the same map, relates
to transcription seven times more strongly.

`Pf_chromatin accessibility` therefore stays empty, and its verdict is `missing` in the sense that
matters — a usable measurement has not been published for this organism, even though an ATAC-seq
experiment has.

## Sixteenth pass: model confidence, fetched one protein at a time

**Toxoplasma 114 of 119, Plasmodium 25 of 103, combined 139 of 222.**

`Pf_fold confidence / disorder` is filled from AlphaFold DB, **5,098 of 5,720 genes**. The bulk
proteome archive is not at the documented path for this organism, so the summaries were fetched per
protein from the API — 5,306 requests — which turned out to be better than the archive would have
been: the API returns the mean pLDDT *and* the fraction of each model in each confidence band,
without downloading or parsing a single structure.

**The fractions matter as much as the mean.** A protein that is half well-folded and half disordered
scores the same mean as one that is uniformly mediocre, and those are not the same protein. The slot
asks about disorder as well as confidence, so both ship.

**The accession problem.** A gene can carry several UniProt entries — 876 do, mostly the variant
surface families where each field isolate's allele has its own. The fetch tries them in order and
`alphafold_accession` records which one supplied the model, so a number can be traced to the
structure it came from rather than to a gene that has eight.

### Validated on an ordering, and one number that looks wrong and is not

Proteins carrying a recognised InterPro domain model at median pLDDT **73.5** against **56.0** for
those without (p = 2e-159) — a domain is a thing that folds.

The correlation with protein length is **−0.555**, which would be alarming in most proteomes and is
correct in this one: *P. falciparum* is famous for long low-complexity asparagine insertions, which
are exactly what AlphaFold models with no confidence. Recorded here because the next person to check
it will have the same moment of doubt.

## Seventeenth pass: the myristoylome, and three states rather than two

**Toxoplasma 114 of 119, Plasmodium 26 of 103, combined 140 of 222.**

`Pf_N-myristoylation` filled from PMID 34695132. The evidence for a substrate is **not** being pulled
down — background comes down too — but coming down **less** when the transferase is inhibited, so the
loader requires significance *and* a negative difference. A positive difference under a blocked
transferase would be a protein that came down more without it, which is not what a substrate does.

**Three states, not two.** 16 substrates; 593 proteins assayed and not substrates; 5,111 genes never
in the pulldown, which stay missing. Collapsing the last two would tell the map that the whole
proteome had been tested for myristoylation by one experiment that saw 609 proteins.

Sparse because the biology is — Plasmodium has roughly thirty predicted NMT substrates — and the list
validates itself: **GAP45, ARO, CDPK1, Rab-5B, ARF1 and ISP3** are the canonical N-myristoylated
families in apicomplexans, and all are present.

## Where the Plasmodium arm stands

**0 to 26 this session**, from nothing: a node table, its own graph, and ten datasets. Filled slots
span sequence, orthology, domains, strain variation, fitness, seven life-stage transcriptomes,
translation, export, model confidence and disorder, phosphorylation, palmitoylation, myristoylation,
chromatin perturbation, and three relation layers.

One slot was built and refused — `chromatin accessibility`, where *var* genes came out **more**
accessible than the genome at large — and one column was renamed rather than shipped under the slot
it was fetched for, when the TMT proteome turned out to be compositional.

The remaining 77 need acquisition of the same kind: a deposit or supplement at a time, each with a
check that would fail if the numbers were not what the slot says.

## Eighteenth pass: the acetylome, and two columns held to two standards

**Toxoplasma 114 of 119, Plasmodium 27 of 103, combined 141 of 222.**

`Pf_acetylation` filled from PMID 26813983: 1,145 genes, 2,163 localised sites.

The judgement here is that **identifying an acetylated peptide and localising the acetyl group to a
particular lysine are different claims**, so the two columns are built to different standards. The
flag uses every identification; the count uses only sites with an Ascore of 0.75 or better, because
a site count is meaningless if you do not know which lysine it is on.

That mattered: the sheet is titled **"Final Ac-K List"** and is **not** pre-filtered on localisation
— Ascores run all the way down to 0. Taking its length as a site count would have been wrong by
about a quarter. The word "Final" refers to the identification list, not to site confidence, and
this is now the fifth time a source's own label needed reading past rather than taking at face value.

Self-validating: fourteen histones appear, and the most heavily acetylated proteins are the PHD
finger proteins, the MYST acetyltransferase and the coactivator ADA2 — the acetylation machinery
itself, which is what any acetylome should be led by.

## Nineteenth pass: the strain-accession problem reaches the Plasmodium arm

**Toxoplasma 114 of 119, Plasmodium 28 of 103, combined 142 of 222.**

`Pf_lactylation` filled — and the interesting part is not the PTM. The study reports against the
**NF54** annotation rather than 3D7, so joining on the accession string would have dropped all 186
genes without a word. That is exactly the failure the Toxoplasma identity layer exists to prevent
(`TGGT1` deposits against a `TGME49` table), met here for the first time on this arm.

### Orthology as identity, and the check that licenses it

`plasmodium.strain_map` resolves NF54 to 3D7 through orthogroups holding **exactly one gene on each
side**, so the paralogous surface families are dropped rather than guessed at — guessing there would
attach a measurement to the wrong family member, which is worse than losing it.

Orthology is a claim about ancestry and this needs a claim about identity, so the map is checked
against protein length: **96.8% of the 4,310 pairs have exactly the same length, 99.2% within 5%**.
That is what it should look like when one line was cloned from the other, and a test fails if a
future PlasmoDB release breaks it. Pairs differing by more than half are dropped; pairs with no
length are kept, because unknown is not a contradiction.

A first attempt validated the map on product-description text instead and got 63%, which looked
alarming until the disagreements turned out to be the same protein worded differently
("SufB protein" against "iron-sulfur cluster assembly protein SufB"). Length is the better check
precisely because it does not depend on annotation prose.

144 of 186 genes resolve. Site counts use a 0.75 localisation cut and the flag does not — the same
split as acetylation, for the same reason.

**The map is reusable.** It is not lactylation-specific, and the next Plasmodium source reported
against NF54, 7G8, Dd2 or any other strain can go through the same door.

## Twentieth pass: an m6A refusal that turned into an isoform fill

**Toxoplasma 114 of 119, Plasmodium 29 of 103, combined 143 of 222.**

Opened PMID 40316999 for its m6A data and refused it: the *falciparum* arm of that table is a
**43-gene intersection with *P. vivax***, not a methylome, and a column called "RNA modification"
built on it would represent the Pf epitranscriptome as forty-three genes when it runs to thousands.
`Pf_RNA modification` stays empty.

The same supplement carries a **SQANTI long-read isoform table** — 2,498 transcript models over 1,857
genes — which fills `Pf_splicing / isoform use` instead. `full-splice_match` is the reference
transcript recovered and is deliberately not counted as novel; only `novel_in_catalog`,
`novel_not_in_catalog` and `fusion` are, which is what the Toxoplasma column of the same name counts.
238 models are not in the annotation.

Absence stays missing rather than reading as 1: a gene with no long-read model was not sequenced
deeply enough to say, which is a different statement from having a single transcript. The obvious
correlation holds and is worth stating rather than hiding — more expressed genes yield more models
(rho +0.36) — which is detection depth, and is why the count is not read as isoform diversity.

Worth noting as a pattern: this is the second time a source opened for one slot filled a different
one. The ESCRT round did the same thing, and both times it happened because the supplement was read
rather than the abstract.

## Twenty-first pass: two derived slots, and a construction shared rather than copied

**Toxoplasma 114 of 119, Plasmodium 31 of 103, combined 145 of 222.**

`Pf_transcription · maximum observed across stages` and `Pf_life-cycle stage label (derived)`, both
computed from the nine stage columns and both declaring `derived_from` so leakage closure excludes
them together with the measurements underneath.

The stage call **reuses `cellcycle.stage_enrichment`** rather than reimplementing it — the function
gained a `stages` parameter and the Toxoplasma default is untouched. That is the same decision as
copying the graph constructions: if the two arms' stage labels are ever compared, a difference should
mean the biology differs and not that one arm z-scored and the other did not. A test asserts the
Toxoplasma default survives being passed a different map.

Only **310 of 5,720** genes are labelled, because the shared rule leaves a gene unlabelled unless one
stage leads the next by half a z-unit. That is the point: a label that is really a coin toss looks
like a measurement in every table it reaches.

The class counts need the same caveat the Toxoplasma arm carries, and for the same structural reason.
Ookinete takes 181 of the 310 — not because it uses more genes, but because ring, trophozoite and
schizont are highly correlated with one another and rarely win by a margin, while the mosquito stages
are separable. The margin rule is working; the interpretation is what needs care.

Two smaller things this pass also fixed: `kind` must name the assay underneath a derivation and not
the derivation itself (the entry first said `derived`, which the registry test correctly refused —
it is RNA-seq, and `derived_from` is what says it is computed), and `plasmodium` needed numpy, which
only surfaced because the first derived column used it.

## Twenty-second pass: antibody epitopes, and a better key than the Toxoplasma arm has

**Toxoplasma 114 of 119, Plasmodium 32 of 103, combined 146 of 222.**

`Pf_seroreactivity / antigenicity` from IEDB directly: 14,610 assay records over 444 antigens,
reducing to **7,366 distinct epitopes across 434 genes**.

**Distinct sequences, not assay records.** MSP1 alone carries 1,739 epitopes, so counting records
would rank antigens by how many groups have studied them rather than by how much of the protein
antibodies recognise — the same reasoning `iedb` records for the Toxoplasma arm.

Where the two arms differ is the key, and the Plasmodium one is better through no merit of mine.
IEDB's Toxoplasma antigen names are verbatim ToxoDB product descriptions, so that arm has to match
descriptions and loses eleven antigens to generic names. The falciparum names carry the **UniProt
accession**, which resolves 434 of 444. The 209 PlasmoDB accessions naming more than one gene are
dropped rather than assigned: an epitope belongs to a protein, and giving it to whichever paralogue
sorted first would be inventing the answer.

Validated against the history of the field rather than a number: **MSP1 is the top antigen and CSP —
the RTS,S vaccine antigen — is present.** Those are the two most studied antigens in the organism.

Absent is absent and not zero, because IEDB records what somebody tested.

## Twenty-third pass: the other half of IEDB, and why it is a separate slot

**Toxoplasma 114 of 119, Plasmodium 33 of 103, combined 147 of 222.**

`Pf_T-cell epitope content` from `tcell_search`: 6,134 assay records over 44 antigens, **1,542
distinct epitopes**. Same UniProt keying and same distinct-sequence counting as the antibody half.

The two halves are kept apart, and the numbers are the argument for it: **434 antigens carry an
antibody epitope and only 44 carry a T-cell one.** Pooling them, or filling either slot with the
other's number, would answer one question with the other — which is exactly the mistake `iedb`
records for the Toxoplasma arm, where ToxoDB's undifferentiated epitope count could not answer a
question about antibodies.

The loader reads whichever halves are on disk, so one fetch failing costs its own column and not the
other.

## The Plasmodium arm at 33

Thirteen datasets and two computed layers, built this session from a standing start. The pattern that
emerged is worth stating compactly, because it is what the remaining 70 slots need:

1. **Find the source, then read past its labels.** Five sheets or fields this session said something
   other than what they contained.
2. **Key it correctly.** Product descriptions, UniProt accessions and strain annotations each needed
   a different door, and one of them needed a new one built.
3. **Validate on an ordering that must hold whatever the numbers are** — kinases phosphorylated more
   than average, domains folded better than their absence, MSP1 the top antigen, CSP present.
4. **Encode missingness in as many states as the experiment produced**, which was two states
   sometimes and three others.
5. **Refuse when a check comes back backwards**, which happened twice.

## Twenty-fourth pass: febrile stress, and a null check that is not a failure

**Toxoplasma 114 of 119, Plasmodium 34 of 103, combined 148 of 222.**

`Pf_transcription · under stress / conversion` from the febrile-temperature series: wild type and two
mutants at 37 °C and at the 41 °C of a malarial fever.

Shipped as **conditions** rather than as a 41-versus-37 contrast — the same restraint as the Sir2
entry, for a different reason. There the check on the contrast **contradicted itself**; here it came
out **null**. Heat shock proteins move by a median log2 of +0.08 against −0.07 for everything else
(p = 0.2), so a fever does not measurably induce them.

That null is consistent with what is known rather than evidence of a fault: this organism's chaperones
are constitutively high rather than stress-induced, and the genes that *do* rise — Maurer's cleft
two-TM proteins and stevor, at four to five log2 — match published fever-driven surface remodelling.
But a null result on the one available prediction is not a validation, and a derived column would
imply it had passed one.

Worth separating three outcomes this campaign has now produced, because they want different actions:

| check result | example | action |
|---|---|---|
| **backwards** | *var* genes more accessible than the genome | refuse the data |
| **contradictory** | *var* up in two stages and down in a third under Sir2 KO | ship conditions, not the contrast |
| **null** | heat shock proteins unmoved by fever | ship conditions, and say the prediction was null |

Only the first is a reason to reject a source. The other two are reasons to ship it one step further
back and write down what was and was not shown.

## Twenty-fifth pass: two findings about the CATALOG, not the data

**No change to the totals: Toxoplasma 114 of 119, Plasmodium 34 of 103, combined 148 of 222.**

Searching for the remaining Plasmodium fitness sources turned up two problems in the slot list itself.
Both are recorded here rather than acted on, for a reason given at the end.

### Four Plasmodium slots are unanswerable by construction

`Pf_fitness · in vivo peritoneum`, `· in vivo lung`, `· in vivo liver` and `· in vivo spleen` carry
the contexts **mouse peritoneum, mouse lung, mouse liver, mouse spleen**. Those are the Toxoplasma
in-vivo CRISPR screen sites, mirrored one-for-one onto *P. falciparum* — **a human parasite that does
not infect mice.**

These four are not empty for want of data. The question cannot be asked. `_pf_mirror` already knows
this class of problem exists — its `stage_words` regex deliberately catches `macrophage|hff|bmdm` as
host contexts needing conversion rather than mirroring, with the comment "HFF is the fibroblast
Toxoplasma is cultured in and Plasmodium does not grow in one". The mouse organs slipped through
because they read as generic rather than as host contexts.

*P. falciparum* in vivo fitness IS measurable — humanised SCID mice, controlled human malaria
infection — so the right correction is to re-context, not to delete. For rodent malaria the organ
contexts are real, which is the second finding's business.

### Instruction 39's transfer slots were never built

Instruction 39 specifies them explicitly: "`Pf_fitness · transferred from Pb` is a slot, marked as a
transfer, `derived_from` the *berghei* column it came from", with the leakage rule that
`target_family` closure must span species. **There are zero transfer slots in the catalog.**

That matters now rather than in the abstract, because the data is sitting there: PlasmoDB serves
`GenesByPhenotype_pberANKA_phenotype_Bushnell_functional_profiling` — the PlasmoGEM *berghei*
knockout growth phenotypes — and the *knowlesi* piggyBac screens beside it. None can fill a
*falciparum* slot directly, and all three could fill transfer slots that do not exist yet.

### Why both are recorded rather than fixed

Re-contexting the four mouse-organ slots would **reduce the Plasmodium denominator from 103**, which
would raise the filled percentage without filling anything. Adding instruction 39's transfer slots
would **raise** it. Doing the first without the second, at a point where the count is what is being
asked about, is not a call to make unilaterally — so the honest move is to write both down together
and let them be decided as one change.

The next pass should do both at once: convert the four mouse-organ contexts to the in-vivo settings
*falciparum* is actually studied in, and add the transfer slots instruction 39 designed, then fill at
least one from PlasmoGEM. That is a net increase in questions asked, and it is the shape instruction
39 settled with the user.

## Twenty-sixth pass: both catalog faults fixed together, denominator unchanged

**Toxoplasma 114 of 119, Plasmodium 34 of 103, combined 148 of 222 — the same numbers as before the
change, which is the point.**

The previous pass recorded two faults and deliberately did not act, because fixing one without the
other would have moved the count in one direction for no reason. Both are now fixed as one change and
the totals are untouched.

**`mouse` joins `TOXO_ONLY_CONTEXTS`.** The four mirrored mouse-organ slots are gone: *P. falciparum*
is a human parasite and does not infect mice, so `fitness · in vivo peritoneum / lung / liver /
spleen` were unanswerable by construction rather than empty. The word excludes exactly those four --
the three brain slots were already caught by `brain` and the screen-specific ones by
`STUDY_SPECIFIC`.

**In their place, four slots that can be asked.** One `fitness · in vivo` with the setting falciparum
is actually studied in (humanised mouse, controlled human malaria infection), which is one setting
rather than four organs; and the three cross-species transfers instruction 39 specified and the
catalog never had:

* `Pf_fitness · transferred from Pb`
* `Pf_fitness · liver stage transferred from Pb`
* `Pf_fitness · transmission transferred from Pb`

Those exist because the data exists and cannot fill a *falciparum* slot directly. PlasmoGEM's
*berghei* knockout phenotypes and the *knowlesi* piggyBac screens are the best evidence there is for
falciparum liver-stage and transmission fitness, and transferring them through orthology is a claim
that has to be **visible** rather than folded into the measured slot.

Each transfer shares its `target_family` with the measured slot it stands in for, which is instruction
39's leakage rule in its own words: *"target_family closure must span species. Transfer berghei fitness
onto falciparum, hold out falciparum fitness, and 'recover' it, and you have measured orthology, not
biology."* Two tests now hold it -- one asserts no transfer slot is alone in its family and that each
declares itself orthology-derived, the other that no Plasmodium slot asks about a mouse organ.

−4 + 1 + 3 = 0. The catalog got more correct and the score did not move.

### The source for the transmission transfer slot, found and not yet used

`Pf_fitness · transmission transferred from Pb` now exists and its data is located: **PMID 37708854 /
PMC7618085**, a genome-scale *berghei* screen for gametocyte-to-sporozoite transition. Supplementary
`mmc2.xlsx` sheet `Table S1A` carries **3,428 PBANKA genes** with a log2 ratio and a p-value, reached
through the nested `EMS208536-supplement-Supplementary_Material.zip`.

Three things stand between that and a filled slot, and they are the work, not the finding:

1. **What the ratio means.** The column is headed `WT/2.33 avg Log2`. `2.33` is a line or pool name
   and the direction has to be established from the paper before the column can be named, exactly as
   the ToxoDB `fold_change_avg` direction and the vesicle EV-minus-cell transform had to be.
2. **A PBANKA -> PF3D7 map.** The NF54 map cannot be reused. That one is an IDENTITY map validated on
   protein length because 3D7 was cloned from NF54; *berghei* is a different species and its
   orthologues will differ in length, so length is the wrong check here. The right one is orthogroups
   holding exactly one gene on each side, with the SAME conservative drop of multi-gene groups, and a
   different validation -- probably that known-conserved essential families map and the variant
   surface families do not.
3. **The transfer has to stay visible.** The slot already declares itself orthology-derived and shares
   `target_family` with `Pf_fitness · transmission`, so the leakage rule is in place and tested. The
   column must be named for what it is -- transferred *berghei* fitness -- and never merged into the
   measured slot.

`Pf_fitness · transferred from Pb` and `· liver stage transferred from Pb` still need their own
sources; PlasmoDB's `GenesByPhenotype_pberANKA_phenotype_Bushnell_functional_profiling` is the
blood-stage one and sits behind an EDA `filter` parameter rather than a plain attribute, which is why
it was not fetched here.

### Correction to the entry above: that source is the wrong KIND of measurement

Checked before use, and refused. `WT/2.33 avg Log2` is **differential gene EXPRESSION** between two
*berghei* lines — ANKA 2.34, which produces gametocytes, against ANKA 2.33, which does not — and the
positive tail is genes enriched in gametocytes. It is a transcriptome comparison, **not a knockout
fitness screen.**

So it cannot fill `Pf_fitness · transmission transferred from Pb`, which is a fitness slot. The
paper's title — "Identification of genes required for *Plasmodium* gametocyte-to-sporozoite
transition" — reads exactly like a screen, and the previous entry here described it as one on the
strength of a per-gene log2 with a p-value. That was wrong and is corrected rather than deleted,
because the mistake is the instructive part: **a log2 with a p-value tells you nothing about what was
measured.** Every source this campaign refused had one.

It also cannot be redirected to a transcription slot. Gametocyte-enriched expression is
`Pf_transcription · gametocyte`, which the Su seven-stage series already fills, and using this beside
it would be one question answered twice.

**What the transfer slots actually need is a knockout screen scored in the mosquito**, which is what
PlasmoGEM's *berghei* transmission data is. That remains behind PlasmoDB's EDA `filter` parameter.
The three transfer slots stay empty and correctly so; the machinery to hold them honestly is built and
tested, and the next attempt should start by confirming that a candidate table's rows are MUTANTS
rather than transcripts.

## Twenty-seventh pass: measured contacts, with the host end thrown away

**Toxoplasma 114 of 119, Plasmodium 35 of 103, combined 149 of 222.**

`Pf_interaction · crosslink MS` from PMID 41966402: **73 parasite-parasite pairs over 102 genes**,
built as an `xlms` edge layer named to match the Toxoplasma one that answers the same slot.

Two filters, neither optional. The experiment crosslinked **parasite inside erythrocyte**, so a third
of the 106 protein pairs have a human protein at one or both ends — spectrin, band 3, protein 4.2.
Those are real contacts and they are a **host bridge**, not a parasite-parasite edge, so they are
dropped from this layer rather than indexed against a table that has no row for them. And a pair whose
two ends resolve to one gene is a homomeric crosslink: evidence the protein self-associates, not an
edge between two genes, and drawing it would put a zero-length line in the graph.

Validated on complexes that have to be there rather than on a count: **EXP2, PTEX150 and HSP101 —
three subunits of the PTEX translocon — crosslink to one another**, prohibitin 1 to prohibitin 2, and
RAP1 to RAP2. A contact map that missed those would not be measuring contacts, and the test asserts
the PTEX pair specifically.

Note what this source also contains and what was not taken: sheet `(F)` is 360 PPIs at 5% FDR and
`mmc5` is complex clusters, which is `Pf_complex membership` — a second slot, from the same study,
needing its own reading. Left for the next pass rather than half-done.

## Twenty-eighth pass: complexes from the same study, as a different question

**Toxoplasma 114 of 119, Plasmodium 36 of 103, combined 150 of 222.**

`Pf_complex membership` from the study whose crosslinks filled the pair slot last pass — the same
source, a different question. That one asks which pairs touch; this asks which assembly a protein
sits in. **128 genes in 42 complexes.**

`complex_spans_host` is the informative column and is kept rather than dropped. Seven of the 47
clusters contain human proteins as well as parasite ones, which is not contamination: the experiment
crosslinked parasite inside erythrocyte, so a complex reaching into the host is a finding. But a
parasite gene in one of those has partners this table cannot name, and a reader taking `complex_size`
at face value would over-count its parasite neighbours. Only parasite members get a row; the host
members belong to a bridge table.

The documented Plasmodium column count in HANDOFF drifted by four during this pass and the
species-aware column test caught it immediately — which is what that test was extended for two passes
ago, so it has already paid for itself.

## Twenty-ninth pass: a bug in last pass's own layer, found by reading the next slot

**Toxoplasma 114 of 119, Plasmodium 36 of 103, combined 150 of 222 — and the crosslink layer went
from 73 edges to 79, because six of them should never have been dropped.**

Setting up `Pf_interaction · with host proteins` meant looking at the host-parasite crosslinks that
the `xlms` layer had deliberately excluded. One of the "host" proteins was `sp|Q6ZMA7|Pfs16`.

**Pfs16 is a *Plasmodium* gene.** Q6ZMA7 resolves to `PF3D7_0406200`.

The layer identified which end was host by matching `PF3D7_\w+` in the mapping field. The source
writes some rows with a UniProt symbol instead of the accession, so those rows failed the pattern and
were classified as host — which **dropped six real parasite-parasite contacts** and would have put a
parasite protein into a host bridge as though it were human. Resolving through the accession index
recovers all 79.

This is precisely the failure the Toxoplasma identity layer exists to prevent, and it happened anyway
because **a regex on an accession field looks like resolution and is not.** The lesson instruction 41
already recorded for deposits keyed on `TGGT1_` applies to any field that carries identifiers in more
than one notation, which is most of them. A test now fails if the count drops back below 79.

Worth noting how it surfaced: not from a validation of the layer itself — the PTEX check passed on 73
edges just as it does on 79 — but from starting the *next* slot and finding the discarded pile had a
parasite protein in it. The thrown-away half of a filter is worth reading.

## Thirtieth pass: the discarded half becomes a bridge, and the species guard reaches its fourth place

**Toxoplasma 114 of 119, Plasmodium 37 of 103, combined 151 of 222.**

`Pf_interaction · with host proteins` is filled from exactly the crosslinks the `xlms` layer throws
away: **10 parasite-to-human pairs, 7 human proteins.** It is a bridge and not an edge for the reason
instruction 39 gives — the pair's two ends live in different tables and a human protein has no index
in the parasite one.

Second bridge in the project, after the Toxoplasma host IP-MS one, and building it needed a
**species-aware bridge lookup**: both arms key their bridge `host`, because both cross to a human
protein, so the name cannot say whose contacts these are and only the parasite end can. That is the
**fourth** place the species guard has had to go — `declared_columns`, `resolve`, `is_filled`'s edge
branch, and now its bridge branch — and the cause is identical each time: a vocabulary deliberately
shared between the arms so they can be compared.

The rule is now worth stating as a rule rather than four incidents: **whenever the two arms share a
name, the thing that distinguishes them must be the data, not the name.** Columns, graph layers and
bridges have all needed it. The next shared vocabulary will need it too.

Validated on an interaction that is in the textbooks: **MESA crosslinks to erythrocyte ankyrin**, and
the rest of the human side is stomatin, calpain, actin and spectrin beta — the membrane skeleton,
which is what an exported parasite protein should be touching.

And the pass before this one is why the count is 10 and not 15: reading the discarded pile through the
accession index rather than a pattern moved five pairs back to the parasite side where they belonged.

## Thirty-first pass: adding three columns deleted three, silently

**Toxoplasma 114 of 119, Plasmodium 38 of 103, combined 152 of 222.**

`Pf_transcription · noncoding and antisense transcription` is filled from the antisense partners of
the sense columns already in the table — same runs, so directly comparable. They sit **beside** their
sense partners rather than being reduced to a ratio, because the denominator is what makes a ratio
interpretable and a reader should see both.

Validated on what antisense has to be: median **1.58 against 12.70** for sense in the same samples,
an eight-fold minority strand. The highest antisense sits on U6 spliceosomal RNA, SRP RNA and rRNA
fragments — structured non-coding RNAs where strand assignment is genuinely ambiguous, which is a
caveat rather than a fault and is worth knowing before reading the column.

### The fault worth more than the slot

Fetching the antisense columns made the slot count **fall from 151 to 149** while a slot was being
added. The matcher used plain substring containment, and

> `sense - asexual blood stages` is a SUBSTRING of `antisense - asexual blood stages`

so each sense entry suddenly matched two headers, failed its one-match test, and vanished. **Adding
three columns deleted three**, with no error — the loader's own guard against ambiguity is what did
the deleting.

Two things now prevent it. The matcher requires the label to start the header or follow a
non-alphanumeric character, and a test asserts that **no declared sample is absent from the built
table** — which is the check that would have caught this instantly and which every loader with a
declared column list should have.

Noticing it at all depended on watching the total go the wrong way. A pass that only checked "did the
new slot fill?" would have shipped the loss.

## Thirty-second pass: the first measurement the two arms can be compared on

**Toxoplasma 114 of 119, Plasmodium 39 of 103, combined 153 of 222.**

`Pf_codon usage / translation efficiency`, computed through the **same code as the Toxoplasma arm**
rather than reimplemented — `codons.codon_usage` gained a table and node-table parameter, the way
`cellcycle.stage_enrichment` gained `stages`. ENC and GC3 are definitions and the CAI reference set is
"the ribosomal proteins" in both arms, so two implementations could only differ by being wrong in one
of them.

Sharing it buys something no other column in this project has yet had: **a measurement the two arms
can be compared on, where the comparison is itself the validation.**

| | Plasmodium | Toxoplasma |
|---|---|---|
| GC3 median | **0.150** | 0.583 |
| ENC median | **37.6** | 53.9 |

*P. falciparum* has the most AT-rich genome of any eukaryote, so extreme codon bias — a low ENC — is
what has to appear, and it does. A test fails if the two arms ever converge, which would mean either a
genome was mixed up or they stopped computing the same quantity. That test is only meaningful because
the construction is shared; had each arm computed its own ENC, a difference would have been
uninterpretable.

Practical note: PlasmoDB's sequence report returned 422, 400 and 500 to three different request
shapes. The static release FASTA is what works, and 5,389 transcripts reduce to 5,318 genes.

## Thirty-third pass: two EC fields that are two kinds of evidence

**Toxoplasma 114 of 119, Plasmodium 40 of 103, combined 154 of 222.**

`Pf_enzyme classification`. PlasmoDB serves **two** EC fields, and the distinction is the whole
decision: one is curated for this organism, the other inferred from the gene's OrthoMCL group.

They are separate columns and `has_ec` counts only the curated one. Merged they would read as 1,584
genes with an EC and give no way to tell which **335 were never annotated in this organism at all** —
inference standing exactly where annotation should. The slot lists the curated column first and its
policy is `one`, so the leading candidate wins and the derived field is available to be chosen
deliberately rather than by default.

For scale: 1,220 curated here against 1,313 in Toxoplasma, on a proteome two-thirds the size.

This is the same distinction as `localization · measured` versus `· transferred`, and as the
`fitness · transferred from Pb` slots added earlier — a value inferred from orthology is not the same
claim as one measured in the organism, and the map's job is to keep saying which is which.

## Thirty-fourth pass: host degree, and a third state the classifier was missing

**Toxoplasma 114 of 119, Plasmodium 41 of 103, combined 155 of 222.**

`Pf_host interaction degree`, derived from the host bridge: **117 genes seen in the crosslink data, 10
with a host partner.**

Not a `fillna(0)`. The Toxoplasma column of the same name IS 0 everywhere without a curated host
target, and that is right there — its source is a curated table covering the literature. This source
is **one experiment**, so a gene it never detected has not been shown to lack host partners. The 117
genes seen carry a count (zero included, because being crosslinked only to parasite proteins is a real
observation) and the other 5,603 stay missing.

### The third state

Writing the test for that exposed a flaw in the classifier all three crosslink consumers share. It
asked "is this a parasite gene I know?" and treated **no** as host. But there are **three** answers:

* a known parasite gene,
* a **parasite** protein the node table does not carry — a deprecated accession, a gene dropped from
  the annotation,
* genuinely host.

Reading the middle case as host inflated host degree and would have put a parasite protein into a host
bridge. Reading it as parasite would index a row that does not exist. Either way the pair is unusable,
so `_classify` now reports it as parasite-with-no-gene and every consumer skips it.

**The shipped numbers were unaffected** — every accession in this file is in the table, so 79 edges and
10 bridge pairs are the same before and after. This is a fix for the next file, not a correction of
this one, and it is worth saying which: the previous two passes found faults that HAD cost real data,
and conflating the two kinds would overstate this one.

## Thirty-fifth pass: translation efficiency computed and refused

**No change to the totals: Toxoplasma 114 of 119, Plasmodium 41 of 103, combined 155 of 222.**

`Pf_translation efficiency · asexual blood stage` has everything it needs in the table already —
polysomal and steady-state arms of one experiment at three IDC stages — so the ratio is one line. It
was computed, checked twice, and **not shipped**, because the two checks disagree.

First, a scale error worth recording because it nearly went unnoticed: the arms are LINEAR (0 to
12,988), not log, so a difference is not a ratio. `log2((poly+1)/(steady+1))` is the quantity.

**The check that passes.** Ribosomal proteins are the most heavily translated things in any cell, and
their TE is **+0.805 against −0.741** for everything else (p = 1e-25). That is exactly right.

**The check that fails.** TE correlates **negatively** with codon adaptation — rho −0.10 to −0.19
across the three stages — where the textbook expectation is positive: heavily translated genes carry
optimal codons. It is worse than merely negative, because CAI's reference set here **is** the ribosomal
proteins, so the genes with the highest TE are the genes that define high CAI, and the global
correlation still runs the other way.

Something dominates that relationship which I cannot name, and in this genome there are candidates —
whether codon bias predicts expression in *P. falciparum* at all is contested, and CAI in a genome
this AT-rich may track base composition more than translational demand. But "there are candidate
explanations" is not the same as knowing which, and a column called *translation efficiency* asserts
that the number measures translation efficiency.

Under the rule this campaign settled, two checks disagreeing is the **contradictory** case, whose
answer is to ship the conditions rather than the derived contrast. Both conditions are already shipped
and already answer two other slots, so there is nothing left to ship one step back — which makes this
the case where the rule resolves to a refusal.

Twelfth source built or computed and then refused. The slot keeps its verdict of missing, and what it
needs is a translation-efficiency measurement whose relationship to codon usage someone has already
had to explain.

## Thirty-sixth pass: `drug sensitivity` re-framed, and the one path left

**No change: Toxoplasma 114 of 119.**

The "suspect the sentence" test applied to `drug sensitivity`, since it worked for `lipid composition`
and for `resistance-conferring mutation`. The sentence recorded against this slot is *"no genome-wide
chemogenomic screen has been published for Toxoplasma"*, and every sweep tested exactly that.

But the slot asks whether **disrupting a gene changes survival under a compound**, and a genome-wide
screen is only one instrument that answers it. Individual knockout studies answer it one gene at a
time, and there are many.

Searching for those turned up one promising hit — *"A combined genetic and chemical approach for
identifying novel antifungal compounds"*, open access, a systematic platform screening knockout mutants
against a 2,704-compound library. It is ***Fusarium graminearum***. A fungus. Caught by reading the
abstract rather than the title, which is the fourth wrong-organism candidate this campaign after
*Dictyostelium*, *Cryptosporidium* and *Theileria*.

**What the reframing does leave is a curation job with a statable bar**, and it is the same shape as
the one that filled `resistance-conferring mutation`:

> a gene is in if a knockout or knockdown shows a MEASURED shift in sensitivity to a NAMED compound,
> read from primary text, with the direction and the compound recorded per row.

Candidates seen while searching: equilibrative nucleoside transporters against purine analogues
(PMC12607627), the ZFT iron/zinc transporter (PMC12875612), TgGSK3 (PMC12589562), Aurora kinases
(PMC12707362) — all open access. That would be perhaps ten to twenty genes, sparse because the biology
is, and it needs the same per-row product cross-check that makes a curated table safe.

It is not started. What is now recorded is that the sentence was wrong in the same way the lipid
sentence was — "no screen exists" is not "the question cannot be answered" — and that the remaining
path is curation rather than search, so the next pass need not sweep the catalogues a seventh time.

## Thirty-seventh pass: `drug sensitivity` filled — the sentence really was the problem

**Toxoplasma 115 of 119. Plasmodium 41 of 103. Combined 156 of 222.**

Last pass re-framed this slot and stopped at "the remaining path is curation". This pass took it.

**PMID 41025776** deletes three equilibrative nucleoside transporters and measures each against two
toxic nucleoside analogues. That is precisely what the slot asks — does disrupting a gene change
survival under a compound — and it had been invisible for six sweeps because every sweep searched for
a *screen*.

| gene | | Ara-A | 5-FU | background |
|---|---|---|---|---|
| TGME49_244440 | TgAT1 | resistant | resistant | parental |
| TGME49_233130 | TgENT3 | resistant | sensitive | ΔTgAT1 |
| TGME49_500147 | TgENT2 | unchanged | unchanged | parental |

Three decisions in that table:

* **`unchanged` rows are kept.** A transporter deleted with no effect on analogue sensitivity is a
  result. Dropping those two rows would leave the column looking like a list of hits, which is how a
  curated table starts lying.
* **TgENT3's rows say `ΔTgAT1`.** They were measured in a double knockout, so they are
  genetic-interaction results. Reading a double mutant's phenotype off one of its genes is its own
  error, and `background` is the column that stops it.
* **Target engagement is excluded and tested for.** "An inhibitor of this protein kills the parasite"
  is a different slot; a test asserts no row's evidence reads that way. Two of the four candidates
  named last pass — TgGSK3 and the Aurora kinases — are exactly that and were not used.

The per-row product check earned its keep immediately: the annotation calls TGME49_244440 **"adenosine
transporter AT1"**, which is independent confirmation the accession is TgAT1 and not a transposition.

### What this changes about the four that remain

The claim "five slots need experiments nobody has run" was wrong, and it was wrong for six passes. One
of the five needed *reading*, not an experiment. The four left are:

`transcription · in IFN-gamma macrophage`, `translation · per cell-cycle phase`, `protein turnover`,
`fitness · in vivo gut`.

Those are still `missing`, and two of them were closed by catalogue-complete sweeps rather than by
phrasing. But the drug-sensitivity case is a standing warning against trusting that: **a sweep tests
the sentence it was given.** For each of the four, the sentence to suspect is written in the slot table's
`searched` field, and the question to ask is what instrument OTHER than the one swept for could answer
the slot.

## Thirty-eighth pass: `fitness · in vivo gut` filled, and the catalogue that was never complete

**Toxoplasma 116 of 119. Plasmodium 41 of 103. Combined 157 of 222.**

Last pass ended with a warning against my own sweeps: *a sweep tests the sentence it was given.* This
pass turned that on the four remaining slots and it cost three of the four their verdicts — not
because the data appeared, but because the sweeps had been narrower than the claims I wrote from them.

### The catalogue that was complete for the wrong set

`protein turnover` was closed on "all 201 Toxoplasma deposits in PRIDE, enumerated not sampled". That
is catalogue-complete **for PRIDE**, which is not the same statement as complete for proteomics.
ProteomeXchange aggregates five repositories, and Toxoplasma has **233** deposits across them — 192
PRIDE, 29 iProX, 6 MassIVE, 5 jPOST, 1 PeptideAtlas. Forty-one deposits my sweep never saw.

None measures turnover. The negative survived, but it had been resting on a smaller catalogue than I
claimed for six passes.

### The trap in the seven numbers

Then the literature instead of the archives: 467 open-access Toxoplasma papers mentioning
cycloheximide, half-life, protein stability or turnover, full text fetched and grepped for a half-life
attached to a gene name. Seven statements matched.

**All seven are ExPASy ProtParam predictions.** Six say exactly `30 h`; the seventh says `> 20 h
(yeast)`. ProtParam returns a constant keyed on the N-terminal residue, and every hit is an in-silico
vaccine-design paper reciting it next to a molecular weight and an extinction coefficient.

Ingesting them would have filled `protein turnover` for eight genes with a column that is **one amino
acid in disguise** — and it would have looked like the campaign's cheapest win. This is the label trap
in its purest form: not a mislabelled dataset, a *predicted* quantity wearing a measured quantity's
name. Refused, and the reason is now written into the slot's verdict so the next sweep doesn't
rediscover them as a find.

### The sixth label lie

Widening GEO the same way — no organism filter, because a dual RNA-seq series tagged host-only is
invisible to one — surfaced **GSE204926**, titled *"Toxoplasma IWS1 determines fitness in
interferon-γ-activated host cells and mice"*, organism tagged `Toxoplasma gondii`. Exactly the
missing slot, by its title.

Its samples are `WT parasite_1..3` and `IWS1 KO parasite_1..3`, and the design says *"freshly isolated
tachyzoites"*. It is a genotype contrast in extracellular parasites. The IFN-γ in the title is what
the gene is *for*, not what the experiment *varied*. Refused — the sixth time a source's own label
named a measurement it had not made.

### Closing a slot on set membership instead of a phrase

`translation · per cell-cycle phase` was closed by reading nine ribo-seq series. Better: intersect two
catalogues. Every BioProject matching ribosome/polysome/translatome (7) against every BioProject
matching synchronised/cell-cycle/sorted (29). **The intersection is empty.** No submission is both.
That closes the slot on a set operation rather than on a judgement about phrasing, and it cannot rot
the way a sentence can.

### The slot that fell

`fitness · in vivo gut` was closed on "no pooled screen through the enteroepithelial stages" — true,
and it will stay true, because the sexual cycle runs only in a felid and nobody has put a barcoded
library through a cat. But the slot asks whether disrupting a gene costs the parasite oocysts, and a
pooled screen is one instrument that answers it. **Feeding one knockout to a cat and counting what it
sheds is another,** and four labs have done it.

| gene | | oocyst yield | sporulation | background |
|---|---|---|---|---|
| TGME49_227100 | Grx5 | reduced, 2.5e7 → 8e6 | 80% → 30% | parental, Pru |
| TGME49_287510 | AAH1 | reduced, ~1e3 vs 1e6–1e7 | not measured | parental, ME49 |
| TGME49_212740 | AAH2 | reduced, ~10-fold | 75–80% → ~60% | parental, ME49 |
| TGME49_212740 | AAH2 | reduced | not measured | **ΔAAH1** |
| TGME49_285940 | HAP2 | reduced, <1,600 vs 11e6 | **abolished, 0%** | parental, CZ H3 |
| TGME49_276850/60/70/80 | LEA850–880 | **unchanged, 30 vs 34 million** | unchanged | ΔLEA cluster |

Decisions worth keeping:

* **The four `unchanged` rows are the strongest in the table.** All four LEA genes were deleted
  *together*, and oocyst yield did not move. A joint deletion that changes nothing proves more than a
  single knockout could: redundancy cannot be hiding the phenotype. The authors also report the
  sporulation dip (~65% vs ~85%) as inside the 60–90% range they see normally, so it is recorded as
  unchanged rather than mined as an effect.
* **Yield and sporulation are separate columns**, because HAP2 sheds a few mis-shapen oocysts that
  *never* sporulate while Grx5 sheds a third as many that sporulate poorly. One column would erase
  the difference between a fertilisation block and a maturation defect.
* **Magnitudes stay in `evidence`, not in a numeric column.** Cats, strains and inocula differ across
  four papers; a float here would invent a precision the experiments do not have.
* **Two papers were found and refused** for failing the bar: one says the cat experiment "should be
  carried out", one says oocysts were seen but "the numbers were not quantified". A test asserts
  neither phrasing can reappear as a row.
* **The HAP2 caveat is carried, not dropped.** That line has a second mutation in an intron of
  TGME49_223060 which the authors disclose and argue against. A caveat the source volunteers is not
  mine to silently delete.

### The identity check that nearly failed

The AAH paper gives `TgME49_212740` for AAH2 and **never gives AAH1's accession**. The node table's
product strings were no help either: I searched `hydroxylase` and got nothing, because ToxoDB spells
both products *"aromatic amino acid hydro**l**ase"*. The paper also says the two genes are "very
closely related and located on chromosome V", and the accessions I had — 212740 and 287510 — are
nothing like adjacent, unlike the tandem LEA cluster 276850/60/70/80.

ToxoDB's record endpoint returned 500 for every attribute, so the check had to come from the local
identity index, which carries a `gene_name` field: **AAH1 = TGME49_287510, AAH2 = TGME49_212740**.
Authoritative, and independent of the product string that hid them. This is the standing rule paying
out again — *resolve through the identity index, never a regex and never a guess.*

### Rebuilding it exposed two silent-loss bugs and a five-tree layout

The rebuild needed `toxonet/data/interim/nodes.parquet`, which lives outside the repo, so I pointed
`STARPLAST_DATA` at the tree that has it. The build **succeeded, exit 0**, wrote its parquet, reported
`8140 nodes, 13 edge types` — and carried 305 columns instead of 390. Eighty-eight datasets had
vanished. Nothing errored, and the three new columns were present and correct, so any check that only
asked *"did my slot land?"* would have passed while the table lost 22% of itself.

Four more attempts landed on 352, 334, 334 and 391 before the layout was right. `BASE` is the parent
of the dataset root, and loaders join it five different ways, so the build needs ONE directory holding
all of:

`starplast/` · `datasets/` · `toxonet/` · `toxo_stage_atlas/` · `starpath_*` · `.claude/skills/…/corpus/`

They were split across two parents. Symlinked into place, and the required paths are now checked
explicitly rather than discovered by a column count going down.

Then the last two columns refused to build, and each was a real bug:

* **`proteomics.load_all` was called without `resolve`.** Two of the quarantine tables it reaches are
  keyed on `TGGT1_`; without the identity layer they join zero rows. That is 7,384 and 2,348 values,
  and the build still exited 0. One-line fix, and the restored coverage matches the old numbers
  *exactly*, which is what says the fix is right rather than merely non-empty.
* **`add_strain_accessions` contradicted its own module docstring.** The documented rule is that a
  cross-strain accession resolves *when its suffix exists in the node table*. The implementation only
  registered accessions enumerated in the strain TSV — and that file is incomplete: `TGGT1_212960`,
  `_251570`, `_297960` and `_310430` are absent from its 8,637 rows while all four ME49 genes are in
  the node table, so four genes of the splitCas9 screen silently failed to join. Now completed by
  suffix, which can only add: each suffix maps to exactly one ME49 gene.

Those two fixes gave **14 columns more coverage than the table had before**: the GSE99395 ribo-seq arms
went 7,437 → 7,506 genes each, `cdpk1_thiophospho_peptides` 361 → 367, the differentiation reporter
235 → 237. A slot-filling pass found bugs in the identity layer because it insisted the rebuild lose
nothing.

**The check that caught all of it**, and the one to keep running: diff the rebuilt table against
`HEAD` on *columns gained, columns lost, and per-column coverage in both directions*. Column presence
alone would have missed the 4-gene screen regression and the 1,306-value `best_model_agreement` wipe,
because both columns still existed. A build that exits 0 is not evidence that a build was correct.

### Two tests that were passing for the wrong reason

`test_a_missing_study_file_reports_and_returns_empty` and `test_a_missing_primary_matrix_reports_rather_than_raising`
set `STARPLAST_DATA` to an empty directory and then assert the loader reports a missing file. But
`paths.find` walks *every* root it can name, and one is derived from the package location, not the
environment — so neither test ever guaranteed the absence it is named for. They passed because the
repo's dataset tree happened not to contain those two files, and they broke the moment it did. Both
now patch `paths.dataset_roots` to name exactly one root. Setting an env var is not isolation when the
resolver has other candidates.

## Thirty-ninth pass: the queue was partly a bug, and the first measured translation on this arm

**Toxoplasma 116 of 143. Plasmodium 44 of 128. Combined 160 of 271** — per unit, gene 142/198,
pair 15/19, metabolite 3/6, host_gene 0/48.

### The work queue this session opened with was partly an artefact

The handoff said seventeen Plasmodium slots were past discovery because they already carried
resolved candidate studies. Reading the candidates rather than counting them says otherwise, in two
separate ways.

**Seven of them were Toxoplasma papers**, published on Plasmodium slots by the generator.
`NEW_SHARED` holds the questions both parasites have, and BOTH of its per-slot fields are about
Toxoplasma: a pattern names a column of the Toxoplasma node table, and a citation was found while
filling the Toxoplasma copy of the question. The patterns were already stripped on the way to the Pf
arm — `codon_` flipped the Pf codon-usage slot the day three Toxoplasma sequence columns arrived —
and the candidates were not. So `Pf_resistance-conferring mutation` cited two TgMAPK1 papers,
`Pf_invasion and egress phenotype` a Toxoplasma splitCas9 screen, `Pf_complex membership` a
Toxoplasma crosslinking interactome, `Pf_seroreactivity` the Toxoplasma IEDB query,
`Pf_host ESCRT recruitment` an unpublished Toxoplasma imaging screen, and
`Pf_essentiality in a second background` a Toxoplasma differentiation-reporter series.

That is the **fifth** place this same rule has had to be enforced, after `declared_columns`,
`resolve`, `is_filled`'s edge branch and its bridge branch. Where the two arms share a NAME, only the
DATA may say which parasite something is about — and a citation's organism is data too.

**The rest carry `automatically proposed; verify assay and parasite-gene shape`, and the warning is
earned.** `Pf_RNA-binding protein targets` proposes two *Anopheles* antibacterial-immunity papers;
`Pf_invasion and egress phenotype` a *Chromera velia* motility paper; `Pf_essentiality in a second
background` a leishmaniasis diagnostic; `Pf_transcription · liver stage` a human hepatocyte circadian
study; `Pf_target engagement` an essential-oil antioxidant screen. Discovery on those slots is not
done. It was never started.

**What IS past discovery is a different file.** `instructions/open/41_candidates_v2.json` holds 76 Pf
keys of GEO- and PRIDE-indexed candidates WITH accessions, and the atlas does not merge it — it
merges `31_candidates.json`, the PubMed-proposal file. 39 of the empty Pf slots appear there. That is
the queue worth working, and it still needs each candidate read rather than counted: its
`S-nitrosylation` and `chromatin accessibility` lists are RNA-seq studies, because the proposer
matched an axis term rather than an assay.

Two tests now, because one of them cannot see the whole fault. The structural one rebuilds the Pf
catalog with the candidate file out of the way and refuses any candidate that was not written for the
Plasmodium arm; it catches all six slots. The text one reads the shipped catalog for a whole-word
Toxoplasma name and catches exactly one, because the ESCRT citation names no organism at all.

### `Pf_translation · per cell-cycle phase`, filled by ribosome profiling

**GSE58402** (PMID 25493618), five points of the intraerythrocytic cycle, both arms of the
experiment: ribosome footprints and matched mRNA, as RPKM per gene. 3,501 genes, 61%, and uneven
across the cycle — 2,182 at the ring, 1,174 at the merozoite. The slot was previously answerable
only through polysome-associated RNA, which is what is ON ribosomes rather than how much ribosome is
on it; this is the first MEASURED translation on this arm.

The mRNA arm joins `Pf_transcription · per cell-cycle phase` as a second dataset beside the
polysomal study's steady-state arm — the same question, a second instrument — and never as a second
slot. One experiment's two conditions do not become two questions.

**Stage labels are the deposit's own, so they were checked against a study sharing no sample.** Each
of ring, early trophozoite, late trophozoite and schizont correlates highest with its own stage in
the independent PlasmoDB seven-stage series (rho 0.65–0.77). The merozoite has no counterpart there
and lands on the ring, which is the neighbouring point of the cycle rather than a contradiction — a
free merozoite is a ring that has not invaded yet — and the test asserts the four, not the five.

**The measurement behaves**: ribosomal proteins carry far more footprint than everything else at
every stage (median log1p 5.6–7.9 against 3.5–4.1, p ≤ 2e-18).

Strain **W2**, not 3D7. The gene set is 3D7's, so nothing is mis-joined, but the surface-antigen
families are where a strain difference would show and that is where to distrust this column.

### The second arm gets an identity layer, because the first source that needed one found zero genes

The deposit is keyed on **pre-2012 accessions** — `PFE0630c`, `PF13_0222` — and the Plasmodium arm
had no identity layer at all. A string join found **0 of its 3,629 ids**. This is precisely the
failure the Toxoplasma layer exists to prevent (the 2019 in-vivo screen, 0 of 8,140 unresolved and
168 resolved), met on the second arm for the first time.

`fetch_names` now fetches PlasmoDB's identity report through the same WDK service and the same code
path, into its own committed file — two identifier spaces in one index is the merge this project
refuses everywhere else. 9,106 previous ids resolve; **66 are claimed by two current genes each**,
which is what a gene model being SPLIT looks like from the other side, and those are withdrawn rather
than assigned to whichever row came first. A current accession outranks another gene's history, since
most files mix both forms.

Three things are refused rather than resolved, about 100 ids of 3,629: the ambiguous ones above, the
deposit's `-a`/`-b` split entries (RPKM is already length-normalised, so neither summing nor
averaging two segments means anything), and two source ids landing on one current gene inside one
file, which is a merge with the same problem.

### Translation efficiency, computed a second time and refused a second time

The thirty-fifth pass computed TE from the polysomal and steady-state arms and refused it: ribosomal
proteins came out right, but TE correlated NEGATIVELY with codon adaptation where the textbook
expects positive, and the entry closed with *"something dominates that relationship which I cannot
name"*.

Ribosome profiling reproduces it — different instrument, different strain, a decade earlier:
rho(TE, CAI) is −0.01 to −0.12 across the five stages, and the RPF LEVEL itself is −0.15 against CAI
at every stage, so it is not the ratio's fault. The ribosomal-protein check also inverts at the
schizont (TE −0.65 against −0.12 for everything else). Two checks disagreeing is the **contradictory**
case, whose answer is to ship the conditions rather than the contrast. Both conditions ship. TE does
not.

**But the replication says where to look, and this one is worth keeping.** `codon_cai_ribosomal` is
built against each arm's own ribosomal proteins by one shared construction. In Toxoplasma the
reference set scores at the top of its own index (0.771 against 0.714, p = 4e-22). In *P. falciparum*
**it does not separate at all** (0.710 against 0.717, p = 0.37), while ENC does, weakly (36.5 against
37.6). So the quantity failing to behave in the TE check is the CODON INDEX, not the footprints —
an AT-rich genome compressing codon-usage signal, not a bug in either arm. That is only visible
because the two arms share one implementation; two would have made it noise. A test now pins both
directions.

### The rebuild, and the check that has to travel with it

`scripts/build_plasmodium.py` is new, and it exists because `build_all` made the TABLE reproducible
while leaving the INVOCATION to whatever was typed at a prompt — including the diff, which is the
part that matters. It builds, diffs against the shipped cache on columns gained, columns LOST and
per-column coverage in BOTH directions, and **refuses to write** if anything went backwards;
`--allow-loss "why"` is the override and it requires a sentence, so that a loss someone accepted and
a loss nobody noticed cannot look the same in the history.

This build: **+10 columns, −0, 5,720 genes unchanged, no column's coverage moved in either
direction.** That is what a correct additive build looks like, and it is the first one here that
said so by itself.

The graph was diffed too, and it should be: every edge layer is byte-for-byte the same number of
edges (`orthogroup` 1,741, `domain` 24,123, `coexpression` 63,158, `xlms` 79) because those are built
from declared columns rather than from every numeric one — and **`xyz` moved**, because two of the
ten new columns clear the >50% coverage floor `embed_features` applies to this arm and the layout
therefore has 73 features instead of 71. A new measurement re-laying the map is the design working,
not a regression, but it is the kind of change that should never be discovered later.

### One thing to carry into the next fetch

`sources.geo_sample_files` is new: many deposits keep their per-gene tables under each SAMPLE and
only a `_RAW.tar` at series level, and `geo_supplementary` reads the series directory. This one reads
the deposit's own `filelist.txt`, takes the files whose names end in a declared suffix, and derives
only the directory a named file lives in — a derivation that is wrong 404s immediately, which is the
failure mode to prefer over the one that silently skipped nineteen *Cryptosporidium* organisms.
Skipped names are logged and counted. 139 MB of coverage tracks stayed where they were.

### `Pf_secretome / excreted`, filled — and the coverage number it nearly published

Second fill of the pass, from the same queue: **PXD006925 / PMID 28944300**, extracellular vesicles
purified from a Kenyan clinical isolate. The deposit is raw-only — 24 RAW files and no RESULT — so
the per-protein numbers come from the paper's supplement, which is the usual shape here.

The sheet worth reading is the paper's own compilation: the union of two independent EV preparations
with a membership column each, so *how many studies saw this* is a fact in the file rather than a
join someone has to get right. 184 proteins, 53 of them in both. Its other columns are seroreactivity
and antibody-array results from unrelated studies — claims about immunity, not about vesicles — and
are deliberately not read; they belong to other slots. The sheet's trailing rows are its reference
list, and citations contain accessions, which is why resolution goes through the identity index
rather than a `PF3D7_` regex.

The column is named `ev_studies`, not `is_secreted`: vesicle proteomics is the instrument this
question has data for, and the second name would assert a route the measurement does not establish.

**Checked in both directions, and only one of them is a reason to trust it.** 20.1% of the 184 carry
a signal peptide against 10.2% of the proteome (p = 8e-05); exported proteins run 6.0% against 3.3%,
same direction and not significant at this size; RESA, KAHRP, MSP1 and Ag332 are all present. It is
**also** six times more expressed than the rest of the proteome (median blood-stage expression 71.3
against 12.2, p = 8e-31). That is what mass spectrometry on a vesicle preparation returns, it is the
same confound as hyperLOPIT tracking abundance on the other arm, and a test now pins both directions
so the caveat cannot quietly leave the prose.

**The number this nearly published.** Shipped first with a companion boolean completed as False for
every other gene — the convention the other mass-spectrometry columns on this arm follow — the atlas
graded the slot **A at 100% coverage** for an experiment that identified 184 proteins. Those columns
pool proteome-wide assays, where "never observed" is an answer; this is one preparation from one
isolate, and the genes it did not report are the six-times-less-expressed ones, so a False there is
a detection limit written down as a negative result. The flag was dropped, the count alone ships, and
the slot reads **C at 3.2%**, beside the Toxoplasma secretome slot at C and 2.0%.

Worth someone's attention, and not changed here because it would alter shipped columns: the same
flag-and-count convention is why `Pf_phosphorylation`, `Pf_acetylation`, `Pf_lactylation` and
`Pf_palmitoylation` all read **100% coverage** while their counts cover 2,503, 1,145, 144 and 503
genes. For a pooled proteome-wide re-analysis that is defensible — the yes/no question really was
asked of every gene — but the atlas prints one number and it is the flag's.

**The refusal guard fired on its first real use**, and it was right to: dropping the boolean lost a
column, `scripts/build_plasmodium.py` refused to write, and the build went through only under
`--allow-loss "the False-filled EV flag was absence rendered as measurement; the count alone is
shipped"`. That sentence is now in the shell history and in this file, which is the difference the
override exists to make.

**Where the atlas stands after both fills: Toxoplasma 116 of 143, Plasmodium 45 of 128, combined 161
of 271** — gene 143/198, pair 15/19, metabolite 3/6, host_gene 0/48.

### What the next pass should take, in order

1. **`Pf_phosphorylation · kinase-substrate`** — PXD005207 (PfCDPK1) and PXD009465 (PfPK7) are both
   raw-only, so it is the papers' supplements again; the Toxoplasma arm answers the same question
   with `cdpk1_thiophospho_peptides`, so the shape is known.
2. **`Pf_interaction degree · IP-MS`** — PXD008219 (Kelch13/Eps15/clathrin), PXD006155 (EPIC),
   PXD008208 (MSRP6). Pair-indexed, which is the axis this arm is weakest on.
3. **`Pf_glycosylation · asexual blood stage`** — PXD033470, C-mannosylation of TSR-domain adhesins;
   a PARTIAL deposit, so again the supplement.
4. **`Pf_cell-cycle timing label`** — the one derived label the arm is missing, and the sources
   exist: a 48-point IDC series, or the single-cell atlas the v2 list proposes. The construction is
   already shared with the Toxoplasma arm (`cellcycle.stage_enrichment`), which is the reason to do
   it rather than invent a second one.

All four are in `41_candidates_v2.json` with accessions. None of the four is in the seventeen-slot
queue this session opened with.
