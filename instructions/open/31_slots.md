# Slots: one place per question, and a choice about what fills it

Proposed 2026-08-13. Two problems this solves, and they are different problems.

**Dominance by column count.** Eighteen raw RNA columns and three summary ones answer *one* question
between them — how much is this gene transcribed — while `n_xlink_partners` answers another on its
own. The current `BLOCKS` already normalise each block to equal variance before weighting, so this
is half solved: what is missing is that a block is defined by a column-name regex, which is a fact
about how a file was written rather than about what it measures.

**Choice between datasets that answer the same question.** There are three transcriptomes in the
cache, and nothing anywhere says they are alternatives. A slot makes that explicit: several
candidates, one decision, recorded.

## What a slot is

    slot        the question, in the map's own terms  ("how much is this gene translated")
    candidates  the datasets that answer it, each with its coverage and its assay
    policy      what to do when more than one is chosen
    provenance  which dataset supplied which gene, kept per gene, not per slot

`policy` is the part worth arguing about. Four, and each is right somewhere:

| policy | what it does | when it is right |
|---|---|---|
| `one` | use the chosen dataset, ignore the rest | assays that are not comparable |
| `average` | rank-normalise each, then mean | replicates of the same measurement |
| `fill` | take the best-covered, fill its gaps from the next, record the source | the orthoLOPIT pattern: measured first, transferred second |
| `separate` | keep both, share the slot's variance between them | when disagreement is the signal |

`average` rank-normalises first because this project already learned that lesson: the seven CRISPR
screens differ 64-fold in spread, and averaging them raw is an average of one screen. `fill` writes
a `<slot>_source` column, exactly as `compartment_source` does today, so an inference never reads as
a measurement.

**The guard gets better, not worse.** `excluded_for` currently reasons about columns and datasets;
with slots it can exclude by slot, which is what shared provenance and shared quantity were
approximating on 2026-08-13.

---

## The taxonomy

Coverage measured on the shipped cache (8,140 genes). Grades: **A** genome-wide and measured,
**B** partial or targeted, **C** indirect or derived, **—** nothing.

### Baseline abundance

| slot | filled by | coverage | grade | what would improve it |
|---|---|---|---|---|
| transcription, baseline | `expr_*`, `rna108740_*` | 8,131 (99.9%) | **A** | — |
| translation, baseline | none | 0 | **—** | ribosome profiling: [29228904](https://pubmed.ncbi.nlm.nih.gov/29228904/) `GSE99395`, [38782906](https://pubmed.ncbi.nlm.nih.gov/38782906/) `GSE245775` |
| protein abundance, baseline | `proteome_*` (undocumented), `protein_ibaq_log2` | 3,005 (36.9%) | **B** | [40716488](https://pubmed.ncbi.nlm.nih.gov/40716488/) `PXD063409` cytoskeleton; [25867681](https://pubmed.ncbi.nlm.nih.gov/25867681/) `PXD000297` proteogenomics |
| phosphorylation, baseline | `n_phosphosites`, `phospho_*` | 1,603–8,140 | **B** | site-level: [33363051](https://pubmed.ncbi.nlm.nih.gov/33363051/) `PXD020655`, [31508380](https://pubmed.ncbi.nlm.nih.gov/31508380/) `PXD007777` |
| other acyl marks (acetyl, lactyl, nitrosyl) | none | 0 | **—** | [37562054](https://pubmed.ncbi.nlm.nih.gov/37562054/) `PXD040368`, [36216028](https://pubmed.ncbi.nlm.nih.gov/36216028/) `PXD022700`, [37959749](https://pubmed.ncbi.nlm.nih.gov/37959749/) `PXD046083` |
| ubiquitination / SUMOylation | none | 0 | **—** | [40348811](https://pubmed.ncbi.nlm.nih.gov/40348811/) `PXD045018`, [40590555](https://pubmed.ncbi.nlm.nih.gov/40590555/) `PXD054719` |
| glycosylation | none | 0 | **—** | [39912628](https://pubmed.ncbi.nlm.nih.gov/39912628/) `PXD056853` (O-fucose) |
| palmitoylation | none | 0 | **—** | nothing deposited found in 24,679 abstracts — a real gap in the field |
| protein turnover / stability | none | 0 | **—** | thermal profiling [35976251](https://pubmed.ncbi.nlm.nih.gov/35976251/) `PXD033642` |
| RNA stability / half-life | none | 0 | **—** | [39899594](https://pubmed.ncbi.nlm.nih.gov/39899594/) `PRJEB67890` (iron, post-transcriptional) |

### Time and stage

| slot | filled by | coverage | grade | what would improve it |
|---|---|---|---|---|
| transcription per cell-cycle stage | none — only the phase LABEL (`cellcycle_phase`) | 873 labelled | **C** | the source of that label carries per-phase expression: `xue_singlecell` re-ingested as values, not a call |
| translation per cell-cycle stage | none | 0 | **—** | no Toxoplasma cell-cycle Ribo-seq found |
| transcription per life-cycle stage | `expr_tachy/cyst/sporulated`, `rna206344_*` | 7,974 (98.0%) | **A** | bradyzoite subtypes [41580398](https://pubmed.ncbi.nlm.nih.gov/41580398/) `GSE311669`; merozoite [25757795](https://pubmed.ncbi.nlm.nih.gov/25757795/) `PRJEB7935`; sexual stages [38093015](https://pubmed.ncbi.nlm.nih.gov/38093015/) `GSE222819` |
| translation per life-cycle stage | none | 0 | **—** | [29228904](https://pubmed.ncbi.nlm.nih.gov/29228904/) `GSE99395` covers tachyzoite→bradyzoite |
| stress response | `stress_*` (undocumented) | 7,880 (96.8%) | **A** | already in the cache; no block can see it |
| single-cell heterogeneity | none | 0 | **—** | [40630957](https://pubmed.ncbi.nlm.nih.gov/40630957/) `GSE295224`; [41278634](https://pubmed.ncbi.nlm.nih.gov/41278634/) brain bradyzoites |

### Fitness

| slot | filled by | coverage | grade | what would improve it |
|---|---|---|---|---|
| in vitro fitness, HFF | `fit_invitro_hff` | 7,325 (90.0%) | **A** | — |
| in vitro fitness, macrophage (naive) | `fit_naive_bmdm` | 7,402 | **A** | — |
| in vitro fitness, macrophage (IFN-γ) | `fit_ifng` | 7,402 | **A** | [36916910](https://pubmed.ncbi.nlm.nih.gov/36916910/) as a second measurement |
| in vivo fitness, lung / spleen / liver / peritoneum | `fit_invivo_*` | 7,395 (90.8%) | **A** | — |
| in vivo fitness, brain / chronic | `invivo_*` (undocumented) | 7,663 (94.1%) | **A** | in the cache; no block can see it |
| fitness under oxidative stress | none | 0 | **—** | [34163449](https://pubmed.ncbi.nlm.nih.gov/34163449/) `PRJNA707360` |
| fitness of hyperLOPIT-unassigned proteins | none | 0 | **—** | [39082802](https://pubmed.ncbi.nlm.nih.gov/39082802/) `GSE253884` — aimed at this map's largest gap |
| genetic interaction (synthetic lethality) | `crispr_gra17_*` | 7,553 | **B** | more Δ backgrounds; only GRA17 exists |
| stage-conversion phenotype | none | 0 | **—** | thin: [40941387](https://pubmed.ncbi.nlm.nih.gov/40941387/) `PRJNA1305388` |
| drug / chemogenomic sensitivity | none | 0 | **—** | 11 papers, none deposited |

### Regulation

| slot | filled by | coverage | grade | what would improve it |
|---|---|---|---|---|
| transcription-factor regulation (one column per TF) | `morc_*` only (undocumented) | 7,841 (96.3%) | **C** | ApiAP2 ChIP: [29788176](https://pubmed.ncbi.nlm.nih.gov/29788176/) `GSE106864`; AP2X-7 [40920096](https://pubmed.ncbi.nlm.nih.gov/40920096/) `PRJNA1281009`; AP2 cascade [39774584](https://pubmed.ncbi.nlm.nih.gov/39774584/) `PRJNA1134678`; BDP1 [37350586](https://pubmed.ncbi.nlm.nih.gov/37350586/) `GSE228853`; SWI/SNF [41193469](https://pubmed.ncbi.nlm.nih.gov/41193469/) `GSE287537`; SNF2L [40593611](https://pubmed.ncbi.nlm.nih.gov/40593611/) `GSE268652`; GCN5a [40879678](https://pubmed.ncbi.nlm.nih.gov/40879678/) `GSE286088` |
| chromatin state (histone marks, variants) | none | 0 | **—** | [35164683](https://pubmed.ncbi.nlm.nih.gov/35164683/) `GSE104347`, [21179246](https://pubmed.ncbi.nlm.nih.gov/21179246/) `GSE22100` |
| chromatin accessibility | none | 0 | **—** | [38093015](https://pubmed.ncbi.nlm.nih.gov/38093015/) `PRJNA921935` (pre-sexual stages) |
| RNA modification (m6A, 5mC) | none | 0 | **—** | [34263725](https://pubmed.ncbi.nlm.nih.gov/34263725/) `GSE168155`, [40830525](https://pubmed.ncbi.nlm.nih.gov/40830525/) `GSE294543` |
| splicing / isoforms | none | 0 | **—** | nanopore isoforms [33688018](https://pubmed.ncbi.nlm.nih.gov/33688018/) `PRJNA606986`; Cdc5 [40263328](https://pubmed.ncbi.nlm.nih.gov/40263328/) |

### Relations between genes (edges, not columns)

| slot | filled by | coverage | grade | what would improve it |
|---|---|---|---|---|
| co-transcription with other genes | `coexpression` edges (49,293) | genome-wide | **A** | recompute per stage rather than pooled |
| co-translation with other genes | none | 0 | **—** | needs the Ribo-seq above first |
| co-fitness across screens | `cofitness` edges (88,997) | genome-wide | **A** | — |
| interaction, crosslink MS | `xlms` edges (2,842 pairs / 1,630 genes) | 20.0% | **B** | the only measured-contact layer; more XL-MS is scarce |
| interaction, IP-MS (parasite–parasite) | `ip_ms` edges (**64 pairs**) | 0.6% | **C** | **55 studies already downloaded and not curated** |
| interaction, proximity labelling (BioID) | none as edges | 0 | **—** | **44 studies already downloaded**; plus [36214684](https://pubmed.ncbi.nlm.nih.gov/36214684/) `PXD032102` |
| structural similarity | `struct` edges (11,684 pairs) | 28.7% | **A** | — |
| interaction with host proteins (per host protein) | none | 0 | **—** | MYR1-dependent set [32075880](https://pubmed.ncbi.nlm.nih.gov/32075880/) `PXD016383`; exportome TurboID [38747635](https://pubmed.ncbi.nlm.nih.gov/38747635/) |

### Where the protein is, and what it is exposed to

| slot | filled by | coverage | grade | what would improve it |
|---|---|---|---|---|
| localization, measured | `compartment` (hyperLOPIT) | 3,827 (47.0%) | **B** | nothing better exists for Toxoplasma |
| localization, transferred | `ortholopit_*`, `compartment_best` | 8,140 with source | **A** for the pattern | this is the `fill` policy, already implemented once |
| membrane topology | `dtm_class` (GLOB 6,053 / TM 1,142 / SP 676 / SP+TM 266), `n_tm`, `has_signal_peptide` | 8,140 | **A** | — |
| exposure to the host cytosol | none | 0 | **—** | derivable now from topology + PVM proteomes; measured by exportome TurboID [38747635](https://pubmed.ncbi.nlm.nih.gov/38747635/), MYR1-dependent effectors [29615509](https://pubmed.ncbi.nlm.nih.gov/29615509/) `GSE109830`, PVM translocation [31366709](https://pubmed.ncbi.nlm.nih.gov/31366709/) `GSE122786` |
| secretome / exported | none as a slot | 0 | **—** | [41137792](https://pubmed.ncbi.nlm.nih.gov/41137792/) `PXD028969` |
| host transcriptional effect per effector | `hosttx_T2`, `hosttx_padj` | 252 (3.1%) | **C** | the dual perturb-seq matrix itself: `GSE229505` |

### Sequence, structure and conservation

| slot | filled by | coverage | grade |
|---|---|---|---|
| domain content | `n_interpro`, `has_domain` | 8,140 | **A** |
| fold confidence / disorder | `mean_plddt` | 6,480 (79.6%) | **A** |
| conservation breadth | `has_pf_ortholog`, `has_cp_ortholog`, `lineage_specific` | 8,140 | **A** |
| paralogy | `paralog_number` | 8,140 | **A** |
| sequence basics (length, TM, pI) | `length`, `n_tm`, `tm_kd_*` | 8,140 | **A** |
| strain variation | none | 0 | **—** (nothing deposited found) |

### Not a measurement, and must never be a feature for inference

| slot | filled by | note |
|---|---|---|
| literature attention | `n_publications`, `n_papers_*`, `attention_depth`, `lit_tier` | The negative control. It belongs in the taxonomy so that it can be **excluded on purpose** rather than forgotten — which is exactly how it leaked into the control's own embedding until 2026-08-13. |

---

## What this says about priorities

1. **Four slots are filled by data already in the cache that no block can reach**: stress response,
   in vivo brain/chronic fitness, TF regulation (MORC), and the 3,005-gene proteome. No downloads.
2. **Translation is the largest genuinely empty axis** — baseline, per stage, and co-translation all
   hang off two Ribo-seq datasets (`GSE99395`, `GSE245775`).
3. **The interaction slots are the worst-served relative to what is already downloaded**: 64 IP-MS
   pairs in the map against 55 studies of IP-MS and 44 of BioID sitting on disk.
4. **"One column per transcription factor" is the slot that most needs the mechanism**, because it
   is the one that would otherwise dominate: seven ChIP datasets over a few dozen factors is a
   hundred-odd columns answering one question, which is precisely what per-slot normalisation is for.
