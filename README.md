# starplast

A 3D browser for the *Toxoplasma gondii* knowledge map. 8,140 genes are points in space, twelve kinds of
relation are toggleable edges, and clicking a gene shows everything that is actually known about it —
including when the answer is "nothing".

The single idea it is built around: **most of this proteome has never been studied, and a tool that does
not say so will mislead you.** 5,574 of 8,140 genes are named in no paper at all, and of those that are
named, most appear only in a screen's hit table rather than in any title or abstract. So the interface
distinguishes measurement from inference from absence everywhere, and never lets one read as another.

Why the map is a UMAP rather than a force-directed layout, how attention correction works, and every
other design decision: **[HANDOFF.md](HANDOFF.md)**. Methods prose for publication:
**[MATERIALS_AND_METHODS.md](MATERIALS_AND_METHODS.md)**.

## Install and run

```bash
pip install starplast            # the program AND the GPU stack, where this platform has wheels
pip install starplast-cpu        # the program alone, no CUDA
pip install -e .                 # this checkout, editable  (add ".[gpu]" for the CUDA stack)
starplast-install-gpu            # add GPU support later; picks the CUDA set from the driver
starplast
```

Needs a display and OpenGL. Nothing else: the built cache ships inside the package, so the application
runs offline, with no dataset, no configuration and no build step.

Rebuilding the cache from the published sources is a separate job and needs the raw data:

```bash
python -m starplast.paths            # where data is being resolved from, and what is missing
python -m starplast.fetch_names      # one-off: ToxoDB symbols, previous and strain accessions
python -m starplast.build_graph      # ~5 min; writes the cache
pytest tests/ -q                     # headless, no network
```

If a file cannot be found, `python -m starplast.paths` prints which locations were searched and which
datasets are absent. Two environment variables override the search:

| variable | what it points at |
|---|---|
| `STARPLAST_CACHE` | the built cache the application reads |
| `STARPLAST_DATA` | the raw published datasets, needed only to rebuild |

The one deliberate exception to running offline is **structure coordinates**. AlphaFold models and
crosslink complexes are gigabytes, so they are fetched when you click a gene and cached under
`~/.cache/starplast`. Everything needed to *decide* something is already on disk; only the picture is
downloaded.

## What you see

| dark | light |
|---|---|
| ![the map, dark theme](docs/screenshots/map_dark.png) | ![the map, light theme](docs/screenshots/map_light.png) |

A gene is selected in both: the halo marks it, and its edges fan out to the genes it is related to.
Grey is *unknown*, never a category and never zero — 5,574 genes are named in no paper at all.

The level of detail follows the data rather than invented tiers. Coarse spatial structures found in
the embedding, then orthogroup centroids, then genes:

| galaxy | orthogroup | gene |
|---|---|---|
| ![galaxy tier](docs/screenshots/lod_0_compartment.png) | ![orthogroup tier](docs/screenshots/lod_1_orthogroup.png) | ![gene tier](docs/screenshots/lod_2_gene.png) |

The coarsest tier is computed from the embedding itself, not from compartment. Compartment centroids
were tried first and collapsed into a blob in the middle of the screen, because a compartment's genes
are spread across the whole map and the mean of scattered points is the middle — HANDOFF decision 4.

### Controls

Settings live in the menu bar; the panel keeps only what is used continuously.

| | |
|---|---|
| left-drag / scroll | rotate / zoom (Navigate mode), or draw a gate (Select mode) |
| left-click a gene | select it; fills the evidence panel |
| **View ▸ Left mouse button** | Navigate or Select; rotation free or about one axis; gate as a 2D lasso or a 3D brush |
| right-click the map | spin, level of detail, point size, export, reset view |
| search box | gene ID or product text, then flies to it |
| filter by | any categorical column — compartment, cell-cycle phase, attention depth, and others |
| double-click a value | fly to that class's centroid |
| **View** | level of detail, color by, point size, spin, theme, preferences |
| **Preferences ▸ lighting** | soft or GPU volumetric ray-traced shadows; mouse/selection target; neutral, cool, or warm light; flat discs or OpenGL GGX glossy/metallic spheres |
| **Edges** | the thirteen relation types, "draw all active edges", attention correction, and why they are never combined |
| **File** | export the image, the visible genes as CSV, or the active graph as GraphML |
| **Tools** | console, running jobs, the walk gallery, and an assistant that is told what is on screen |

Edges draw only for the selected gene unless *Edges ▸ Draw all active edges* is on, which draws the
strongest 20,000 per type.

Every results table in the analysis panel behaves the same way: **click a row to see the map it is
about** — a walk or search row is rebuilt exactly, on the same genes, with the clustering its score
counted; a clustering row is applied to the map on screen; an Inference or Validation row colors the
map by the clustering it was scored against. **Right-click a table to save it as CSV** (the whole
result, not the screenful shown) or to copy the selected rows.

In **Select** mode a drag gates a set of genes: a 2D lasso takes everything behind it, which is what
"grab that visual cluster" means, while a 3D brush takes a ball in world space around the gene you
press on, which is what you want when the cloud is deep and a lasso would also catch the far side.
The gated set's composition appears in the evidence panel — led by how many of its genes carry no
label, since those are the candidates a gate is drawn to find — and *File ▸ Export gated selection*
writes it out. Constraining rotation to one axis makes a view reproducible; a free orbit never
returns to the same one twice.

A hyperparameter walk fills the **gallery** along the bottom as it runs: one thumbnail per
configuration, appearing the moment it is computed, either as a grid to compare or as a single view
stepped through with a slider. Clicking one opens that embedding in the main view, where it behaves
like any other map. Genes outside the walk's subsample have no position in it, so they are hidden
rather than drawn at the origin.

## Datasets included

Every dataset the map is built from, with what it measures and where it came from. This table is
generated from `starplast.datasets.REGISTRY`; a test fails if it drifts from the registry.

### DNA — genetic perturbation and DNA-level readouts

| Dataset | Type of data | Coverage | Reference |
|---|---|---|---|
| Arrayed splitCas9 imaging screen | What a parasite looks like when a gene is off: egress, actin, apicoplast, replication | 319 genes screened, 99 with a phenotype, 35 at egress | PMID [35538310](https://pubmed.ncbi.nlm.nih.gov/35538310/); `Nat Microbiol 41564-2022-01114 Supplementary Tables 2 and 3B` |
| Differentiation reporter CRISPR screen (COMPUTED) | Guide enrichment in reporter-positive parasites against the bulk population | 235 genes | `GSE132237` |
| Fitness in the reporter strain (COMPUTED) | Guide depletion over eight passages of ordinary growth | 262 genes | `GSE132237` |
| GRA12 strains and mouse subspecies | Median L2FC in vitro and in vivo, DISCO score; two screens | 236 / 232 | GRA12 is a common virulence factor across Toxoplasma gondii strains and mouse subspecies; PMID [40240328](https://pubmed.ncbi.nlm.nih.gov/40240328/) |
| GRA17 synthetic-lethal screen | RH and RH-delta-gra17 phenotype by passage; MAGeCK p-values | 7,553 (genome-wide) | Genome-wide CRISPR screen identifies genes synthetically lethal with GRA17, a nutrient channel encoding gene in Toxoplasma; PMID [37498952](https://pubmed.ncbi.nlm.nih.gov/37498952/) |
| HDAC3 occupancy (CUT&TAG, COMPUTED) | Mean HDAC3 CUT&TAG coverage over the promoter, relative to the genome mean | 8,140 genes (100%) | `GSE277553` |
| Histone H4 acetylation (ChIP-chip, via ToxoDB) | Genome-wide H4 K5/K8/K12/K16 acetylation score within 1 kb of the gene | 7,515 genes (92%) | `ToxoDB Hakimi/Ali H4 acetylation` |
| Host ESCRT recruitment screen (UNPUBLISHED) | Per-gene effect on host TSG101 recruitment to the vacuole, by two models | 13 and 8 genes | Olafsson EB et al., A pooled image-based CRISPR screen identifies EAF1 as a T. gondii modulator of ESCRT subversion. bioRxiv 2026 (under submission); `spaCR screen, bioRxiv 10.64898/2026.07.08.737057` |
| Host-transcription effector screen | Hotelling T2 plus full per-effector host-response signature | 252 screened / 22 full signatures | High-throughput identification of Toxoplasma gondii effector proteins that target host cell transcription; PMID [37827122](https://pubmed.ncbi.nlm.nih.gov/37827122/) |
| In vitro CRISPR fitness (HFF) | Competitive growth in fibroblasts | 7,325 (90.0%) | A Genome-wide CRISPR Screen in Toxoplasma Identifies Essential Apicomplexan Genes (Sidik et al. 2016); PMID [27594426](https://pubmed.ncbi.nlm.nih.gov/27594426/) |
| In vivo CRISPR composite scores | Peritoneum, lung, liver, spleen composite scores | 7,395 (90.8%) | PMID [31481656](https://pubmed.ncbi.nlm.nih.gov/31481656/); `ToxoDB tgonGt1CrisprFunc*` |
| In vivo CRISPR platform | Mean log fold-change across replicates | 168 | A CRISPR platform for targeted in vivo screens identifies Toxoplasma gondii virulence factors in mice; PMID [31481656](https://pubmed.ncbi.nlm.nih.gov/31481656/) |
| In-vivo fitness of hyperLOPIT-unassigned proteins | Two targeted libraries tested during mouse infection | measured at build time | Tachibana Y et al., CRISPR screens identify genes essential for in vivo virulence among proteins of hyperLOPIT-unassigned localization. mBio 2024; PMID [39082802](https://pubmed.ncbi.nlm.nih.gov/39082802/); `GSE253884;GSE253885` |
| Macrophage CRISPR screens | Naive BMDM and IFN-gamma survival | 7,402 (90.9%) | Wang Y et al., Genome-wide screens identify Toxoplasma gondii determinants of parasite fitness in IFN-gamma-activated murine macrophages. Nat Commun 2020;11:5258; PMID [33067458](https://pubmed.ncbi.nlm.nih.gov/33067458/) |
| Oxidative-stress CRISPR screen | Screening score per gene under oxidative challenge | 7,384 genes (91%) | PMID [34163449](https://pubmed.ncbi.nlm.nih.gov/34163449/); `PMC8216390 Data Sheet 1` |
| Promoter accessibility (ATAC-seq, COMPUTED) | Mean ATAC coverage over the promoter, relative to the genome mean | 7,988 genes (98%) | `GSE313048` |
| Young 2019 in vivo screen | In vivo fitness | 115 | Young J et al., A CRISPR platform for targeted in vivo screens identifies Toxoplasma gondii virulence factors in mice. Nat Commun 2019;10:3963; PMID [31481656](https://pubmed.ncbi.nlm.nih.gov/31481656/) |

### Transcription — RNA abundance

| Dataset | Type of data | Coverage | Reference |
|---|---|---|---|
| Alkaline-stress differentiation transcriptome | Unstressed tachyzoites and alkaline-stressed bradyzoites | 7,880 (96.8%) | Waldman BS et al., Identification of a Master Regulator of Differentiation in Toxoplasma. Cell 2020;180:359-372.e16; PMID [31955846](https://pubmed.ncbi.nlm.nih.gov/31955846/); `GSE132248` |
| Antisense transcription (via ToxoDB) | Percentile of antisense signal at this gene, across the life cycle | 8,140 genes (100%) | `ToxoDB full life-cycle transcriptome, Antisense` |
| BFD2-bound transcriptome (RIP-seq, COMPUTED) | Enrichment of each transcript in the BFD2 immunoprecipitation | 7,463 genes (92%) | `GSE223620` |
| Bradyzoite restriction-checkpoint transcriptome | Cyclin perturbations in tachyzoite and bradyzoite conditions | measured at build time | `GSE200962` |
| CPSF4 RNA-processing perturbation transcriptome | RNA response at 7, 24 and 48 hours after CPSF4 depletion | measured at build time | Farhat DC et al., A plant-like mechanism coupling m6A reading to polyadenylation safeguards transcriptome integrity. eLife 2021;10:e68312; PMID [34263725](https://pubmed.ncbi.nlm.nih.gov/34263725/); `GSE168155` |
| Enteroepithelial stage transcriptome (via ToxoDB) | Expression in the feline enteroepithelial stages against tachyzoites | 7,739 genes (95%) | `ToxoDB Ramakrishnan enteroepithelial` |
| Expression in infected macrophages (via ToxoDB) | Expression percentile in ME49-infected murine macrophages | 8,140 genes (100%) | `ToxoDB Saeij 29 strains` |
| Feline merozoite transcriptome | Merozoite expression with matched tachyzoite comparators | measured at build time | Behnke MS et al., Toxoplasma gondii merozoite gene expression analysis with comparison to the life cycle. BMC Genomics 2014;15:350; PMID [24885521](https://pubmed.ncbi.nlm.nih.gov/24885521/); `GSE51780` |
| In vivo brain-stage transcriptome | Tachyzoites, acute/chronic whole brain, and purified bradyzoites | 7,663 (94.1%) | Garfoot AL et al., Proteomic and transcriptomic analyses of early and late-chronic Toxoplasma gondii infection shows novel and stage specific transcripts. BMC Genomics 2019;20:859; PMID [31726967](https://pubmed.ncbi.nlm.nih.gov/31726967/) |
| Life-cycle stage enrichment (DERIVED) | Which stage a gene's own expression is highest in | 1,911 of 8,140 genes called | *citation not yet confirmed* |
| MORC depletion and BFD1 perturbation transcriptome | MORC knockdown, BFD1 knockout and BFD1 stabilization series | 7,841 (96.3%) | `PXD058095` |
| Novel transcript models (Nanopore, via ToxoDB) | How many novel TALON transcript models long reads support for this gene | 798 genes (10%) | `ToxoDB Stuart/Ralph nanopore` |
| Oocyst sporulation series | Unsporulated / sporulating / sporulated, 2 replicates (6 columns) | 7,974 (98.0%) | `GSE206344` |
| Plasmodium falciparum life-stage and polysomal RNA | Transcript abundance across seven life stages, and what is on ribosomes | 5,720 P. falciparum genes | `PlasmoDB: Su seven stages, Bunnik polysomal IDC, Gomez-Diaz mosquito stages` |
| Plasmodium long-read transcript models | Transcript models per gene, and how many the annotation does not contain | 1,857 genes, 2,498 models, 238 novel | PMID [40316999](https://pubmed.ncbi.nlm.nih.gov/40316999/); `Malar J 05376 Supplementary Data 2` |
| Plasmodium peak expression and stage label (DERIVED) | Maximum expression across stages, and which stage a gene belongs to | 5,720 genes for the maximum, 310 labelled | *citation not yet confirmed* |
| Plasmodium transcription at febrile temperature | Wild type and two mutants at 37 C and at the 41 C of a malarial fever | 5,791 genes | `PlasmoDB Pfal3D7 Febrile temps RNA-Seq` |
| Plasmodium transcription under Sir2 knockout | Wild type and sir2a / sir2b knockout at ring, trophozoite and schizont | 5,615 genes | `PlasmoDB Sir2 KO Marray` |
| Primary brain-cell parasite differentiation time course | Parasite base mean and log2 fold-change at days 1, 2, 4, 7 and 14 | measured at build time | Mouveaux T et al., Primary brain cell infection by Toxoplasma gondii reveals spontaneous bradyzoite differentiation and modification of neuron biology; PMID [34610266](https://pubmed.ncbi.nlm.nih.gov/34610266/); `GSE168465` |
| Pru tachyzoite / 72-hour bradyzoite stage array | Matched tachyzoite and alkaline-induced bradyzoite expression | 7,253 genes | `GSE22258` |
| Sexual development in the cat (single-cell atlas) | Enrichment at 8 days post-infection, when gametogony happens | 4,463 genes (55%) | PMID [41929010](https://pubmed.ncbi.nlm.nih.gov/41929010/); `PMC13042011 supplementary media-2` |
| Single-parasite transcriptional atlas (cell cycle) | Measured cell-cycle phase per gene, and pseudotime cluster | 873 genes phased, 7,499 clustered | Xue Y et al. eLife 2020;9:e54129; PMID [32065584](https://pubmed.ncbi.nlm.nih.gov/32065584/) |
| Stage transcriptome | Tachyzoite, day 3/5/7, in vivo tissue cyst (12 columns) | 7,739 (95.1%) | `GSE108740` |
| Synchronized tachyzoite cell-cycle transcriptome | Two replicates across blocked, asynchronous and hourly release states | measured at build time | Behnke MS et al., Coordinated progression through two subtranscriptomes underlies the tachyzoite cycle of Toxoplasma gondii. PLoS ONE 2010;5:e12354; PMID [20865045](https://pubmed.ncbi.nlm.nih.gov/20865045/); `GSE19092` |
| m6A methylome (MeRIP peaks) | How many m6A peaks the authors called on this gene in tachyzoites | 837 genes (10%) | PMID [34324585](https://pubmed.ncbi.nlm.nih.gov/34324585/); `PLoS Pathogens 1009335 Table S3A` |
| mRNA stability after actinomycin D | Proportion of transcript remaining after five hours of transcription block | 412 genes | PMID [39899594](https://pubmed.ncbi.nlm.nih.gov/39899594/); `PLoS Pathogens 1012857 Table S12` |

### Translation — protein abundance

| Dataset | Type of data | Coverage | Reference |
|---|---|---|---|
| AP2XII-1/AP2XI-2 perturbation total proteome | Replicate abundance and log2 fold-change during pre-sexual conversion | 3,005 (36.9%) | Antunes AV et al., In vitro production of cat-restricted Toxoplasma pre-sexual stages. Nature 2024;625:366-376; PMID [38093015](https://pubmed.ncbi.nlm.nih.gov/38093015/); `PXD039400, PXD042658` |
| Co-translation layer (COMPUTED) | Gene pairs whose ribosome footprints covary | 6,231 edges over 7,437 genes | *citation not yet confirmed* |
| Differentiation ribosome profiling (eIF1.2) | RPF and RNA counts, and their ratio, in tachyzoites and pre-bradyzoites | 7,880 genes (97%) | PMID [38782906](https://pubmed.ncbi.nlm.nih.gov/38782906/); `GSE245775` |
| Host-context parasite ribosome profiling | Parasite ribosome footprints, RNA and translation efficiency in two HFF states | measured at build time | Holmes MJ et al., Simultaneous Ribosome Profiling of Human Host Cells Infected with Toxoplasma gondii. mSphere 2019;4:e00292-19; PMID [31167946](https://pubmed.ncbi.nlm.nih.gov/31167946/); `GSE129869` |
| Intracellular/extracellular ribosome profiling | Ribosome footprints, matched RNA and relative translation efficiency | measured at build time | Hassan MA et al., Comparative ribosome profiling uncovers a dominant role for translational control in Toxoplasma gondii. BMC Genomics 2017;18:961; PMID [29228904](https://pubmed.ncbi.nlm.nih.gov/29228904/); `GSE99395` |
| Oocyst developmental-stage iTRAQ proteome | iTRAQ abundance ratios across oocyst developmental stages | 2,079 (25.5%) | Possenti A et al., Proteomic Differences between Developmental Stages of Toxoplasma gondii Revealed by iTRAQ-Based Quantitative Proteomics. Front Microbiol 2017;8:1732; PMID [28626452](https://pubmed.ncbi.nlm.nih.gov/28626452/); `PXD003765` |
| Pru proteome and IP abundance | Median log2 iBAQ across replicates | 748 (9.2%) | `PXD043808, PXD065585` |

### Post-translation — properties of the folded protein

| Dataset | Type of data | Coverage | Reference |
|---|---|---|---|
| BioID/TurboID supplement membership corpus | Number of downloaded proximity-labeling studies whose supplement names each gene | measured at build time | *citation not yet confirmed* |
| C. parvum hyperLOPIT | Donor labels for orthoLOPIT transfer | 1,107 usable | Guerin et al. 2023 |
| CDPK1 substrates (thiophosphate labelling) | Thiophosphorylated peptides per gene from analog-sensitive CDPK1 | 361 genes | PMID [37933960](https://pubmed.ncbi.nlm.nih.gov/37933960/); `eLife 85654 supplementary file 6` |
| Calcium thermal-shift proteome (mineCETSA) | How far a protein's melting curve moves when calcium is added | 2,348 proteins | PMID [35976251](https://pubmed.ncbi.nlm.nih.gov/35976251/); `PMC9436416 Supplementary file 3` |
| Crosslinking MS interactome | How many proteins this one crosslinks to | 494 proteins | PMID [40874616](https://pubmed.ncbi.nlm.nih.gov/40874616/); `mBio 02159-25 supplementary file s0004` |
| Cyst wall interactome | Strongest bait signal and how many baits saw the protein | 56 proteins | PMID [32019789](https://pubmed.ncbi.nlm.nih.gov/32019789/); `PMC7002340 Data Set S1` |
| Foldseek structural similarity | TM-align over Toxoplasma AlphaFold models, TM >= 0.7 | 11,684 pairs / 2,338 genes | *citation not yet confirmed* |
| Host proteins at the vacuole | How enriched a host protein is at the parasitophorous vacuole | 12 host proteins | PMID [34898650](https://pubmed.ncbi.nlm.nih.gov/34898650/); `PLoS Pathogens 1010138 supplementary table` |
| IP-MS of tagged baits | Replicated pulldown vs untagged control | 64 pairs / 48 genes | `PXD043808, PXD065585` |
| IP-MS supplement membership corpus | Number of downloaded pulldown studies whose supplement names each gene | measured at build time | *citation not yet confirmed* |
| Lysine acetylome (GCN5b) | Acetylation sites reported per gene | 3,921 genes measured | `PXD079431` |
| Lysine lactylome | Lactylation sites reported per gene | 515 genes measured | `PXD031526` |
| MYR1 host interactome (bridge) | Host proteins co-immunoprecipitating with the parasite protein MYR1 | 219 host proteins, 1 parasite gene | PMID [32075880](https://pubmed.ncbi.nlm.nih.gov/32075880/); `PXD016383` |
| Monomethylarginine proteome (via ToxoDB) | Monomethylarginine sites reported per gene | 368 genes | `ToxoDB Yakubu monomethylarginine` |
| N-myristoylated proteome | The authors' confidence that this protein is myristoylated, 3 high to 1 low | 65 substrates | PMID [32618271](https://pubmed.ncbi.nlm.nih.gov/32618271/); `eLife 57861 supplementary file 4` |
| O-fucosylated glycoproteins (AAL pulldown) | Peptide identifications in the AAL lectin pulldown, per gene | 394 genes | `PXD004426` |
| Oocyst-versus-tachyzoite phosphoproteome | Measured-site counts and strongest up/down phosphosite ratios | 1,603 (19.7%) | Wang Z-X et al., Comparative Phosphoproteomic Analysis of Sporulated Oocysts and Tachyzoites of Toxoplasma gondii Reveals Stage-Specific Patterns. Molecules 2022;27:1109; PMID [35164288](https://pubmed.ncbi.nlm.nih.gov/35164288/); `PXD017032` |
| P. falciparum LOPIT | Donor labels for orthoLOPIT transfer | 1,646 usable | Chisholm SA et al., The spatial proteome of the Plasmodium falciparum schizont. Nat Commun 2026;17:6192 -- CONFIRM against the file on disk; PMID [42218142](https://pubmed.ncbi.nlm.nih.gov/42218142/) |
| PVM proximity labelling | Whether the study placed this protein at the parasitophorous vacuole membrane | 1,274 genes (73 positive) | PMID [34749525](https://pubmed.ncbi.nlm.nih.gov/34749525/); `mBio 00260-21 Data Set S1` |
| Phosphosite counts | Count of phosphosites per protein, no positions | 1,175 (14.4%) | Treeck M et al. 2011 -- CONFIRM against the file on disk |
| Plasmodium N-myristoylome (NMT-inhibitor sensitive) | Proteins whose click-chemistry capture drops when N-myristoyltransferase is blocked | 16 substrates of 609 assayed | PMID [34695132](https://pubmed.ncbi.nlm.nih.gov/34695132/); `PLoS Biol 3001408 S11` |
| Plasmodium lysine acetylome | Acetylated lysines per gene, and whether the gene was seen acetylated at all | 1,145 genes, 2,163 localised sites | PMID [26813983](https://pubmed.ncbi.nlm.nih.gov/26813983/); `Sci Rep 19722 S2` |
| Plasmodium lysine lactylome (resolved from NF54) | Lactylated lysines per gene, reported against NF54 and resolved to 3D7 | 144 genes | PMID [41417877](https://pubmed.ncbi.nlm.nih.gov/41417877/); `PLoS Genet 1011991 S1` |
| Plasmodium palmitome (observed only) | Proteins observed S-palmitoylated, with the motif prediction deliberately excluded | 503 proteins | PMID [36250062](https://pubmed.ncbi.nlm.nih.gov/36250062/); `Front Cell Infect Microbiol Table 3` |
| Plasmodium phosphosites (re-analysis of all public data) | Distinct phosphorylated residues per gene, pooled across every public study | 16,318 sites over 2,503 genes | `PXD046874` |
| Protein melting temperature (mineCETSA) | Where this protein's melting curve sits, in degrees | 3,120 proteins (38%) | PMID [35976251](https://pubmed.ncbi.nlm.nih.gov/35976251/); `eLife 80336 supplementary file 3` |
| Proximity labelling | Proximity partners reported per gene | 1,734 genes measured | `PXD059579` |
| Proximity-labeling corpus | 42 BioID/TurboID/APEX studies with a tagged Toxoplasma protein | 28 studies with data, 127 files | *citation not yet confirmed* |
| Pulldown corpus | 55 IP-MS / co-IP studies with a tagged Toxoplasma protein | 29 studies with data, 140 files | *citation not yet confirmed* |
| S-nitrosylation (iodoTMT) | S-nitrosylation sites reported per gene | 660 genes measured | `PXD046083` |
| S-palmitoylome (Foe 2015, via ToxoDB) | 17-ODYA enrichment per gene, against hydroxylamine and against palmitate | 470 and 488 genes | PMID [26468752](https://pubmed.ncbi.nlm.nih.gov/26468752/); `ToxoDB Foe palmitome` |
| Secreted-fraction partition | How a secreted protein splits between the soluble and vesicular fractions | 165 proteins | PMID [40874616](https://pubmed.ncbi.nlm.nih.gov/40874616/); `ToxoDB Ramirez-Flores vesicles` |
| StarPath crosslink MS | Measured physical proximity; residue-level crosslinks and Chai-1 complexes | 2,842 pairs / 1,630 genes | Mapping a Toxoplasma gondii interactome by crosslinking mass spectrometry and machine learning (2025); PMID [40874616](https://pubmed.ncbi.nlm.nih.gov/40874616/) |
| T. gondii hyperLOPIT | Subcellular compartment, MAP and MCMC, with posteriors | 3,827 (47.0%) | A Comprehensive Subcellular Atlas of the Toxoplasma Proteome via hyperLOPIT (Barylyuk et al. 2020); PMID [33053376](https://pubmed.ncbi.nlm.nih.gov/33053376/) |
| Ubiquitination / SUMOylation (GlyGly) | GlyGly sites reported per gene | 128 genes measured | `PXD042937` |

### Reference — not a study result

| Dataset | Type of data | Coverage | Reference |
|---|---|---|---|
| AlphaFold DB | Per-gene mean pLDDT; coordinates fetched on demand | 6,480 (79.6%) | Varadi et al. 2024 NAR (database); Jumper et al. 2021 Nature (method); `UP000001529 (taxid 508771), AlphaFold DB` |
| Antibody epitopes (IEDB) | Distinct antibody epitope sequences per gene | 34 genes, 222 distinct epitopes | `IEDB bcell_search` |
| Codon usage bias (COMPUTED) | Effective number of codons, GC3, and codon adaptation index | 8,140 genes (100%) | `ToxoDB ME49` |
| Drug sensitivity per gene (CURATED) | Knockouts with a measured shift in sensitivity to a named compound | 3 genes, 2 compounds | PMID [41025776](https://pubmed.ncbi.nlm.nih.gov/41025776/); `mBio, PMID 41025776` |
| Enteric / sexual-cycle fitness per gene (CURATED) | Gene disruptions carried through the feline stage with oocyst output measured | 8 genes, 4 studies | PMID [28288194](https://pubmed.ncbi.nlm.nih.gov/28288194/); `PMIDs 28288194, 30728393, 36809045 and PMC12942651` |
| Enzyme classification (ToxoDB) | EC number per gene, and whether it has one | 1,313 enzymes of 8,140 genes | `ToxoDB ME49` |
| IEDB epitopes mapped to genes (via ToxoDB) | How many IEDB epitopes ToxoDB maps to this gene | 221 genes | `ToxoDB / IEDB` |
| InterPro domains | Domain identity and count | 8,140 | *citation not yet confirmed* |
| Membrane lipid composition of parasite vesicles | Lipid species abundance, and its proportion against the host cell | 194 lipid species | PMID [41716462](https://pubmed.ncbi.nlm.nih.gov/41716462/); `Front Cell Infect Microbiol 1745625 Tables 1-3` |
| Metabolome and isotope labelling under iron deprivation | Steady-state metabolite levels and the fraction labelled from glucose or glutamine | 1,102 metabolites | PMID [41925342](https://pubmed.ncbi.nlm.nih.gov/41925342/); `mBio 03788-25 Tables S3 and S5` |
| OrthoMCL orthogroups | Orthogroup assignment and cross-species bridge | 16,793 groups | `OrthoMCL release 6.21` |
| Plasmodium T-cell epitopes (IEDB) | Distinct T-cell epitope sequences per gene | 44 antigens, 1,542 distinct epitopes | `IEDB tcell_search` |
| Plasmodium antibody epitopes (IEDB) | Distinct antibody epitope sequences per gene | 434 antigens, 7,366 distinct epitopes | `IEDB bcell_search` |
| Plasmodium codon usage (COMPUTED) | Effective number of codons, GC3, and CAI against the ribosomal proteins | 5,318 genes | `PlasmoDB-68 Pfalciparum3D7 AnnotatedCDSs` |
| Plasmodium complexes from crosslinking MS | Which crosslink-derived complex a gene belongs to, and whether it reaches the host | 128 genes in 42 complexes | PMID [41966402](https://pubmed.ncbi.nlm.nih.gov/41966402/); `Cell Rep mmc5 Clusters` |
| Plasmodium crosslinking MS contacts | Protein pairs joined by a measured crosslink | 79 parasite-parasite pairs | PMID [41966402](https://pubmed.ncbi.nlm.nih.gov/41966402/); `Cell Rep mmc1 sheet D` |
| Plasmodium enzyme classification (PlasmoDB) | EC number per gene, curated and orthology-derived kept apart | 1,220 curated, 335 more from orthology | `PlasmoDB ec_numbers and ec_numbers_derived` |
| Plasmodium export prediction (ExportPred) | Predicted export to the erythrocyte, as an ordinal confidence tier | 440 genes called at some threshold, 191 at the default | `PlasmoDB GenesByExportPrediction` |
| Plasmodium falciparum 3D7 gene attributes | The second species: sequence, orthology, domains, strain SNPs and piggyBac fitness | 5,720 P. falciparum genes | `PlasmoDB GenesByTaxon attributesTabular` |
| Plasmodium host interaction degree (COMPUTED) | How many host proteins a gene was crosslinked to, where it was looked at | 117 genes seen, 10 with a host partner | *citation not yet confirmed* |
| Plasmodium model confidence and disorder (AlphaFold DB) | Mean pLDDT per protein, and the fraction of it at each confidence band | 5,098 of 5,720 genes | `AlphaFold DB API, per UniProt accession` |
| Plasmodium relation layers (COMPUTED) | Gene pairs sharing an orthogroup or a domain, and pairs whose stages covary | 1,741 + 24,123 + 63,158 pairs | *citation not yet confirmed* |
| Plasmodium to human contacts (crosslinking MS) | Parasite protein to erythrocyte protein, measured as a crosslink | 10 pairs, 10 parasite genes, 7 human proteins | PMID [41966402](https://pubmed.ncbi.nlm.nih.gov/41966402/); `Cell Rep mmc1 sheet D` |
| PubMed Central open-access full texts | Sectioned JATS XML | 6,667 articles | *citation not yet confirmed* |
| PubMed abstracts | Titles and abstracts for co-mention and attention | 33,924 records | *citation not yet confirmed* |
| Strain variation (ToxoDB HTS SNPs) | SNPs per gene across every sequenced strain, split by effect | 8,140 genes (100%) | `ToxoDB ME49` |
| ToxoDB gene identity | Symbols, previous IDs, product descriptions | 8,843 ME49 genes | `ToxoDB ME49` |
| Validated resistance-conferring mutations (CURATED) | Mutations shown to CAUSE drug resistance by putting them back into a clean background | 1 gene, 3 substitutions, 3 compounds | PMID [24533298](https://pubmed.ncbi.nlm.nih.gov/24533298/); `Int J Parasitol Drugs Drug Resist, PMIDs 24533298 and 25941623` |
## Licence and citing

Cite the original studies, not this table. The repository redistributes derived facts and the identity
tables, never the source articles' files.

The cell drawing under the compartment list is **SwissBioPics** artwork by Philippe Le Mercier (SIB
Swiss Institute of Bioinformatics), used under [CC BY 4.0](http://creativecommons.org/licenses/by/4.0).
The program shows the outlines only, filled from its own palette; the attribution travels with the
widget, whose tooltip carries it, because the artwork's own credit block is not drawn.
