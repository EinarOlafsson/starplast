# Datasets left to fold in — the Toxoplasma map

Written 2026-08-13, from three sources: an audit of the shipped cache, the raw dataset tree on this
machine, and a survey of **24,679 PubMed records** with *Toxoplasma* or *gondii* in the title or
abstract (every one NCBI holds; the query has to be sliced by year because PubMed will not return
more than 9,999 records for one query, however it is paged).

Of those, 913 carry a signal that they produced a dataset. 181 of them are open access and name an
accession in their full text, which is the only way to tell a dataset from a figure without reading
the paper — abstracts almost never name a deposit. 105 of those measure the PARASITE rather than the
host, which is the only shape this map can hold: a value per Toxoplasma gene.

The *Plasmodium* map is shelved by decision; see the note at the end.

---

## 1. Already in the cache, and no map can see it

**82 of the 149 numeric columns in the shipped cache match no feature block.** They are built,
committed, and unreachable from every embedding the program can make:

| family | columns | genes measured | what it is |
|---|---|---|---|
| `morc_` | 18 | 7,841 | MORC depletion series |
| `stress_` | 10 | 7,880 | stress / bradyzoite conversion |
| `invivo_` | 14 | 7,663 | in vivo screen replicates |
| `proteome_` | 15 | 3,005 | a proteome **four times larger** than the 748-protein one the registry documents |
| `oocyst_` | 8 | 2,079 | oocyst iTRAQ |
| `phospho_` | 5 | 1,603 | phosphosite up/down ratios |

None of them has a registry entry either, so `search.excluded_for` cannot reason about them — and
provenance gaps are exactly where the circularity guard leaked on 2026-08-13. **Nothing to download.
This is block regexes and registry entries, and it is the largest gain per hour available.**

## 2. Downloaded, not ingested

- **Dual perturb-seq (PMID 37827122, `GSE229505`).** In the cache as two columns for 252 genes:
  `hosttx_T2` and `hosttx_padj`. What is missing is the per-effector host transcriptional
  signature — the only host-side readout this project could have, and the only way to ask which
  effectors do the same thing to the host cell.
- **The BioID corpus**: 44 studies, 127 files, on disk under `datasets/post_translation/BioID/`.
- **The IP-MS corpus**: 55 studies, 140 files. The map currently carries **64** `ip_ms` pairs.

## 3. The survey: what to add, grouped by the axis it opens

Each row is a study that measured Toxoplasma genes and deposited the data. Accessions are as named
in the paper's own full text.

### Regulation — the axis the map has none of

| PMID | year | accession | what it would add |
|---|---|---|---|
| [35164683](https://pubmed.ncbi.nlm.nih.gov/35164683/) | 2022 | GSE104347, GSE87834 | histone variants genome-wide; H2A.X/H2A.Z occupancy per gene |
| [29788176](https://pubmed.ncbi.nlm.nih.gov/29788176/) | 2018 | GSE106864, GSE109086 | ApiAP2 cooperative binding at virulence genes; ChIP-seq |
| [39774584](https://pubmed.ncbi.nlm.nih.gov/39774584/) | 2024 | PRJNA1134678 | the ApiAP2 cascade that builds daughter cells; RNA-seq per factor |
| [40920096](https://pubmed.ncbi.nlm.nih.gov/40920096/) | 2025 | PRJNA1281009 | TgAP2X-7, a cell-cycle transcription factor |
| [41193469](https://pubmed.ncbi.nlm.nih.gov/41193469/) | 2025 | GSE287537, PXD060237, PXD064041 | SWI/SNF remodeler: division and expression, with proteome |
| [40879678](https://pubmed.ncbi.nlm.nih.gov/40879678/) | 2025 | GSE286088, GSE286090 | GCN5a at telomeres, priming latency |
| [21179246](https://pubmed.ncbi.nlm.nih.gov/21179246/) | 2010 | GSE22100 | GCN5-A under alkaline stress, the classic bradyzoite trigger |

### Translation — orthogonal to every RNA column already here

| PMID | year | accession | what it would add |
|---|---|---|---|
| [29228904](https://pubmed.ncbi.nlm.nih.gov/29228904/) | 2017 | GSE99395 | ribosome profiling across differentiation: translational efficiency per gene |
| [41925342](https://pubmed.ncbi.nlm.nih.gov/41925342/) | 2026 | PRJEB67890, PRJEB83013, PXD066828 | translational and metabolic remodelling under iron starvation |

### Post-translational marks beyond a phosphosite count

| PMID | year | accession | what it would add |
|---|---|---|---|
| [37562054](https://pubmed.ncbi.nlm.nih.gov/37562054/) | 2023 | PXD040368 | lysine acetylation across life stages |
| [37959749](https://pubmed.ncbi.nlm.nih.gov/37959749/) | 2023 | PXD046083 | cysteine S-nitrosylation, proteome-wide |
| [36216028](https://pubmed.ncbi.nlm.nih.gov/36216028/) | 2023 | PRJNA791485, PXD022700 | protein lactylation, with its own RNA-seq |
| [33363051](https://pubmed.ncbi.nlm.nih.gov/33363051/) | 2020 | PXD020655 | iTRAQ phosphoproteome of tachyzoites, site level |
| [31508380](https://pubmed.ncbi.nlm.nih.gov/31508380/) | 2019 | PXD007777 | phosphoproteome, type I vs type II |
| [35976251](https://pubmed.ncbi.nlm.nih.gov/35976251/) | 2022 | PXD033642, PXD033650, PXD033713, PXD033765 | thermal and temporal proteome profiling; PP1 substrates |
| [40901993](https://pubmed.ncbi.nlm.nih.gov/40901993/) | 2025 | PXD064226 | PP2A-2 during cytokinesis |
| [37933960](https://pubmed.ncbi.nlm.nih.gov/37933960/) | 2023 | PXD019677, PXD039426, PXD039431, PXD039432, PXD039434, PXD04408, PXD044080, PXD044081 | CDPK1 substrate map |

### Proteome — the current one is 748 proteins

| PMID | year | accession | what it would add |
|---|---|---|---|
| [40716488](https://pubmed.ncbi.nlm.nih.gov/40716488/) | 2025 | PXD063409 | cytoskeleton proteome |
| [28626452](https://pubmed.ncbi.nlm.nih.gov/28626452/) | 2017 | PXD003765 | iTRAQ across developmental stages |
| [25867681](https://pubmed.ncbi.nlm.nih.gov/25867681/) | 2015 | PXD000297, PXD000298 | proteogenomics, Toxoplasma and Neospora |

### Stages the map cannot currently see

| PMID | year | accession | what it would add |
|---|---|---|---|
| [41580398](https://pubmed.ncbi.nlm.nih.gov/41580398/) | 2026 | GSE311669 | bradyzoite subtypes |
| [25757795](https://pubmed.ncbi.nlm.nih.gov/25757795/) | 2015 | PRJEB7935 | merozoite expansion, distinct from tachyzoite |
| [29996885](https://pubmed.ncbi.nlm.nih.gov/29996885/) | 2018 | PRJNA475428 | stage-conversion transcriptomes |
| [35213559](https://pubmed.ncbi.nlm.nih.gov/35213559/) | 2022 | GSE191153 | clonal tachyzoite transcriptional heterogeneity |
| [42020723](https://pubmed.ncbi.nlm.nih.gov/42020723/) | 2026 | *not named in the abstract or the OA text* | single-cell atlas of sexual development in the feline gut (2026) |
| [41278634](https://pubmed.ncbi.nlm.nih.gov/41278634/) | 2025 | *not named in the abstract or the OA text* | single-cell brain-derived bradyzoites, AP2XI-6 cycle (2025) |

### Fitness under conditions the seven screens do not cover

| PMID | year | accession | what it would add |
|---|---|---|---|
| [39082802](https://pubmed.ncbi.nlm.nih.gov/39082802/) | 2024 | GSE253884, GSE253885 | in vivo virulence among hyperLOPIT-UNASSIGNED proteins - aimed at the same gap this map has |
| [34163449](https://pubmed.ncbi.nlm.nih.gov/34163449/) | 2021 | PRJNA707360 | genome-wide screen for oxidative-stress defence |
| [36916910](https://pubmed.ncbi.nlm.nih.gov/36916910/) | 2023 | GSE3920 | fitness under interferon gamma |

### RNA modification

| PMID | year | accession | what it would add |
|---|---|---|---|
| [34263725](https://pubmed.ncbi.nlm.nih.gov/34263725/) | 2021 | GSE168155, PRJNA705300, PXD024326 | m6A reading coupled to polyadenylation |
| [40830525](https://pubmed.ncbi.nlm.nih.gov/40830525/) | 2025 | GSE294543 | 5-methylcytosine methylome across lineages |

## What was deliberately left out

- **Host-side studies** — mouse astrocytes, decidual immune cells, cat organ proteomes, serum
  miRNA. Real data, no value per parasite gene. The exception is the dual perturb-seq matrix above,
  where the host response is indexed BY a parasite gene.
- **Metabolomics** (87 candidates) — measured per metabolite, not per gene. It would enter as a
  phenotype of a mutant set, not as a column.
- **Clinical, seroprevalence and veterinary survey work**, which is most of the corpus.
- **Single-protein studies.** A study of one protein is a fact; the map needs a column.

## The caveat that applies to every row above

Each of these adds columns measured on a subset, and **missingness is the map's oldest artefact**: a
gene absent from most assays is absent from most of the feature matrix, which is what made
`attention_depth` look like a recoverable target. Every dataset added here should be checked with
`embedding.missingness_leak` before its block is used in a search, and the block should carry its
coverage in the interface the way the existing ones do.

## Plasmodium, for the record

Shelved on 2026-08-13 in favour of finishing Toxoplasma. When it comes back: only **35.1%** of
Toxoplasma genes share an orthogroup with a *P. falciparum* gene and **27.3%** are one-to-one, and
the genes that lack an ortholog are the lineage-specific effectors this project exists to study — so
Plasmodium data belongs in a map of its own, not as columns on this one. 2.9 GB of PlasmoDB genomes,
proteins and GFF are already on disk.
