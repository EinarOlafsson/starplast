# starplast — Materials and Methods

Reference for the methods section of the PLOS ONE manuscript. Every number here was read from the built
cache (`data/nodes.parquet`, `data/graph.npz`, `data/mentions.parquet`, `data/crosslink_models.parquet`)
at commit `29b0bd9`, not quoted from earlier drafts.

> **Citations marked ⚠ must be confirmed before submission.** They are datasets whose *content* is
> verified on disk but whose originating publication is not recorded in the pipeline. Do not cite them
> from this file alone.

---

## 1. Data sources

### 1.1 Gene universe and annotation

| Data | Provides | Coverage | Source | On disk |
|---|---|---|---|---|
| Node table | Gene universe, protein features | 8,227 → **8,140** after de-duplication on `gene_id` | upstream `toxonet` build | `toxonet/data/interim/nodes.parquet` |
| Gene products | Product descriptions | 8,140 | OrthoMCL / ToxoDB ME49 | `datasets/orthomcl_toxoplasma_gondii_ME49.csv` |
| InterPro | Domain identity and count | 8,140 (`n_interpro`) | InterPro ⚠ version | `datasets/interpro_tgon.csv` |
| Gene identity | Symbols, previous IDs, GT1/VEG accessions | 8,843 ME49 / 8,637 GT1 / 8,563 VEG | **ToxoDB REST API**, retrieved 2026-08-11 | `data/toxodb_identity.tsv`, `data/toxodb_strain_{gt1,veg}.tsv` |
| AlphaFold confidence | `mean_plddt` | 6,480 (79.6%) | AlphaFold DB ⚠ version | in node table |

### 1.2 Localisation

| Data | Provides | Coverage | Source | On disk |
|---|---|---|---|---|
| *T. gondii* hyperLOPIT | `lopit_map`, `lopit_mcmc`, posteriors; 26 compartments | **3,827 / 8,140 (47.0%)** | ⚠ Barylyuk et al., *T. gondii* hyperLOPIT | `datasets/lopit_toxoplasma_gondii_ME49.csv` |
| *P. falciparum* LOPIT | Donor labels for transfer | 3,000 rows, 1,646 usable after vocabulary mapping | ⚠ | `datasets/lopit_plasmodium_falciparum_3D7.csv` |
| *C. parvum* hyperLOPIT | Donor labels for transfer | 1,742 assigned, 1,107 usable after mapping | ⚠ Guérin et al. 2023 | `datasets/lopit_cryptosporidium_parvum_MEASURED_Guerin2023.csv` |
| Vocabulary dictionary | Maps 87 species-specific compartment terms → **12 unified categories** | tgon 26, cpar 22, tbru 20, pfal 18 terms | this project | `datasets/lopit_vocabulary_dictionary.csv` |
| Orthogroup bridge | Tg↔Pf↔Cp↔Tb gene correspondence | 16,793 orthogroups | OrthoMCL | `datasets/MASTER_parasite_wide_by_orthogroup.csv` |

### 1.3 Transcriptomics

| Data | Provides | Coverage | Source |
|---|---|---|---|
| Stage series | Tachyzoite, day 3/5/7, in vivo tissue cyst; 12 columns | 7,739 (95.1%) | **GEO GSE108740** |
| Oocyst sporulation | Unsporulated / sporulating / sporulated, 2 replicates each; 6 columns | 7,974 (98.0%) | **GEO GSE206344** |

All 18 raw FPKM columns ship. Derived summaries: `expr_tachy`, `expr_cyst`, `expr_max`, `expr_sporulated`
(log₂(mean FPKM + 1) over the relevant columns).

### 1.4 Genetic screens

| Screen | Provides | Coverage | Source |
|---|---|---|---|
| In vitro CRISPR (HFF) | `fit_invitro_hff` | 7,325 (90.0%) | ⚠ *confirm* — genome-wide fitness screen |
| In vivo CRISPR: peritoneum, lung, liver, spleen | `fit_invivo_{PE,lung,liver,spleen}` | 7,395 (90.8%) | ToxoDB attributes `tgonGt1CrisprFunc*`; composite scores correspond to **PMID 31481656** |
| Naive BMDM, IFN-γ | `fit_naive_bmdm`, `fit_ifng` | 7,402 (90.9%) | ⚠ *confirm* — macrophage screens |
| Young 2019 in vivo | `fit_invivo_young2019` | **115** | ⚠ *confirm* |
| **GRA17 synthetic-lethal** | RHΔ*gra17* phenotype, Δ−WT differential, candidate flag | **7,553** (genome-wide) | **PMID 37498952 / PMC10409377** |
| **GRA12, screen 1 and 2** | Median L2FC in vitro / in vivo, DISCO score | 236 / 232 | **PMID 40240328 / PMC12003902** |
| **In vivo CRISPR platform** | Mean log fold-change across replicates | 168 | **PMID 31481656 / PMC6722137** |
| **Host-transcription effectors** | Hotelling *T*² statistic, adjusted *p* | 252 | **PMID 37827122 / PMC12033024** |

The four named screens were retrieved from publisher supplementary material (PMC `/bin/` paths return
404; PLOS `article/file?id=…&type=supplementary`, Springer `static-content.springer.com/esm/…` and the
Europe PMC REST `supplementaryFiles` endpoint all resolve). Stored under
`datasets/DNA/CRISPR_screen/<PMID>/` with `META.json` and `SOURCES.md`.

**The screens share neither a sign convention nor a scale.** Naive-BMDM and IFN-γ are inverted relative to
the others and standard deviations differ ~15-fold. Rank-normalise and verify direction before pooling.

### 1.5 Protein-level measurements

| Data | Provides | Coverage | Source |
|---|---|---|---|
| Phosphosites | `n_phosphosites` (count only, no positions) | **1,175 (14.4%)** | ⚠ *confirm* |
| Protein abundance | median log₂ iBAQ across replicates | **748 (9.2%)** | **PXD043808**, **PXD065585** (Pru) |

**This is not a proteome.** Both are immunoprecipitation experiments (424 and 594 proteins); the figure is
enrichment, not whole-cell abundance, and must not be reported as proteome-wide coverage.

### 1.6 Interactions and structure

| Data | Provides | Coverage | Source |
|---|---|---|---|
| Crosslink MS | `xlms` edges, weight = crosslink count | 2,842 pairs / 1,630 genes | **StarPath** DSS XL-MS ⚠ cite the StarPath resource |
| Crosslink complexes | 4 Chai-1 models per pair, residue positions, satisfaction | 2,843 pairs; 2,439 with local CIFs | StarPath + `starpath_crosslink_mining` |
| IP-MS | `ip_ms` edges | 64 pairs / 48 genes | **PXD043808**, **PXD065585** |
| Structural similarity | `struct` edges, Foldseek TM-align **TM ≥ 0.7** | 11,684 pairs / 2,338 genes | Foldseek over **6,900** *T. gondii* AlphaFold models |
| Proximity-labelling corpus | 42 BioID/TurboID/APEX studies, 28 with data | 127 files | catalogued from PubMed; `datasets/post_translation/BioID/<PMID>/` |
| Pulldown corpus | 55 IP-MS/co-IP studies, 29 with data | 140 files | `datasets/post_translation/IPMS/<PMID>/` |

> The proximity/pulldown corpora are **downloaded and indexed but not yet parsed into edges**. Do not
> describe them as integrated.

**StarPath uses RH88 accessions whose numbering does not correspond to ME49** (`TGRH88_016370` is
`TGME49_210408`, not `TGME49_216370`). Mapping uses the alias column shipped with the export.

### 1.7 Literature

| Data | Provides | Coverage | Source |
|---|---|---|---|
| Abstracts | Title + abstract | **33,924** records | PubMed, *Toxoplasma* query |
| Open-access full texts | Sectioned JATS XML | **6,667** articles | PubMed Central OA |

Both are machine-local; the committed cache is what the application ships.

---

## 2. What each component does

The pipeline is four resolution layers feeding one builder; the application reads only the built cache.

### 2.1 `identity.py` — gene identity

Resolves every string the literature uses for a gene to one canonical ME49 accession. Match kinds, in
descending confidence: `accession` (current), `accession_prev` (pre-2012), `accession_strain` (GT1/VEG),
`symbol`, `alias` (Tg-prefixed).

- **Cross-strain mapping is by numeric suffix**, validated rather than assumed: among ME49/GT1 pairs
  sharing a suffix and an OrthoMCL release, **99.53%** share an orthogroup (VEG: **99.46%**).
- **Ambiguity is recorded, never resolved by guessing.** A string claimed by two genes at the same
  confidence tier is withdrawn and retained for reporting (**153** strings). Across tiers the more
  specific identifier wins, so a gene's current accession is never withdrawn because another gene once
  carried it.
- **Precision guards**, each pinned by a regression test: digit-free symbols require an upper-case surface
  form (`HOOK`, `CLAMP`, `CLIP`, `SPARK`, `REMIND` are all valid symbols and ordinary English words);
  tokenisation is Unicode-aware (an ASCII class truncates French *Santé* to `Sant`, which is a symbol);
  an accession's prefix is not re-read as a symbol (`TGGT1_…` is not the gene symbolled `GT1`); and
  *Toxoplasma* strain designations are blocked.

Final index: 19,368 resolvable strings plus 15,192 strain accessions.

### 2.2 `corpus.py` — documents

Presents abstracts and full texts as one `Document` stream of typed `Section`s (`title`, `abstract`,
`body`, `caption`). **Reference lists are excluded**: a bibliography names every gene its cited papers
named, which would credit a paper with the content of its citation graph. Figure/table captions are
emitted once (an unguarded body scan double-counts them, since JATS wraps caption text in `<p>`).

### 2.3 `literature.py` — mentions

Writes `data/mentions.parquet`: one row per gene × document × section × match kind, with source and
confidence tier. Every literature figure derives from this table, so nothing requires re-scanning ~40,000
documents.

**Co-mention units.** For an abstract the unit is the record; for a full text it is the **paragraph**,
because two genes named in one paragraph plausibly stand in a relation while two genes named anywhere in a
10,000-word paper mostly do not. Units naming more than **12** genes are excluded as lists or hit tables.

**Attention correction** (on by default). For two genes appearing in *n*₁ and *n*₂ of *N* units:

    expected = n₁ · n₂ / N
    residual = log₂( (observed + 0.5) / (expected + 0.5) )

*N* is the number of **units** — independence is defined over units — and *n*ᵢ are counted over the same
population co-mention is counted over, i.e. after list-like units are excluded.

**Depth of attention** is read off document structure, not assigned as weights:

| tier | definition | genes |
|---|---|---|
| `focal` | named in a paper **title** | 286 |
| `substantive` | named in an **abstract** | 464 |
| `incidental` | named only in **body or caption** | 1,816 |
| — | named nowhere | 5,574 |

### 2.4 `localisation.py` — hyperLOPIT and orthoLOPIT

- MAP and MCMC assignments are kept as separate columns. They **disagree for 980 of 3,827 genes (26%)**;
  `lopit_methods_agree` exposes this.
- **orthoLOPIT**: for a gene with no native call, measured labels of its *P. falciparum* and *C. parvum*
  orthologs are consulted. All terms map through the vocabulary dictionary to 12 unified categories first;
  a label transfers only when **every donor agrees**. *T. brucei* is excluded (not apicomplexan; its
  compartment set lacks the apical secretory organelles). Adds **126** genes → **3,953 / 8,140 (48.6%)**.
- `compartment` remains measured; `compartment_best` carries the fallback and `compartment_source`
  records provenance (`hyperLOPIT` / `orthoLOPIT` / `unknown`).

*Note for the manuscript:* `MASTER_parasite_wide_by_orthogroup.csv` ships a `cpar_lopit_native` column
that is empty on every row despite the measured Crypto data existing; this pipeline joins from source.

### 2.5 `screens.py`, `interactions.py`, `structures.py`

`screens.py` normalises the four published screens and the iBAQ tables to one row per gene. **Accessions
are resolved through `identity.py`**: the 2019 in vivo screen cites pre-2012 accessions for every gene and
contributed 0 of 8,140 rows until routed through the identity layer, then 168.

`interactions.py` lifts `xlms`, `ip_ms` and `struct` from the upstream edge table (already ME49-keyed) and
assembles `crosslink_models.parquet` — residue positions, model files, and whether the predicted pose
places the crosslinked residues within reach.

`structures.py` resolves coordinates **on demand** (local mirror, else AlphaFold DB, cached in
`~/.cache/starplast`). Coordinates are the only thing not shipped.

### 2.6 `build_graph.py` — the built cache

**Embedding.** 3D UMAP (`n_neighbors=25`, `min_dist=0.25`, Euclidean, `random_state=0`) over 43 columns:
16 numeric (3 expression, 7 fitness, pLDDT, paralogs, InterPro count, phosphosites, has-domain,
lineage-specific), z-scored after median imputation, concatenated with 27 hyperLOPIT one-hot columns × 0.5.

**Edge types** (12; never merged):

| edge | n | definition |
|---|---|---|
| `comention` | 435 | ≥2 shared abstracts |
| `comention_ft` | 7,733 | ≥2 shared full-text paragraphs |
| `xlms` | 2,842 | DSS crosslink MS; weight = crosslink count |
| `ip_ms` | 64 | replicated pulldown vs untagged control |
| `struct` | 11,684 | Foldseek TM ≥ 0.7 |
| `orthogroup` | 3,452 | shared orthogroup, groups ≤ 60 |
| `coexpression` | 49,293 | top-25 neighbours, *r* ≥ 0.95 |
| `cofitness` | 88,997 | top-25 neighbours, *r* ≥ 0.90 |
| `compartment` | 118,712 | shared hyperLOPIT class, groups ≤ 250 |
| `domain` | 10,399 | shared InterPro domain, 2–60 members |
| `structural_hole` | 255 | **derived** — see 3.1 |
| `unwritten_interaction` | 2,673 | **derived** — see 3.2 |

### 2.7 `app.py` — the browser

PyQt6 + OpenGL. Level of detail follows the data (hyperLOPIT compartment → orthogroup/module → gene →
that gene's evidence). Six colour modes. Points use translucent blending with depth testing — additive
blending saturates 8,140 overlapping points to white and destroys the colour encoding. Edge opacity scales
with edge weight, without which the attention toggle is visually inert.

---

## 3. Derived measures

### 3.1 Structural holes (n = 255, over 457 genes)

A gene pair that **co-expresses across the stage series AND co-behaves across the seven CRISPR screens,
yet appears together in no abstract and no open-access paragraph.**

Two confound controls are load-bearing and should be stated in the manuscript:

1. **Homology cannot be one of the two legs.** `orthogroup` and `domain` are one evidence family, not two,
   because paralogs almost always share domains. Counting them separately made 53 of the first 66
   candidates pure paralogy; allowing homology to pair with expression admitted 291 further pairs of which
   **76% were same-orthogroup**. Requiring both independent *phenotype* measurements leaves 3 paralogs
   among 255 pairs.
2. **`compartment` is excluded entirely** — 118,712 edges is too unspecific, and hyperLOPIT assignment
   tracks abundance, so it would preferentially link well-expressed (hence well-studied) genes.

Six pairs have both endpoints well studied.

### 3.2 Unwritten interactions (n = 2,673, over 1,602 genes)

Pairs **measured** to interact (`xlms` ∪ `ip_ms`) that appear together in no abstract and no open-access
paragraph: **2,673 of 2,906 measured pairs (92%)**, of which **147** join two genes that are each
individually well studied. A stronger claim than a structural hole — the interaction was observed, not
predicted. The union of the two source types is explicit and labelled; both remain separately selectable.

### 3.3 Literature coverage

| | genes | of proteome |
|---|---|---|
| named anywhere | 2,566 | 31.5% |
| …only in passing (body/caption) | 1,816 | 71% of coverage |
| **named in a title or abstract** | **750** | **9.2%** |
| named in a title | 286 | 3.5% |
| named nowhere | 5,574 | 68.5% |

Across 5,690 distinct documents. Coverage rose from 601 to 2,566 on adding full texts and the identity
layer, but **attention rose only 601 → 750**: nearly everything full texts add is a gene in a hit table.

---

## 4. Limitations to disclose

1. **Median imputation leaks study effort into node position.** Missing values are imputed to the column
   median before embedding. Genes lacking fitness data sit **0.38 map-radii** from those that have it
   (pLDDT 0.21; phosphosites 0.22). `n_phosphosites` is missing for 85.6% of genes, so that column is
   effectively a binary "appeared in a phosphoproteomics experiment" indicator.
2. **hyperLOPIT contributes only 1.1% of the embedding's feature variance** (0.17 of 16.17). It is listed
   as a design input but does not materially affect position.
3. **UMAP preserves local neighbourhoods, not global distances.** "These genes are neighbours" is
   interpretable; "these clusters are far apart" is not. Both derived edge types are computed from edges,
   never from embedding distance, and are unaffected by 1–3.
4. **The predicted complexes mostly fail.** Of 2,397 scored crosslink pairs, **60% place no crosslink
   within reach**; median interface ipTM is **0.16**; only **162** are both confident (ipTM ≥ 0.5) and
   crosslink-consistent. Many failures are dense-granule proteins, which are disordered. The crosslink is
   the measurement; the model is a hypothesis about the pose.
5. **Coverage is uneven and "not measured" is not "no effect".** Targeted libraries leave most genes
   untested (GRA12: 236 of 8,140). The interface renders these as "—".
6. **The two literature sources are different populations.** Abstracts cover the field; full texts are
   only what publishers deposited open access. Full-text coverage cannot be quoted as coverage of
   *Toxoplasma* research.
7. **The full-text corpus grew during development** (5,493 → 6,667 files). `build_graph` logs the count it
   used; that log is authoritative over any written figure.
8. **The GRA12 screens are not replicates** (in-vivo L2FC *r* = 0.41) and are kept as separate columns.

---

## 5. Reproducibility

```bash
pip install -e .
python -m starplast.fetch_names     # ToxoDB identity tables (network)
python -m starplast.build_graph     # rebuilds data/ (~5 min)
pytest tests/ -q                    # 51 tests, headless, no network
starplast
```

Deterministic given fixed inputs: UMAP uses `random_state=0`; identity, screen and literature layers are
pure functions of their input files. The application requires no network and no dataset; the committed
cache (11 MB) carries all 95 columns for 8,140 genes.

Software: Python 3.10, PyQt6 6.7.1, pyqtgraph 0.13.7, umap-learn, scikit-learn, pandas, NumPy, pyarrow.
⚠ Record exact versions with `pip freeze` at submission.
