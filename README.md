# starplast

A 3D browser for the *Toxoplasma gondii* knowledge map. 8,140 genes as points in space, twelve kinds of
relation as toggleable edges, and one panel per gene showing everything that is actually known about it —
including when the answer is "nothing".

```bash
pip install -e .
python -m starplast.fetch_names     # one-off: ToxoDB symbols, previous IDs, strain accessions
python -m starplast.build_graph     # one-off, ~4 min: builds data/graph.npz
starplast

pytest tests/ -q                    # 51 tests, headless, no network
```

Needs a display (PyQt6 + OpenGL). **The app is standalone**: the committed cache (11 MB) carries every
measurement for all 8,140 genes — 85 columns covering 8 CRISPR fitness screens, 4 published screens, the
full 18-column transcriptomic series, protein abundance, sequence and structure metadata — so `starplast`
needs no dataset, no network and no build step. `build_graph` is only for regenerating the cache.

The one deliberate exception is **coordinates**. 6,538 AlphaFold models and 12,265 crosslink complexes are
gigabytes, so `structures.py` resolves them on demand: a local mirror if present, else AlphaFold DB,
cached under `~/.cache/starplast`. Everything you need to *decide* is offline; only the picture is fetched.

## What makes it different from a network viewer

**Position means something — with two caveats worth knowing.** Node coordinates are a 3D UMAP embedding
of a 43-feature matrix: stage expression, seven CRISPR fitness screens, hyperLOPIT compartment, paralog
number, domain count, phosphosites, mean AlphaFold pLDDT. Two genes near each other are biologically
similar, and a force-directed layout would look similar while meaning nothing.

The caveats: missing values are median-imputed, so "was this gene measured" is faintly visible as geometry
(genes lacking fitness data sit 0.38 map-radii from those that have it); and the compartment block, at
×0.5 over 27 one-hot columns, contributes only **1.1%** of the matrix variance rather than the balanced
share the design intends. Both are open issues. Note that `structural_hole` and `unwritten_interaction`
are computed from **edges**, never from embedding distance, so neither is affected.

**Twelve edge types, never merged.**

| edge | n | source |
|---|---|---|
| `xlms` | 2,842 | **StarPath DSS crosslink MS** — measured physical proximity in a cell lysate |
| `ip_ms` | 64 | replicated IP-MS of a tagged bait, against an untagged control |
| `struct` | 11,684 | **Foldseek TM-align ≥ 0.7** over 6,900 Toxoplasma AlphaFold models |
| `comention` | 435 | 33,924 PubMed abstracts (2 shared abstracts minimum) |
| `comention_ft` | 7,733 | 6,667 open-access full texts, **per paragraph** (2 shared paragraphs minimum) |
| `orthogroup` | 3,452 | OrthoMCL |
| `coexpression` | 49,293 | GSE108740 stage series, top-25 neighbours at r ≥ 0.95 |
| `compartment` | 118,712 | hyperLOPIT |
| `cofitness` | 88,997 | 7 CRISPR screens, top-25 at r ≥ 0.90 |
| `domain` | 10,399 | shared InterPro domain |
| `structural_hole` | 255 | **derived** — see below |
| `unwritten_interaction` | 2,673 | **derived** — measured to bind, never written about |

The first three are qualitatively different from everything else: they are **measurements, not
correlations and not text**. A crosslink says two residues were covalently joined in a living-cell lysate,
so those proteins were within the crosslinker's reach. Structural similarity needs no orthology, so it
reaches lineage-specific effectors that homology edges cannot see. Crucially, none of the three is
attention-biased — they are the first edge types here that can connect a gene nobody has written about to
one everybody has.

They answer different questions and disagree with each other; merging them into one "relatedness" score
would be the single easiest way to make this tool lie. The two co-mention types stay separate for the same
reason: one is drawn from every abstract in the field, the other only from papers a publisher deposited
open access. Full-text co-mention is counted per *paragraph*, because two genes named in one paragraph
plausibly stand in a relation while two genes named anywhere in a 10,000-word paper mostly do not.

**Attention correction is on by default.** Raw co-mention edges reproduce the literature's popularity
contest — ROP18 and GRA16 become hubs because they are studied, not because they are central. The default
view shows log2(observed / expected) given each gene's own publication count, so an edge means *more
co-mention than these two genes' individual fame predicts*. Switch it off and the status bar says
`RAW co-mention (attention-biased)`.

It is what makes the full-text layer worth having. Ranked by raw count, full-text co-mention returns the
famous pairs — ROP18/ROP5, SAG1/GRA6. Ranked by the corrected residual it returns HDAC3/MORC, MIC1/MIC4
and AP2XII-1/AP2XI-2: actual complexes, surfaced because they co-occur far more than their individual fame
predicts.

**Two numbers matter, and conflating them is the mistake this app exists to prevent.**

| | genes | of proteome |
|---|---|---|
| named **anywhere** in the readable literature | 2,566 | 31.5% |
| …of those, named only **in passing** (body or caption) | 1,816 | 71% of coverage |
| named in a **title or abstract** — the honest attention figure | **750** | **9.2%** |
| named in a **title** — the paper is about it | 286 | 3.5% |
| named **nowhere** | 5,574 | 68.5% |

Coverage rose from 601 to 2,566 when full texts and a real identity layer were added. **Attention barely
moved: 601 → 750.** Almost everything the full texts add is a gene sitting in a screen's hit table, named
once and never discussed. Reporting 2,566 as "genes the field has studied" would repeat exactly the error
the attention correction exists to prevent, so the app tiers every gene by *where* it is named — `focal`
(title), `substantive` (abstract), `incidental` (body/caption only) — and a gene reached only through hit
tables is labelled **"Listed, not studied"** in its panel. The tiers are read off document structure, not
assigned as weights.

Genes below 7 abstracts get an explicit warning: absence of evidence there is absence of attention. Read
the two sources separately — abstracts cover the whole field, full texts are only the openly deposited
subset, so full-text coverage answers a different question and cannot be quoted as if it covered
*Toxoplasma* research generally.

**Structural holes: where the biology connects two genes and the literature never has.** This is what the
map is for — not what the field says, but where its own data says there is a crossing nobody has made.

A hole is a pair that **co-expresses across the stage series AND co-behaves across the seven CRISPR
screens, yet appears together in no abstract and no open-access paragraph**. 255 such pairs over 457
genes; 6 have both endpoints already well studied, which makes those the sharpest.

Two confound controls do the real work, and both were found by looking at what the first version returned:

- **Homology cannot be one of the two legs.** `orthogroup` and `domain` are one fact, not two — paralogs
  almost always share domains. Counting them separately made 53 of the first 66 candidates pure paralogy.
  Allowing homology to pair with expression admitted 291 more, **76% of them same-orthogroup**, because
  paralogs co-express *because* they are paralogs. Requiring both independent phenotype measurements
  leaves 3 paralogs in 255 pairs.
- **`compartment` is excluded entirely.** Sharing one of 27 hyperLOPIT classes is real co-localisation but
  hopelessly unspecific at 118,712 edges, and hyperLOPIT assignment tracks abundance — so it would
  preferentially link well-expressed genes, which are the ones already well studied.

A hole is a *derived* relation — the absence of a co-mention across a pair the measurements agree about —
and a hypothesis generator, not evidence. The sharpest six:

| | |
|---|---|
| LMF1 — ATPTG9 | mitochondrion–IMC tether · ATP synthase subunit |
| BFD2 — AC11 | bradyzoite formation deficient · apical cap protein |
| RON13 — ISP2 | rhoptry neck kinase · IMC sub-compartment protein |
| GRA44 — CLIP | dense granule protein · CLAMP-linked invasion protein |
| SRS51 — GRA12D | SAG-related sequence · dense granule protein |
| TGME49_224620 — BAG1 | hypothetical protein · bradyzoite antigen |

**92% of measured interactions have never been written about.** Of 2,906 protein pairs shown by
crosslink MS or IP-MS to physically interact, **2,673 appear together in no abstract and no open-access
paragraph** — and 147 of those join two genes that are *each* individually well studied. That is the
`unwritten_interaction` layer, and it is a stronger claim than a structural hole: the interaction was
observed, not predicted. It is an explicit, labelled merge of `xlms` and `ip_ms`; both stay separately
toggleable.

**How the binding happens — with the models, and with their failure rate stated.** StarPath ships four
Chai-1 predicted complexes per interaction (12,265 CIF files on local disk). `crosslink_models.parquet`
joins each gene pair to its crosslinked residue positions, its model files, and whether the predicted
pose actually puts those residues within reach. Selecting a gene lists its partners, the residues, and
the verdict.

The verdict is usually negative, and the app says so:

| of 2,397 scored pairs | |
|---|---|
| model explains **no** crosslink | 1,447 (60%) |
| explains ≥ 1 | 950 (40%) |
| median interface ipTM | **0.16** |
| **both confident (ipTM ≥ 0.5) and crosslink-consistent** | **162** |

So the crosslink is the measurement; the model is a guess at the pose, and only ~162 pairs have a picture
worth trusting. Those are flagged `model usable` in the panel. Many of the failures are dense-granule
proteins, which are largely disordered — exactly where structure prediction should be expected to fail.

**Genes are resolved through an identity layer, not a symbol table.** The literature does not use one
identifier for a gene. `starplast/identity.py` resolves all of them to one canonical ME49 accession:

| form | example | genes reached |
|---|---|---|
| current accession | `TGME49_208830` | 1,210 |
| previous accession | `TGME49_008830` (pre-2012, still cited) | 369 |
| strain accession | `TGGT1_208830`, `TGVEG_208830`, mapped by numeric suffix | 1,150 |
| symbol | `GRA16`, `GRA-16` | 1,127 |
| Tg-prefixed alias | `TgGRA16` | 744 |

Only the first row existed before. The other four are why coverage moved.

The suffix mapping is verified rather than assumed: ME49/GT1 pairs sharing a suffix and an OrthoMCL
release share an orthogroup 99.53% of the time (VEG: 99.46%). Strings claimed by two genes are withdrawn
and recorded, never guessed — 153 of them.

Matching precision is guarded, because a false match does not crash anything, it just hands a gene
attention it never had. Digit-free symbols must appear in upper case (`HOOK`, `CLAMP`, `CLIP`, `SPARK` and
`REMIND` are all real symbols and all real English words); tokenisation is Unicode-aware (an ASCII-only
class carved `Sant` out of French `Santé` and handed one gene 23 abstracts); an accession's own prefix is
not re-read as a symbol (`TGGT1_209030` is not a mention of the gene symbolled `GT1`); and *Toxoplasma*
strain designations are blocked outright. Each of those is pinned by a test.

## Finding structures that predict something they were never told

The point of the apparatus. `starplast/search.py` walks combinations of datasets and hyperparameters and
scores each by how well it recovers a label that was **excluded from the embedding**:

```python
from starplast.search import search
R, per_label = search(nodes, target="lopit_unified", sample_size=3000)
```

The target and everything that restates it are removed from every map first — for `lopit_unified` that is
six columns, including `compartment`, `lopit_map` and `lopit_mcmc`. Precision and recall are kept
separately, because a cluster that is 100% apicoplast but holds 5% of apicoplast proteins is useless for
inference while one holding 97% at 90% precision is the goal.

**What 164 runs actually found, and why the guard matters.** With hyperLOPIT leaked into the embedding,
the battery reports it separating clusters at V = 0.96 — a spectacular result, and an artefact. Held out
properly, the best structure over 164 combinations reaches **mean F1 0.362**:

| label | n | precision | recall | F1 | from |
|---|---|---|---|---|---|
| apical secretory | 95 | 0.68 | 0.53 | **0.59** | expression + fitness |
| nucleus | 398 | 0.46 | 0.62 | 0.53 | expression + fitness + interactions |
| cytosol | 132 | 0.27 | 0.81 | 0.41 | expression + fitness |
| ER | 88 | 0.34 | 0.51 | 0.41 | published screens + protein features |
| mitochondrion | 167 | 0.25 | 0.33 | 0.29 | fitness + screens + interactions |

So localisation is **only weakly predictable** from expression and fitness — apical secretory best, which
fits, since secretory proteins carry distinctive stage-expression profiles. That is a smaller claim than
the circular version, and it is the true one. Every run records its full recipe, seed and scores, so a hit
can be rebuilt exactly.

## Interpretation rules the UI enforces

- **"unassigned" is grey, not a 27th compartment.** hyperLOPIT assignment tracks protein abundance, so a
  missing call means unknown, never absent.
- **Fitness scores are competitive growth, not essentiality.** The in vitro screen is predictable from
  protein features (R² = 0.45); the five in vivo screens are not (−0.11 to +0.10). The panel says so.
- **Missing values are rendered grey, never mapped onto the colour scale.** Missingness here is
  informative — large secreted proteins are exactly the ones AlphaFold DB skips.
- **Colour maps match the data.** Sequential for ordered quantities, diverging only where values
  genuinely straddle a midpoint, categorical for unordered classes. A diverging map on a positive
  quantity invents a midpoint; a sequential map on a residual hides its sign.
- **Depth of attention is categorical, not a ramp.** `focal` / `substantive` / `incidental` are read off
  where a paper names a gene; shading them along a gradient would imply a measured quantity. Genes named
  nowhere stay grey with everything else that is unknown rather than zero.
- **Draw caps are stated.** Above 20,000 drawn edges the status bar reports what was dropped.
- **Edge opacity encodes edge weight.** At one flat alpha the strongest and weakest edges look alike, and
  the attention toggle — which reorders exactly that quantity — becomes visually almost a no-op.
- **Points occlude, they do not sum.** The scatter uses translucent blending with depth testing. With
  additive blending, 8,140 overlapping points saturated to white and every colour mode rendered as one
  featureless blob while the array-level tests all still passed.

## Controls

| | |
|---|---|
| left-drag / scroll | rotate / zoom |
| left-click a node | select; fills the evidence panel |
| search box | gene ID or product substring, then flies to it |
| double-click a compartment | fly to that compartment's centroid |
| level of detail | compartment (galaxy) → orthogroup (system) → gene (planet) |
| edge checkboxes | per-type; edges draw for the selected gene unless "draw all" is on |

## How the literature layer is put together

Four small modules instead of one scan, so every figure can be re-derived without re-reading 40,000
documents:

| module | job |
|---|---|
| `identity.py` | every string the literature uses for a gene → one canonical accession |
| `corpus.py` | abstracts JSONL and PMC JATS XML → one `Document` stream, sectioned, references excluded |
| `literature.py` | → `data/mentions.parquet`, the tidy table everything else derives from |
| `build_graph.py` | composes those into edges, node columns and the UMAP embedding |

`data/mentions.parquet` is the auditable intermediate: one row per gene × document × section ×
match kind, carrying the confidence tier. Coverage, publication counts and both co-mention layers are all
recomputable from it.

## What ships, and what it covers

Coverage varies enormously by assay, and the app shows "—" for not-measured rather than implying zero:

| data | columns | genes covered |
|---|---|---|
| CRISPR fitness (in vitro, 4× in vivo, naive BMDM, IFN-γ) | 7 | ~90% |
| CRISPR, Young 2019 in vivo | 1 | 115 |
| **GRA17 synthetic-lethal, genome-wide** (PMC10409377) | 3 | **7,553** |
| GRA12 strains/subspecies, 2 screens (PMC12003902) | 6 | 236 / 232 |
| in vivo CRISPR platform (PMC6722137) | 1 | 168 |
| host-transcription effectors (PMC12033024) | 2 | 252 |
| transcriptomics: stage series + **oocyst sporulation** | 18 + 4 | ~95% |
| protein abundance, log2 iBAQ (Pru IP) | 1 | 748 |
| phosphosite count | 1 | 14.4% |
| AlphaFold pLDDT | 1 | 79.6% |
| sequence, length, TM, signal peptide, LOPIT posterior, Pfam, InterPro | 12 | ~100% |
| crosslink / IP-MS / structural-similarity partners | 4 | edge-derived |

Two honest gaps. The **GRA12 screens are not replicates** — their in-vivo L2FCs correlate at r = 0.41 —
so they are kept as separate columns. And **there is no deep proteome in this tree**: the only mass-spec
abundance is two Pru immunoprecipitation experiments (424 and 594 proteins), which is enrichment, not
coverage.

## Data sources

Built from `/mnt/firecuda2/Claude/toxoplasma_projects` (home) — see `HANDOFF.md` for the full table and the
work-machine paths. `data/toxodb_identity.tsv` and `data/toxodb_strain_{gt1,veg}.tsv` are committed and
regenerable with `python -m starplast.fetch_names`. The open-access full texts live on a machine-local
disk; where that disk is absent the build falls back to abstracts only and says so.

## Status

v0 + v1 complete: precomputed embedding, GL scatter, picking, evidence panel, edge toggles, attention
correction, level-of-detail. v1.1 adds the identity layer, the full-text layer and a real test suite.
v2 (species switching via cross-species orthology) and v3 (continuous star-map zoom) are deliberately
deferred — see `HANDOFF.md`.
