# starplast — point a session here

A 3D knowledge-map browser for *Toxoplasma gondii* (v0/v1), built to be the front end for the information
map. **Everything a fresh session needs is in this file.** Created 2026-08-10.

---

## What this is

A PyQt6 desktop app. Every gene is a node in a 3D space; edges are relations drawn from real datasets. You
click and zoom from compartment level down to a single gene's evidence. Position is **not** an arbitrary
force-directed layout — it is a UMAP embedding of a multimodal feature matrix, so proximity means biological
similarity.

Terminal entry point: `starplast`

## Run it

```bash
cd /mnt/firecuda2/Claude/toxoplasma_projects/starplast
pip install -e .            # installs the console_scripts entry point
python -m starplast.fetch_names   # one-off: ToxoDB identity tables (needs network)
python -m starplast.build_graph   # one-off: builds data/graph.npz (~4 min)
starplast                   # launch
pytest tests/ -q            # 48 tests, headless, no network, ~1 s
```

If the GL widget fails on a headless machine, that is expected — this needs a display. The test suite is
headless and does not.

## Design decisions, and why (do not silently reverse these)

**1. Position = UMAP of a feature matrix, not force-directed layout.**
Force-directed position is aesthetic and arbitrary; two adjacent nodes mean nothing. UMAP over
(expression × 7 fitness screens × compartment × orthology breadth × domain content × disorder) gives
positions where proximity is interpretable. Precomputed and cached — never laid out live.

**2. Twelve edge types, toggleable, never merged silently.**
"The knowledge map" is not one graph. Each edge type answers a different question:

| edge | source | note |
|---|---|---|
| `comention` | 33,924 PubMed abstracts | **attention-biased — see decision 3** |
| `comention_ft` | 6,667 open-access full texts, per paragraph | attention-biased; **OA subset, not the field** |
| `orthogroup` | OrthoMCL (16,794 groups) | clean; the cross-species bridge later |
| `coexpression` | GSE108740 stage series | quantitative |
| `compartment` | hyperLOPIT (26 compartments) | clean |
| `cofitness` | the 7 CRISPR screens | real and underused by the field |
| `domain` | InterPro shared domains | weakest; included for completeness |
| `xlms` | StarPath DSS crosslink MS | **measured physical proximity**; not attention-biased |
| `ip_ms` | replicated IP-MS vs untagged control | few, but the most direct binding evidence here |
| `struct` | Foldseek TM >= 0.7 over 6,900 AF models | needs no orthology, so it reaches lineage-specific effectors |
| `structural_hole` | **derived** from the above | biology links them, literature does not -- see decision 3c |
| `unwritten_interaction` | **derived** | measured to bind, never written about -- see decision 3d |

**2b. Gene identity is its own layer, and cross-strain accessions map by numeric suffix.** (Added v1.1.)
The literature calls one gene `TGME49_208830`, `TGME49_008830`, `TGGT1_208830`, `GRA16` and `TgGRA16`.
`identity.py` resolves all of these to one canonical accession, tagging each match with a confidence kind
so coverage can be reported by tier. `TGGT1_`/`TGVEG_` → `TGME49_` by numeric suffix is a VEuPathDB naming
convention, verified not assumed: 99.53% orthogroup agreement for GT1, 99.46% for VEG, within release.
Strings claimed by two genes at the same tier are withdrawn and recorded (153), never guessed.

**3. Attention correction is not optional decoration — it is a correctness feature.**
Raw co-mention edges reproduce the literature's popularity contest: GRA16 and ROP18 become bright hubs,
and the 78 genes named in ≤6 abstracts stay dark, regardless of biology. The app therefore ships an
**attention-corrected mode** that shows co-mention *residual from expected given each gene's publication
count*, and it is **on by default**. Turning it off shows raw co-mention and the UI says so. Without this
the app would be confidently misleading.

The full-text layer is the clearest demonstration that this works. Ranked by raw count it returns
ROP18/ROP5 and SAG1/GRA6 — fame. Ranked by corrected residual it returns HDAC3/MORC, MIC1/MIC4 and
AP2XII-1/AP2XI-2 — real complexes. Two details of the arithmetic matter and were both wrong at one point:
the expectation's denominator is the number of **units** (independence is defined over units, not over the
sum of per-gene counts), and per-gene unit counts are taken over the **same** population co-mention is
counted over, i.e. after list-like units are excluded.

**3b. What the map draws must be visible, and that is testable.** (Added v1.1.) The scatter blended
additively, so 8,140 overlapping points summed to white and every colour mode rendered as one blob --
while every array-level test passed, because they checked the colour array and never the render. Points
now occlude (translucent + depth test) and edge alpha scales with edge weight, without which the attention
toggle is visually almost a no-op. Both are pinned by tests that read the GL state. **Render the app and
look at it before believing a display claim.**

**3c. Structural holes are the point of the app, and their confound controls are not optional.**
(Added v1.1.) A hole is a pair that **co-expresses across the stage series AND co-behaves across the seven
CRISPR screens, yet appears in no abstract and no open-access paragraph**: the data says these belong
together and the field has never said it. 255 pairs over 457 genes, 6 with both endpoints already studied.

Two rules keep it from being a paralogy detector, and both were found by looking at output, not by
reasoning:

- **Homology (`orthogroup`/`domain`) is one family and can never be one of the two legs.** Counted as two,
  53 of the first 66 candidates were pure paralogy. Allowed to pair with expression, it admitted 291 more,
  **76% same-orthogroup** — paralogs co-express *because* they are paralogs, which is one fact twice. The
  strict rule leaves 3 paralogs in 255 pairs.
- **`compartment` is excluded entirely** — 118,712 edges is far too unspecific, and hyperLOPIT assignment
  tracks abundance, so it would preferentially link the well-expressed genes that are already well studied.

Relaxing either rule puts one protein family at the top of every ranking. A hole is *derived* — the
absence of a co-mention across a pair the measurements agree about — and is a hypothesis generator, not
evidence. Say so wherever it is used.

**3d. Measured binding is not co-mention, and 92% of it is unwritten.** (Added v1.2.) `xlms`, `ip_ms`
and `struct` come from work already done in this tree (`toxonet/data/interim/edges_v3.parquet`, already
keyed by `TGME49_`, so they bypass the identity layer). Of 2,906 measured interacting pairs, **2,673 are
in no abstract and no open-access paragraph**, and 147 of those join two genes that are each well studied.
That is `unwritten_interaction`: a stronger claim than a hole, because the interaction was observed rather
than predicted.

**The models come with their failure rate attached.** 12,265 Chai-1 CIFs sit in `starpath_dump/cifs`, and
`crosslink_models.parquet` joins each pair to its crosslinked residues, its model files and whether the
pose puts those residues in reach. **60% of scored models explain no crosslink at all; median interface
ipTM is 0.16; only 162 pairs are both confident and crosslink-consistent.** The crosslink is the
measurement, the model is a guess at the pose -- never present the picture as the evidence. Many failures
are dense-granule proteins, which are disordered, so this is expected rather than alarming.

> StarPath identifies proteins by **RH88** accession and that numbering does **not** correspond to ME49:
> `TGRH88_016370` is `TGME49_210408`, not `TGME49_216370`. Use the alias column the export ships. The
> suffix rule that works for GT1 and VEG (decision 2b) is wrong here and would silently mis-assign.

**4. Level of detail is data-driven, not invented tiers.**
galaxy = hyperLOPIT compartment (26) → solar system = orthogroup / co-expression module → planet = gene →
surface = that gene's evidence (papers, domains, screens, phenotypes).

## Data sources — all already on this disk

| file | what it gives |
|---|---|
| `../toxonet/data/interim/nodes.parquet` | 8,227 genes; **7 CRISPR fitness screens**, mean pLDDT, paralog number, orthogroup, domains, stage expression, phosphosites |
| `../datasets/lopit_toxoplasma_gondii_ME49.csv` | hyperLOPIT compartment + posterior |
| `../datasets/MASTER_parasite_wide_by_orthogroup.csv` | orthogroups across Tg/Pf/Cp/Tb |
| `../datasets/interpro_tgon.csv` | InterPro domains |
| `../datasets/orthomcl_toxoplasma_gondii_ME49.csv` | gene products (names) |
| `../datasets/stagetranscriptome_GSE108740_FPKM.xlsx` | tachyzoite, day 3/5/7, **in vivo tissue cyst** |
| `../.claude/skills/toxoplasma-scientist/corpus/pubmed_toxoplasma.jsonl` | **33,924 abstracts** with MeSH |
| `/mnt/wd4tb/skill_corpora/toxoplasma-scientist/*.xml` | open-access full texts, **6,667 at last build** (machine-local) |
| `data/toxodb_identity.tsv`, `data/toxodb_strain_{gt1,veg}.tsv` | symbols, previous IDs, GT1/VEG accessions (committed; `fetch_names.py`) |
| `../toxonet/data/interim/edges_v3.parquet` | **xlms / ip_ms / struct edges, already keyed by TGME49_** |
| `../starpath_{interactions.csv,crosslinks.json}` | crosslink pairs, residue positions, RH88->ME49 aliases |
| `../starpath_crosslink_mining/crosslink_satisfaction.csv` | do the Chai-1 models explain the crosslinks |
| `../starpath_dump/cifs/` | **12,265 Chai-1 complex CIFs, 6.7 GB** (in the synced tree) |
| `/mnt/wd4tb/skill_corpora/*/veupathdb/` | 431 genome/protein/GFF files, 9.3 GB, incl. coccidian relatives |

## Context you need to not repeat mistakes

Read `.claude/skills/toxoplasma-scientist/SKILL.md` before interpreting anything. In particular:

- **The CRISPR "fitness score" is competitive growth in fibroblasts**, not essentiality. Protein features
  predict it at R² = 0.453 and predict the five other screens at **−0.105 to +0.102** — so cofitness edges
  from different screens are not interchangeable.
- **hyperLOPIT assignment tracks protein abundance.** A missing compartment is not evidence of absence, so
  "unassigned" must be rendered as unknown, not as a 27th compartment.
- **Domain annotation partly records study effort**, not conserved architecture. Domain edges are the
  weakest type for that reason.
- **Never compare a gene set against the whole proteome** — 86% of adequately powered motif claims fail a
  matched background, and set median protein length predicts apparent enrichment at ρ ≈ +0.8.

## Start a session on this project by pasting this

```
Read /mnt/firecuda2/Claude/toxoplasma_projects/starplast/HANDOFF.md and continue starplast.
v0+v1+v1.1+v1.2 are done and pushed to github.com/EinarOlafsson/starplast (private); 48 tests pass
headless. Do not re-derive the design decisions in that file.
Next: <state what you want — e.g. "v2 species switching", or "rank the structural holes">.
```

Fill the `Next:` line in before sending — leaving the placeholder just costs a round trip.

## State at handoff — verified 2026-08-11 (v1.2)

**Built and working.** `identity.py`, `corpus.py`, `literature.py`, `build_graph.py`, `fetch_names.py`,
`interactions.py`, `app.py`, `tests/`, `pyproject.toml`, `README.md`. **48 tests pass headless** (`pytest tests/ -q`), covering
identity resolution, every precision guard, JATS parsing, the mentions table, the attention arithmetic,
the attention-depth tiering, and the app itself offscreen: 8,140 nodes, all 12 edge types, all 3 LOD levels, all 6 colour modes, picking,
search (`GRA16` → TGME49_208830), edge toggles, attention toggle.

**Numbers as built** (do not quote the older estimates):

| | |
|---|---|
| genes | 8,140 (node table, deduplicated) |
| compartments | 27 including `unassigned` |
| co-mention edges, abstracts | 435 (≥2 shared abstracts) |
| co-mention edges, full text | 7,733 (≥2 shared paragraphs) |
| structural holes | 255 pairs over 457 genes (6 with both endpoints studied) |
| measured binding | xlms 2,842 / ip_ms 64 / struct 11,684 |
| **unwritten interactions** | **2,673 (92% of all measured pairs); 147 both-studied** |
| crosslink models | 2,843 pairs, 2,439 with local CIFs, **only 162 trustworthy** |
| orthogroup / coexpression / compartment / cofitness / domain | 3,452 / 49,293 / 118,712 / 88,997 / 10,399 |
| genes named in 33,924 abstracts | 745 (9.2%) |
| genes named in 6,667 OA full texts | 2,546 (31.3%) |
| genes named in **either** | **2,566 (31.5%)** — was 601 |
| genes named in **neither** | 5,574 (68.5%) |
| distinct papers naming ≥1 gene | 4,579 |
| cache size | 2.4 MB app cache + 2.6 MB ToxoDB identity tables, all committed |

**Depth of attention — the number to quote, and the one that reframes the rest:**

| tier | meaning | genes |
|---|---|---|
| `focal` | named in a paper's **title** — the paper is about it | 286 |
| `substantive` | named in an **abstract** | 464 |
| `incidental` | named only in a **body or caption** | 1,816 |
| — | named nowhere | 5,574 |

**Coverage went 601 → 2,566. Attention went 601 → 750.** 71% of the coverage gain is `incidental` — genes
sitting in a screen's hit table, named once and never discussed. Quoting 2,566 as "genes the field has
studied" would repeat precisely the error the attention correction exists to prevent. The app labels such
genes **"Listed, not studied"**, and `attention_depth` is a colour mode.

**601 → 2,566 is better identity, not more literature.** The old code read only current `TGME49_`
accessions, so old `TGME49_0xxxxx` ids and the `TGGT1_`/`TGVEG_` strain ids papers use interchangeably
resolved to nothing and were dropped by an intersection with the node table — silently. Genes reached:
current accession 1,210, strain 1,150, symbol 1,127, Tg-alias 744, previous accession 369. Only the first
existed before.

**The caveats that must travel with 2,566:**

- It is still a **lower bound**. A gene with no symbol, never cited by accession, is invisible.
- **The two sources are different populations.** Abstracts cover the field; full texts are only what
  publishers deposited open access (6,667 papers against 33,924 abstracts). Full-text coverage answers a
  different question and must not be quoted as coverage of *Toxoplasma* research.
- **The full-text corpus grew under the build** — 5,493 → 6,667 files during this session — so every
  full-text number is a snapshot. `build_graph` logs the file count it actually saw; trust that over this
  table.

**Already tried, do not redo: co-mention restricted to focal papers.** `comention` is *already* that
layer. Rebuilding it over every readable paper's title and abstract, rather than the abstract corpus
alone, adds 22 papers and 20 edges, because the PubMed corpus already contains 1,109 of the 1,131
open-access papers that name a gene focally. Measured, not assumed. The corollary is worth keeping: the
1,816 `incidental` genes participate in `comention_ft` and nothing else, so that layer is the only edge
type connecting genes nobody has written about focally.

**Not built:** v2 (species switching, cross-species orthology edges) and v3 (continuous star-map zoom).
Deliberately deferred — v3 is most of the effort and least of the value. Note that v2 got cheaper: the
identity layer already carries GT1 and VEG accessions.

**GitHub:** `git@github.com:EinarOlafsson/starplast.git`, **private**, branch `main`, first commit
`61812bd`.

> **Hazard:** this working copy sits inside the Syncthing tree, so `.git` replicates to the work machine.
> `cellect/` already does this, so it is established practice here — but commit from **one machine at a
> time**, or Syncthing will fork objects in `.git` and produce `*.sync-conflict-*` inside it.

## The queue this app serves

`../new project ideas/PROPOSALS_FROM_LITERATURE.md` holds 20 proposals mined from the whole literature.
Eight are queued as analyses: B1 (exportome residue grammar), B2 (motif-audit survivors), B3
(cross-coccidian interface, 37 genomes), B4 (*Hammondia* — what *Toxoplasma* can do that it cannot),
B5 (strain-specific clinical outcomes), B6 (audit published effector claims), B7 (attention map),
B8 (cross-screen discordance). B9 dropped as overstated.

**The information map** — entities and relations extracted from the 33,924 abstracts, then inference over
the map's *structural holes* rather than over authors' stated gaps — is the governing piece of work and is
what this app displays. It was in progress when context ran out. B7 is one of its six hole types, so
building the map gets B7 nearly free.

**Known outstanding bug:** the *Cryptosporidium* VEuPathDB fetch silently downloaded nothing (CryptoDB
organism listing threw an HTTPError, was logged, and the job reported completion). Same failure class as
two earlier bugs — a job reporting success having done nothing. Needs a retry with the fixed
`abbrev()` in `.claude/skills/download_veupathdb.py`.
