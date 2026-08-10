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
python -m starplast.build_graph   # one-off: builds data/graph.npz (~2 min)
starplast                   # launch
```

If the GL widget fails on a headless machine, that is expected — this needs a display.

## Design decisions, and why (do not silently reverse these)

**1. Position = UMAP of a feature matrix, not force-directed layout.**
Force-directed position is aesthetic and arbitrary; two adjacent nodes mean nothing. UMAP over
(expression × 7 fitness screens × compartment × orthology breadth × domain content × disorder) gives
positions where proximity is interpretable. Precomputed and cached — never laid out live.

**2. Six edge types, toggleable, never merged silently.**
"The knowledge map" is not one graph. Each edge type answers a different question:

| edge | source | note |
|---|---|---|
| `comention` | 33,924 PubMed abstracts | **attention-biased — see decision 3** |
| `orthogroup` | OrthoMCL (16,794 groups) | clean; the cross-species bridge later |
| `coexpression` | GSE108740 stage series | quantitative |
| `compartment` | hyperLOPIT (26 compartments) | clean |
| `cofitness` | the 7 CRISPR screens | real and underused by the field |
| `domain` | InterPro shared domains | weakest; included for completeness |

**3. Attention correction is not optional decoration — it is a correctness feature.**
Raw co-mention edges reproduce the literature's popularity contest: GRA16 and ROP18 become bright hubs,
and the 78 genes named in ≤6 abstracts stay dark, regardless of biology. The app therefore ships an
**attention-corrected mode** that shows co-mention *residual from expected given each gene's publication
count*, and it is **on by default**. Turning it off shows raw co-mention and the UI says so. Without this
the app would be confidently misleading.

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
| `/mnt/wd4tb/skill_corpora/toxoplasma-scientist/*.xml` | ~6,800 open-access full texts (machine-local) |
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

## State at handoff

**Built:** `starplast/build_graph.py` (data prep), `starplast/app.py` (v0+v1 UI), `pyproject.toml`.

**Not built yet:** v2 (species switching, cross-species orthology edges) and v3 (continuous star-map zoom
polish). Both deliberately deferred — v3 is most of the effort and least of the value.

**GitHub:** repo `starplast`, private, on account `EinarOlafsson`.

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
