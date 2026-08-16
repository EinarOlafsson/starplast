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
cd /mnt/firecuda2/Claude/repo/starplast        # the working copy on this machine
pip install -e .            # installs the console_scripts entry point
pip install -e ".[gpu]"     # optional: cuML and CuPy, for CUDA 12 -- see below
python -m starplast.fetch_names   # one-off: ToxoDB identity tables (needs network)
python -m starplast.build_graph   # one-off: rebuilds starplast/data/ (~5 min)
starplast                   # launch
pytest tests/ -q            # 2,285 tests, headless, no network, ~4 min
pytest tests/ -q -m slow    # the real build and the pdoc pass, ~2 min
```

If the GL widget fails on a headless machine, that is expected — this needs a display. The test suite is
headless and does not.

**GPU acceleration is optional and off by default.** `Preferences ▸ compute` reports what it found
and what that will buy. cuML moves UMAP and HDBSCAN themselves; CuPy or torch move the array work,
which measured about 1.5x on large distance matrices and *negative* on rank scaling -- that path was
deleted rather than shipped. Install with `pip install starplast-gpu`, or `pip install -e ".[gpu]"`
in a checkout; both resolve to the same wheels, since the metapackage's only dependency is the
extra. `starplast-install-gpu` picks the CUDA set from the driver and shows the command before it
runs anything -- **the newest set the driver could run is not the right answer**, because drivers
are backward compatible and cu12 wheels run on a CUDA 13 driver while cu13 wheels do not run on a
CUDA 12 one; it takes what torch in the environment was built against, else the lowest runnable set.
The switch defaults ON where a backend is importable, since installing two gigabytes of CUDA wheels
is a deliberate act, and every run records which library built it (`backend` column, and in the
saved recipe). **A map built by cuML's UMAP is a different map of the same data**, not the same map faster,
and the log says so whenever it happens.

**Two interpreters, and the second is the one that matters.** The suite runs in
`~/anaconda3/envs/spacr` (pandas 2.3.3); the user runs `~/anaconda3/envs/starplast` (pandas 3.0.5),
where the editable install points at this checkout. Three bugs have shipped that were invisible on
pandas 2 and total on pandas 3 — a read-only array written in place, and `.astype(str)` keeping NA
where it used to give `"nan"`. Re-check anything touching a dataframe under the second interpreter.

There is a second, **stale** checkout at `../toxoplasma_projects/starplast` from before the move.
Nothing reads it; do not commit into it.

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
additively, so 8,140 overlapping points summed to white and every color mode rendered as one blob --
while every array-level test passed, because they checked the color array and never the render. Points
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

**3e. Standalone means every measurement ships; coordinates are the one exception.** (Added v1.3.) The
cache is 30 MB and carries 358 columns for all 8,140 genes, and lives INSIDE the package
(`starplast/data/`) so a wheel carries it and `paths.py` resolves it with no configuration. An earlier `keep` allowlist silently shipped
3 of 18 RNA columns and 7 of 8 fitness screens; the build now ships every column that survives, with an
explicit drop list. Structures resolve on demand (`structures.py`) because 6,538 AlphaFold models plus
12,265 crosslink CIFs are gigabytes and git is the wrong place for them.

**Screen accessions must go through the identity layer.** Papers cite whatever id was current when they
were written. The 2019 in vivo screen uses pre-2012 `TGME49_0xxxxx` ids for *every* gene: unresolved it
contributed 0 of 8,140, resolved it contributes 168. Any future dataset join has this failure mode.

**Coverage is wildly uneven and "not measured" is not "no effect".** GRA17 is genome-wide (7,553); GRA12
is 236 targeted genes across two screens that are NOT replicates (in-vivo L2FC r = 0.41, kept separate);
the host-transcription screen is 252. There is **no deep proteome in this tree** — the only abundance is
two Pru IP experiments totalling 748 proteins, which is enrichment, not coverage. Do not present it as
proteome-wide.

**3f. The search is the point; the circularity guard is what makes it valid.** (Added v1.4;
substantially revised 2026-08-12; **the guard was found leaking three more ways on 2026-08-13 — read
the entry below this one before quoting any number in this one.**)
`search.py` walks dataset combinations x hyperparameters and scores each by recovery of a label
*excluded from the embedding*. With hyperLOPIT leaked in, the battery reports it separating clusters at
V = 0.96. Held out, the number is far smaller — and it took a negative control to find out how much
smaller.

**A negative control changed how every other number reads.** Scoring `attention_depth` — how much a
gene has been *studied* — first gave mean F1 0.654, above both measured targets. The per-label table
said why: essentially all of it came from one class, the never-named genes, at F1 0.77 over 1,439 of
them. The map was not ranking attention. It was separating genes that have been **measured** from genes
that have not, which is a property of missingness — a gene absent from most assays is absent from most
of the feature matrix. `compartment` was carried by the same artefact, its best label being
`unassigned`, and hyperLOPIT assignment tracks abundance so `unassigned` is largely "too scarce to
call".

So `score_recovery` now excludes labels meaning "not measured" — unassigned, unknown, empty, nan — and
all four targets were re-run under that rule, 328 runs each:

    stage_enriched_derived   0.709  →  0.709   positive control, unchanged
    cellcycle_phase          0.489  →  0.489   measured, unchanged
    attention_depth          0.654  →  0.339   NEGATIVE control
    compartment              0.484  →  0.228   measured

The two targets with no absence class did not move at all, which is the rule behaving as it should.
The negative control now sits below both measured targets and **passes**: the map is not substantially
organised by study effort among genes that have been studied.

**But localization (0.228) is now below the negative control (0.339).** Stated plainly: the map recovers
how much a gene has been studied better than it recovers where the protein is. That caveat must travel
with any localization claim. Cell cycle at 0.489 is the only target that clears the floor, and it does
so comfortably — which reverses the assumption the project started from. Location was the target that
had data; cell cycle is the one the map organises.

**Then all four were re-run on the full proteome, and the subsample turns out to have been flattering.**
Every number falls, the two measured targets furthest — 0.709→0.675, 0.489→**0.398**, 0.339→**0.322**,
0.228→**0.193**. The ranking survives, in the expected order, which is worth something. The margin does
not: cell cycle clears the negative control by 0.076 at full scale rather than 0.150.

**Read `n_labels_recovered`, not `mean_f1`.** At full scale the positive control scores 0.675 — the
highest of the four — from a two-cluster solution in which all three stage labels best-match the *same*
cluster, with per-label recalls of 1.000, 1.000, 1.000. Precision is then just prevalence: oocyst is 72%
of the labelled set, scores F1 0.838 alone, and carries the mean. Mean F1 over labels is inflated by a
trivial partition whenever one class dominates; `n_labels_recovered` is not. On that column:

    stage_enriched_derived   1 of 3    POSITIVE control — the majority class, nothing else
    cellcycle_phase          0 of 5    measured
    attention_depth          0 of 3    NEGATIVE control
    compartment              1 of 24   measured

So **no claim of the form "the map recovers X" is supported at full scale.** The subsample's cell-cycle
result does not survive 8,140 genes. See `results/full_proteome_2026_08_12/`, which ships the run script.

**Recovery is not prediction.** `search.predictions()` was run on the winning cell-cycle structure and
produced nothing: its purest cluster is 52% one phase against an 80% bar, and no cluster in the top ten
configurations exceeds that. F1 0.489 is carried by recall; there is not enough precision to place an
uncharacterised gene. See `results/predictions_2026_08_12/`.

Precision and recall are never blended: a cluster that is 100% apicoplast holding 5% of apicoplast
proteins is useless for inference.

**Every run stores its full recipe, seed and scores.** A hit nobody can rebuild is not a result.

**3f-2. The guard leaked three more ways, and every number in 3f predates the fix.** (Added
v0.18.0, 2026-08-13.) Found by running a real search in the real window and reading what won.

Recovering `compartment`, the winning configuration was built on the `localization` block:
`lopit_prob_map`, `lopit_prob_mcmc`, `lopit_methods_agree` — hyperLOPIT's own posteriors and its
methods-agreement flag. The label being recovered is hyperLOPIT's assignment. **On the shipped cache
those three columns recover it at mean F1 0.259, above interactions (0.192), protein features (0.171)
and every expression and fitness block.** Three columns beating eighteen RNA columns and eight CRISPR
screens is not the map finding biology.

Neither existing mechanism could see it. Measured association is 0.29 and 0.43 — far under the 0.8
threshold, because a posterior does not restate *which* compartment a protein is in. Declared
derivation sees nothing either: the label is not computed from the posterior. They are the same
experiment's **other outputs**, which is a third thing, and `excluded_for` now excludes them by
reading provenance from the registry.

The same question asked of the other targets found two more:

- **The negative control was embedding its own inputs.** `attention_depth` is `np.select` over
  `n_papers_focal / substantive / incidental`, and none of the three was declared or excluded (they
  score 0.25–0.43 against the tiering, because a count is not a restatement of a tier while
  determining it completely). A control that sees its own inputs scores too high, and every measured
  target then looks worse than it is by exactly that much — including the headline in 3f that
  localisation is recovered worse than study effort.
- **The same quantity, measured another way.** `ortholopit_label` (localization transferred from
  *P. falciparum* and *C. parvum* orthologs) sits at 0.72 against `compartment`; `lit_tier` at 0.71
  against `attention_depth`. Both just under the threshold, which is what a near-copy does.
  `search.SAME_QUANTITY` now groups columns by what they are an estimate OF, written out rather than
  pattern-matched.

So the guard is now three tests, and each was a bug before it was a rule: measured association for
the undeclared copy, declared derivation for the joint function, shared provenance and shared
quantity for everything the same work produced by another route. `validate.circularity_error` applies
the same three, because the Validation tab was refusing only the label column by name and would have
put a validated-looking number on a map built from that label's own experiment.

**Where the leak was reachable, and where it was not.** The published full-proteome battery
(`search.search` with `block_sets=None`) sweeps a hard-coded six-block base —
`expression_summary`, `expression_raw`, `fitness_screens`, `published_screens`, `protein_features`,
`interactions` — which contains neither `localization` nor `literature`. So no published row could
ever have used the leaking columns, and the re-run under the fixed rules confirms it for all four
targets: **0.675, 0.398, 0.322 and 0.193, unchanged to the last digit, winning block combination and
per-label best F1 included.** See
`results/full_proteome_2026_08_13_provenance/`.

The 3,000-gene subsample battery in `results/search_2026_08_12_corrected/` used the same six-block
base — checked, not assumed: 0 of its 328 rows name `localization`. **So no published number in this
project is affected by any of the three leaks.**

The interface is the other case. `AnalysisPanel.run_search` builds its combinations from **every**
block that has columns in the table, `localization` included, so a search run from the Search tab
could and did select it: on a 600-gene subsample the localization block won outright at mean F1
0.259 against 0.192 for the best measurement block. That is where this was found, and it is the
configuration a user would have believed.

**3g. A walk emits its configurations; a gene with no position in a map is not drawn in it.**
(Added v0.4.0, 2026-08-12.) `tuning.walk_umap_iter` yields one `WalkStep` per configuration as it is
computed — scores, coordinates, and a boolean mask over the node table — and `walk_umap` is that
collected and ranked. The table and the gallery therefore fill one row and one thumbnail at a time,
which is the difference between watching a 288-configuration sweep and waiting half an hour for it.
The subsample is taken in table order for the same reason the mask exists: drawn by `rng.choice` it
came back in random order, so row *i* of the embedding was an arbitrary gene and every attempt to say
which gene a point was named a different one.

The display half is a decision, not a detail. A walk map covers a subsample, so most of the proteome
has **no position in it**, and those genes were being left at the origin — 7,340 of 8,140 in a lump in
the middle of the map, drawn at low alpha, pickable, counted in the status bar and included in the
class filter. That is absence rendered as a measurement, in the application whose first rule is that
it must not be. `Window.placed` now records what the displayed embedding covers; unplaced genes are
drawn at alpha 0, excluded from `visible_mask` and unpickable. Hiding alone is not enough — an
invisible point is still the nearest point to a click where it sits.

**4. Level of detail is data-driven, not invented tiers.**
galaxy = coarse spatial structure of the embedding → solar system = orthogroup / co-expression module →
planet = gene → surface = that gene's evidence (papers, domains, screens, phenotypes).

The galaxy tier **was** hyperLOPIT compartment, and that was a mistake worth recording because it looked
like a rendering fault. Drawing one centroid per compartment produced a knot of dots in the middle of the
screen — correctly. A compartment's genes are scattered across the whole map, so their mean lands near the
middle of it: measured on the shipped embedding, the median compartment centroid sits **0.17** of the map
radius from the centre while the median compartment's own members spread **0.37** of the radius. The
centroids were an honest average of a quantity with no spatial meaning, drawn as though it had one.

This is the same fact the held-out search reports from the other direction — `compartment` is the
worst-recovered target in the map, below the study-effort control — which is why it could not have been
fixed by adjusting the drawing. A tier keyed to compartment could only ever have shown a central blob.

The tier is now computed from the embedding itself (`lod.py`): connected components over an occupancy
grid, no clustering library, deterministic, numbered largest-first so a landmark keeps its name. On the
shipped map that is 5 structures covering 8,123 of 8,140 genes, with centroids 0.40–0.97 of the radius
from the centre and spreads of 0.04–0.15 — each one now closer to its own members than to the middle,
which is the property the compartment version failed. The smallest is 76 genes, of which **25 carry a compartment
at all, and 24 of those are `PM - peripheral 1`** — 96% of the labelled quarter, not of the structure,
and the difference is the whole first rule of this project: 51 of those 76 genes have no localization,
which is not the same as being peripheral. Independently, `PM - peripheral 1` is the single
best-recovered compartment label in the full-proteome search (F1 0.500, measured before the
circularity fix in 3f-2). Two methods finding the same structure is the reason to believe it is there;
neither of them says what the other 51 genes are.

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
| `starplast/data/toxodb_identity.tsv`, `starplast/data/toxodb_strain_{gt1,veg}.tsv` | symbols, previous IDs, GT1/VEG accessions (committed; `fetch_names.py`) |
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

## Where the work list lives

`instructions/` tracks unfinished work across sessions: `START_HERE.md` orients a cold start,
`INDEX.md` is the status table, `open/` holds one file per unfinished task, and `done/` records what
was finished and how it was verified. Read `START_HERE.md` first — it carries the traps that cost
real time to find. Tasks 30–34 landed together in v0.31.0.

`skills/` holds reusable techniques worked out here, also installed under `.claude/skills/`.

## Start a session on this project by pasting this

```
Read /mnt/firecuda2/Claude/repo/starplast/instructions/START_HERE.md and HANDOFF.md, then
continue starplast. The working tree is v0.31.0; 2,285 tests pass headless with every module at
100% coverage. Check `git status` before assuming it has been published. Do not re-derive the design
decisions in that file.
Next: <state what you want — e.g. "v2 species switching", "search a new target", or
"curate hit lists from the 65 PDF-only interaction studies">.
```

Fill the `Next:` line in before sending — leaving the placeholder just costs a round trip.

## State of the application — verified 2026-08-14 (v0.31.0)

**2,285 tests pass headless** (`pytest tests/ -q`, ~4 min) and **every module is at 100% coverage**
(9,983 statements). No `pragma: no cover` anywhere: a Qt-thread body is covered by calling it
directly, and a branch that genuinely cannot run is deleted. Two functions were deleted in the last
pass on that rule, and writing one of the missing tests found a real defect in `objectives.adjusted`.

What the window does now, beyond the map itself:

* **A sweep is watched, not waited for.** `search.search` and `tuning.walk_umap_iter` emit one step
  per configuration — scores, coordinates, labels, and a mask of the genes it covers — so the table
  and the **gallery** (`gallery.py`, grid and scroll) fill one row and one thumbnail at a time. Click
  a thumbnail and it becomes the central map, where a gene is clickable like any other.
* **Clusterings are kept** (`runs.py`) with their recipe and their gene mask, named by timestamp and
  renameable, and they colour the map from the same one place every colour comes from.
* **Annotations** (`annotations.py`) cannot be saved without the validated precision, recall and
  enrichment for their category. That refusal is the feature.
* **Results save and load** (`results.py`): one table is a CSV whose first line names which table it
  is; every tab at once is a zip of those plus a manifest. A loaded row is as clickable as a computed
  one — clicking it rebuilds its map — and a file whose kind does not match the tab is refused.
* **A user's own table imports** (`importer.py`, `File ▸ Import data…`) from CSV/TSV/Excel/parquet,
  with the quantification guessed from the RANGE and the reason shown, every preprocessing choice
  offered rather than assumed, and the whole record kept. Columns are prefixed and joined in memory;
  the cache on disk is never written.
* **The compartment list has a parasite beside it** (`celldiagram.py`), drawn in outline on the
  panel's own ground, with exactly one compartment filled — the selected one, in its list colour —
  and clickable both ways.
* **Logging** (`logging_util.py`) is opt-in with per-level console control, and jobs
  (`jobs.py`) can be stopped, inspected and have their traceback copied.

The identifier spelling is American throughout as of 2026-08-12 (`color_of`, `normalize`,
`localization.py`); `embedding.RENAMED_BLOCKS` translates the block name inside recipes saved before
that, because a block name lives in every stored recipe and an unrecognised one rebuilds a different
map in silence.

**The data layer beneath it, verified 2026-08-11 (v1.2).** `identity.py`, `corpus.py`,
`literature.py`, `build_graph.py`, `fetch_names.py`, `interactions.py`: identity resolution, every
precision guard, JATS parsing, the mentions table, the attention arithmetic, the attention-depth
tiering, and the app offscreen over 8,140 nodes, all 12 edge types, all 3 LOD levels, every colour
mode, picking, search (`GRA16` → TGME49_208830), edge toggles, attention toggle.

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
| cache size | 17 MB inside the package (`starplast/data/`), including the ToxoDB identity tables, all committed |

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
genes **"Listed, not studied"**, and `attention_depth` is a color mode.

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

**The *Cryptosporidium* VEuPathDB fetch is fixed and has been run — 2026-08-13.** It had downloaded
nothing and reported completion: the CryptoDB organism listing 404'd because the webapp path is
`cryptodb`, not the `crypto` that dropping "db" from the project name gives, and the error was logged
while the job printed `DONE files=0`.

Two more silent skips turned up on the retry, both of the same shape — work that did nothing and
reported a number as though it had. `abbrev()` left punctuation in the species token, so
`Cryptosporidium sp. chipmunk LX-2015` became `Csp.chipmunkLX2015` and 404'd on all four of its
files; and the derived name is simply wrong for some organisms, because VEuPathDB keeps hyphens in
`CbaileyiTAMU-09Q1` and `CparvumIOWA-ATCC` and drops "genotype I" from `Cspchipmunk37763`. Deriving a
name and hoping is now replaced by **resolving it against the release index** (`release_dirs`,
`resolve`), matched on letters and digits alone; anything that still does not resolve is printed by
name, counted, and downgrades the run to PARTIAL rather than vanishing.

Result: 19 of 19 *Cryptosporidium* organisms resolved, 61 files, 354 MB. Thirteen carry the full set;
the other six publish a genome and nothing else, which was verified against the site rather than
assumed — those assemblies are unannotated upstream.

The same fix applies to every other skill that script serves, and **the other corpora have not been
re-run**: any organism whose name needed the index to resolve was silently skipped there too.

---

# Session 2026-08-13/14 — lighting, discovery, persistence, headless

Versions 0.24.0 → 0.28.0. Everything below is pushed to `main`. Read this section before touching
`discovery`, `metrics`, `optimize`, `interpret`, `searches`, `discover`, `rays` or `sprite`.

## 1. What shipped, in order

**Display (0.24.1 – 0.26.0).** Four reported faults, each a different cause, all measured before
and after:

* The grid was shaded through the point-cloud shader, whose normal is the direction out of the
  cloud's centre. The grid sits *below* that centre, so its inferred normal pointed down and every
  overhead light gave it exactly zero. It has a real normal now (up), and its **alpha** moves with
  the light as well as its colour — a colour-only change was below the threshold of visible.
* Corner and pointer lights were placed in **world** coordinates, so "top left" meant top left of
  the data and drifted the moment the map was orbited. They are placed in the camera frame now.
* `setMouseTracking(True)` was never called, so Qt delivered move events **only while a button was
  held**. The pointer light — the default source — moved only during a drag, which is also what
  orbits the camera. This was the whole of "I can't see the mouse light".
* Finishes were invisible because a Phong lobe of 48 is ~4° wide and lands on one or two points of
  a sparse scatter. Fixed by **contrast, not highlights**: matt sits on a raised ambient floor and
  the shiny finishes on a lowered one (`lighting.FINISHES[*]["ambient"]`, floor `MIN_AMBIENT`).
  Median point now moves 15–32 of 255 between finishes; before, 1–3.

**Sphere sprites (0.25.0).** pyqtgraph draws a scatter as point sprites and its texture is
`pData[:] = 255` with only the alpha shaped into a disc — every pixel inside a ball is the same
colour, so *no* per-point shading can make one look spherical. There is no fragment shader, but the
fixed-function stage MODULATES texture × vertex colour, so `sprite.texture()` builds a shaded unit
sphere per finish and uploads it (`sprite.ShadedScatter`, flushed at paint because a texture needs a
current GL context). "2D" is now a finish, not the only option.

**Rays (0.26.0).** `rays.py` — coarse occupancy grid, transmittance sampled along the segment,
first-hit marching for bounces. Four settings under **rays**: `none / shadows / emitter / bounce`.
Not ray tracing and says so in its own docstring: it is volume sampling of a density grid, the
trick a volume renderer uses for smoke, which is the honest model for points with no surfaces.

**Lighting simplification (0.32.0, supersedes the 0.24–0.26 UI).** The four old ray choices, corner/
orbit sources, light count/speed, and five look-alike finishes are no longer exposed. Preferences
now has `off / soft / ray traced`, exactly three interaction targets, three color moods, and
`flat / glossy 3D / metallic 3D`. The later explicit request for ray tracing is implemented as one
honestly named volumetric shadow-ray mode: one segment per gene through the density grid. It is not
Vulkan/path tracing and the UI says so. Bounces and emitter ornaments are gone. Fixed full-map GL
comparisons and llvmpipe frame timings live in `results/lighting_2026_08_14/`; regenerate them with
`scripts/benchmark_lighting.py`.

**GPU PBR renderer (0.33.0, supersedes the 0.32 point surface implementation).** Glossy and metallic
are now GLSL sphere impostors, not CPU-lit textures: per-fragment hemisphere normals, GGX BRDF,
procedural studio reflection, Fresnel, and curved `gl_FragDepth`. The 48³ density grid uploads once
as an `R32F` 3D texture and the vertex shader marches 24 shadow samples per gene/light. A static
ray-traced scene produced one framebuffer hash over 12 paints. Mouse light positions are low-pass
eased because the gene under a 2D cursor is discrete; this removes the shadow-field teleport at the
boundary between overlapping points. On llvmpipe, GPU rays take 7.7–7.8 ms for all 8,140 genes versus
15.3 ms for the NumPy fallback. The shader bridge covers pyqtgraph 0.13's fixed-function scatter and
0.14's VBO renderer; never restore a lookup of the private `pointSprite` shader name, which 0.14
removed. Evidence: `results/pbr_lighting_2026_08_14/`.

**Discovery stack (0.27.0)** — five new modules, all at 100% coverage:

| module | what it is |
|---|---|
| `discovery` | guilt-by-association and layer-disagreement findings, discrete and continuous, BH-corrected, circularity-marked |
| `metrics` | AUPRC / AUROC / lift, ARI, AMI, silhouette, kNN purity vs its own chance level, trustworthiness, permutation floor |
| `optimize` | hill climbing with restarts and caching over the joint embedding × clustering space; 7 objectives |
| `interpret` | ranks findings by strength × reach × novelty and writes them as claims with caveats |
| `searches` | a whole run saved and reloadable (0.28.0) |

Plus t-SNE and PCA (`embedding.METHODS`), kmeans and agglomerative (`clustering.ALGORITHMS`), and
HDBSCAN's `cluster_selection_epsilon` / `cluster_selection_method`. New **7 · Discover** tab.

**Persistence (0.28.0).** Every finished climb writes itself to `~/.cache/starplast/searches/<name>/`
as `manifest.json` + `configs.csv` + `findings.json` + `labels.npz`. Coordinates are NOT stored (a
recipe plus a seed rebuilds them); labels ARE (clustering libraries do not reproduce across
versions). The manifest fingerprints the node table by hashing gene ids **in order**; a mismatch
prints across the top of the reading rather than silently renaming every gene implicated.

**Headless (0.28.0+).** `starplast-discover`, `starplast/discover.py`. No Qt anywhere in the import
path — there is a test that asserts this in a fresh interpreter. Each `--task` is saved as it
finishes and skipped if already saved, so an interrupted batch resumes by re-running the same
command.

```bash
starplast-discover --task guilt:compartment_best --budget 100 --restarts 3
starplast-discover --task disagreement:compartment_best:fit_invitro_hff --exclude literature,localization
starplast-discover --list
starplast-discover --read bigA_00_guilt_compartment_best
```

## 2. Numbers worth carrying forward

* **kmeans beats HDBSCAN on this map.** First large run: kmeans at 40–60 clusters scored 158–266;
  HDBSCAN scored 143 while calling **52% of genes noise**. HDBSCAN at the old default
  (`min_cluster_size=25`) returns **2 clusters over 8,140 genes at 0% noise** — which is why every
  finding from it was an artefact.
* **Winning maps do not contain expression.** Both localisation climbs chose `interactions` +
  `protein_features` (+ `fitness_screens` for the 27-way target). Nobody has explained this yet.
* Best localisation run: 42 findings, mean AUPRC 0.276 at **5.1× prevalence**, kNN lift 1.9.
* `lopit_unified` run: 97 findings but **224 clusters** — see §4.

## 3. What was tried and changed on measurement, not taste

* **Fixed-count shadow sampling was wrong.** Twelve samples over a long ray steps *over* a thin
  occluder and reports a clear line through a wall. Sample count now follows ray length in cells.
* **Per-ray step, not global.** Absorption scaled by the longest ray's step charged a gene two cells
  from the light the same optical depth as one across the map: the pool around the pointer vanished
  (1.02× the rest). Per-ray: 1.53×, with 59% of the map genuinely occluded.
* **Lissajous as position, not direction.** Lights ranged 0.2–1.7× the intended radius and spent 4%
  of frames *inside* the cloud. On a shell now; brightness swing 21% → 9%.
* **Leave-one-out biases AUROC below 0.5** — measured 0.36 on four balanced classes in clusters of
  60. Do not read 0.5 as the floor; `metrics.chance_level` measures it by permutation, and it must
  shuffle *both* sides (scoring shuffled labels against true membership measures the wrong thing).
* **Interest score saturates.** A q of 1e-40 about two well-published genes outranked a q of 1e-8
  about twenty-five until `CERTAIN` capped it.
* **The three guards in `discovery`** (`MIN_LIFT`, `MIN_PURITY`, `MAX_SHARE`) exist because the
  first real run reported "2,085 unlabelled genes are cytosol" from a cluster covering half the
  proteome at 1.5× background, q = 1e-32. Significance is a statement about sample size.
* **`_split` standardised by within-cluster spread** made a cluster whose values were all −3.0±0.05
  look like a dramatic subdivision. It uses the layer's spread across the whole map now.

## 4. Known limitations — measured, not pending implementation

1. **Yield predominantly rewards fragmentation.** The crossed-factor diagnostic is now implemented.
   On the four saved winners, `compartment_best × cellcycle_phase` explains 0/40, 0/40, 0/31 and
   18/97 sibling clusters; only the final run clears the permutation null (18.6%, q=0.035). Fine
   clustering can resolve real conjunctions, and `discovery.conjunction` reports those explicitly,
   but it does not explain most of the 60–573-cluster solutions. Treat raw yield as inflated until it
   is normalized by reach rather than finding count.
2. **One cluster, several categories of one layer.** Those remain competing alternatives. Categories
   from different layers are now conjunction candidates instead, so the wording no longer blurs the
   two cases.
5. **The first large run was killed** (exit 137, user force-quit) after 2 of 10 tasks. That is why
   the per-task process split and the resume-by-skip exist. Not a bug, but the reason for the shape.
6. **Cluster labels are stored, coordinates are not.** If `rebuild` ever stops being deterministic,
   a reloaded search's map and its labels will disagree and nothing will notice.

## 5. Future scope, not an open numbered task

* **The malaria map remains separate.** The authoritative 24-slot Pf catalogue and its cached
  candidates now exist, and processed source files were acquired under
  `datasets/plasmodium_acquisition_2026_08_14/`. There is deliberately no Pf node table or
  combined-organism view: the Toxoplasma application must not attach malaria measurements to
  Toxoplasma genes through sparse orthology. Building that second map is a future product task.
* The ten-task headless run was stopped at 7 of 10 by request, to free the GPU for that work. Results land in
  `~/.cache/starplast/searches/bigA_*` … `bigE_*`; read them with `starplast-discover --read NAME`.

## 6. Conventions this session confirmed

* Verify a UI claim by **driving the real widgets** and measuring what reaches the renderer, in
  8-bit steps at the **median**. ~12/255 is roughly the visible floor. "Changed by >2/255" is noise
  and passed twice while nothing on screen changed.
* Every module stays at 100% coverage, no `pragma`. Tests assert **orderings** where possible
  (perfect > corrupted > shuffled), because those survive the numbers being tuned.
* Bump the version for feature work; no `Co-Authored-By` trailer.

## Standing goal, set 2026-08-15: finish every open instruction

Five are open. Suggested order, cheapest-unblocking-first:

1. **34 leftovers / coverage back to 100%** — it is 99% now (51 lines, worst `holdout_cv` at 94%).
   Do this first: everything below lands on top of it, and the project rule has held since 28.
2. **39 — species tables and host bridges.** Structural, and 40 and 41 both depend on it. Start with
   the `Toxo_` -> `Tg_` rename while it is still only a string change.
3. **40 — the slot tree window.** Small, and it is the thing that makes 41 auditable while it runs.
4. **41 — fill the slots.** The long one. Fix the Plasmodium candidate query FIRST (see the file);
   the current 16 off-target citations come from `malaria` being used as a standalone query term.
5. **37 (lighting material lab)** and **38 (pan-apicomplexan archive)** — 38 overlaps 41's fetching;
   read both before starting either, and fold 38's archive into 41's procedure rather than building
   two fetchers.

The atlas of all 239 slots, filled and empty, is published and regenerates from
`scripts/generate_slot_table.py` plus `starplast/data/slots.json`.
