# Next session: start here

Updated **2026-10-08** at the **0.54.0** release. Read this page first: where things stand, the rules,
the code map, and **what is left to do, in priority order (section 4)**. `HANDOFF.md` holds the
long-standing design decisions, `instructions/` is the task ledger (`instructions/INDEX.md`), and
`CHANGELOG.md` lists every change by release.

**Current direction (user, 2026-10-07):** explore published datasets for selected organisms and
their hosts; start from a gene, class, localization or label and inspect every applicable strategy,
precomputed results, agreements/conflicts and a ground-truth-tested meta-inference. Scorecards
must make accuracy, capacity and the underlying test evidence accessible at every level.
The controlling plan is **[instruction 64](instructions/open/64_information_space_and_inference_atlas.md)**:
40 bounded action cards with dependencies and acceptance tests. **13/40 complete (32.5%)**:
64.01 provides the entity/query schema and provenance-preserving exact alias resolver
(`starplast/query.py`, 472 focused checks passed). 64.02 provides the offline inventory
(`starplast/inventory.py`, 525 focused checks passed): all 162 sources reconcile to 180 rows
including 16 scoped refusals. Evidence and gaps are in
`results/information_inventory_2026_10_07/README.md`. 64.03 adds typed provenance for 628 addresses, a measured mapping-loss audit and exact
reproduction of all 15,437 installed mouse macrophage TPM values. 64.04 records 529 targets,
646 benchmark candidates and explicit validation gaps for all 39 strategies × four species.
Zero independent biological benchmarks admitted; 440 output grades remain unresolved.
64.05 freezes nested roles and training-only source exclusions; two real candidate
cohorts and all ten deliberate leakage refusals verified (532 checks passed, two existing skips).
64.06 provides shared training-only baselines and synthetic/null controls (426 checks passed).
64.07 declares all 39 strategy capabilities, 40 component-test roles and 624 synthetic query routes (855 final checks passed). 64.08 adds immutable artifact roles and cache invalidation (840 checks passed). 64.09 reconciles 1,038,372 rows in 4,112 cohorts, preserving both accuracy denominators (835 checks passed). Continue at **64.10**, categorical frozen row adapters, and remaining host expression/gene-space work under 65.04/64.21/64.22. See
`results/ground_truth_registry_2026_10_07_v2/README.md` and the provenance artifact.
**Reporting:** whenever an action completes, replace its percentage with a green tick and repost
the entire progress table (40 instruction-64 rows, four instruction-65 rows, four instruction-66 rows, watchdog 67.01 and two instruction-68 rows: 51 total).
Keep the stored tracker and completion evidence current.

**New user priority (2026-10-07):** [instruction 65](instructions/open/65_dataset_selection_audit.md)
audits all dataset choices against literature alternatives. For comparable admitted datasets,
prefer citations per year since first publication, with modest recency and comprehensiveness
bonuses. Use `starplast/dataset_selection.py`, retain the factors and provider/snapshot, and
do not guess missing metadata or promote unsuitable popular papers. Host storage is present
(7 Hs/3 Mm sources, 20,989/15,590 protein rows; 12/50 host slot views filled), while independent
host gene spaces/inference packs remain pending under 64.21/64.22/64.28/64.29. Instruction 65
adds four bounded audit/promotion rows to the progress report.

**Instruction 65.01–65.03 complete:** all 162 sources accounted for; 95 recorded
PMID/DOI identities resolved; all 297 slots searched, 2,032 distinct publications
retained. Six diagnostic quantity/context traps triaged, no sources promoted.
**576 checks passed; one optional pdoc documentation module skipped locally.**
Read `results/dataset_selection_2026_10_07/README.md`; executed notebooks and
checksummed primary responses retain the evidence. **65.04 remains open**:
origin-paper/deposit lineage, 231 truncated searches, measured assay coverage and
validated replacements. First pages by total citations/date can miss the highest
annual citation-rate paper. Prioritize host rhoptry journal-version equivalence,
2017 quantitative RBC complement, GTEx v10 lineage, then the all-source review queue.
Do not mistake metadata resolution or pending abstract hits for admitted biology.

**Source:** `/media/carruthers/mnt3/claude/repo/starplast` (the package is `starplast/` inside it),
branch `nightly`. `main` is the release branch at 0.54.0. Nightly also contains the October 7
literature-verifier audit, the instruction-64 action plan, the entity/query contract and inventory;
runtime inference algorithms, parasite tables/graphs and shipped claims retain their released baseline.
Nightly corrects host ambiguity projections and numeric error metrics; frozen historical
calibration/benchmark records retain their original input/code lineage.

---

**Active source recovery:** [instruction 66](instructions/open/66_source_recovery_and_host_candidates.md).
Use the archive root `/media/carruthers/mnt3/claude/toxoplasma_projects/datasets` itself.
GTEx v10 has a verified relocated-file receipt. `scripts/recover_pmc_sources.py` uses the
current public PMC cloud API: 12 named missing inputs recovered, 57 processed source
bindings available in `results/source_recovery_pmc_2026_10_07_v2/`. Ten attempted sources
remain unresolved, including one declared PMID/PMCID mismatch. Preserve these refusals.
The 2017 RBC processed MaxQuant archive and two DOI-matched publisher supplements are
now external with verified checksums; RBC source review (66.02) is complete: retain the distinct installed fraction/surface
evidence, keep absolute abundance as a separate candidate under 65.04. Rhoptry
version comparison (66.03) is complete: all 20,010 genes × three fields exactly
match the journal. Canonical citation updated; raw v2 provenance and measurements
retained. The audit found 39 ambiguous symbol-to-protein projections assigned by
the old first-accession policy: 66.04 now withholds these projections with preserved original gene evidence.
An additional 43 processed supplementary inputs were recovered for unmatched
legacy derivations; they require transform/table association review, not renaming. Original
source data stays external; no runtime dataset replacement is implied by acquisition.
The proposed worker-stagger/retry feature was cancelled as a wrong-session request.

## 1. Where things stand

| | |
|---|---|
| Version | **0.54.0** on PyPI (2026-10-05), verified: wheel + sdist on PyPI, GitHub release with both |
| Branches | Runtime/data baseline `1249c55`; nightly adds the handoff, October 7 literature audit, instruction-64 action plan, query contract and evidence inventory |
| Tables | *T. gondii* `nodes.parquet` 8,140 genes; *P. falciparum* `pf_nodes.parquet` 5,720; host proteins 36,579 (Hs 20,989 + Mm 15,590) |
| Strategies | 39 in 9 families, each with method, techniques, scorecard task, calibration grade, card, explainer and worked examples |
| Calibration | 7,640 self-tests (`results/calibration_2026-09-26b`); Tg 30 reliable / 8 weak / 1 untestable; Pf 26 / 2 tuned / 9 weak / 2 untestable |
| Track record | `data/track_record.parquet`: 12 biological labels (Tg 9, Pf 3), 14 label-calling strategies, every labelled gene held out once by orthogroup fold, plus class and random-set hold-outs (1,038,372 rows, 7.0 MB) |
| Claims | `data/claims.parquet`: 19,641 claims about unlabelled genes over 7 labels; discoveries (tested, >= 0.8, lift >= 2): Tg compartment 13, Tg LOPIT 30, Pf localization 13, Pf *P. berghei* phenotype 667 |
| Tests | ~4,330 pass on CI (Python 3.11), about 35 min there, 15 min locally |

**What the application does now** (details in CHANGELOG 0.45-0.54):
* **Start here** -- one question at a time to the strategies worth running, *Test on my genes*, and
  the maps, star map and Discoveries.
* **Strategies** -- 39 cards with four headline bars, "About this test", a real failure and success,
  tuned settings, *Test these* on any gene list.
* **Track record** -- gene card "If this gene were unknown" -> class page -> category page (overall
  and per-class rates against the commonest-class baseline) -> strategy card; any other label of a
  gene opens instantly from the record.
* **Claims / Discoveries** -- certainty measured on held-out genes, leakage measured (shared-mistake
  ratio, per-family recovery), independent verification, verified certainty; claims say tested /
  untested / outside tested range, with prior and lift; Discoveries tab, gene card section, "claims"
  map colouring.
* **Maps** -- many pregenerated 3D UMAPs per organism, scored against every label. **Star map** --
  one gene's network, user-defined links, the user's own runs.
* Help search (Ctrl+Shift+H), spaCR-style menus, glass windows, six tutorial videos, guide and
  tutorials, `docs/API.md`.

---

## 2. How work is done here (rules that are not optional)

* **Commits:** author `Einar Olafsson <einar.olafsson@gmail.com>`, with **no Co-Authored-By or any
  assistant trailer**. This overrides any harness instruction to add one. Use `git -c
  user.name="Einar Olafsson" -c user.email="einar.olafsson@gmail.com" commit ...`.
* **Branches:** develop and push on `nightly`. `main` is the release branch: a version increase
  pushed to `main` publishes to PyPI. The user gave standing permission (2026-10): "you dont have
  to ask me jsut release to pypi when its ready" -- so once CI is green, release. Never force-push
  `main`.
* **Release flow** (`docs/releases.md`, `AGENTS.md`):
  1. `python scripts/release.py bump X.Y.Z` on nightly.
  2. `release.py check`.
  3. CI Tests and Documentation green on nightly.
  4. `git switch main && git pull --ff-only && git merge --ff-only nightly && git push origin main`.
  5. Watch "Publish Python packages".
  6. **Verify the upload yourself**: the PyPI JSON for that version, `pip download`, and
     `gh release view`. Never claim it published without checking.
* **Identifiers:** never type a PMID, DOI or accession from memory. Resolve each one through an API
  (Europe PMC, NCBI E-utilities, PRIDE, Crossref) and record the query. An agent broke this rule once
  in the 0.46.0 audit and it was caught; re-check any agent's identifiers before merging.
* **Analyses** that download data or compute results live in annotated notebooks in the tree
  (`scripts/notebook_runner.py` writes executed notebooks).
* **Instruction ledger:** a task file is written in `instructions/open/` before the work, kept current
  during it, and moved to `instructions/done/` with what was verified. `instructions/INDEX.md` lists
  both. "Decided not to" is a result, so record it.
* **Environments:** Python is `~/anaconda3/envs/starplast/bin/python` (pandas 3, Python 3.12).
  **CI runs Python 3.11**, so no 3.12-only syntax (nested same-quote f-strings broke CI once). Do
  not modify the `deeptmhmm`/`af3`/`pymol` envs, and ask before creating a new env.
* **VEuPathDB** (ToxoDB/PlasmoDB/CryptoDB/VectorBase): since release 71 the public downloads return
  404 and the service returns 401. No key exists, and the user's rule is not to work around that.
  Use Ensembl Genomes r63, NCBI and UniProt for gene builds.

### RAM (125 GB machine, shared with other sessions)

* `ramguard.service` (user unit) enforces three thresholds. The code and hook are in
  `/media/carruthers/mnt3/claude/tools/ramguard/`.

  | At | It does |
  |---|---|
  | 80 GB | warns |
  | 100 GB | freezes the largest job tree, and the PreToolUse hook denies heavy commands |
  | 110 GB | kills |

  Don't build another guard.
* Run heavy jobs under a cap, e.g. `systemd-run --user --scope -p MemoryMax=40G ...`.
* Keep long jobs at about 80 GB total or less; systemd-oomd kills VS Code scopes under pressure.
* **Shared RAM ledger:** `/media/carruthers/mnt3/claude/ram_leases.md`. Before a heavy job, add a
  line (`session | job | peak GB | start | expected end`) and check current use plus leases stays
  under ~80 GB; remove it when done. The autotrade session (repo-02) uses it too and will pause its
  queue on request. The RAM hook blocks any command containing "python" at >= 100 GB.

### Tests

```bash
cd /media/carruthers/mnt3/claude/repo/starplast
QT_QPA_PLATFORM=offscreen ~/anaconda3/envs/starplast/bin/python -m pytest -q -p no:cacheprovider
```

* Under systemd the PATH lacks `~/anaconda3/bin`; add it, or 3 proteomics tests fail on `bsdtar`.
* `tests/conftest.py` creates **one QApplication for the whole session** and closes each module's
  windows. If the application is recreated, Qt's style globals outlive it, and the glass style made
  `setStyleSheet` segfault. Do not reintroduce per-module applications.

### Calibration sweep (after any data, slot or strategy change)

```bash
SNAP=/media/carruthers/mnt3/claude/starplast_calib_snapshot
rsync -a --delete starplast scripts pyproject.toml $SNAP/          # FREEZE the code: workers re-import
cd $SNAP && python scripts/calibrate_strategies.py --count          # size it
systemd-run --user --unit=starplast-calib-X --collect -p MemoryMax=62G --working-directory=$SNAP \
  -E PATH=$HOME/anaconda3/bin:/usr/bin:/bin ~/anaconda3/envs/starplast/bin/python \
  scripts/calibrate_strategies.py --out <repo>/results/calibration_<new-dir>   # NEW dir: an old one resumes
cd <repo> && python scripts/calibrate_strategies.py --publish --out results/calibration_<new-dir>
```

* It takes about 1.5–4 hours, depending on load, and is memory-scheduled.
* `--publish` writes `starplast/data/strategy_calibration.json`, the README calibration and scorecard
  blocks, `docs/calibration.md` and `docs/scorecards.md`.
* Then rebuild the tutorials: `QT_QPA_PLATFORM=offscreen python scripts/build_tutorials.py`
  (about 10 minutes).
* Re-measure the shipped verdicts with `scripts/strategy_selftests.py [--organism Pf] [--only k1,k2]`.

### Adding data

1. Put the raw files in the dataset root (`STARPLAST_DATA=/media/carruthers/mnt3/claude/toxoplasma_projects/datasets`),
   in `<level>/<kind>/<PMID|accession>/` with `URLS.txt` and `SHA256SUMS.txt`.
2. Write a loader and a `Deposit` in `starplast/deposits.py`, and add a `Dataset` in
   `starplast/datasets.py`.
3. Add a section to `scripts/derive_deposits.py`; it re-executes the notebook, which should reproduce
   the paper's headline result.
4. Merge with `scripts/add_deposits.py`; it refuses to lose values.
5. Place the columns in slots with `scripts/generate_slot_table.py`.
6. Add the leakage family in `starplast/search.py` (`SAME_QUANTITY`).
7. Run `scripts/leakage_audit.py` and require 0 gaps and 0 residual leaks.
8. Rebuild the layouts with `scripts/rebuild_layouts.py`.
9. Re-run the calibration.

---

## 3. Map of the code

| Module | What it is |
|---|---|
| `starplast/claims.py` | **Claims** (0.54): `shared_mistakes`, `family_recovery`, `certainty_model`, `verification`, `verified_model`, `candidates`, `recipe`, `generate`, `applicability`; shipped readers `shipped`, `recipes`, `gene_html`, `claim_html`, `discoveries` |
| `starplast/discoveries_panel.py` | The Discoveries tab |
| `starplast/track_record.py` | **Track record** (0.52-0.53): `evaluate`, `evaluate_sets`, `my_list`, `alone`, `summary`, `labels`, `default_target`, `gene_html`/`class_html`/`target_html`/`list_html`, `record_phrase`, `beats_baseline`, `_baselines` |
| `scripts/build_track_record.py`, `scripts/build_claims.py` | Rebuild the shipped record (`--workers 4`, ~1 h, ~1 GB/worker) and then the claims (minutes). Claims read the record: rebuild in that order. |
| `scripts/audit_claim_verifiers.py` | Audit corrected literature verifiers against the shipped recipes without promoting claims; executed measurements and decisions in `results/claim_verifiers_2026_10_07/` |
| `scripts/notebook_track_record.py`, `scripts/notebook_claims.py` | Write the executed notebooks `notebooks/track_record_2026_10_03.ipynb`, `notebooks/claims_2026_10_04.ipynb` |
| `starplast/strategies.py` | Strategy/TestResult/Context framework, bans (`banned`, `banned_layers`), `knn_vote`, `propagate` |
| `starplast/strategy_catalog.py`, `strategy_graph.py`, `strategy_learning.py` | Strategies 01-32, 33-34, 35-39 |
| `starplast/scorecard.py`, `techniques.py`, `calibration.py` | 6 tasks / 51 metrics, 40 techniques, grades and intervals |
| `starplast/strategy_panel.py`, `strategy_card.py`, `strategy_explainers.py` | The Strategies tab and its cards |
| `starplast/guided.py`, `guided_panel.py` | Start here (question tree without Qt; the panel) |
| `starplast/umap_gallery*.py`, `star_map.py`, `star_edges.py` | Maps tab, star map and its link store |
| `starplast/organisms.py` | The organism registry (one `Space` per species) |
| `starplast/search.py` | Leakage closure: `excluded_for`, `excluded_layers`, `LAYER_SOURCES` |
| `starplast/app.py` | The window; evidence-panel links route through `_detail_link` (`starplast://gene|class|target|record|alone|claim/...`) |
| `starplast/help_index.py`, `help_search.py`, `glass.py`, `theme.py`, `ambient.py` | Help search, glass, menus, background (`keep_behind` fixes the blank panel) |

---

## 4. Current execution plan

Follow [instruction 64](instructions/open/64_information_space_and_inference_atlas.md) and its
linked action cards. It covers both clarified goals: evidence-space exploration and a gene/class/label
inference atlas, backed by appropriate ground-truth tests, reusable precomputation, dependence-aware
meta-analysis and scorecards at every level. 64.01–64.09 are complete; continue remaining adapters/source admission and
respect each card's dependencies. Completed cards record fixtures, validation and limitations.
**Latest user priority (2026-10-08): function across Discoveries.** Browser coverage
68.01 is complete; **68.02 is 25%** with pinned domain names, a frozen EC recovery
pilot and actual functional scorecards. Continue independent activity/source admission,
Plasmodium/domain targets, further adapters and calibrated deployment using the
existing 64 contracts. Localization
was the first developed prediction path; the old menu/default filters hid broader
annotation and untested-label coverage. Preserve all existing host/audit work.
The census leaves missingness causes unknown, identifies 16 rows lacking slot context and
retains unattributed host bridge records; do not turn these gaps into negative measurements.

The user's standing preferences still apply: **quality before quantity**, compact information with
click-through detail, and measured leakage. The full progress table lives in instruction 64 and is
reposted after each item completes, using a green tick in place of its percentage.

### Supporting backlog from earlier instructions

The list below preserves previous work and scientific gaps. Its numbering is historical; the
implementation order and bounded acceptance conditions now come from instruction 64.

### A. Claims (instruction 63)
1. **Wider independent coverage.** Only 6-22% of unlabelled Tg genes are reached by any independent
   check (physical partners, shared fold). Find or build more verifiers and measure each one's
   shared-mistake ratio; next candidates: orthology transfer from Pf (cross-species), the Lourido
   density screen, newer LOPIT. **Literature audit completed October 7:** no Tg candidate passes
   both independence and calibration; one Pf phenotype candidate passes the pooled screen but
   adds only two in-range claims, with one agreement and one disagreement. No promotion. Read
   `results/claim_verifiers_2026_10_07/README.md`; do not assume excluded literature is independent.
2. **Step 3: frozen claims tested prospectively.** Freeze claims with version and hash; when a dataset
   is added, score the frozen claims against it automatically and show the result. The strongest test.
3. **Genes outside the tested range** (42% of unlabelled Tg genes on compartment, mostly genes LOPIT
   could not detect). Measure whether certainty holds for them in any way (e.g. other-stage proteomes),
   or keep saying "not measured".
4. **Combine claims across labels** where consistent (compartment vs LOPIT unified agree for the same
   gene) and measure whether agreement across labels raises precision.
5. Labels with no claims: Tg cell-cycle phase has no independent verifier (229 calibrated, untested
   claims); the binary screens are uncalibrated. Measure whether per-class calibration rescues them.

### B. Track record (instruction 62)
6. Check stage 2's acceptance: the pooled record against the scorecard numbers per strategy.
7. Show per-class recall on the class page (the category page has it).

### C. Usefulness (instruction 60, items still open)
8. Methods paragraph + resolved citation list + BibTeX for any run (item 2).
9. Result provenance (organism, version, table fingerprint, seed, columns) and session save/load (item 3).
10. Figure export (PDF/SVG) from every panel (item 4).
11. A gene digest: coverage, extreme measurements as percentiles, which strategies named it (item 5).
12. Per-target calibration in the UI, and weight Start here's ranking by it (item 6).
13. Missingness / dominance / noise warning on every result (item 7).
14. **"Drop in my own screen" end to end** -- a real bug: imported columns do not reach the
    strategies' contexts (item 8).

### D. The 100 biological questions (instruction 59)
15. Q01-Q50 were executed on branch `worktree-agent-a5f1423763f196a27` (b0d8932, WIP, unreviewed; the
    agent stopped at a spend limit). Review, merge, write `docs/questions.md`, and link questions from
    Start here.

### E. Organism spaces (instruction 53; the larger plan)
16. Build the Hs and Mm spaces from `<STARPLAST_DATA>/spaces/Hs`, `/Mm`; publish organism packs
    (`docs/space-packs.md`); a Space menu (parasites / hosts / vectors).
17. Cryptosporidium (Cp) and *P. berghei* (Pb) spaces; *Anopheles* (Ag, As), rat (Rn), cat (Fc) --
    mind the id-space traps in 53.
18. Leakage families and calibration per space; then cross-space bridges (orthology, host-pathogen PPI,
    vector-parasite) and pathogen-host-vector inference.

### F. Data (instructions 38, 41)
19. Fill remaining slots (41): fetch, verify, key and load the discovered candidates. 38's archive
    machinery is built.
20. The next Tg and host data audit.

### G. Housekeeping
21. A fresh calibration sweep after the 0.50-0.54 data and strategy changes (section 2 has the recipe).
22. Lower the organism-literal ratchet as files move to the registry.
23. Retry active learning with calibrated (conformal) uncertainty; plain uncertainty sampling was no
    better than random (`done/54`).
24. Twelve `worktree-agent-*` branches remain locally; all their work is merged except
    `worktree-agent-a5f1423763f196a27` (item 15). Delete the rest after a final check.

### For the user (not for an agent)
* Look at the glass windows on the real desktop (drag, resize, compositor were verified offscreen only).
* Delete `/media/carruthers/mnt3/claude/starplast_video_scratch/` (permission was denied to the agent)
  and the junk `<STARPLAST_DATA>/spaces/As/abundance/PXD001647/` (Hydra/Daphnia files).

---


## 5. Traps already met (read before touching these areas)

* **Leakage is measured, never assumed** (user rule, 2026-10-04). A signal peptide is location
  information with no shared dataset. `claims.shared_mistakes` and `claims.family_recovery` are the
  measurements; a dataset-overlap argument is a first screen only.
* **The track record's runner must apply the same bans as the strategies.** It once walked label
  diffusion's default `layer=coexpression` for the derived stage label and "recovered" it at 91-100%;
  `track_record._permitted` now swaps a banned layer. `_derived` labels are not recorded at all.
* **Categorical columns:** groupby without `observed=True` builds the full cross product (a gene card
  took 1.5 s); `value_counts` lists zero-count categories (class pages once named confusions that
  never happened). Convert to str or pass `observed=True`.
* **Base rates pose as knowledge.** Uncalibrated recipes make no claims, and every claim carries its
  class prior and lift: the binary screens would otherwise offer thousands of "confident" claims that
  were the 86-89% base rate of "no phenotype".
* **A strategy that absorbs every evidence type (stacking) cannot be independently tested**; recipes
  pick the generator that independent evidence CAN test (`claims.candidates`).
* **Window width:** a new dock's control rows set the window minimum. Put rows in a horizontal scroll
  strip (`discoveries_panel._strip`, `star_map.py`) or `test_display` fails at 800 px.
* **Every public function needs a docstring and every colour mode a menu explanation** -- CI tests
  enforce both, and new files must not add "Tg"/"Pf" literals (use `organisms.TOXOPLASMA`).

* **pandas deep-copies `Series.attrs`** on every operation, so per-class scores ride in a holder whose
  `__deepcopy__` returns itself (`strategies._ClassScores`).
* **Conformal per-class thresholds** are infinite for a class with fewer than ceil(1/α)−1 calibration
  genes, and such classes fall back to the overall threshold. Coverage alone is no test, which is why
  the verdict is set efficiency.
* **Strategy families must hold consecutive numbers**, and a test enforces it. A new strategy that
  does not extend the last family needs a new family.
* **README prose is capped at 200 lines**, not counting the generated dataset, calibration and
  scorecard blocks.
* **ToxoDB PTM `modification_site` is the peptide start**, not the residue.
* **Cellpose**: always use the full path to the `cpsam_v2` weights. Bare names silently load the wrong
  checkpoint. This applies to spaCR work too.
* Things found by the acquisition agent:
  * HPA is CC BY 4.0 now.
  * The DepMap portal is behind Cloudflare; 24Q4 came from figshare+.
  * The cat has three id spaces, *A. gambiae* NCBI moved to AGAMI1_ ids, and rat STRING uses old
    ENSRNOP ids.
  * gene2pubmed barely covers the parasites.
* **Parallel agents** worked well in isolated git worktrees with non-overlapping files. Merge
  conflicts were only ever in `CHANGELOG.md`, where you keep both sides. Re-verify an agent's
  identifiers and look at its screenshots before merging.

---

## 6. Records worth reading

| Record | Covers |
|---|---|
| `instructions/done/50_data_audit_2026_09.md` | The first data audit, and a citation that was wrong |
| `instructions/done/51_strategy_calibration.md` | How calibration works |
| `instructions/done/52_data_audit_2026_09_b.md` | The second audit: six Pf datasets, five refusals, the accession error |
| `instructions/done/54_scorecards_and_advanced_models.md` | Scorecards, strategies 35–39, active learning not shipped |
| `instructions/open/53_organism_spaces.md` | The organism-spaces plan |
| `instructions/open/62_holdout_track_record.md` | Track record design, stages, status, traps |
| `instructions/open/63_claims_generate_and_verify.md` | Claims design and every measurement behind it |
| `notebooks/track_record_2026_10_03.ipynb`, `notebooks/claims_2026_10_04.ipynb` | The executed record of both |
| `docs/scorecards.md`, `docs/calibration.md`, `docs/strategies.md` | Every number, generated |
| `docs/tutorial/index.html` | Six tutorials and the complete guide |

## October 8 source recovery completion

66.01 is complete: all 162 sources accounted for, 69 processed bindings, 11 container manifests and 538 candidate files, with explicit association, transform and access gaps. Read `results/source_recovery_final_2026_10_07_v2/README.md`. CSPA exact installed-field reproduction verified; its non-detection semantics corrected without numerical change. Frozen benchmark readers support explicit checksum-verified relocation. **717 relevant checks passed.** Instruction 66 is 4/4 complete; continue 64.10/64.11 frozen benchmark adapters. Numerical admission/replacements remain 65.04. Earlier code identities remain archived with their immutable artifacts.

64.07 complete: `starplast/capabilities.py` and `results/capability_contracts_2026_10_08_v2/`; distinct pair/class/evidence/gene outputs and explicit unavailable host adapters. No new biological performance claims. Next 64.08 artifact roles/identity/invalidation, then 64.09 shared aggregation.

64.08 complete: `starplast/artifacts.py`, `results/artifact_contracts_2026_10_08_v2/`; six task types and five storage roles, ten exact synthetic round trips. Frozen split guards reject misassigned fit/test entities and known benchmark genes as unknown deployment. No estimator fitted. Next 64.09 shared row-level aggregation/reconciliation.

64.09 complete: `starplast/record_scorecards.py`, `results/record_reconciliation_2026_10_08_v3/`; all legacy outcomes match standard metrics on identical cohorts. Original truth/exclusion/nested-fit lineage remains unresolved, not retroactively admitted. Settings/seeds/modes/named sets remain separate; no independent-sample inflation or component/per-gene probability claims.

64.10 is 35%: first `starplast/label_records.py` feature-kNN adapter (`8b59643`), with native
training rank parity and frozen train distributions for held-out values. The
prediction-grade compartment pilot retains all 560 test genes and full native
class scores; `results/label_knn_pilot_2026_10_08_v5/` is the rank-transform run.
Earlier z-score/serialization diagnostics retain original code. Scorecard
identities now preserve one-ULP differences. Continue remaining categorical
adapters/full outer coverage and independent truth acquisition; host-symbol
correction is complete under 66.04. No runtime strategy, installed values or calibration
has changed. The missing categorical ortholog-transfer adapter is also verified:
`starplast/transfer_records.py`, `results/ortholog_transfer_pilot_2026_10_08/`.
Same 560 test genes: 123 calls/437 abstentions, 65 correct/58 wrong; all-hidden
surrogate agreement 0.116071 versus called-only 0.528455. Native method has no
class-score/support output; no scores invented. Both target/source are stored
spatial predictions with context/independence gaps. 464 checks passed and exact
projection/call/card/baseline replay verified (commit `445e5f1`).
`starplast/conformal_records.py` adds kNN conformal sets with verified base model
and separate calibration labels, exact native parity on frozen roles. Shared
cards retain set coverage/size/singleton/empty share/native efficiency and class
metrics. `results/conformal_label_pilot_2026_10_08_v2/`: empirical prediction-grade
coverage 0.973214, mean set size 24.455357/26, efficiency 0.061786, zero singleton
calls. All-training-classes control has coverage 1, size 26, efficiency 0. 908
checks passed; optional pdoc skipped. No biological/exchangeability guarantee.
Continue missing `holdout_search`/`multiplex_modules`, logistic conformal variant,
full outer coverage and independent truth. The former two need explicit label
benchmark-task declarations alongside their native cluster tests; graph/map
transductive construction and learned representation fitting must follow 64.05.

GT-SPATIAL-01 truth review: `results/spatial_truth_source_review_2026_10_08/`
retrieves/verifies seven primary spatial companions, complete marker fields and
Plasmodium Figure 3 outcomes. Toxoplasma 718 final markers include its 62 new IFA
results; those cannot independently test the final classifier. Nine selected
Plasmodium IFA outcomes map exactly and are outside both marker fields; six
other attempted targets remain unknown. Coarse taxonomy, assay/context, source
dependence and target-selection gaps remain. 426 relevant checks passed; exact
complete-row/media/mapping replay. No biological admission or runtime changes.
Continue host GTEx/FANTOM source lineage (65.04) and remaining categorical
adapters independently of unavailable content.

Host expression review **GT-HOST-TX-01**:
`results/host_expression_source_review_2026_10_08/` downloads primary GTEx v11
summary/LCM README and original Atlas `tpmss.tsv`, with verified associations and
provider/acquisition hashes. Atlas parsed values exactly match the archived table.
v11 remains a candidate with no new citation-age clock. Full gene rows preserve
PAR_Y, which native split-at-dot collapses. Current replay has three missing/two
extra GTEx projections and Atlas mapping/cutoff gaps. The documented 0.5 filter
per region reproduces every overlapping legacy value; small mapping/coverage gaps
remain. All rows/gaps/refusal diagnostics retained. 534 relevant checks passed;
counts/values/hashes/policy refusal replayed. 65.04 is 4%; 66.04 is complete.
Continue **GT-HOST-TX-02**, host mapping correction and original gene evidence.
No installed values or calibration changed.

Historical 66.04 first builder milestone (then 25%): `host.uniprot_index` now withholds ambiguous
symbols; Ensembl mappings/schema unchanged. `deposits.k562_rhoptry_evidence`
retains all original source rows before protein projection.
`results/host_symbol_mapping_2026_10_08/`: 20,010 exact gene records, 39 ambiguous
and 1,271 unmapped retained. The original builder compared 18,700 unambiguous
candidate protein rows using pandas' default tolerance; this did not prove
binary identity. 574 relevant checks passed, two existing skips. The later final
migration below completed withdrawal/gene preservation and strict precision
review, preserving pre-existing source/legacy encoding differences explicitly.

## Latest controlled completion: 66.04, 2026-10-08

[66.04](instructions/done/66_04_host_symbol_projections.md) completes instruction 66
(4/4): 39 ambiguous rhoptry protein projections withheld (117 cells), every other
installed human value exactly preserved, all 20,010 original gene records shipped.
Verified withdrawal ledger, exact external backups and executed migration/replay
are in `results/host_symbol_mapping_2026_10_08/migration/`. 575 follow-up checks
passed, two existing skips; two initial failures resolved. Strict precision review
corrects earlier tolerance-based exactness wording and retains pre-existing source/
legacy encoding gaps. Mouse/parasite feature/graph/track-record/claim hashes unchanged.
Total completed actions: 16/48; instruction 64 remains 9/40, categorical records
35%, literature-wide admission/replacement audit 4%. Continue frozen benchmark
adapters; host expression lineage and host gene spaces remain separate open work.

Latest 64.10 partition: native logistic conformal pilot completed and exact native
replayed in `results/conformal_logistic_pilot_2026_10_08/`. Fixed C=0.5/alpha=0.1,
train-only rank/model state; separate 569 calibration/560 test genes. Stored-label
set coverage 0.905357, mean size 6.975/26, zero singleton calls. 1,088 relevant
checks passed. Categorical action now 35%; no new completion tick, independent
biology/full outer coverage/missing adapters remain pending.

64.11 metric prerequisite now 10%: numeric cards retain MAE/RMSE for constant
baselines and small answered cohorts, with matched answered-row baseline skill.
`results/numeric_metric_review_2026_10_08/` contains executed analytic controls.
544 final checks pass; original historical calibration/data remain unchanged.
Next bounded numeric adapter: training-only native ridge on the existing
`fit_invitro_hff` candidate, original direct-experiment grade with unresolved
source/units/admission gaps; freeze scope before expensive fitting.

Latest numeric partition: 64.11 now 25%, native ridge on the existing HFF
fitness candidate, canonical `results/ridge_value_pilot_2026_10_08_v2/`.
Training-only exclusions/ranks/model; 4,015 train/1,106 tune/1,103 calibration/
1,101 test genes, 331 features. All training matrices/native predictions replay
exactly; MAE 1.269871, RMSE 1.590402, Spearman 0.701225, matched training-mean
MAE/MSE skill 0.335439/0.474868, coverage 1. Source/units/mapping/independent
biological admission remain gaps. Constant-float baseline decile diagnostic
preserved; revised cards retain identical values/models and correct code hashes.
1,123 final checks passed. Continue other numeric/categorical adapters and
biological truth; no additional full action complete (still 16/48).

## Active continuation/watchdog, 2026-10-08

The user explicitly requests a watchdog to keep the list moving and log stops.
Native goal continuation is active for unfinished instructions 64/65/66.
`starplast-work-watchdog.timer` is enabled/active, polling once/minute under
128 MB. `.starplast-watchdog/events.jsonl` retains private stop/progress events;
read it and state.json at resume, get actual goal status, record current item and
reason through `scripts/work_watchdog.py record`. Before a turn ends, log its
actual reason; progress reports are automatic-continuation boundaries, not
completion. Four earlier runtime invalid_prompt/content-restriction errors were
observed and logged, alongside eight other historical turn boundaries. Do not
bypass restrictions or invent a cause after an abrupt loss of runtime activity.
403 checks pass; completed [67.01](instructions/done/67_01_work_continuation_watchdog.md).
Total now 17/49 complete; 64.10 stays 35%, 64.11 25%, 65.04 4%. Repost all 49
rows whenever a complete item earns a tick. Details: [watchdog](docs/work-watchdog.md).

Latest numeric interval partition: **64.11 now 35%**, verified train-only ridge
reference plus separate 1,103 calibration genes and all 1,101 test genes in
`results/value_conformal_pilot_2026_10_08/`. Native point predictions, bounds,
absolute-error quantile and row/baseline cards replay exactly. Fixed nominal
coverage 0.9; stored-value coverage 0.921889, half-width 2.762314, finite
availability 1, mean width 5.524628 (2.517286 eligible-truth SD). Explicit
finite/unbounded/unavailable interval statuses and glossary metrics prevent
coverage being mistaken for capacity. Point-only cards gain null interval fields;
historical snapshots retain original values/code. 1,009 relevant checks pass,
10 initial fixture module-reference failures resolved. No biological admission
or exchangeability guarantee; units remain unresolved. Full model state is
verified upstream, not duplicated as a standalone deployment bundle.
No full action completed (17/49); 64.10 remains 35%, 65.04 4%.
Continue controlled remaining numeric adapters/variants and biological truth;
full outer coverage and shared user-facing scorecard component remain open.

Latest host partition: **64.21 human gene space now 30%**. Executed
`results/human_gene_space_foundation_2026_10_08/` builds 58,988 actual canonical
source-gene rows from verified GTEx v10, preserving all 59,033 original source
records and all 118,066 selected source cells exactly. The 45 qualified PAR_Y
records keep their 90 cells separately; they never overwrite ordinary X-gene
measurements. Reviewed UniProt crossrefs retain all 22,320 pairs, with 39,630
unmapped, 19,208 one-protein and 150 multi-protein genes; 1,713 source genes have
a shared protein. Counts/source positions/annotation versions are separate
metadata and never biological features. No protein feature is projected.
Per-column metadata plus typed measurement provenance retain exact source
hashes, median-TPM units, cultured-adult-fibroblast versus HFF context gap, reference
scope and redistribution/benchmark gaps. 618 checks passed; two wrapper source-role
errors before output writes recorded in `results/human_gene_space_preflight_2026_10_08/`.
All inspected installed host/parasite table/graph/track-record hashes unchanged.
No host space registered or distributable pack claimed. Continue controlled
reference/license/graph/independent-opening gates for this source set, then mouse
gene foundation and remaining sources. GTEx v11 still pending 65.04 comparison;
no claim that v10 is the best available replacement. Still 17/49 complete.

Latest Discoveries partition: **68.01 complete; 18/51 full actions complete**.
`starplast/discovery_labels.py` plus the actual Discoveries panel now browse 41
Toxoplasma and 30 Plasmodium labels, with function first. InterPro/Pfam/EC classes,
phenotypes, stages, structural labels and annotation flags no longer require a
generated claim to appear. Search terms/IDs, select overlapping class members,
open genes and existing legacy label/class scorecards. Claim filters/export/map
behaviour preserved; verifier-column rendering and empty-recipe lookup corrected.
Executed `results/discoveries_label_coverage_2026_10_08_v4/` independently replays
every variable and all 21,787/18,879 functional memberships exactly. Source tables,
claims, recipes and 1,038,372 held-out rows have unchanged hashes. **413 final
checks pass, one optional pdoc skip**; earlier broader app/navigation/help checks
pass with one corrected new fixture assertion. Three audit diagnostics are kept:
EC replacement codes mentioned in descriptions are not extra assigned classes,
empty InterPro delimiters and internal EC semicolons preserve source semantics.
No new functional claims or held-out rows: **68.02 is the next priority**, with
declared curation/domain-prediction/orthology lineage, current ontology and
unknown-negative semantics before fitting or reporting biological accuracy.
Instruction 64 categorical/numeric actions remain 35%; human gene space 30%;
dataset promotion/admission 4%. No native-goal completion or release implied.

Latest functional partition: **68.02 now 25%; still 18/51 whole actions complete**.
Official ENZYME/InterPro/Pfam nomenclature is pinned. All 34,983 domain memberships
and original descriptions preserved; 7,926 of 8,043 IDs have current names, 112
InterPro IDs unresolved and five Pfam IDs withdrawn with no reassignment. Current
names do not establish original assignment releases or independently measured activity.
Strict EC resolution gives 1,226/1,050 complete curated-source profiles in the two
parasites, plus 1,055 separate Plasmodium orthology-derived profiles.

Frozen FN-EC-01 seed23/k15/min-share0.3 fits training-only numeric kNN: 679 train,
175 tune,190 calibration,182 test; tuning/calibration unused. 349 inputs explicitly
withhold EC/domain/homology/attention summaries; graph empty. Exact recovery71/182
(39.0%), coverage171/182 (94.0%), correct71/171 among calls (41.5%); majority66/182
(36.3%). All native scores/model/card outputs replay exactly, including11abstentions
and an unseen full profile. This is weak source-profile recovery; **zero independently
admitted activity benchmarks, calibrated probabilities or deployment claims**.
Canonical pilot `results/functional_ec_knn_pilot_2026_10_08_v4/`, identity
`349317afd8f24b860d94ae1e11280a6130bf451700b02b9d2145fab85fc99162`.
Source reviewv2/domain metadatav4/UI bundlev4 and
`results/functional_discoveries_audit_2026_10_08/` retain executed evidence.
Actual Discoveries offers label/strategy/profile/seven major-class/control cards,
gene outcomes and selected-domain provenance. Both organisms' kinase searches
reach InterPro/Pfam; Plasmodium benchmark availability remains explicitly zero.
Original node/claim/recipe/record hashes unchanged. 653 distinct focused checks
verified (652 broad pass, then19 functional UI/reader checks after correcting
one HTML-tooltip fixture). Earlier Qt tuple lookup and method-docstring failures fixed.

Next bounded function work: protected-group Plasmodium EC recovery; complete
InterPro/Pfam profiles; source admission and independent activity truth; additional
adapters/calibration/deployment. Domain-content and domain-edge families differ,
and search.LAYER_SOURCES omits domain: explicitly exclude domain graph edges and
derived maps/representations with sentinel tests before graph benchmarks. Attention
aliases need a general exclusion follow-up; already explicitly withheld here.
Three worker agents were explicitly authorized by the user and delivered metadata,
UI scorecards and the ordinary precompute framework. Keep concurrent small jobs
within the combined2GB lease and root Qt/app checks serial; defer GPU/>=8GB jobs
during the plaque priority. Two actual worker invalid_prompt stops were logged;
root continued independent tasks without retrying restricted content. Partial
`scripts/review_functional_activity_truth.py`, `results/functional_activity_truth_2026_10_08/`
and `tests/test_functional_contract_audit.py` remain untouched/unreviewed/uncommitted,
not admitted evidence. Private watchdog logs retain actual causes; timer remains active.

Latest precompute partition: **64.25 now35%**, serial typed job graph, complete
parent content IDs, input snapshots, validated/corruption-refusing resume,
interruption recovery, failures/blocked descendants, atomic journal/writer lock
and cooperative attempt/time/currentRSS/disk/package budgets. Executed
`results/precompute_ec_views_2026_10_08/` retains all182 real frozen EC test rows,
recomputes exact profile/seven-member cards, stops after one partition, resumes
only two unfinished builders and replays with zero builders. Stages0.350/0.957/
0.328s; process lifetime peak187.3MiB, typed artifacts166,244bytes, recorded evidence
688,151bytes before addedREADME.83 pure framework/artifact/split checks pass;
34 runner checks also pass within the final broad app run using deterministic
fixtureRSS. Remaining production feature/map/model/deployment builders and adapter
coverage keep this item open. No new tick (18/51). Functional work pushed as682a399;
no original runtime strategy, dataset promotion, native-goal completion or release.

Latest integration partition: **68.02 30%, 64.17 50%; still 18/51 whole actions**.
Three authorized workers delivered complete recorded-domain profiles, the exact
functional coverage matrix and shared scorecard presentation; they continue bounded
profile-target, production-feature and class-card software partitions. Domain census
retains all34,983 memberships, including withdrawn/missing-current identifiers.
Coverage has880 addresses/11installedlabels/oneverifiedrecovery artifact and zero
independent biology/calibrated deployment. Discoveries exposes current-organism
240/320address coverage with exact recorded-result routes. Exact installed-node
hash/table binding prevents imported/altered contexts from borrowing archived accuracy.
Functional source exclusions close registered/declared derived inputs and reject
cycles, unresolved provenance and missing encoded questions. Existing artifacts
are not retroactively rewritten.

Shared immutable cards expose definitions/full record/export, source/split/control
gaps and readable evidence counts; actual human source status retains58,988canonical
genes without host inference/registration. Task-aware set metadata preserves native
algorithm schemas. Profile/major class populations remain distinct. The final
executed UI source replay is recorded in
`results/functional_coverage_ui_2026_10_08_v3/`; earlier prototypes remain immutable.
FN-EC-02 stopped with an actual runtime content-access restriction; do not retry or
reconstruct that task. Restricted partial activity-truth files remain unreviewed.
Watchdog logs the actual stop and ongoing independent work.

GitHubCI for b5b854b failed before regressions because default pytest collected
archived immutable test copies with duplicate module names. Commit043a2ad limits
default discovery to live tests and is verified pushednightly; normal collection
now succeeds without archive paths. Preserve archivedcopies. The replacement CI
run is pending; do not claim full-suite success or release readiness.

Final integration validation: **738 focused checks pass**, package-version check
passes, V3 actual UI/source audit passes, all51manifest-covered outputfiles verified.
Prepared-profile target software has25pure checks; numeric intermediate feature
builder/runner has59pure checks under400MB. The latter pins source/split/exclusion
and imported guardcode, preserves native train ranks/frozen queryECDF, and rejects
complex/attention inputs. No actual production feature pilot yet;64.25stays35%.
Root combined test RSS fixtures are deterministic; actual external limits remain.
Continue a bounded executed feature-operator pilot and domain prepared-target
capacity census, then remaining categorical/numeric/source admission partitions.

Latest bounded audit partition (2026-10-08): three worker jobs are terminal.
Prepared complete-Pfam targets are canonical in
`results/functional_profile_targets_2026_10_08_v3/`; earlier serialization failures
remain preserved. Source/target/split/capacity replay retains all 8,140 genes and
333 unsupported test profiles without fitting. The numeric operator pilot in
`results/precompute_features_2026_10_08/` exactly replays all 349 columns and 1,226
ordered EC entities, including child-only resume and zero-builder replay.
Its 401.6 MiB observed process peak differs from external 400 MiB cgroup accounting;
do not describe process peak as below 400 MiB. The six-task presentation audit in
`results/scorecard_task_views_2026_10_08/` passes 539 software checks, not biology.
Output receipts (20/33/20 files), input/code hashes and three script parsers pass.
**64.25 is now 45% (4–6 h); 64.17 stays 50%; 68.02 stays 30%; still 18/51 complete.**
User now requests minimum tokens and smallest tasks first: integrate ready evidence
and choose bounded existing acceptance gaps before expanding models or sweeps.
Restricted activity-truth/FN-EC-02 work remains unreviewed and must not be retried.


Latest acceptance: **64.17 complete; 19/51 full actions, 10/40 instruction-64 actions**.
Dedicated calibration navigation and explicit unavailable date/version indicators
preserve raw source snapshots and native metric schemas. The canonical desktop
packet is `results/desktop_scorecards_2026_10_08_v2/`: 1,573 checks across 14
actual Qt cards, 40 verified receipts, 0.68 seconds/150.0 MiB process peak.
The first attempt preserves an incidental JSON metric-order diagnostic, corrected
by comparing full metric records by identity with exact values/definitions.
Fresh Discoveries UI/source replay is canonical in
`results/functional_coverage_ui_2026_10_08_v4/`; all 51 nested receipts/input hashes
pass, as do 90 focused view/browser/host/class/Discoveries tests and version check.
This completes reusable presentation only; no biological truth, calibrated
functional deployment, host installation or scientific adapter admission added.
Choose the next small browser/navigation acceptance gap while larger biological
source and model work remains open. Restricted partials stay untouched.

Latest acceptance: **64.18 complete; 20/51 full actions, 11/40 instruction-64 actions**.
Tools → Datasets shares the Slots window. All 162 sources reconcile to 180 scoped
records, with explicit host/entity qualification, original stored values, source
cards, bounded paging and separate session imports. Processed availability is
separate from raw-assay provenance; stored slots do not imply experimental grade.
Canonical `results/dataset_browser_2026_10_08_v8/`: 3,900 notebook checks,
27 current input hashes, 24 verified receipts, 98.05 seconds/1,082.8 MiB peak.
505 focused regressions pass, one existing GL-context skip; six prior CI failures
and native interval checks pass (24). Full CI for old 9d7752b failed; its six
specific regressions are fixed, replacement nightly CI must still run.
Two native Qt test failures remain recorded with unproven causes; after explicit
new slot-test iterator cleanup, the same combined scope passes. A worker's 400 MiB
explanation audit was OOM-killed, then passed serially. Preserve all diagnostics.
No biological truth, strategy fits, host installation or source promotions added.
Next smallest bounded action: 64.19 gene evidence entry, reusing the query,
provenance, scorecard and dataset-browser contracts. Human foundation stays 30%,
mouse 0%; 65.04 stays 4%, 68.02 stays 30%. Restricted partials remain untouched.

Latest acceptance: **64.19 complete; 21/51 full actions, 12/40 instruction-64 actions**.
Exact organism-qualified canonical and recorded alias lookup exposes collisions,
including index targets absent from the current table. Product matches require
an explicit candidate click. Unavailable alias targets clear stale selection.
The gene card opens condensed/all-column evidence grouped by biological question,
with exact original values, source cards and dataset-address navigation. Label
and class routes preserve complete legacy ledger results; changed tables cannot
borrow that accuracy. Replacing the dialog releases its predecessor.
Canonical `results/gene_evidence_2026_10_08_v5/`: 1,246 executed checks over
442 Tg/167 Pf columns, GRA16/AMA1 aliases and two recorded label/class routes;
22 current input hashes, 20 output receipts, 102.68 seconds/804.1 MiB peak.
473 relevant regressions pass with one existing GL-context skip; final widget
width checks pass (4). Initial audit diagnostics preserve a wrong Pf column and
an unscored AMA1 class; no annotation or scorecard was fabricated to pass.
Next small action64.20: labels/classes as entry points, reusing Discoveries,
query/provenance/dataset/gene contracts. Scientific coverage, hosts and dataset
promotions remain open. Preserve restricted partials and nightly/full-CI gates.

Latest acceptance: **64.20 complete; 22/51 full actions, 13/40 instruction-64 actions**.
Tools → Labels and classes opens all 41 Tg/30 Pf native labels, function first,
with membership/unknown counts and separate records for all 39 strategies.
Original complete protein profiles, overlapping classes and known False values
remain intact. Class precision/recall/confusion and label overall/macro metrics
retain exact legacy settings/seed/cohort; missing biology/hierarchy stays explicit.
Actual class→gene→label→source navigation passes for both organisms. Ambiguous
source questions clear stale cards and records. Changed contexts refuse accuracy.
Final immutable V4 Tg/Pf packets: 418+459 checks, 18 input hashes/20 receipts each;
178.64/181.66 s, 1849.4/1373.2 MiB peaks, serial under the same 1900 MiB cap.
494 relevant checks pass; version/whitespace pass. Preserve first categorical
failure, successful V2 and combined V3 OOM diagnostics. No biological truth,
source promotion, strategy fitting, deployment or host installation added.
Latest small partition: 68.02 training-only complete-Pfam baseline controls on
verified prepared-target V3. The exact control partition is accepted: 16 focused checks, 35 input hashes and
22 receipts verify; the executed 400 MiB-capped notebook finishes in 6.676 s.
All 635 test genes/333 unsupported profiles remain; majority/prevalence recover
19/2 profiles. This is annotation-control arithmetic, not biological accuracy.
68.02 stays 30%; classifier fitting/source admission/deployment remain open.
Full goal remains active; nightly CI/release gates and restricted partials remain.

Latest documentary gate: human V5 terms/reference review is accepted, with 25
foundation/source receipts verified twice, five authoritative captures, five
GTEx V10→GENCODE39/GRCh38 checks and 19 output receipts. Conditional GTEx public
redistribution and UniProt CC BY4 documented; future pinned v10 must show older
version/source/date attribution. GENCODE39 is GRCh38.p13. V11 linkage, pack notices,
reference-gene mapping/graph/opening, publication and biological admission remain
open. Human stays30%, mouse0%,65.04 stays4%; no source/foundation promoted.
Preserve V1–V4 diagnostics/captures; actual V5 review111.2MiB/0.44s under400MiB.
Functional attention-alias closure64 checks pushed37f6770. Complete-Pfam native
kNN pilot is frozen at k15/share0.3 with train-only source-closed numeric inputs.
Initial categorical parity failure is NaN/None representation; only missingness
normalization changed,13 synthetic checks pass. Root verifies actual full matrix
and native outputs under800MiB before canonical replay. Two signal9 worker
partials remain diagnostic, without asserting an unproven OOM cause. No new tick.

Latest functional acceptance: canonical complete-Pfam kNN V4 at
`results/functional_profile_knn_2026_10_08_v4/`; fixed k15/share0.3, 348 inputs,
2,381 training/635 test genes, 333 unsupported test profiles. Exact native ranks,
votes/support/scores, refit, serialization and typed artifact replay pass.
12 correct/eight wrong/615 abstentions: all-eligible recovery1.89%, coverage3.15%,
correctness among calls60%; majority control recovers19/635 (2.99%). This weak
annotation recovery is not admitted biology, calibrated confidence or deployment.
47 input/47 output/12 artifact receipts independently verified in the acceptance
notebook;124 focused software checks pass. Runtime46.61s/975.92MiB process peak
serial under1900MiB, all workers terminal. Earlier unit/population metadata refusal
and authoritative800MiB systemd OOM preserved; only categorical missingness,
metadata and retained-copy lifetime changed. No relaxed numeric checks/settings.
Whole actions22/51,64.10 stays35%,68.02 stays30%. Next small partition: generalize
the existing EC-only functional result reader/packager to complete domain profiles,
then expose this exact Pfam benchmark and class/gene cards in Discoveries. Independent
truth admission, other mechanisms, human/mouse packs and dataset promotions remain.

Latest UI integration accepted: Discoveries now ships exact EC+complete-Pfam
benchmarks, with namespace-specific profiles,479 Pfam profile cards,583 domain
membership cards and all635 held-out genes/333 unsupported profiles. Recorded
presence/complement never claims biological negatives. Native control scopes and
full split identity/role counts retained. Candidate/merge/source/UI receipts pass;
118 relevant regressions pass. Actual UI V3 audit10.17s/1,067,978,752byte process
peak under1900MiB, all workers terminal. Preserve V1 authoritative800MiB OOM and
V2 incorrect export-envelope key diagnostic. Shipped bundle externally pinned;
no fit, biological admission or deployment. 68.02 now35%, whole actions22/51.
Next smallest bounded function work: protected-group Plasmodium EC preparation
using existing source/exclusion/controls contracts, preserving missingness and
source-grade gaps. Restricted activity-truth partials remain untouched; native
goal/watchdog active, host spaces and dataset promotions still open.
