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
40 bounded action cards with dependencies and acceptance tests. **9/40 complete (22.5%)**:
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
the entire progress table (40 instruction-64 rows, four instruction-65 rows and four instruction-66 rows).
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
strategies, calibration, data and shipped claims remain at the released baseline.

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
the old first-accession policy: correct these across host sources under 66.04.
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
meta-analysis and scorecards at every level. 64.01–64.06 are complete; start at **64.07** and
respect each card's dependencies. Completed cards record fixtures, validation and limitations.
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

66.01 is complete: all 162 sources accounted for, 69 processed bindings, 11 container manifests and 538 candidate files, with explicit association, transform and access gaps. Read `results/source_recovery_final_2026_10_07_v2/README.md`. CSPA exact installed-field reproduction verified; its non-detection semantics corrected without numerical change. Frozen benchmark readers support explicit checksum-verified relocation. **717 relevant checks passed.** Instruction 66 is 3/4 complete; continue 64.07 strategy capability contracts and 66.04 host mapping correction. Numerical admission/replacements remain 65.04. Earlier code identities remain archived with their immutable artifacts.

64.07 complete: `starplast/capabilities.py` and `results/capability_contracts_2026_10_08_v2/`; distinct pair/class/evidence/gene outputs and explicit unavailable host adapters. No new biological performance claims. Next 64.08 artifact roles/identity/invalidation, then 64.09 shared aggregation.

64.08 complete: `starplast/artifacts.py`, `results/artifact_contracts_2026_10_08_v2/`; six task types and five storage roles, ten exact synthetic round trips. Frozen split guards reject misassigned fit/test entities and known benchmark genes as unknown deployment. No estimator fitted. Next 64.09 shared row-level aggregation/reconciliation.

64.09 complete: `starplast/record_scorecards.py`, `results/record_reconciliation_2026_10_08_v3/`; all legacy outcomes match standard metrics on identical cohorts. Original truth/exclusion/nested-fit lineage remains unresolved, not retroactively admitted. Settings/seeds/modes/named sets remain separate; no independent-sample inflation or component/per-gene probability claims.

64.10 is 30%: first `starplast/label_records.py` feature-kNN adapter (`8b59643`), with native
training rank parity and frozen train distributions for held-out values. The
prediction-grade compartment pilot retains all 560 test genes and full native
class scores; `results/label_knn_pilot_2026_10_08_v5/` is the rank-transform run.
Earlier z-score/serialization diagnostics retain original code. Scorecard
identities now preserve one-ULP differences. Continue remaining categorical
adapters/full outer coverage and independent truth acquisition; host-symbol
correction remains 66.04. No runtime strategy, installed values or calibration
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
counts/values/hashes/policy refusal replayed. 65.04 is 4%; 66.04 remains open.
Continue **GT-HOST-TX-02**, host mapping correction and original gene evidence.
No installed values or calibration changed.

66.04 first builder partition (25%): `host.uniprot_index` now withholds ambiguous
symbols; Ensembl mappings/schema unchanged. `deposits.k562_rhoptry_evidence`
retains all original source rows before protein projection.
`results/host_symbol_mapping_2026_10_08/`: 20,010 exact gene records, 39 ambiguous
and 1,271 unmapped retained; all 18,700 unambiguous candidate protein rows match
the installed three scores exactly. 574 relevant checks passed, two existing
skips, source replay verified. Installed wrong protein projections remain
pending migration, explicit withdrawal/gene-preservation provenance and affected
non-loss/layout/calibration/benchmark checks. Continue that bounded partition.

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
