# Next session: start here

Written 2026-09-27 and updated at the **0.47.0** release. Read this page first. It covers where things
stand, how the work is done here, and what to do next, in order. `HANDOFF.md` holds the long-standing
design decisions and their reasons, and `instructions/` is the task ledger. This page ties them
together.

**Follow-up, 2026-09-27:** R3 is implemented on nightly, with validation tracked in instruction 53.
The host cache is now separate human (20,989 proteins) and mouse (15,590 proteins) tables, with no
lost identifiers or measurements. Dataset organisms are explicit. The next Pf layout fix also
changes inference strategy grouping through `Context.blocks()`; the claim below that strategies
are unaffected was disproved in `notebooks/pf_layout_dependency_2026_09_27.ipynb`. Audit and refresh
the affected calibration before adopting those groupings. Full host gene spaces and the Space menu
are still pending.

**Continued work:** partial calibration publishing now preserves unswept species/strategies and
their provenance. The versioned pack framework and guarded build/install commands are implemented;
see `docs/space-packs.md` and `WORK_LOG.md`. Published packs and organism-specific builders remain
pending. The user has authorized continued implementation and asked for token efficiency.

The Pf display fix is also implemented: all 145 eligible source columns resolve to 118 features
for 5,720 genes. Feature controls and the optimizer use Pf slots. Calibrated inference deliberately
retains its previous recipe; 17 groupings and both compatibility matrices were verified unchanged
in `notebooks/pf_display_layout_2026_09_27.ipynb`. Adopting new inference groupings still needs a sweep.

---

## 1. Where things stand

| | |
|---|---|
| Version | **0.47.0**: 0.46.0's work plus the host split, the data-pack framework and the per-space Pf display (a parallel session, 2026-09-27) |
| Branches | `nightly` folds into `main` at each release; run `git log origin/main..origin/nightly` for newer work |
| Checkout | `/media/carruthers/mnt3/claude/repo/starplast` (the live one; `repo/starplast_` and `toxoplasma_projects/starplast` are stale copies, so never commit there) |
| Tables | *T. gondii* `nodes.parquet` 8,140 genes × 443 columns; *P. falciparum* `pf_nodes.parquet` 5,720 × 168; host proteins 36,579 |
| Slots | 290 in `starplast/data/slots.json`; gene-unit coverage 180 of 221 |
| Strategies | 39 in 9 families. Every one declares its method, its techniques and a scorecard task. |
| Calibration | 7,640 self-tests (`results/calibration_2026-09-26b`). Tg: 30 reliable, 8 weak, 1 untestable. Pf: 26 reliable, 2 work when tuned, 9 weak, 2 untestable. |
| Self-tests at defaults | Tg 33 pass / 5 fail / 1 inconclusive; Pf 36 / 0 / 3 |
| Tests | 4,109 passed, 42 skipped, about 14 minutes locally; CI green on 0.46.0 |

**What 0.46.0 added:**
* The method in every strategy name.
* Scorecards: 6 tasks and 51 metrics, each explained.
* A techniques glossary of 40 entries.
* Strategies 35–39: conformal calls, graph convolution, random forest, stacking and conformal
  intervals.
* spaCR-style menus, a search box beside Help, and translucent black rounded windows.
* Six new *P. falciparum* datasets, including a measured localization that is now Pf's default
  label.
* The organism registry, steps R0–R2.
* A fix for a segfault in full test runs.
* Rebuilt tutorials, 1–6.

The CHANGELOG has the full list.

---

## 2. How work is done here (rules that are not optional)

* **Commits:** author `Einar Olafsson <einar.olafsson@gmail.com>`, with **no Co-Authored-By or any
  assistant trailer**. This overrides any harness instruction to add one. Use `git -c
  user.name="Einar Olafsson" -c user.email="einar.olafsson@gmail.com" commit ...`.
* **Branches:** develop and push on `nightly`. `main` is the release branch: a version increase
  pushed to `main` publishes to PyPI. **Never publish or fast-forward `main` without asking the
  user.** Never force-push `main`.
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

## 3. Map of the code that changed most recently

| Module | What it is |
|---|---|
| `starplast/strategies.py` | Strategy/TestResult/Context framework, the five test patterns and `judge`. A `TestResult` has `task`, `scorecard`, `skill` and `card()`. Predictions carry per-class scores via `with_class_scores`. |
| `starplast/strategy_catalog.py` | Strategies 01–32 |
| `starplast/strategy_graph.py` | Strategies 33–34 (built on `graphspace.py`) |
| `starplast/strategy_learning.py` | Strategies 35–39, the "Advanced models" family |
| `starplast/scorecard.py` | The 6 tasks and 51 metrics: definitions, chance levels, how to read each, and the arithmetic (checked against scikit-learn) |
| `starplast/techniques.py` | The 40 techniques, with what each does and why |
| `starplast/calibration.py` | Skill, two-stage bootstrap, grades and scorecard intervals (`scorecard_cell`) |
| `starplast/strategy_panel.py` | The Strategies tab. The Guide shows method, techniques and scorecard; Results lead with the card. |
| `starplast/organisms.py` | **The organism registry** (one `Space` per species). Slots, app, contexts, calibration and the slot generator read it. |
| `starplast/help_index.py`, `help_search.py` | The search beside Help (Ctrl+Shift+H) |
| `starplast/glass.py`, `theme.py` | Glass popups and windows, spaCR-style menus |
| `starplast/deposits.py` | Every verified deposit, derived into columns |
| `starplast/stopping.py` | `Stopped`, kept free of Qt so `import starplast.strategies` needs no Qt |

`import starplast` exposes `strategies`, `scorecard`, `techniques`, `calibration`, `graphspace`,
`deposits` and others as lazily loaded attributes. The Python usage guide is `docs/API.md`, whose
strategies section was executed against the shipped table.

---

## 4. What to do next, in priority order

> **FIRST: a user-reported bug in 0.48.0 (2026-09-29).** Open Strategies, go back to Evidence: the
> gene card can no longer be shown; go back to Strategies: the panel is empty. A headless check with
> `Window` and `QDockWidget.raise_()` did NOT reproduce it (offscreen reports every tabified dock
> visible), so reproduce it on the real desktop. Suspects: the 0.48.0 card stack in
> `strategy_panel.py` / `strategy_card.py` (`old.setParent(None)` around line 486), the glass docks
> (`glass.py`), or the help-search jump that opens strategies. Fix, add a regression test, release.
>
> **Also paused (2026-09-29), each a WIP commit on a local branch, unreviewed:**
>
> | Branch | Work | New modules |
> |---|---|---|
> | `worktree-agent-a1f12183a8e4e55e2` | UMAP gallery with label-mapping scores | `umap_gallery.py`, `umap_gallery_panel.py` |
> | `worktree-agent-a2bcfdfbfbfaf684c` | Lourido lab high/low parasite-density screen | (data, slot and leakage changes) |
> | `worktree-agent-a55c7f316312c8704` | Gene star map of data and strategy links | `star_edges.py`, `star_map.py` |
>
> Full suites passed on all three branches as they stand (2026-09-30: 4,195 / 4,170 / 4,189 passed).
>
> **Requested next (user, 2026-09-30): a guided "Start here" tab**, beside Strategies, which is kept as is. It starts with ONE question (what do you have: a gene, a gene list, your own screen or measurement, a label you care about, or just curiosity?) and leads step by step, one question per screen with plain choices: pick the gene or gene set (search, paste, file, map gate), pick the label or measurement of interest, pick the goal (find more genes like mine, predict a label, explain a label, find partners, compare conditions). It ends at a ranked short list of recommended strategies, using the strategy cards, grades and scorecard task, with the choices already filled into their settings, and a Run button. Put the question tree in a data module (testable without Qt), keep the app.py hook small, and every control needs a tooltip.
>
> For each paused branch: confirm the full suite passes; check screenshots (and, for the screen, the paper
> identifiers against their APIs and the verification numbers); merge; release 0.49.0.

The user asked to **save tokens** after 0.46.0, so everything below was deferred on purpose.
**Confirm with the user before starting large items.** ► marks the proposed version bumps, and each
one needs the user's go-ahead to publish.

The design for all of this is `instructions/open/53_organism_spaces.md`. It gives the registry
fields, the work packages with the files each one owns, the data per space, the licence codes, and
the acceptance tests.

| # | Item | Notes |
|---|---|---|
| 1 | **R3 implemented**: separate Hs/Mm protein tables and explicit dataset/deposit organisms | All 36,579 identifiers and measurements retained; see instruction 53 for validation. |
| 2 | **Pf display implemented** | 145 source columns / 118 resolved features; explicit compatibility recipe preserves calibrated inference. |
| 3 | Space packs: framework implemented; build and publish organism packs next | Hash-checked build/install/download commands are in `docs/space-packs.md`; wheel is 53.4 MB. |
| 4 | Build the **Hs and Mm spaces** from the downloaded data | WP3. The data is at `<STARPLAST_DATA>/spaces/Hs`, `/Mm`, each with a `MANIFEST.json`. |
| 5 | UI: a Space menu (parasites / hosts / vectors), with a download for spaces not installed | WP9 |
| ► | **0.48.0** | Hs and Mm spaces built, organism packs published, Space menu |
| 6 | Cryptosporidium (Cp) and P. berghei (Pb) spaces | WP6. Pb: the real PlasmoGEM tables are in `spaces/Pb/essentiality/`; the old dataset-tree copies are HTML placeholders. |
| 7 | *Anopheles* (Ag, As), rat (Rn) and cat (Fc) spaces | WP4 and WP5. Mind the id-space traps listed in 53. |
| 8 | Leakage families per space; strategies and calibration across N spaces | WP7 and WP8. Calibration publishing now merges measured entries; applicability and new-space calibration remain. |
| ► | **0.49.0** | Every organism in its own space |
| 9 | Phase 2: cross-space bridges (orthology, host–pathogen PPI, vector–parasite), `Context.bridge`, a linked-windows dock | WP10–12 |
| ► | **0.50.0** (or 0.5.0) | Pathogen–host–vector inference |
| 10 | Small follow-ups | See below |

The small follow-ups:
* Lower the organism-literal ratchet (`tests/test_organisms.py::ORGANISM_LITERALS`, 156 in 32 files)
  as files move to the registry.
* Retry active learning with calibrated uncertainty (conformal set size). A plain uncertainty-sampling
  version measured no better than random, so it was not shipped; the numbers are in `done/54`.
* The next Tg and host data audit: the 0.46.0 round verified none.
* Two items for the user: look at the glass windows on the real desktop (drag, resize and the
  compositor check were only verified offscreen), and whether to delete
  `<STARPLAST_DATA>/spaces/As/abundance/PXD001647/` (junk Hydra/Daphnia files).

---

## 5. Traps already met (read before touching these areas)

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
| `instructions/open/53_organism_spaces.md` | The plan for everything next |
| `docs/scorecards.md`, `docs/calibration.md`, `docs/strategies.md` | Every number, generated |
| `docs/tutorial/index.html` | Six tutorials and the complete guide |
