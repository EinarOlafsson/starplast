# Start here

Read this, then `INDEX.md` for the status table, then the numbered file for whatever you pick up.
`../HANDOFF.md` carries the project's standing decisions and is the source of truth for why anything
is the way it is.

## First five minutes

```bash
cd /mnt/firecuda2/Claude/repo/starplast     # this machine
git pull && pip install -e .
QT_QPA_PLATFORM=offscreen PYQTGRAPH_QT_LIB=PyQt6 python -m pytest tests/ -q
```

Expect ~1,855 passing, 4 skipped, and roughly four minutes. Anything red is from the last session,
not from you.

## The list is clear

Every task in `INDEX.md` has landed (2026-08-13, v0.17.1). What a fresh session should know before
adding to it:

- **`search.search` streams.** `on_run` gets a `RunStep` per embedding — its best clustering under
  the objective in force, its per-category scores, its coordinates and its labels. `walk_umap_iter`
  does the same for a plain hyperparameter walk. Drive those rather than re-running anything.
- **`Window.placed`** is the mask of genes the displayed embedding covers; unplaced genes are hidden,
  unpickable and excluded from `visible_mask`. **`Window.gated`** is the set someone drew a gate
  around, which is what annotation takes as input.
- **The color-by panel is the one place color is chosen** — columns, kept clusterings (`runs.py`)
  and binned quantities. `Window.category_values()` is what everything reads.
- **Annotations refuse to be saved without a validated precision**, live in their own file, and draw
  in a color used for nothing else. That refusal is the point of the feature, not an obstacle in it.
- **Results save and load** (`results.py`): one table is a CSV whose first line names which table it
  is, every tab at once is a zip of those plus a manifest, and a loaded row is as clickable as a
  computed one. A file whose kind does not match the tab is refused rather than loaded.
- **A user's own table imports** through `importer.py` and `File ▸ Import data…`, with the
  preprocessing offered rather than assumed and every choice recorded. Imported columns are prefixed
  and joined in memory; the cache on disk is never written.
- **One thing is known and deliberately not on the list**: the *Cryptosporidium* fetch bug recorded
  in `HANDOFF.md`.

## Two environments, and the second one is the one that matters

    /home/olafsson/anaconda3/envs/spacr       pandas 2.3.3   <- tests run here
    /home/olafsson/anaconda3/envs/starplast   pandas 3.0.5   <- the user runs here

**Three bugs have now shipped that were invisible on pandas 2 and total on pandas 3.** After any
change touching dataframes, run it under the second interpreter too:

```bash
/home/olafsson/anaconda3/envs/starplast/bin/python -c "..."
```

The pattern each time: an array arrives read-only from a library and is then written in place
(`embedding.embed`, `build_graph.embed`, `embedding.build_matrix`), or `.astype(str)` keeps NA where
it used to produce "nan". `app.as_text` and `embedding.as_text` exist for the second.

The user also runs a second checkout on another machine
(`carruthers@carruthers:/media/carruthers/mnt3/claude/repo/starplast`). **Unpushed work is invisible
there**, so push before asking them to test anything.

## What this application is for

Letting a user look at many structures, colored by many held-out variables and by clusters, and
judge whether a cluster means anything. Everything else serves that. When a change makes it harder
to see a structure or easier to believe a cluster, it is wrong however well it works.

## The standing constraints

- **Measurement, inference, absence and annotation never read as one another.** Grey means unknown;
  never zero, never a category.
- **No candidate without the number that says how much to believe it.** On the full proteome, one row
  out of 7,710 (configuration x compartment) reaches F1 0.5 -- `PM - peripheral 1`, at exactly 0.500.
  A cluster that looks pure is a lead, not a result.
- **Every objective has a degenerate maximizer.** Singletons win any purity score; one giant cluster
  wins any recall score. `objectives.py` enforces the floors and `EXPLANATION` has the measured
  table. It is shown in the app under Help.
- **Task 17 must not ship before task 20**, which it now does not have to: 20 landed. What survives
  is the requirement — a saved annotation carries the validated precision and recall for its
  category, the cluster's composition, and the enrichment over prevalence.

## Conventions worth knowing before you write anything

- **American spelling everywhere**, text and identifiers both: `color_of`, `color_mode`,
  `normalize`, `localization.py`. The rename landed on 2026-08-12 (`done/25b`). A block name lives
  inside every saved recipe, so `embedding.RENAMED_BLOCKS` translates old spellings on the way in --
  any future rename of a block needs an entry there or it silently rebuilds a different map.
- Tooltips say WHY a control exists and what choosing badly costs, not what it is called. Bounded
  controls explain their bounds. A test fails if any control has none.
- Every public module, class, function and method has a docstring. A test fails otherwise.
- **Results tables go through `AnalysisPanel.results_table`**, which gives them a row action (click a
  row, see the map it is about) and a right-click menu (save as CSV, copy rows). A table wired by
  hand is a table that silently lacks both.
- **Coverage is 100% on every module** (7,389 statements, 1,849 tests) and stays there. The way to
  cover a Qt-thread body is to call it directly, never a pragma; genuinely unreachable branches get
  deleted. Two functions were deleted rather than covered in the last pass, and writing one of the
  missing tests found a real defect in `objectives.adjusted`.
- **The suite must not write into the user's own state.** `STARPLAST_STATE` points kept runs,
  annotations and saved embeddings at a temporary directory for the run (`conftest.py`), the way
  QSettings is already redirected. Both isolations exist because both were violated, and each was
  found by looking at a real machine rather than by a failing test.
- Commit messages explain the reasoning and admit what was got wrong. Write them to a file and use
  `git commit -F` -- backticks in `-m` have twice executed shell commands here.
- Bump the version for feature work (currently 0.17.3).

## Where things are

    starplast/app.py             the window, menus, drawing, exports
    starplast/analysis_panel.py  the six analysis tabs
    starplast/objectives.py      what "good structure" means, and the guards
    starplast/validate.py        putting an error rate on an annotation
    starplast/jobs.py            background work, stopping, state
    starplast/lod.py             the three levels of detail
    starplast/gallery.py         thumbnails of a sweep, grid and scroll
    starplast/runs.py            kept clusterings, with their recipe and their mask
    starplast/results.py         saving and loading what the tabs computed
    starplast/importer.py        a user's own table, with the preprocessing offered
    starplast/celldiagram.py     the parasite under the compartment list
    instructions/INDEX.md        status table
    results/                     every published number, with the script that made it
