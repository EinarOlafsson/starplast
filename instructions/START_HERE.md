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

Expect ~1,501 passing, 4 skipped, and roughly four minutes. Anything red is from the last session,
not from you.

## The one thing to do first

**Task 23** — the left panel becomes "color by". Then 22 (the cell diagram), 25 (American spelling)
and 24 (the logo drafts).

**Tasks 16 and 17 are done** (2026-08-12, v0.9.0 and v0.10.0). `search.search` streams a `RunStep`
per embedding — its best clustering under the objective in force, its per-category scores, its
coordinates and its labels — and the Search tab fills a per-configuration-and-category table as it
runs, with `search.frontier` marking what nothing beats on both mean and best F1.
`annotations.AnnotationStore` is a CSV of its own that **refuses** a row with no validated precision
or a precision under 10%, is never written into the node table, and draws in a colour used for
nothing else.

**Task 18 is done** (2026-08-12, v0.8.0): the left button has a Navigate mode (free orbit or
constrained to one axis, so a view can be returned to) and a Select mode (2D lasso in screen space,
3D brush in world space). A gate fills `Window.gated`; the evidence panel shows the set's composition
and `File ▸ Export gated selection` writes it out.

**Task 19 is done** (2026-08-12, v0.7.0): `logging_util.get_logger(__name__)` anywhere, levels
DEBUG–ERROR, an opt-in rotating file, and a console level that is separate from the file's and is
never silent — a warning nobody enabled a log to see is a silent failure with extra steps. Console
lines carry a `[LEVEL]` prefix, which is what the pane colours and filters by. **Log the outcome of
anything that can fail quietly**, with the URL or the recipe attached; that is the whole point of
it.

**Task 20 is done** (2026-08-12, v0.5.0), and two of its findings outlive it:

- The circularity guard never fired. It compared a category VALUE against a list of column names,
  and the panel passed BLOCK names, so the tab that exists to refuse circular questions would have
  scored a target the map was built on. It now tests the label COLUMN, before the job is submitted
  as well as inside it. `AnalysisPanel.used_columns()` is what to pass anywhere else this matters.
- Every score on screen now carries what it has to be read against: `prevalence` and `lift` in the
  per-category table, `overall_frac_category` and `enrichment` on every candidate. On the real table
  the best cluster is 10.2% nucleus-chromatin against a map that is 9.4% — enrichment 1.09.

**Task 15 is done** (2026-08-12, v0.4.0), so the contract it was blocking is now this:
`tuning.walk_umap_iter` yields one `WalkStep` per configuration as it is computed — the scores, the
coordinates, and a boolean mask over the node table saying which genes those coordinates are for.
`walk_umap` is that collected and ranked, with an `on_step` callback, so existing callers are
unchanged. The walk table and the gallery dock each fill one row and one thumbnail at a time.
**Task 16 is unblocked, and should drive that iterator rather than re-running the walk.**

Two consequences anything touching the map needs to know:

- `Window.placed` is the mask of genes the displayed embedding has coordinates for. A walk map covers
  a subsample, so most genes have no position in it: they are hidden, unpickable, and excluded from
  `visible_mask`, which is what points, centroids, edges and the gene count all read. They used to be
  left at the origin, drawn and clickable.
- `embedding.normalise` puts every map at the same extent, and is what to use for anything that
  produces coordinates — otherwise the camera has to be re-framed for each one.

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

- **American spelling in user-facing text.** Identifiers are still British in places
  (`colour_of`, `colour_mode`, `localisation.py`) -- see task 23; do not rename them piecemeal.
- Tooltips say WHY a control exists and what choosing badly costs, not what it is called. Bounded
  controls explain their bounds. A test fails if any control has none.
- Every public module, class, function and method has a docstring. A test fails otherwise.
- **Results tables go through `AnalysisPanel.results_table`**, which gives them a row action (click a
  row, see the map it is about) and a right-click menu (save as CSV, copy rows). A table wired by
  hand is a table that silently lacks both.
- Coverage is 100% on the modules that have it, and the way to cover a Qt-thread body is to call it
  directly, never a pragma. Genuinely unreachable branches get deleted.
- Commit messages explain the reasoning and admit what was got wrong. Write them to a file and use
  `git commit -F` -- backticks in `-m` have twice executed shell commands here.
- Bump the version for feature work (currently 0.10.0).

## Where things are

    starplast/app.py             the window, menus, drawing, exports
    starplast/analysis_panel.py  the six analysis tabs
    starplast/objectives.py      what "good structure" means, and the guards
    starplast/validate.py        putting an error rate on an annotation
    starplast/jobs.py            background work, stopping, state
    starplast/lod.py             the three levels of detail
    instructions/INDEX.md        status table
    results/                     every published number, with the script that made it
