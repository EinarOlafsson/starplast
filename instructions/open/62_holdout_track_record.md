# 62 · The track record: hold out genes, and show where every strategy is right and wrong

**Status: open (designed 2026-10-03).** The user's request, verbatim:

> "have tests for each strategy and analysis type where we hold out genes and the test is to see if
> the strategy recapitulates the known label or information on that gene or set of genes. Wherever
> possible tests should hold out as many genes together or individually and the scorecards should
> clearly show how often the strategy got the information right [and] when it got the information
> wrong. It should be possible to click your way to see this on every possible level, starting at
> the gene level then on for each label in each information class. This is a massive amount of work,
> but if structured correctly would be extremely useful to have!"

And the standing principle, to apply everywhere, not only here:

> "the information shown to the user should be as condensed and intuitive and short as possible, but
> always with the ability to click your way to deeper and longer descriptions."

## What exists, and what is missing

| Exists | Missing |
|---|---|
| Self-tests hide 25% of a label once and score the recovery | Every labelled gene held out, so each has an outcome |
| Scorecards: 6 tasks, 51 metrics, per run | Per-GENE right/wrong, kept and browsable |
| Calibration over settings × targets × seeds, with intervals | Per-class and per-gene breakdown surfaced anywhere |
| `TestResult.details`: per-class counts for one run | A record that survives the run and can be clicked into |
| `per_target` skill in the shipped calibration | It is read by nothing (`guided.py` uses only `grade`) |

So the arithmetic exists; what is missing is **coverage** (every gene, not a quarter) and **a kept,
browsable record** at every level.

## The design

### 1. Coverage: k-fold by orthogroup, not one 25% split

Replace the single split with **k folds over the labelled genes, grouped by orthogroup**, so every
labelled gene is held out exactly once per (strategy, target, seed) and no gene is ever predicted by
a model that saw it or its paralog. Five folds by default; the existing `Context.groups()` already
gives the orthogroup, and `hide()` already hides whole groups.

Two hold-out modes, both offered:

* **Together (the default):** the k-fold above. Honest, complete, and one pass costs about what five
  of today's self-tests cost.
* **Alone (on demand):** leave ONE gene out and ask every applicable strategy about it. This is the
  question a biologist actually asks — "could the software have told me this gene is a rhoptry
  protein?" — and it is only affordable per gene, on request, cached afterwards.

### 2. The ledger: one row per held-out gene

`starplast/track_record.py` writes a parquet table, one row per (organism, strategy, target, setting
key, seed, fold, gene):

| field | meaning |
|---|---|
| `gene` | row index into the organism's node table |
| `strategy`, `target`, `setting_key`, `seed`, `fold` | what produced the row |
| `truth`, `prediction` | label codes (int16); `prediction` is null when the strategy abstained |
| `correct` | bool, null when it abstained |
| `support` | the strategy's own score for its call (vote share, probability, field strength) |
| `mode` | "together" or "alone" |

For value strategies the same shape with `truth_value`, `predicted_value`, `abs_error` and
`percentile_error` instead of the label columns.

About 12 bytes a row. One organism × one target × ~15 label strategies ≈ 122k rows ≈ 1.5 MB.

**What ships:** the default target of each organism (`compartment`, `lopit_pf_location`) across every
applicable strategy. **What is computed on demand** into `paths.user_cache_dir()`: any other target,
and every "alone" run. The wheel is 65 MB of a 100 MB limit, so nothing large may ship; if the
shipped ledger ever exceeds ~5 MB it moves to a downloadable pack (the pack framework exists).

### 3. Aggregation: the same numbers at four levels

`track_record.summary(level, ...)` returns one tidy frame per level, each row carrying `n`, `right`,
`wrong`, `abstained`, the rate, and its **Wilson 95% interval** (already in `calibration.wilson`),
plus the chance level for that cell:

1. **Gene** — "14 strategies tried this gene's compartment; 9 right, 3 wrong, 2 abstained."
2. **Class** (one value of one label, e.g. *rhoptries 1*) — per strategy and pooled: recall, the
   classes it is confused WITH (the confusion row, which is the "when it got it wrong" the request
   asks for), and how learnable the class is at all.
3. **Information class** (one label column, e.g. *compartment*, *fitness class*) — per strategy, and
   the best strategy for that column.
4. **Strategy** — what the scorecard already shows, now backed by the per-gene rows, plus "where it
   fails": the classes and genes it gets wrong most.

Rules that keep it honest: never a bare percentage without its interval and `n`; say "too few to
judge" below `n = 5` rather than printing 100%; mark a cell *circular* when the target's closure
could not be fully removed (the gallery already computes this); and keep abstentions visible, since a
strategy that answers 10% of genes perfectly is not better than one that answers all of them well.

### 4. The UI: condensed first, click deeper, at every level

A new **Track record** view, and entry points from where the question arises:

* **Gene card (evidence panel):** one line — "Known compartment: dense granules. Strategies recover
  it 9 times in 14." Click → the per-strategy table for that gene; click a strategy → that run's
  scorecard; click the class → the class level.
* **Strategy card:** one line under the four bars — "Right on 62% of held-out genes; weakest on
  *apical 1* (2 of 11)." Click → the classes it fails, then the genes.
* **A browser** (its own tab or a page of the Strategies tab): information classes → classes →
  genes, each row condensed to a rate with its interval, expanding on click. Sortable by "most
  learnable", "least learnable", "most disagreed about".
* Every number has the same hover sentence it has elsewhere (`scorecard.explain`).

### 5. Cost, and what is affordable

From the calibration timings: a cheap label strategy (kNN, logistic, network vote) runs in seconds on
the shipped table; the expensive ones (multiplex modules ~96 s, map walks ~38 s) are minutes per
fold. So:

* **Full k-fold, shipped:** the ~15 cheap label strategies × default target × 5 folds ≈ 20 minutes
  per organism.
* **The expensive strategies:** keep the existing single split, and say so in the UI rather than
  pretending the coverage is complete.
* **"Alone" runs:** one gene × the cheap strategies ≈ seconds, computed when asked.

## The work, in stages

| Stage | What | Acceptance |
|---|---|---|
| 1 | `track_record.py`: k-fold runner, ledger schema, writer/reader, tests on the planted table | Every labelled gene appears exactly once per strategy/target/seed; no gene predicted by a model that saw its orthogroup |
| 2 | `summary(level=...)` for the four levels, with Wilson intervals, chance levels, abstentions and the confusion row | Numbers agree with the existing scorecard when pooled to the strategy level |
| 3 | Build and ship the default-target ledger for both organisms; notebook; size check | Ledger ≤ 5 MB; wheel stays under 75 MB |
| 4 | UI: gene card line + drill-down; strategy card line; the browser | Every level reachable by clicking; every cell has `n` and an interval; tooltips everywhere |
| 5 | "Alone" mode on demand, cached in the user cache | One gene, every cheap strategy, in seconds; cached result reused |
| 6 | Wire into Start here ("this label is recovered well/badly for your genes") and the questions page | — |

### Status (2026-10-03)

* Stages 1-3 done: `starplast/track_record.py`, `scripts/build_track_record.py`,
  `starplast/data/track_record.parquet` (190,624 rows, 1.4 MB), set hold-outs (each class, and
  random sets of 1-500), `notebooks/track_record_2026_10_03.ipynb`. Stage 2's "agrees with the
  scorecard when pooled" is NOT yet checked.
* Stage 4 done except a standalone browser: gene card -> class page -> category page -> back, and the
  strategy card line links to the category page. These HTML pages in the evidence panel ARE the
  browser for now; a separate widget is only worth building if they prove too small.
* Stage 5 partly: `track_record.my_list(ctx, genes, target)` hides a user list together (Python
  only; no UI, no cache yet). `track_record.alone(ctx, gene, target)` hides one gene with its
  orthogroup for ANY label, five fast strategies, ~9 s on Tg, cached per release/organism/label/gene
  in the user cache. The gene card links each of the gene's other labels to it (`starplast://alone/
  <row>/<target>`), run on the job runner; the answer replaces the panel only if the user is
  still on that gene. The user list runs from Start here (*Test on my genes*) -> grid in the
  evidence panel. Stage 5 done.
* Stage 6 partly: Start here recommendations append `track_record.record_phrase` (rate + the
  commonest-class baseline) when the record covers the label in hand. Not yet: the questions page.
* FOLLOW-UP: layer_propagation is graded "reliable" by calibration (its per-layer test) but is
  BELOW the commonest-class baseline on held-out compartment (15% vs 18%). `guided.recommend`
  now lists below-baseline strategies last (`track_record.beats_baseline`). Still open: the
  calibration verdict for label diffusion should be re-examined.
* Found on the way: categorical `value_counts` named zero-count classes as confusions (fixed, with
  a regression test). Whole-class hold-outs score `right == 0` by construction; only `together`
  and `placed_at` mean anything there, and the notebook says so.

## Traps

* **Never pool across settings** when reporting a gene's record: a gene called right at one setting
  and wrong at another is not "50% right", it is a gene whose answer depends on the setting, and that
  is the interesting fact. Report per setting, default setting first.
* **Abstention is not a wrong answer**, and must not be averaged into one.
* **A class with 3 genes** cannot have a meaningful rate; show the count and say so.
* The existing self-tests must keep working unchanged — this is an addition, not a replacement.
