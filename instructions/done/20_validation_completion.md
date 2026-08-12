# Validation: refit, candidates, orthogonal evidence, per-category scores — DONE 2026-08-12

The tab existed and scored categories; what was missing was everything that makes the score usable
or trustworthy. Four things were asked for. All four are in, and the first of them turned out to be
a repair rather than an addition.

## The guard never fired

`masked_recovery` refused when `category in used_columns` -- a category VALUE against a list of
column names, which is never true. The application made it worse: the tab passed
`columns_for(...).keys()`, the embedding's BLOCK names, so the check was comparing "dense granules"
against "expression_summary". The tab that exists to refuse circular questions would have scored
`compartment` against a map built on `compartment` and reported the number without comment.

Now `circularity_error` tests the label COLUMN against the real column names (categoricals
included), the panel passes those, and the refusal is checked **before** the job is submitted as
well as inside it -- waiting for a job to fail in order to be told the question cannot be asked is
the wrong shape for that answer. It is shown as an explanation with what to do next, not as a red
failed job with a traceback.

Rendered and looked at, which caught the second half: the note said "Not scored" while the previous
run's table of precisions sat underneath it. A refusal now clears the scores and the candidate list,
because a table under an explanation reads as the explanation's result.

## refit was a flag that changed only the sentence

`refit=True` set a field that made `summary()` say "re-embedded per fold" while nothing was
re-embedded. It now requires a `rebuild` callable and uses that fold's labels; asking for it without
one is an error. The panel supplies one that rebuilds the embedding with `seed + 1 + fold` and
re-clusters. `refit` is recorded on every row of the result, because the difference between a fresh
map per fold and one fixed map is invisible in the numbers.

What re-fitting actually buys is worth being precise about: the label never feeds the embedding
either way -- the guard above refuses when it does -- so the gain is that the estimate stops being
conditional on one layout. UMAP moves noticeably between seeds at this size.

## Candidates, and what has to travel with them

Clicking a category opens the genes its cluster would have you annotate. Every row carries
`cluster_frac_category`, `cluster_frac_contradicting`, `overall_frac_category` and `enrichment`.
The last two are new: on the real table the best cluster is 10.2% nucleus-chromatin, which sounds
like something until you see that the map is 9.4% nucleus-chromatin -- enrichment 1.09, i.e.
nothing. The prevalence is computed over the same denominator as the cluster fraction, since the two
are meant to be divided by each other; over labelled genes alone it read 1.3x higher.

## Orthogonal evidence

`orthogonal_support` counts, for each candidate, how many genes already carrying the category share
an **orthogroup**, a **Pfam** or an **InterPro** domain with it -- agreement from evidence the map
never saw. A column that fed the embedding is refused and the refusal logged: its agreement would be
the map agreeing with itself. Zero support is reported as zero, never as a blank, because it is the
common answer and it is the one to be careful about.

## Per-category scores in Inference

`clustering.per_category` takes each category's best cluster and reports precision, prevalence,
lift, recall and F1 -- the same "best single cluster per label" convention `search.score_recovery`
uses, so the two cannot disagree about what recovery means. A feature-level Cramer's V of 0.2
describes 27 compartments weakly smeared across every cluster and one compartment falling out
cleanly, and only this table separates them.

Two guards, both from looking at the rendered table rather than from reasoning about it:

- **lift**, because ranked by F1 the top row was a majority class -- a cluster holding 90% of
  everything "recovers" a 90%-prevalent label at F1 0.95 while saying nothing. Lift 1.0 means the
  cluster is no more that category than the map is.
- **`min_in_cluster=5`**, because ranked by lift the top row was then a cluster holding ONE
  apicoplast protein at 5x enrichment: true, meaningless, and indistinguishable at a glance from a
  result. It is the singleton exploit in per-category form.

## Two bugs found on the way

- **The target was not aligned to the clustered genes.** `labels` comes from an embedding that may
  have dropped rows; the tab scored it against the full label column, so with any spec that drops
  genes it compared gene i's cluster with gene j's compartment. The battery already aligned through
  `self.rows`; validation did not.
- **A job's progress note travelled backwards.** `JobRunner._on_progress` copied the queued
  "started" message onto the job whenever Qt delivered it, overwriting whatever the worker had
  reported since -- a walk that had reached configuration 12 reverted to "started". The note is now
  written where it is produced.

`Job` also keeps the exception itself, not only its text, so a caller can tell a deliberate refusal
from a crash by type rather than by parsing a message.

## Verified

- 1,353 tests pass headless. `validate.py`, `clustering.py`, `analysis_panel.py` and `jobs.py` are
  all at **100%** -- including `_Progress`, whose body runs on a Qt-managed thread and is invisible
  to coverage unless called directly.
- The tooltip audit caught the new checkbox with no explanation, which is what it is for.
- Run end to end under the pandas 3 interpreter on the real 8,140-gene table: the refusal fires, 24
  categories score, the candidate list arrives with its numbers, the per-category table fills.
  **This is how the missing `log` parameter was found** -- every panel test had mocked
  `validate_all`, so the arguments the tab passes had never met the function that receives them, and
  the tab raised TypeError the first time it ran for real. There is now a test that runs the job
  against the real module.
- Rendered under Xvfb and looked at: validation, the refusal, and Inference.

**What the numbers say, on the real table**: a 2-cluster HDBSCAN over the full proteome scores
precision 0.022 at recall 0.978 for its best category. Annotating from it would be almost entirely
wrong. That is the expected result for two giant blobs, and it is exactly why the tab exists.

## What is left, and where it belongs

Saving an annotation -- the store, the free-text reasoning, the separate file, the fourth colour --
is **task 17**, which was blocked on this one and is now unblocked. It must still carry these
numbers with every saved row.
