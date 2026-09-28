# 55 · Strategy cards: the same four bars for every strategy, the prose one click away

**Status: complete (2026-09-28, worktree branch off nightly 0.47.0; not merged, tutorials not
rebuilt).**

The user, 2026-09-28:

> it is not very intuitive how to use it. there are now a ton of strategies and a ton of text. I am
> finding it difficult to figure out why or how the strategies worked. I want a minimalistic visual
> style with standardised scorecards showing the same metrics in simple terms up front, then all
> the detailed text one click further. [...] the test for each strategy should come with a text box
> that describes exactly what the strategy does, how it is evaluated, what failing looks like and
> why and what success looks like and why. Best would be if each test actually shows a failed
> result and describes the failure mode and a successful use of the strategy and describes the new
> information gained.

## Plan

1. `scorecard.HEADLINE` / `scorecard.headline()`: four bars per task, in the same positions for
   every strategy -- Better than chance (skill), Reach (coverage or its task's analogue), and two
   task metrics in plain words. Tested.
2. `starplast/strategy_card.py`: painted bars (QPainter), the card, the explainer box and the two
   worked examples. `strategy_panel.py` shows the card first; Guide / Settings / Results move
   behind **Details ▸**.
3. `starplast/strategy_explainers.py`: four plain fields per strategy (what it does, how it is
   evaluated, what failure looks like and why, what success looks like and why), derived from each
   strategy's own explanation and test description. A test requires all four for all 39.
4. `scripts/build_strategy_examples.py` -> `starplast/data/strategy_examples.json`, recorded as an
   executed notebook: one real failure and one real success per strategy per organism, chosen from
   `results/calibration_2026-09-26b/runs.jsonl`; the success is re-run once to list the new calls.
5. `docs/strategy_cards.md`; tutorial captions that name moved controls.

Rules: no algorithm, calibration number or data table changes.

## STATE

* **Headline bars** (`scorecard.HEADLINE`, `REACH`, `headline_bars`). Reach analogues where a
  task has no coverage: clustering = 1 - unclustered share; ranking = recall in the top 10%
  (chance 0.1); set retrieval = genes returned (count); replication = findings made (count).
  Chance ticks come from a fixed level, the card's own metric (replication's null rate), the
  verdict's measured null when the bar IS the verdict metric, or precision/recall over fold
  enrichment for sets. Accuracy and macro F1 have no tick unless measured.
* **Card** (`starplast/strategy_card.py`): painted bars (value, 95% whisker, chance marker), pills
  for method/task/grade, About box (4 fields, one elided line each, click to open), two mini
  scorecards, new-information rows whose gene ids select the gene. Styled through a card-local
  stylesheet by object name, because the app stylesheet beats palettes and fonts. Screenshots
  checked at 620 px dark and 440 px light.
* **Explainers**: 39 x 4 fields, drafted from each strategy's explanation/test description; every
  number in them was checked to occur in the strategy's own text.
* **Examples**: success = best PASS at default, on the strategy's default target where one passed;
  failure = at default, the target that failed most often, median-skill run. Real-data failures:
  Tg 30, Pf 24; noise-table failures: Tg 9, Pf 15. No success (never passed): Tg split_clusters;
  Pf split_clusters, conjunctions. No `run` was skipped (longest 334 s; the 300 s alarm cannot
  interrupt C code). neighbour_space's top 500 pairs are all already measured, so it lists none.
  **Finding:** structural_homology PASSES its self-test on the noise table (0.21 vs 0.16), a
  false alarm the card now says; worth a look.
* Decided not to: show the self-test verdict as a grade when a strategy is uncalibrated (the pill
  says "not calibrated" and the basis line names the verdict).
