# 54 · Scorecards, techniques, and five advanced strategies

**Status: open. The code, the tests and the measurements on the shipped tables are done
(2026-09-26, nightly d69c01b and later). Still to do: the full calibration sweep that fills every
scorecard with intervals.**

The user, 2026-09-26:

* "strategies should have the method they use in their name (UMAP) for example at the end"
* "each strategy should clearly state what techniques it uses and what its score card metrics are.
  these should be as standardized as possible with the same metrics organized in the same way
  whenever possible. metrics should also be explained. (precision, recall, AUPRC, F1, and so on
  there are many more metrics that should be added than these)"
* "add more strategies if you can think of any more"

## STATE

* **Method in the name.** `Strategy.method` is required at registration, and `Strategy.name` is
  "title (method)". The tab shows it and elides from the middle, so the method stays visible. The
  filter matches it.
* **Techniques.** `starplast/techniques.py` defines 40 techniques, each with what it does and why a
  strategy uses it. Every strategy lists its own `techniques`. A test checks that every technique is
  used and every one used is defined.
* **Scorecards.** `starplast/scorecard.py` defines six tasks and 51 metrics. Every metric carries
  its definition, range, chance level and how to read it.

  | Task | Metrics |
  |---|---|
  | label calls | 11: accuracy, coverage, precision of calls, macro precision/recall/F1, weighted F1, kappa, MCC, macro AUROC/AUPRC |
  | ranking | 11: AUROC, AUPRC, lift, prevalence, partial AUROC, R-precision, top-1% precision and enrichment, top-10% recall, best F1, nDCG |
  | set retrieval | 7 |
  | cluster recovery | 8: weighted F1/precision/recall, ARI, NMI, homogeneity, completeness, unclustered share |
  | values | 9: Spearman, Pearson, Kendall, R-squared, NRMSE, MAE, top/bottom decile recall, coverage |
  | replication | 5 |

  * Every metric was checked against scikit-learn or scipy (`tests/test_scorecard.py`).
  * Every self-test, including the 16 custom ones, reports its declared task's full card, in order,
    on the same hidden genes as its verdict. The planted test checks this for all 39 strategies.
  * Label predictors pass per-class scores through `with_class_scores`, so macro AUROC/AUPRC are
    measured, not estimated. Where a strategy has no per-class scores, they are reported missing.
* **The Strategies tab.**
  * The Guide shows the method, each technique explained, and the scorecard table: each metric's
    chance level and reading, with the shipped values and intervals once calibrated.
  * Results lead with the test's card, and each metric explains itself on hover.
* **Calibration.** `calibration.scorecard_cell` bootstraps each metric over targets, then runs.
  `publish` writes the README scorecard tables and `docs/scorecards.md`, which covers both
  organisms, default and tuned settings, and both glossaries.
* **API.** `strategies.metrics()`, `strategies.techniques()`, `Strategy.techniques_table()`,
  `Strategy.scorecard_table()`, `TestResult.card()` and `TestResult.skill`. `docs/API.md` gains
  the strategies section; every example in it was executed against the shipped Tg table.
* **Strategies 35-39** are a ninth family, "Advanced models", in `starplast/strategy_learning.py`.
  All ten self-tests (5 strategies × 2 organisms) pass on the shipped tables.

  | # | Strategy | *T. gondii* | *P. falciparum* |
  |---|---|---|---|
  | 35 | Conformal calls | coverage 0.911 (promised 0.90); set efficiency 0.757 vs null 0.251; 24 classes → mean set 6.6 | 0.924; efficiency 0.394 vs 0.151 |
  | 36 | Graph convolution | 0.491 correct vs chance 0.083; MCC 0.446 | 0.586 vs 0.346 |
  | 37 | Random forest | 0.497; MCC 0.438 | 0.648 |
  | 38 | Stacking | 0.465 | 0.610 |
  | 39 | Conformal values (fitness in HFF) | Spearman 0.726 | 0.493 |

  * Same test, same metrics, Tg compartment: kNN MCC 0.284 → logistic 0.416 → graph convolution
    0.446 / random forest 0.438 (tutorial 4's table).
  * Graph convolution puts 37% of its weight on the network neighbourhood.
* **Tutorials.** Tutorial 4 gains the scorecard and the comparison table. Tutorial 6 (Advanced
  models) is new. The complete guide describes the method, techniques and scorecard, and its
  strategy count now comes from the catalogue (it read "Thirty-two").

## Tested and NOT shipped: choosing the genes to characterise next

The strategy: uncertainty sampling on a logistic classifier, optionally weighted by typicality.
Simulated on the Tg compartment: 25% held out as a test set; the start is 5% or 20% of the rest;
the rest is the pool. The random rows are 5 draws.

| start | budget | random | uncertainty | uncertainty × density |
|---|---|---|---|---|
| 143 | 50 | 0.325 ± 0.011 | 0.332 | 0.328 |
| 143 | 200 | 0.368 ± 0.012 | 0.367 | 0.384 |
| 143 | 500 | 0.405 ± 0.009 | 0.402 | 0.406 |
| 571 | 50 | 0.407 ± 0.003 | 0.420 | 0.400 |
| 571 | 200 | 0.415 ± 0.006 | **0.443** | 0.419 |
| 571 | 500 | 0.426 ± 0.004 | 0.428 | 0.442 |

* At the default (571 / 100, density-weighted) it read 0.415 against random 0.409 ± 0.009: FAIL.
  On Pf it was −0.011 against random.
* The planted table is at 100% accuracy before any label is added, so it cannot pass there either.
* **Verdict: not reliably better than random, so not shipped.**
* Worth retrying only with a model whose uncertainty is calibrated (for example the conformal set
  size from strategy 35), or in a genuinely label-poor slot.

## WHAT TO DO

1. After the data-audit-b and UI branches merge, rerun `scripts/strategy_selftests.py` for all
   strategies on the merged tables.
2. Run the calibration sweep from a frozen snapshot (6,900 runs; see the starplast memory note for
   the recipe). Then `--publish`, which fills the README scorecard tables and `docs/scorecards.md`.
3. Rebuild every tutorial (`QT_QPA_PLATFORM=offscreen python scripts/build_tutorials.py`), then read
   the screenshots.

## HOW TO KNOW IT WORKED

* The README scorecard tables show a mean with a [95% CI] for every strategy (no dashes, except
  metrics a strategy structurally lacks).
* The Guide's "Shipped data" column is filled.
* `tests/test_scorecard.py`, the planted scorecard assertions in `tests/test_strategies.py`, and the
  full suite all pass.

## TRAPS

* pandas deep-copies `Series.attrs` on every operation. The per-class scores therefore ride in a
  holder whose `__deepcopy__` returns itself, or large frames would be copied constantly.
* Per-class (Mondrian) conformal thresholds are infinite for a class with fewer than
  ceil(1/alpha) − 1 calibration genes, so that class would enter every set and nothing could be
  called. Such classes fall back to the overall threshold. The thresholds table lists how many
  labelled genes each class had.
* A conformal set that holds every label always keeps its promise. Coverage alone is not a test,
  which is why efficiency is the verdict.
* Families must hold consecutive numbers (a test enforces it). New strategies that do not extend
  the last family need a new family.
