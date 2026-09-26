# 49 — A Strategies tab: thirty-two ways to infer something, each testing itself

**Status: done 2026-09-25 (0.44.0, on nightly).**

## The request

> I want this to be a useful tool, so I want it to be easy to use it for inference. The user
> should be able to (1) hold out a category, do a hyperparameter search and find a structure that
> maps onto the holdout -- the first idea for this dataset; (2) provide a list of genes and see if
> the hyperparameter walk finds a structure with high precision and recall for those genes and
> one cluster. I want a new section, to the right of Evidence and Analysis, called Strategies:
> at least 30 strategies with tooltips, detailed explanations and walkthroughs -- different ways
> to use this combined data to do inference and generate knowledge. Each strategy should have a
> test function that uses known information and hold-out information to test whether its
> inference mechanism works.

## What exists now

| where | what |
|---|---|
| `starplast/strategies.py` | `Context` (the table, graph, other organism, caches and the leakage guard), `Strategy` / `Param` / `StrategyResult` / `TestResult`, the inference primitives, the five test patterns, `planted_context` |
| `starplast/strategy_catalog.py` | the 32 strategies: runner, tester, tooltip, explanation, walkthrough, test description, parameters |
| `starplast/strategy_panel.py` | the tab: list by family with the shipped verdict, Guide / Settings / Results, Run / Test / Stop, gene-list input (paste, file, example, gate), Show on map |
| `app.Window._strategies` | the dock, tabified after Analysis, wired to `use_embedding`, `use_clusters` and gene selection |
| `scripts/strategy_selftests.py` | measures every strategy on both shipped tables; writes `starplast/data/strategy_selftests.json`, `results/strategies_<date>/` and `docs/strategies.md` |
| `tests/test_strategies.py`, `test_strategy_panel.py`, `test_strategy_edges.py` | catalogue completeness, machinery, planted PASS / null not-PASS for every strategy, the panel, every guard |

Strategies 01 and 02 are the two founding questions, exactly as asked: `holdout_search` and
`geneset_hunt`.

## The five test patterns, and why every verdict means the same thing

Every tester is one of five patterns in `strategies.py`, and every one compares the strategy with
the same procedure run WITHOUT the information it claims to use:

1. **label transfer** -- hide 25% of a label (whole orthogroups at a time, so no gene is recovered
   through a visible paralog), predict it back; null: shuffled visible labels;
2. **set expansion** -- hide 30% of a gene set, rank every gene from the rest; null: random sets;
3. **pairs** -- hide 20% of an edge layer, score hidden edges against random non-edges; null:
   gene identities permuted in the evidence;
4. **values** -- hide 20% of a measurement, rank-correlate; null: the exact chance distribution of
   a rank correlation (see below);
5. **replication** -- make findings on half the genes, check them on the other half; null: the
   second half's evidence scrambled.

PASS = above the null's 95th percentile AND by a stated margin. INCONCLUSIVE = too little could be
hidden to score. On 40 null tables the false-pass rate measured 0-4 in 40 per strategy, which is
what a 95th percentile estimated from 20 draws gives.

`planted_context()` builds an organism where every signal is planted and a `null=True` twin where
every column and edge is dealt out at random. The suite requires every strategy to PASS on the
first and to pass on at most one of three null tables -- one, because a 95th-percentile bar lets 5%
of true nulls through by design and demanding zero would fail on a fair coin.

## Measured on the shipped data (default settings, 2026-09-25)

| | PASS | FAIL | INCONCLUSIVE |
|---|---|---|---|
| *T. gondii* | 26 | 5 | 1 |
| *P. falciparum* | 24 | 5 | 3 |

The strongest on *T. gondii*: structural homology (EC class correct for 0.77 of hidden
annotations against 0.25), crosslink partners (0.52 against 0.07), a classifier on the permitted
measurements (0.43 against 0.08), attention-corrected co-mention (top pairs share a compartment
0.75 against 0.39), link prediction (AUROC 0.81). The founding two pass: the hold-out search's
chosen clusters recover hidden compartment members at F1 0.150 against 0.079 for searches made on
shuffled labels, and the gene-list hunt at 0.123 against 0.043.

The FAILs are reported, not dropped: consensus modules and multiplex communities beat their nulls
by less than the 0.05 margin, cluster-enrichment calls are right 8% of the time, the in vivo
fitness shift is predictable only at rho 0.079, and paralog divergence barely separates
paralogs in different compartments (AUROC 0.54). Each is a measured statement about what this data
does not support.

## Things found on the way, each of which would have made a verdict wrong

* **HDBSCAN's default selection is degenerate on the real proteome.** Excess-of-mass returns two to
  six clusters, the largest holding 45-90% of a 2,500-gene map; no label dominates any of them, so
  every cluster-based call was empty. `leaf` returns dozens of small, purer clusters (median purity
  0.39 against 0.20). Walks now try both; single-map strategies expose the choice.
* **A cluster's commonest label is the wrong test of "a structure that maps onto the holdout".**
  With 26 compartments the plurality label collapses to `nucleus - chromatin` under chance too. The
  test is now the search's own logic: choose each label's best cluster on known genes, score how
  well it holds that label's hidden genes (F1).
* **The prefix fallback in `Context.blocks` leaked.** For tables the catalogue does not describe,
  only the banned column was removed from its group, and `cycle_t2` -- associated with the
  held-out phase at 0.7999, just under the 0.8 threshold -- carried the phase into a map built to be
  blind to it. A group with ANY banned column now goes whole, the rule `search` already applies.
* **Knockout screens predicting knockout screens.** "Predict fitness from other evidence" scored
  rho 0.895 because the other evidence included a Δgra17 screen, a second-background screen and the
  oxidative-stress screen, all filed under other block names. Families now come from the slot
  catalogue's axis; leaving out the target's own kind gives 0.721 -- a real signature, read from
  expression, abundance and sequence.
* **`multiplex_communities` chains.** Joining agreed pairs by connected components put all 480
  genes of the planted table in one community. `join="louvain"` runs community detection on the
  agreement graph; the default is unchanged for existing callers.
* **Three refits on shuffled values are not a null.** On a table with nothing in it they averaged
  -0.10 and let a correlation of 0.10 pass. Rank-correlation tests now use the exact chance
  distribution (sd 1/sqrt(n-1)).
* **Two-means splits noise.** A numeric split finding "replicated" whenever a cluster's values were
  as spread as the map's, because two-means puts plain normal noise 1.6 sd apart. A split now has
  to be wider than random sets of the same size manage.
* **A null that never calls anything is the ideal null, not a missing one.** Precision on zero
  calls is undefined, which made the enrichment test inconclusive exactly when its null was
  behaving perfectly. Counted as zero; precision verdicts need at least ten calls instead.

## Not done, and worth knowing

* **`cellcycle_phase`'s closure removes only itself and its pseudotime.** The GSE19092 cell-cycle
  expression columns it was presumably derived from stay in, so a strategy scored against the
  phase can read it back. No strategy defaults to that label for this reason; the closure itself
  (datasets registry `derived_from`) should be fixed.
* **Split clusters (26) is inconclusive on both arms.** Too few clusters are homogeneous in one
  localization for splits to be significant on half a sample; a table with more genes per cluster,
  or a coarser A label, would give it something to test.
* **Multiplex communities take two minutes on *T. gondii*** -- greedy modularity per layer.
* **The planted tests use a 1,500-gene table for 26, 27 and 28**, which need enough genes per half
  or enough labelled paralog pairs; the default 480 is too small for them to be significant.
