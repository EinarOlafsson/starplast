# Baselines and controls

`starplast.baselines.RECIPES` declares shared baselines, null controls,
denominators and limits for all six scorecard tasks. Helpers provide training-only
majority/prevalence calls, training mean/median predictions, uniform rankings on
a fixed candidate universe, known-positive retrieval and bounded network controls.
Task-specific adapters still need to execute and record these recipes alongside
each strategy. The declarations do not certify existing self-test denominators.

```python
from starplast.baselines import value_baselines, positive_only_retrieval

training_values = truth.loc[list(split.entities("train"))]
predictions = value_baselines(training_values, split)
assert tuple(predictions.index) == split.entities("test")

recovery = positive_only_retrieval(candidate_ranking, known_positives, depth=20)
print(recovery["known_positive_recall"], recovery["precision"])
```

Training truth must cover the complete eligible training cohort, contain no
unknown/nonfinite values and pass the split's fit guard. Missing truth is never
imputed. Majority and class prevalence come only from training labels; hidden
class prevalence cannot influence the baseline. Numerical mean/median likewise
use training values only. Each result records training/evaluation counts, split
identity and the exact ordered universe hash. No alternate denominator is used
when a strategy abstains; later scorecards must distinguish all-eligible accuracy
from accuracy among calls.

Ranking and fixed-size retrieval use the same frozen candidates, query-seed
exclusions and fixed depth as the strategy. The random-ranking helper rejects
seeds that remain in the candidate universe. Its first `k` entries supply uniform
fixed-size retrieval. Candidate order, count, seed and seed exclusions are retained.

Positive-only truth identifies recovery of **known** positives. Unlabelled
candidates may be positive too, so the helper leaves precision, AUROC and AUPRC
unavailable. The observed-positive fraction among returned items is recorded with
that name; it is not presented as full biological precision. Known-positive recall
does not identify recall over all biological positives.

Network controls use undirected, unweighted edge swaps that preserve each node's
degree and contacts with each declared detection stratum. They require explicit
strata for every endpoint, retain the attempt/success counts and never turn rewired
pairs into biological negatives. Weighted or directed data requires a separate
recipe. Sparse/complete graphs can permit no swaps; that outcome remains
`unavailable_no_valid_swaps`. Successful swaps do not establish a uniformly mixed
null, so mixing always stays `not_established` pending a dedicated review.

The executed fixtures in `results/baseline_controls_2026_10_07/` retain training
truth, matched predictions, the split, recipes and checksummed inputs/outputs.
Known planted positives pass a recovery check; reversed scores and a random null
fail that harness. These are synthetic software tests, not independent biological
benchmarks. Real pair negatives, grouped null execution, map/replication adapters
and strategy-specific baseline comparisons remain in subsequent benchmark items.
