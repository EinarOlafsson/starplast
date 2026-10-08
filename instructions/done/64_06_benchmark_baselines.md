# 64.06 · Standardize baselines and negative controls

Status: DONE — ✅, 2026-10-07. Shared recipes/helpers and synthetic controls verified; per-strategy benchmark execution remains in subsequent cards.

Parent: [64 · Information space and inference atlas](../open/64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.04, 64.05 |
| Estimated remaining engineering time | 0 h |
| Existing code to inspect | `starplast/scorecard.py`, `starplast/strategies.py`, `starplast/methods.py`, `starplast/calibration.py` |

## Deliverable

Baseline/control recipes by task, run on the same eligible units and split as each strategy.

## Controlled scope

Shared baseline/control definitions and fixtures only.

Frozen scope: training-only majority/prevalence and mean/median recipes, matched
ranking universes, positive-only scoring limits, degree/detection-stratum
preserving undirected network controls with explicit mixing gaps, and planted/null
fixtures. No runtime strategy changes, fitting on held-out labels, invented pair
negatives or new biological performance claims. Per-strategy execution belongs
to the subsequent task-specific benchmark adapters.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] Label baselines include majority/prevalence; value baselines include train-only central tendency; ranking/retrieval baselines use matched candidate universes.
- [x] Network controls preserve the declared degree/detection constraints; positive-unlabeled tests state what can be scored without verified negatives.
- [x] Planted fixtures detect reversed/null harness failures; split guards reject optimistic held-out fitting; baseline denominators match the tested cohort.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

`starplast/baselines.py` declares recipes for all six tasks. Training-only majority,
prevalence, mean and median helpers enforce the complete eligible training cohort
and exact evaluation IDs. Uniform rankings/fixed-size retrieval preserve the
candidate universe, seed and query-seed exclusions. Positive-only scoring reports
known-positive recovery; precision/AUROC/AUPRC remain unavailable without negatives.

Undirected/unweighted network swaps preserve node degrees and per-node contact
counts with each declared detection stratum. Attempt/success counts and missing
mixing validation are explicit; complete graphs return an unavailable control.
Rewired pairs are not manufactured biological negatives.

Executed synthetic evidence: `results/baseline_controls_2026_10_07/`. Planted
positives recover completely; reversed and random rankings fail the harness.
Training-only numerical/categorical baselines share the same frozen test universe.
Network fixture performs 20 degree/detection-preserving swaps; its mixing remains
unestablished. Complete graph performs zero valid swaps and remains unavailable.
426 relevant checks passed, including baselines, split guards, truth registry,
scorecards, organism and docstrings. Executed verification checks every input/output
hash, matched denominator and Python 3.10 syntax. Commit subject:
**Standardize matched baselines and explicit null-control limits**, on `origin/nightly`.

Scope is shared definitions and fixtures only. No runtime strategy changes,
independent biological admission or new accuracy/capacity claim. Grouped nulls,
map/replication adapters, real pair negatives and per-strategy comparisons remain
under 64.10–64.15; network mixing validation belongs to relationship benchmarking.
Continue with 64.07; source recovery remains the user's concurrent priority.
