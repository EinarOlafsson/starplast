# 64.11 · Keep held-out numeric predictions

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.05, 64.06, 64.07, 64.08, 64.09 |
| Estimated remaining engineering time | 6–10 h |
| Existing code to inspect | `starplast/strategy_catalog.py`, `starplast/strategy_learning.py`, `starplast/scorecard.py`, `starplast/track_record.py` |

## Deliverable

Row-level numeric evaluation records for regression, imputation, condition shifts, numeric ortholog transfer and conformal values.

## Controlled scope

Numeric-task adapter and its declared pilot target per mechanism; no new measurement ingestion.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Mask real observed values before any dependent preprocessing; group/cell/dataset masks state their intended generalization.
- [ ] Cards expose units, MAE/RMSE, skill against the matched train-only baseline, correlation and coverage; prediction intervals expose width and achieved coverage.
- [ ] Derived shifts preserve both source conditions and their exclusions; out-of-range and unobserved values are not scored as truth.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
