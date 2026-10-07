# 64.06 · Standardize baselines and negative controls

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.04, 64.05 |
| Estimated remaining engineering time | 3–5 h |
| Existing code to inspect | `starplast/scorecard.py`, `starplast/strategies.py`, `starplast/methods.py`, `starplast/calibration.py` |

## Deliverable

Baseline/control recipes by task, run on the same eligible units and split as each strategy.

## Controlled scope

Shared baseline/control definitions and fixtures only.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Label baselines include majority/prevalence; value baselines include train-only central tendency; ranking/retrieval baselines use matched candidate universes.
- [ ] Network controls preserve the declared degree/detection constraints; positive-unlabeled tests state what can be scored without verified negatives.
- [ ] Planted positives detect harness failure and null fixtures expose optimistic leakage; all baseline denominators match the tested cohort.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
