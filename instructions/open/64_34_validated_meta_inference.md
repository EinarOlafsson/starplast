# 64.34 · Train and test the meta-inference

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.05, 64.06, 64.08, 64.16, 64.32, 64.33 |
| Estimated remaining engineering time | 8–16 h + compute |
| Existing code to inspect | `starplast/claims.py`, `starplast/strategy_learning.py`, `starplast/scorecard.py`, `starplast/calibration.py` |

## Deliverable

A dependence-aware meta-model per compatible output task/target, initially piloted on one target, then admitted by the same task contract.

## Controlled scope

One target and one output task in the pilot; each additional target/task is a separate measured admission partition.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Train with nested out-of-fold base predictions produced without the outer held-out genes/groups; inner fitting/tuning/calibration cannot see outer truth.
- [ ] An untouched outer test compares meta-inference with the best base selected on inner training, fixed consensus and matched null/baseline.
- [ ] Report per-class/stratum error, calibration, coverage and incremental capacity; superiority is measured on paired/grouped test units.
- [ ] Persist weak/failed cases; an uncalibrated or unsupported meta-model cannot replace a measured method answer.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
