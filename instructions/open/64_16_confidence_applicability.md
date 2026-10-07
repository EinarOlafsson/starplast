# 64.16 · Calibrate confidence and measure where it applies

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.10, 64.11, 64.12, 64.13, 64.14, 64.15 |
| Estimated remaining engineering time | 6–12 h |
| Existing code to inspect | `starplast/claims.py`, `starplast/calibration.py`, `starplast/scorecard.py`, `starplast/strategies.py` |

## Deliverable

Target/class/context/coverage-stratum calibration with explicit applicability and abstention policies.

## Controlled scope

Confidence/applicability layer over registered tasks; no global universal threshold or fabricated per-gene accuracy.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Fit calibrators on held-out calibration partitions and score them on untouched test partitions; report calibration error/proper scores where meaningful.
- [ ] Report common/rare classes, well-measured/understudied genes, mapping/missingness strata and out-of-distribution checks with actual sample counts.
- [ ] Native scores and calibrated probabilities remain distinguishable; no confidence is claimed outside the tested range; numerical promotion rules are recorded per benchmark.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
