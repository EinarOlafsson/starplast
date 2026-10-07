# 64.04 · Register the ground-truth test cases

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.01, 64.03 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `starplast/scorecard.py`, `starplast/calibration.py`, `starplast/track_record.py`, `starplast/datasets.py`, `starplast/search.py` |

## Deliverable

Versioned benchmark entries for each organism/target/task, declaring the truth source, evidence grade, measured population, context, exclusions and available negatives.

## Controlled scope

Register current truth and explicit gaps; data collection for a gap is one deposit-specific follow-up.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Existing biological targets are classified as direct experiment, curation, transfer, prediction, derived quantity or synthetic control.
- [ ] A test entry has a defined evaluation unit and eligibility mask; missing measurements are never manufactured negatives.
- [ ] Each strategy's supported task has a benchmark reference or an explicit unresolved ground-truth gap; synthetic success alone cannot close that gap.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
