# 64.09 · Make all scorecard levels use the same records

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.04, 64.06, 64.08 |
| Estimated remaining engineering time | 4–6 h |
| Existing code to inspect | `starplast/scorecard.py`, `starplast/track_record.py`, `starplast/calibration.py` |

## Deliverable

One aggregation path from individual evaluation records to class, target, method, strategy and organism summaries.

## Controlled scope

Metric/aggregation semantics and reconciliation of existing records; no UI redesign.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] For an identical cohort/settings/split, recomputed metrics match the standard scorecard; existing record-versus-scorecard acceptance is resolved.
- [ ] Cards retain eligible, answered, correct, wrong and abstained counts; all-hidden accuracy and accuracy-among-calls have separate labels.
- [ ] Small samples, uncertain intervals and unavailable metrics are explicit; repeated seeds/folds are not counted as independent biological samples.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
