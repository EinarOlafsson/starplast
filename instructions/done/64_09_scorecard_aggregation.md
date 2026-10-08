# 64.09 · Make all scorecard levels use the same records

Status: DONE — ✅, 2026-10-08. Shared record aggregation and complete legacy denominator reconciliation verified.

Parent: [64 · Information space and inference atlas](../open/64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.04, 64.06, 64.08 |
| Estimated remaining engineering time | 0 h |
| Existing code to inspect | `starplast/scorecard.py`, `starplast/track_record.py`, `starplast/calibration.py` |

## Deliverable

One aggregation path from individual evaluation records to class, target, method, strategy and organism summaries.

## Controlled scope

Metric/aggregation semantics and reconciliation of existing records; no UI redesign.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] For an identical cohort/settings/split, recomputed metrics match the standard scorecard; existing record-versus-scorecard acceptance is resolved.
- [x] Cards retain eligible, answered, correct, wrong and abstained counts; all-hidden accuracy and accuracy-among-calls have separate labels.
- [x] Small samples, uncertain intervals and unavailable metrics are explicit; repeated seeds/folds are not counted as independent biological samples.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

`starplast/record_scorecards.py` provides the shared standard-metric path for six tasks, class/gene/target/strategy/method/organism views, explicit count/metric availability, and cohort separation. Label cards preserve all-hidden accuracy and accuracy among calls separately. Class precision includes false calls from other classes; recall includes abstentions. Positive-only retrieval/ranking never manufactures negatives or precision. Geometry and numerical predictions have no invented categorical correct/wrong count. Outcome rows and metric parameters have explicit content identities.

Every shipped row is reconciled in `results/record_reconciliation_2026_10_08_v3/`: **1,038,372 outcomes, 4,112 distinct cohorts**, zero standard-metric or legacy-precision mismatches. **910 cohorts** have different all-hidden/called rates; 1,852 have abstentions. **13,860 unique organism/gene addresses** are distinct from evaluation count or independent studies. Settings, seeds, modes and named sets retain separate cards. Folds are disjoint parts of a cohort rather than extra samples. Original diagnostics/code retained; notebooks verify all input/output identities and a class false-positive/abstention fixture.

**835 relevant aggregation, scorecard, artifact, capability, docstring and organism checks passed**. Synthetic tests cover all six metric tasks, all-abstention/nullable/duplicate/inconsistent outcomes, unknown-positive semantics, group uncertainty, score ordering, fractional-cluster refusal and changed row/parameter identities. Python 3.10 syntax verified. Missing biological groups leave intervals unavailable; supplied-group bootstrap is descriptive and cannot establish study independence or an unknown gene's probability. Method views retain component accuracy as unavailable until isolated tests.

The legacy file does not pin independent truth/exclusions/nested fit lineage, and this remains an explicit gap in every cohort. Arithmetic is reconciled without retroactive biological admission or numerical prediction/calibration change. The old display remains a compatibility path; 64.17 consumes the shared records. Task-specific frozen row adapters/truth evidence remain 64.10–64.16. Documentation is API-linked and included in the hosted build. The coherent validated change is committed and pushed to nightly.
