# 64.04 · Register the ground-truth test cases

Status: DONE — ✅, 2026-10-07. Candidate registry and explicit biological gaps verified; no independent biological benchmarks admitted.

Parent: [64 · Information space and inference atlas](../open/64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.01, 64.03 |
| Estimated remaining engineering time | 0 h |
| Existing code to inspect | `starplast/scorecard.py`, `starplast/calibration.py`, `starplast/track_record.py`, `starplast/datasets.py`, `starplast/search.py` |

## Deliverable

Versioned benchmark entries for each organism/target/task, declaring the truth source, evidence grade, measured population, context, exclusions and available negatives.

## Controlled scope

Register current truth and explicit gaps; data collection for a gap is one deposit-specific follow-up.

Frozen pilot (2026-10-07): census installed parasite categorical/numeric targets and
host protein measurements, with source-level evidence classifications and immutable
eligibility masks. Cover all 39 strategy tasks in each of the four requested species.
No fitting, downloads, synthetic biological truth, or runtime measurement changes.
Publisher-classified LOPIT labels are predictions from experimental profiles; orthology
transfers and derived labels remain separate. Source lineage and assay-population gaps
block independent biological validation even when an installed value exists.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] Reviewed biological targets distinguish direct experiment, curation, transfer, prediction, derived quantity and synthetic control; unreviewed grades remain explicit gaps.
- [x] A test entry has a defined evaluation unit and eligibility mask; missing measurements are never manufactured negatives.
- [x] Each strategy's supported task has a benchmark reference or an explicit unresolved ground-truth gap; synthetic success alone cannot close that gap.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Schema 1 in `starplast/ground_truth.py` pins truth snapshots, source IDs, entity
universes, ordered eligibility, evidence-grade basis, context, units, exclusions,
negative semantics and admission gaps. Duplicate addresses, cross-species/task
references, altered files, reordered cohorts and invalid masks are refused.

Executed census: `results/ground_truth_registry_2026_10_07_v2/`, 529 target
addresses, 646 candidate entries, 156 strategy/species addresses (39 × four).
50 fields classified direct experiment, 18 derived, 15 predicted and six transfer;
440 remain unreviewed with explicit gaps. Native LOPIT labels are classifier
outputs, cell-cycle inference and fitness composites are derived, and both
berghei fitness arms are transfers. Synthetic fixtures certify software only.

**Zero benchmarks admitted.** Actual assayed populations, source lineage,
independent validation and many measurement units remain unresolved. Candidate
eligibility is stored observational availability, not biological accuracy or
assay capacity. No negative class is admitted automatically. Relationship and
replication strategies have explicit pair/discovery-validation gaps; host entries
are protein candidates, not admitted host gene contexts. The 39 legacy ambiguous
rhoptry mappings block admission pending 66.04. Source-specific collection/review
continues under 65.04/66; split protocols under 64.05. No runtime data, strategies,
existing benchmarks or claims changed.

The initial diagnostic used an overly broad assay-kind shortcut; it was reviewed
and superseded before admission. Its original code/hash and failure explanation
remain alongside the snapshot. Reviewed checks: **436 passed**, including registry,
provenance, scorecard, track-record, docstring and organism checks. Executed
verification validates every census input/output hash and Python 3.10 syntax.
Commit subject: **Register biological truth candidates and explicit validation gaps**,
published to `origin/nightly`. Continue at 64.05.
