# 64.32 · Measure dependence between inference mechanisms

Status: OPEN — 30%, 2026-10-08. First matched source-profile group intervals accepted; broader methods and lineage coverage remain.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.03, 64.08, 64.10, 64.11, 64.12, 64.13, 64.14, 64.15, 64.16 |
| Estimated remaining engineering time | 6–10 h |
| Existing code to inspect | `starplast/claims.py`, `starplast/strategies.py`, `starplast/scorecard.py`, `starplast/calibration.py` |

## Deliverable

Dependence reports by organism/target/context from evidence lineage plus shared-error/residual behavior on matched held-out units.

## Controlled scope

Dependence measurement for current methods; the completed literature audit is retained as existing evidence.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Methods sharing inputs, labels, pretrained representations or fitted outputs declare that overlap, including derived sources.
- [ ] Error dependence is measured on matched test cohorts with group-aware uncertainty; applicability/coverage intersections are reported.
- [ ] A strategy and its own submethods/ensemble bases cannot count as independent votes; insufficient overlap yields an unresolved dependence status.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

### Frozen pilot, 2026-10-08

Compare accepted Pf direct complete EC kNN and random forest on all 152 identical
outer-test genes, retaining four unsupported genes and abstentions. Revalidate
native artifacts, target/context/truth/group/order/training roles before pairing.
Report shared dependency hashes and wrong-call overlap; different algorithms do
not establish source independence. No fitting, confidence calibration, ensemble
selection or biological admission. Group uncertainty and other method/target
partitions remain separate follow-ups. This pilot cannot complete the action.

### Accepted descriptive partition

Canonical [report](../../results/pf_functional_agreement_2026_10_08_v2/report.json)
and [independent acceptance](../../results/pf_functional_agreement_acceptance_2026_10_08_v2/executed.ipynb):
152 genes/145 groups, four unsupported genes retained;143 joint calls.94 agreements
(33 correct,61 same wrong),49 conflicts (seven kNN-only correct,25 forest-only
correct,17 both wrong),nine forest-only calls (one correct). Both wrong on78 genes.
Shared dependencies/training and exact native input/retained feature-column overlaps
are inspectable. No independence, group-aware interval or biological accuracy claim.
61 focused checks pass;101 source outputs/25 artifact receipts/14 current inputs and
13 final packet outputs verified. Report1.186s under800MiB; checks serial under1900MiB,
workers terminal. Initial800MiB test OOM/renderer failures retained in diagnostics.
Commit/push receipt is the nightly history for this partition; no whole-action tick.

### Frozen uncertainty partition, 2026-10-08

Use the accepted paired report SHA940da900…8e98b, all152 genes/four unsupported
genes and145 recorded groups. Follow the existing record-card convention:
500 whole-group resamples with replacement, seed20261008, at least five groups,
2.5/97.5 percentiles with linear interpolation. Both methods use the same draws.
Estimate the five existing paired rates and right-minus-left all-eligible recovery.
Retain numerator/denominator totals, undefined draws, group selections and exact
replay. Empty denominators/small group populations remain unavailable. Intervals
describe source-recovery sampling variation conditional on recorded groups; no
biological independence, calibration, model selection or new fitting is implied.

### Accepted group-uncertainty partition

Canonical V2 packet `results/functional_group_uncertainty_2026_10_08_v2/` preserves
all152 genes/four unsupported/145 groups. Twelve focused checks pass. Independent
acceptance expands every saved group draw to its original gene rows and reproduces
all500 paired numerators/denominators/ratios and six percentile intervals exactly;
10 output receipts/six current inputs verified. No tolerance comparisons.
Agreement source recovery33/94 (35.1%), descriptive interval25.0–45.8%; both-wrong
78/143 (54.5%),46.1–63.0%; right-minus-left all-eligible recovery19/152 (12.5
percentage points),5.3–20.5 points. All500 denominators are defined in this cohort;
empty denominators/fewer than five groups remain unavailable in fixtures. Shared
draws retain method pairing; saved selections/arrays replay exactly. Generic integer
ratio sums refuse malformed/replayed indices and unsupported count capacity.
The shared scorecard exports these intervals without biological independence or
per-gene calibration. This is a separate artifact; app comparison still uses its
previous report. Root800MiB jobs finish under one second; initial11-check fixture
packet retained. No new fit/source/benchmark promotion. Nightly history records
commit/push. Next: expose pinned intervals in the existing comparison view.
