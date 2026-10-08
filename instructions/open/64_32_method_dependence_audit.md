# 64.32 · Measure dependence between inference mechanisms

Status: OPEN — 20%, 2026-10-08. First matched source-profile report accepted; group-aware uncertainty and broader methods remain.

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
