# 64.32 · Measure dependence between inference mechanisms

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

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
