# 64.30 · Show every strategy's answer for a gene

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.17, 64.19, 64.26, 64.27 |
| Estimated remaining engineering time | 6–10 h |
| Existing code to inspect | `starplast/app.py`, `starplast/strategy_panel.py`, `starplast/track_record.py`, `starplast/claims.py` |

## Deliverable

One gene inference profile listing each strategy's applicable outputs, status, supporting evidence, calibrated reliability and conflicting calls.

## Controlled scope

One unified gene profile over the existing evidence page; use available Tg/Pf packs first, then validated host packs.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] All 39 catalogue entries are accounted for; unsupported, abstained, failed, stale and pending answers remain inspectable.
- [ ] A known gene's recovery view uses its held-out prediction, and an unknown gene's claim exposes the matching class/stratum test population.
- [ ] Cached results appear without recomputation; new query-dependent work has visible progress/cancel and attaches its own manifest.
- [ ] One click from each answer reaches its scorecard and deeper method/evidence description.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
