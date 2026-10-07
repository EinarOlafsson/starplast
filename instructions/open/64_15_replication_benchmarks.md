# 64.15 · Test split-half and conjunction findings

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.05, 64.06, 64.07, 64.08, 64.09 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `starplast/strategy_catalog.py`, `starplast/scorecard.py`, `starplast/discovery.py` |

## Deliverable

Finding-level replication records for blind tests, split clusters and conjunctions.

## Controlled scope

Replication-task adapter for the three existing mechanisms.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Discovery and verification populations/groups are disjoint; finding selection and multiple-testing rules are declared before verification.
- [ ] Every finding records both halves, the original question, eligible sample counts, effect/direction and replication outcome.
- [ ] Cards expose replicated findings / tested findings, effect sizes, uncertainty and matched-null rates; split-half replication is identified separately from an independent dataset test.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
