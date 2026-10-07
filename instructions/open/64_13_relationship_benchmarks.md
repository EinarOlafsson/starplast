# 64.13 · Test relationship and link inference

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.05, 64.06, 64.07, 64.08, 64.09 |
| Estimated remaining engineering time | 6–12 h |
| Existing code to inspect | `starplast/strategy_catalog.py`, `starplast/strategy_graph.py`, `starplast/scorecard.py`, `starplast/methods.py` |

## Deliverable

Held-out pair records for link prediction, attention-corrected relations, unwritten links and learned neighbor/network spaces.

## Controlled scope

Relationship-task adapters and benchmark semantics only; acquisition of physical truth is source-specific.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Candidate pairs have endpoint organism IDs, evidence source, known-edge/cold-node regime and matching negative/control rules.
- [ ] A co-mention recovery score is described as literature recovery; inferred physical links have a separate experimental validation requirement.
- [ ] Ranking cards retain prevalence, AUPRC/lift, precision at fixed depth and dependence exclusions; unknown pairs remain unlabeled where negatives cannot be established.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
