# 64.03 · Trace each measurement to its source

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.02 |
| Estimated remaining engineering time | 4–6 h |
| Existing code to inspect | `starplast/datasets.py`, `starplast/deposits.py`, `starplast/identity.py`, `starplast/results.py` |

## Deliverable

A source link for every inventory measurement: paper/deposit, file hash, source species, unit, transformation, mapping and redistribution status.

## Controlled scope

Implement the provenance contract and audit the existing tables; new source acquisition is separately scoped.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Unresolved source records stay visibly unresolved; citation identifiers come from verified records.
- [ ] Gene↔protein and strain/assembly mappings record cardinality, ambiguity and mapping loss.
- [ ] A selected measurement traces through its transforms to the source file; inferred/transferred labels remain distinct from directly measured evidence.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
