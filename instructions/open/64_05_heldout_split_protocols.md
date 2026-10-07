# 64.05 · Freeze the hold-outs and leakage rules

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.04 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `starplast/strategies.py`, `starplast/track_record.py`, `starplast/search.py`, `starplast/methods.py` |

## Deliverable

Reusable split manifests and exclusion manifests for gene groups, masked values, pairs, datasets and biological contexts.

## Controlled scope

Implement split generation and guard checks; no strategy optimization in this item.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Known-gene predictions hold out the gene and the declared homology group; all preprocessing, tuning and calibration obey the same split.
- [ ] Pair tests separate edge completion on known nodes from cold-node/generalization tests and state the regime.
- [ ] Nested selection has disjoint tuning/calibration/final-test partitions; deliberate target/source/derived-layer contamination is caught.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
