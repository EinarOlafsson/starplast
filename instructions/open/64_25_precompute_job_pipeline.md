# 64.25 · Build a resumable precomputation pipeline

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.08, 64.10, 64.11, 64.12, 64.13, 64.14, 64.15, 64.16 |
| Estimated remaining engineering time | 6–10 h |
| Existing code to inspect | `scripts/build_track_record.py`, `scripts/build_claims.py`, `starplast/searches.py`, `starplast/packs.py` |

## Deliverable

A dependency-aware job plan for reusable features/operators/maps/models, hold-out records, deployment outputs and scorecard summaries.

## Controlled scope

Pipeline plus tiny fixtures and one bounded pilot; large sweeps begin only with a measured budget and RAM lease.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Immutable input snapshots and job fingerprints prevent mixed-source runs; interrupted jobs resume by validated completed partitions.
- [ ] A changed dataset/target/settings fingerprint invalidates only declared downstream outputs; cancellation leaves valid earlier partitions.
- [ ] Pilot jobs report runtime/RAM/storage, set a compute budget and package budget, and preserve audit logs and failures.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
