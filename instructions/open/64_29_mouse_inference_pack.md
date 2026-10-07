# 64.29 · Precompute and validate mouse inferences

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.10, 64.11, 64.12, 64.13, 64.14, 64.15, 64.16, 64.22, 64.23, 64.25 |
| Estimated remaining engineering time | 6–12 h + compute |
| Existing code to inspect | `starplast/organisms.py`, `starplast/strategies.py`, `starplast/calibration.py`, `starplast/packs.py` |

## Deliverable

An Mm inference pack with its own measured target/context reliability for a declared initial evidence set.

## Controlled scope

One host and a fixed pilot target set; additional contexts are distinct partitions.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Applicable adapters run on mouse data and tests; human/parasite performance is never reused as measured mouse performance.
- [ ] Tissue, infection/context, missingness and mapping differences remain visible in results and scorecards.
- [ ] Pack integrity, cached/fresh parity and cold-cache organism-switch checks pass for the selected pilot scope.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
