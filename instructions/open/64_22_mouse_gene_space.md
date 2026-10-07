# 64.22 · Build the mouse gene information space

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.01, 64.02, 64.03, 64.04, 64.05, 64.08 |
| Estimated remaining engineering time | 8–16 h |
| Existing code to inspect | `starplast/organisms.py`, `starplast/host.py`, `starplast/spaces/`, `starplast/packs.py`, `starplast/datasets.py` |

## Deliverable

One reproducible Mm gene-space builder using already acquired, verified reference and selected measurement sources.

## Controlled scope

One host species and a named pilot source set; share builder contracts with item 21 while keeping species data separate.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Apply explicit gene/protein mappings and source-specific mouse tissue/condition semantics with measured mapping loss.
- [ ] A notebook reproduces the build and source sanity checks; unique IDs, non-loss and graph-order checks pass.
- [ ] A validated distributable mouse pack opens independently with its own available targets and benchmark gaps.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
