# 64.27 · Precompute the Plasmodium inference atlas

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.07, 64.17, 64.25, 64.26 |
| Estimated remaining engineering time | 4–8 h + compute |
| Existing code to inspect | `starplast/organisms.py`, `starplast/strategies.py`, `starplast/track_record.py`, `starplast/claims.py`, `starplast/packs.py` |

## Deliverable

A Pf atlas using its own targets, graph, available evidence, split/calibration artifacts and fixed settings manifest.

## Controlled scope

One parasite species and a frozen pilot manifest; expand targets in a later bounded partition.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] The same atlas admission/parity checks as Tg pass using Pf's own gene and target spaces.
- [ ] Missing evidence/ground truth stays explicit; Tg metrics, defaults or models cannot silently supply Pf results.
- [ ] Build runtime, RAM, size, coverage and incomplete expensive queries are recorded in the manifest.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
