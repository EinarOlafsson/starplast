# 64.14 · Test map and module recovery

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.05, 64.06, 64.07, 64.08, 64.09 |
| Estimated remaining engineering time | 6–10 h |
| Existing code to inspect | `starplast/strategy_catalog.py`, `starplast/embedding.py`, `starplast/clustering.py`, `starplast/scorecard.py` |

## Deliverable

Persistent cluster/module recovery records for map searches, consensus modules and multiplex modules.

## Controlled scope

Cluster-recovery task and three current mechanisms; no renderer work.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Maps/module selection use permitted training evidence; settings search selects on inner data and evaluates on an untouched outer set.
- [ ] Cards show held-out recovery, fraction placed/clustered, noise, fragmentation/reach controls and class-level confusion.
- [ ] Geometric quality/stability and biological recovery have distinct metrics; withheld labels are not used to choose the reported winner.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
