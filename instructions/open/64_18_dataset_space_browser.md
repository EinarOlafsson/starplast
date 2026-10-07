# 64.18 · Browse the selected organism's datasets

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.02, 64.03, 64.17 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `starplast/slot_tree.py`, `starplast/datasets.py`, `starplast/app.py`, `starplast/organisms.py` |

## Deliverable

A browser over dataset families, biological labels, contexts, coverage and source records for the selected space.

## Controlled scope

Existing datasets/slots in one browser route; no new sources or extra top-level dock by default.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Filters by organism/host, unit, family, measured context and availability reproduce inventory counts.
- [ ] Each source opens its coverage/provenance card and measured entities; unavailable evidence explains its recorded gap.
- [ ] Measured, transferred, predicted and user-imported evidence remain identifiable; source links work from a compact summary.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
