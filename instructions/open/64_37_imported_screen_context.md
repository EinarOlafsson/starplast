# 64.37 · Make user data reach the same inference flow

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.03, 64.04, 64.05, 64.08, 64.23, 64.25, 64.30 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `starplast/importer.py`, `starplast/app.py`, `starplast/strategies.py`, `starplast/guided.py`, `starplast/search.py` |

## Deliverable

Imported measurements rebuild the correct organism context, declare units/evidence family and invalidate affected caches.

## Controlled scope

One end-to-end import fixture for each currently supported parasite/host ID scheme; reuse instruction 60's import repair.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] A real fixture import changes the strategy inputs and its inference/scorecard context in the same session.
- [ ] Imported labels are excluded when used as evaluation truth; missing provenance/leakage declarations produce an explicit unresolved validation state.
- [ ] Shared shipped data remain immutable; save/load retains the exact input hash, organism, transformations and local result artifacts.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
