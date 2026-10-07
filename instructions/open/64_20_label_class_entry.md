# 64.20 · Make labels and protein classes entry points

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.01, 64.02, 64.17, 64.19 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `starplast/track_record.py`, `starplast/slot_tree.py`, `starplast/app.py`, `starplast/guided.py` |

## Deliverable

Search/browse pages for label variables and their values, including localization, curated protein classes and gene sets.

## Controlled scope

Class/label navigation over declared existing annotations; ontology expansion is a separate registered-source task.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] A class page shows measured membership, definition/context, coverage, per-class precision/recall and confusions for each applicable mechanism.
- [ ] A label page shows class balance, overall/macro performance, unmeasured classes and the hierarchy where applicable.
- [ ] Class→gene→label→source navigation preserves organism and context; identical text in different labels or species is not merged.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
