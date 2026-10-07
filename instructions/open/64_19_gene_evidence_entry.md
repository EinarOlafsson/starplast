# 64.19 · Make gene lookup an evidence entry point

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.01, 64.03, 64.17, 64.18 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `starplast/identity.py`, `starplast/app.py`, `starplast/organisms.py`, `starplast/slot_tree.py` |

## Deliverable

A gene page with alias lookup, all measured evidence grouped by biological question, source coverage and links to class/label pages.

## Controlled scope

Lookup and evidence navigation for available spaces; inference profiles arrive in item 30.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Canonical and alias searches select the correct organism/entity and expose ambiguity instead of guessing.
- [ ] Every available measured column for the gene is reachable, with context, units and provenance; missing values keep their meaning.
- [ ] The initial view is condensed and can expand to all evidence; a measured label and its class both have working scorecard routes.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
