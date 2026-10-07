# 64.02 · Inventory the available information space

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.01 |
| Estimated remaining engineering time | 3–5 h |
| Existing code to inspect | `starplast/datasets.py`, `starplast/slots.py`, `starplast/slot_tree.py`, `starplast/organisms.py` |

## Deliverable

A generated inventory of dataset × organism × unit × biological context × evidence family, with availability and coverage.

## Controlled scope

Inventory existing source records and installed tables; acquisition becomes a separate source-specific task.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Every registered source appears once per declared species/unit, including unavailable and rejected sources.
- [ ] Coverage denominators and missingness distinguish unmeasured, unmapped, inaccessible and explicit negative measurements.
- [ ] Inventory row counts reconcile with the registry and installed tables for the existing parasite and host reference data.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
