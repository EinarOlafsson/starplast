# 64.01 · Define the entities and questions

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | None; first contract task |
| Estimated remaining engineering time | 2–4 h |
| Existing code to inspect | `starplast/organisms.py`, `starplast/slots.py`, `starplast/identity.py`, `starplast/track_record.py` |

## Deliverable

One query schema for organism, gene/protein accession, gene list, label, class value, biological context and requested output type.

## Controlled scope

Schema and resolution contract only; one canonical fixture per entity type.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Gene aliases resolve within an explicit organism; ambiguous aliases return choices with provenance.
- [ ] A label such as compartment and a value such as rhoptry are distinct addresses; multilabel/context-specific entities remain representable.
- [ ] One fixture each for a gene, class, label, gene set, continuous trait and pair round-trips without losing its organism or context.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
