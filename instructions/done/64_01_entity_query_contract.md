# 64.01 · Define the entities and questions

Status: DONE — ✅, 2026-10-07. Schema and explicit-organism alias resolution verified.

Parent: [64 · Information space and inference atlas](../open/64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | None; first contract task |
| Estimated remaining engineering time | 0 h |
| Existing code to inspect | `starplast/organisms.py`, `starplast/slots.py`, `starplast/identity.py`, `starplast/track_record.py` |

## Deliverable

One query schema for organism, gene/protein accession, gene list, label, class value, biological context and requested output type.

## Controlled scope

Schema and resolution contract only; one canonical fixture per entity type.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] Gene aliases resolve within an explicit organism; ambiguous aliases return choices with provenance.
- [x] A label such as compartment and a value such as rhoptry are distinct addresses; multilabel/context-specific entities remain representable.
- [x] One fixture each for a gene, class, label, gene set, continuous trait and pair round-trips without losing its organism or context.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Frozen pilot: seven synthetic query fixtures (gene, host protein, class, label, gene set,
continuous trait and organism-qualified pair). No data acquisition, model fitting or biological
accuracy measurement is needed for this transport/resolution contract.

Implementation: `starplast/query.py`; regression and malformed-input checks in `tests/test_query.py`.
Validation: **472 passed** (2026-10-07, Python 3.12, offscreen Qt, CPU, 4 GB cap), running
`tests/test_query.py`, `tests/test_identity_and_corpus.py`, `tests/test_organisms.py` and
`tests/test_docstrings.py`. Whitespace checks passed. Seven canonical synthetic fixtures cover
all query kinds; additional cases cover multilabels, pair direction, multiple mapping sources,
duplicate/missing/unknown JSON fields, bad types, schema versions and foreign-organism inputs.
The initial 469-check run passed; review added protection against assigning known parasite
accessions to host addresses and three regression cases, then the final 472-check run passed.

The API usage is documented in `docs/API.md`. This card's commit, **Define organism-qualified
queries and provenance-preserving alias resolution**, is pushed to `origin/nightly`.

Limits: the schema is infrastructure, not a strategy execution or biological accuracy claim.
Host gene/protein accession membership still requires a validated source/pack. Generic alias
records are caller-supplied verified mappings; the index adapter only resolves aliases present
in its supplied index. Legacy literature matching remains unchanged. Collisions in the legacy
index lack their original mapping kind, so the adapter exposes them as `ambiguous_alias` with
honest index-artifact provenance, not invented row provenance. Ground-truth and UI integration
remain in their dependent cards. Runtime strategies, shipped claims and calibration are unchanged.
