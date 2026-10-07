# 64.23 · Select and install available organism spaces

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.18, 64.21, 64.22 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `starplast/organisms.py`, `starplast/packs.py`, `starplast/paths.py`, `starplast/app.py` |

## Deliverable

A space selector and pack catalogue for parasites and hosts, with explicit available/installable/unavailable states.

## Controlled scope

Integrate the existing two parasite spaces and validated Hs/Mm packs; additional species get separate builder tasks.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Selecting a space updates evidence, targets, strategies, scorecards and caches atomically without retaining another species' context.
- [ ] Published pack versions/hashes, interrupted install and offline operation follow the existing pack contract.
- [ ] Other instruction-53 species use the same admission checklist one species at a time; an unfinished space is never shown as calibrated or ready.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
