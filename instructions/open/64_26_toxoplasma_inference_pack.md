# 64.26 · Precompute the Toxoplasma inference atlas

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.07, 64.17, 64.25 |
| Estimated remaining engineering time | 4–8 h + compute |
| Existing code to inspect | `starplast/organisms.py`, `starplast/strategies.py`, `starplast/track_record.py`, `starplast/claims.py`, `starplast/packs.py` |

## Deliverable

A Tg atlas over an explicitly frozen initial target/settings manifest and the reusable outputs each applicable strategy can precompute.

## Controlled scope

One organism and a declared pilot target/settings set; record the remaining target coverage rather than calling the atlas complete.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Each of the 39 strategies has an output or a specific per-query status; gene-label, numeric, relationship and module outputs preserve their roles.
- [ ] Known genes use held-out records; unknown genes use separately generated deployment predictions; provenance and cards are attached.
- [ ] Canonical class/list queries are cached by query identity; arbitrary new lists/parameters reuse base artifacts then run the remaining bounded work on demand.
- [ ] Random sampled cached outputs match a fresh run of the same frozen configuration within declared tolerances and fit-time conditions.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
