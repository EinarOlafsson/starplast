# 64.08 · Version inference and benchmark artifacts

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.03, 64.05, 64.07 |
| Estimated remaining engineering time | 4–6 h |
| Existing code to inspect | `starplast/strategies.py`, `starplast/results.py`, `starplast/searches.py`, `starplast/packs.py`, `starplast/paths.py` |

## Deliverable

A common manifest for outputs, row-level hold-outs, models and scorecards keyed by organism, table/graph/label hashes, code version, settings, evidence exclusions and seed.

## Controlled scope

Schema, writer/reader and invalidation contract for one small artifact of each output task.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Known-gene held-out outputs and unknown-gene deployment outputs have distinct roles and partition IDs.
- [ ] Changed inputs invalidate dependent artifacts; loaders refuse mismatched gene order, organism, target or model/split identity.
- [ ] Saved artifacts round-trip with provenance, status, support/calibrated confidence and evaluation scope intact.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
