# 64.40 · Verify both goals end to end and release

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.18, 64.19, 64.20, 64.21, 64.22, 64.23, 64.24, 64.26, 64.27, 64.28, 64.29, 64.30, 64.31, 64.33, 64.35, 64.36, 64.37, 64.38, 64.39 |
| Estimated remaining engineering time | 6–10 h + checks |
| Existing code to inspect | `tests/`, `scripts/release.py`, `scripts/build_tutorials.py`, `docs/`, `NEXT_SESSION.md` |

## Deliverable

A release evidence matrix for the two goals, all 39 strategies, supported organisms/hosts, every scorecard level and all unresolved biological gaps.

## Controlled scope

Integration/release verification of admitted scope; unimplemented species and unavailable truth are recorded as remaining work, not completed goals.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] From an installed distribution, users browse a source, look up a gene/class/label, inspect applicable cached inferences and meta-results, and open their grounded cards.
- [ ] Offline/cold-cache/missing-pack/import/organism-switch cases pass; memory/storage/runtime budgets and stale-artifact guards are measured.
- [ ] Every strategy has an appropriate explicit test and its evaluated scope, or remains visibly unvalidated with its truth-acquisition item still open.
- [ ] Fresh calibration and current tutorials/docs accompany changed data/methods; required CI and package checks pass before the standing authorized release flow.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
