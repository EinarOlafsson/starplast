# 64.07 · Declare what each of the 39 strategies answers

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.01, 64.04 |
| Estimated remaining engineering time | 4–6 h |
| Existing code to inspect | `starplast/strategies.py`, `starplast/strategy_catalog.py`, `starplast/strategy_graph.py`, `starplast/strategy_learning.py` |

## Deliverable

An explicit capability map: input entities, output task, required evidence, parameter dependence, benchmark and precomputable components for all current strategies.

## Controlled scope

Catalogue declarations and adapter interface; no additional inference algorithms.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] The map covers all 39 catalogue keys and their underlying techniques; new catalogue entries require a declaration.
- [ ] Every query/strategy pair resolves to an applicable output or a specific reason it cannot answer.
- [ ] Each underlying technique has an explicit validation role: a direct ground-truth prediction test, a measured pipeline contribution/ablation, or known-truth/null controls for transformations and statistical operations; composite accuracy is not attributed to untested components.
- [ ] Maps, rankings, sets, pair scores, numeric estimates and label calls preserve their meanings; adapters cannot invent a gene-label probability from a geometric score.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
