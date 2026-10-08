# 64.07 · Declare what each of the 39 strategies answers

Status: DONE — ✅, 2026-10-08. Catalogue declarations and adapter plans verified; biological component tests remain subsequent benchmark work.

Parent: [64 · Information space and inference atlas](../open/64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.01, 64.04 |
| Estimated remaining engineering time | 0 h |
| Existing code to inspect | `starplast/strategies.py`, `starplast/strategy_catalog.py`, `starplast/strategy_graph.py`, `starplast/strategy_learning.py` |

## Deliverable

An explicit capability map: input entities, output task, required evidence, parameter dependence, benchmark and precomputable components for all current strategies.

## Controlled scope

Catalogue declarations and adapter interface; no additional inference algorithms.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] The map covers all 39 catalogue keys and their underlying techniques; new catalogue entries require a declaration.
- [x] Every query/strategy pair resolves to an applicable output or a specific reason it cannot answer.
- [x] Each underlying technique has an explicit validation role: a direct ground-truth prediction test, a measured pipeline contribution/ablation, or known-truth/null controls for transformations and statistical operations; composite accuracy is not attributed to untested components.
- [x] Maps, rankings, sets, pair scores, numeric estimates and label calls preserve their meanings; adapters cannot invent a gene-label probability from a geometric score.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

`starplast/capabilities.py` declares all 39 catalogue keys and all 40 techniques, with exact query/output/evidence/parameter requirements, benchmark biological units, fit/cache dependencies and test roles. New undeclared catalogue entries/techniques fail closed. The resolver preserves canonical entities, target/seed/context meanings, returns named refusals and never invokes inference or manufactures probabilities. Evidence-family/class summaries and pair rankings retain their actual units; imputation uses transformed percentiles, condition shifts rank residuals, and ortholog output adapts to target type. Candidate benchmark references match effective target, task, organism and unit; independent validation gaps remain explicit.

Frozen artifact: `results/capability_contracts_2026_10_08_v2/`, with executed catalogue and verification notebooks. **624 synthetic query/strategy routes and 78 native plans** verified, including explicit host-adapter refusals. Original diagnostic/code retained. **1,104 broad strategy/contract checks passed**; following final declaration tightening, **855 current capability/query/truth/docstring/organism checks passed**. An initially invalid organism fixture and target-refusal expectation were corrected; no failed test remains. Python 3.10 syntax and all final input/output hashes verified. Documentation is linked through the API and included in the docs build.

This item covers declarations and adapter interface only, without additional algorithms, fitted biological models or numerical data/calibration changes. Technique roles have status `declared_not_measured`; biological validation, direct/ablation/null execution, frozen fitting and benchmark artifacts follow under 64.08–64.16. The coherent validated change is committed and pushed to nightly.
