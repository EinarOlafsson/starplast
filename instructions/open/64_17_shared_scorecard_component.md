# 64.17 · Expose one reusable scorecard component

Status: OPEN — 50%, 2026-10-08. Shared immutable presentation, functional routes and a representative host evidence view delivered; full task/calibration acceptance remains open.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.09, 64.16 |
| Estimated remaining engineering time | 2–4 h |
| Existing code to inspect | `starplast/strategy_card.py`, `starplast/strategy_panel.py`, `starplast/scorecard.py`, `starplast/app.py` |

## Deliverable

A compact card showing result quality, coverage, uncertainty, baseline, benchmark population and freshness, with an expandable full test record.

## Controlled scope

Reusable card and one representative host view; later items wire other entry points.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] The same selected cohort/target/organism/settings produces the same displayed values across views.
- [ ] A score opens its definition, ground-truth source, split, sample sizes, control, failure examples and row-level outcomes.
- [ ] Evidence-quality cards, performance cards and single-gene outcomes have appropriate labels; missing or untested states remain visible.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

`scorecard_view.py` supplies immutable source/metric/detail/scope models, compact and
expanded HTML, definitions and exact supplied-record exports. `scorecard_browser.py`
dispatches validated definitions and registered row/outcome/source links. All six
task types have synthetic presentation checks; evidence/performance/individual
outcome labels remain distinct. Missing source, uncertainty, calibration and
deployment stay unavailable rather than becoming confidence.

Discoveries uses the same component for verified functional label/strategy/control
cards, scoped coverage and the reviewed human source foundation. The host card has
58,988 canonical source genes and preserved mapping ambiguity, context/redistribution
and admission gaps; it has no performance metrics and registers no host gene space.
Complete-profile and overlapping-major-class cards retain separate member-only and
whole-cohort precision populations. Native set metric metadata is task-aware and
does not change Task.metrics or any frozen algorithm output schema.

Remaining acceptance: full task/admission/calibration routes depend on unfinished
64.10–64.16; representative evidence presentation is not a host inference benchmark.
Record canonical UI snapshot/check evidence with the validated nightly commit.
No full-action tick is claimed.

Validation: 738 focused source/target/coverage/presentation/feature-runner/card/UI/
organism/docstring checks pass in the final combined run. The source and native
card/outcome replay is executed in the V3 UI audit. The first broad run found
one missing public-class docstring; two later feature fixtures inherited the Qt
process's RSS budget and were made deterministic (an intermediate fixture-name
typo was also corrected). The separate feature/runner suite passes59checks under
a real400MB external cap. Production memory guards were not relaxed.
