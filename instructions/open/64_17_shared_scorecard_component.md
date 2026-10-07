# 64.17 · Expose one reusable scorecard component

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.09, 64.16 |
| Estimated remaining engineering time | 4–6 h |
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

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
