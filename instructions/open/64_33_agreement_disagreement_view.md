# 64.33 · Show calibrated agreement and conflicts

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.30, 64.31, 64.32 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `starplast/claims.py`, `starplast/app.py`, `starplast/strategy_panel.py`, `starplast/scorecard.py` |

## Deliverable

A gene/class/label comparison view grouping compatible answers and conflicting predictions by evidence family and reliability.

## Controlled scope

Descriptive comparison and matched held-out summaries; the learned meta-model is item 34.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Agreement only compares aligned entity/target/context/units; hierarchical or cross-label compatibility requires an explicit mapping.
- [ ] Cards show how many methods spoke, which families contributed, which disagreed/abstained and observed held-out accuracy conditional on agreement where measurable.
- [ ] Consensus counts alone never manufacture a confidence probability; redundant methods remain inspectable as dependent evidence.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
