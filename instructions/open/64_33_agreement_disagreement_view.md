# 64.33 · Show calibrated agreement and conflicts

Status: OPEN — 20%, 2026-10-08. First descriptive paired scorecard/export accepted; in-app navigation and broader calibrated comparisons remain.

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

### Frozen pilot, 2026-10-08

One descriptive paired source-profile scorecard and inspectable row-level export
for the accepted Pf EC kNN/forest artifacts. Preserve complete profiles and all
152 test genes. Distinguish both abstaining, one answering, agreement correct or
wrong, and conflicting calls. Conditional source recovery uses its stated
denominator; calibrated confidence and biological accuracy remain unavailable.
No learned combination or deployment. App navigation and broader calibrated
comparisons remain follow-ups; this partition cannot complete the action.

### Accepted descriptive partition

[Shared-interface scorecard](../../results/pf_functional_agreement_2026_10_08_v2/scorecard.html)
and exact JSON export retain152 genes and ten distinct outcomes. Conditional source
recovery33/94; correct agreements/all eligible33/152; agreement94/143 joint calls;
shared wrong78/143. Empty denominators are unavailable. Original method settings/cards
remain separate; no confidence, biology, ensemble selection or uncertainty admitted.
Mutated population/truth/groups/support, context, training roles, feature access and
fabricated outcomes/counts are refused.61 checks pass; independent notebook recounts
every count/rate exactly. See64.32 for receipts and failures. This is a separate
export, without a new Discoveries comparison button. Nightly history records the
partition commit/push; no whole-action tick.
