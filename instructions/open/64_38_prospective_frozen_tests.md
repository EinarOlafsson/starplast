# 64.38 · Freeze predictions and score new independent truth

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.03, 64.04, 64.05, 64.08, 64.16, 64.32, 64.34 |
| Estimated remaining engineering time | 6–10 h |
| Existing code to inspect | `starplast/claims.py`, `starplast/datasets.py`, `starplast/deposits.py`, `starplast/scorecard.py`, `scripts/add_deposits.py` |

## Deliverable

Version/hash-frozen claims and meta-inferences plus an evaluator triggered by newly admitted compatible experimental truth.

## Controlled scope

Freeze/evaluate infrastructure and one historical time-ordered or genuinely later dataset pilot with provenance.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Predictions, thresholds and eligible populations are frozen before seeing the new dataset; old/new dataset identity and reuse exclusions are checked.
- [ ] Reports retain tested/correct/wrong/unmapped/unmeasured counts, class/context performance and precision/coverage at frozen thresholds.
- [ ] A later dataset scores archived predictions without fitting them again; its results become visible on the same result/method/meta cards.
- [ ] Until real later truth is available, the feature may pass fixture tests but the biological prospective validation remains open.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
