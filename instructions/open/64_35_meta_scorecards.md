# 64.35 · Expose meta-inference scorecards at every level

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.17, 64.30, 64.31, 64.34 |
| Estimated remaining engineering time | 4–6 h |
| Existing code to inspect | `starplast/scorecard.py`, `starplast/strategy_card.py`, `starplast/claims.py`, `starplast/app.py` |

## Deliverable

Gene-result, class, label, meta-method and organism views of the meta-model's validation and base-method comparisons.

## Controlled scope

UI/readers over item 34's tested artifacts; no extra model training.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Every displayed meta confidence links to the applicable test population, calibration and dependence report.
- [ ] Paired base-versus-meta improvement, disagreement cases and abstention/coverage tradeoffs are reachable from the compact card.
- [ ] Known-gene outer-test outcomes and unknown-gene conditional reliability are labeled correctly; global accuracy does not masquerade as class-specific or gene-specific accuracy.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
