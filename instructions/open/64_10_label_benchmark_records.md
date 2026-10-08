# 64.10 · Complete the categorical-label benchmark records

Status: OPEN — 10%, 2026-10-08. First frozen native-vote adapter implemented; remaining label-capable adapters and full coverage are pending.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.05, 64.06, 64.07, 64.08, 64.09 |
| Estimated remaining engineering time | 6–12 h |
| Existing code to inspect | `starplast/track_record.py`, `starplast/strategy_catalog.py`, `starplast/strategy_learning.py`, `starplast/scorecard.py` |

## Deliverable

Native held-out label outputs, including class-score or prediction-set fields where produced, for all label-capable strategies in the capability map.

## Controlled scope

Categorical-label task only; one new adapter per validated commit, first on a frozen pilot target.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Reuse the existing 14-strategy record; add missing label-capable adapters one at a time with native-prediction parity.
- [ ] Every eligible labelled gene has a held-out call or an explicit abstention; rare and unsupported classes are retained as evaluation gaps.
- [ ] Precision, recall/F1 by class, coverage, baselines and confusion outcomes derive from those rows; conformal calls report coverage and set size/singleton efficiency.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

First bounded change: `starplast/label_records.py` implements feature kNN with
training-only native rank scaling, frozen source exclusions, retained class
scores and explicit abstention. Training ranks match native `Context.matrix`,
including ties/missing values; query ranks use frozen training distributions.
The whole-context rank transform is deliberately withheld under inductive
evaluation. No installed data or runtime strategy is changed.

Pilot: `results/label_knn_pilot_2026_10_08_v5/`, candidate stored Toxoplasma
compartment calls, seed 17, 3,827 eligible genes, train/tune/calibration/test
2,126/572/569/560. Stored calls have prediction grade; recovery is not biological
accuracy. Every test gene and native class-score value is retained. Earlier
diagnostics and their original code remain available. Full-precision score
hashing fixes an independently discovered rounded-JSON identity weakness.

Verification: all 560 calls and 26 native score columns replay exactly; training
vectors exactly match native `Context.matrix`. 196 correct calls, 236 wrong,
128 abstentions; all-hidden recovery 0.35, called-only recovery 0.453704,
majority baseline 0.230357. Artifact identity
`4924dfa57c5fc39a6bd3f92718e031d60d23163b34da869ef102a5c1ec03e350`.
**1,161 checks passed**, one optional pdoc documentation module skipped locally;
445 focused checks also passed before the broader strategy/track-record suite.
This coherent change is committed/pushed to nightly as the first adapter.

Pending: remaining adapters, full outer coverage, biological source admission
and calibration. This item does not earn a completion tick from one pilot.
