# 64.10 · Complete the categorical-label benchmark records

Status: OPEN — 20%, 2026-10-08. Native feature-kNN and the missing ortholog-transfer adapter verified; remaining label-capable adapters and full coverage are pending.

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
This coherent change is committed/pushed to nightly as `8b59643`.

Second bounded change: `starplast/transfer_records.py` adds the previously
unrecorded categorical `ortholog_transfer` path. The native donor orthogroup
aggregation receives no receiver target column; the source-to-label mapping
fits only on receiver training genes. Both endpoints remain prediction-grade.
Source table/mapping identities, semantic review and dependence gaps are explicit.
The existing Plasmodium schizont spatial prediction is a scoped installed-source
pilot, not a new preferred dataset or an admitted biological comparator.

`results/ortholog_transfer_pilot_2026_10_08/`: same frozen 560-gene test cohort,
123 calls, 437 abstentions, 65 correct/58 wrong. Mapping reaches 967/3,827 receiver
genes, 550 training and 126 test genes; three mapped test values lack a learned
source category. All-hidden surrogate recovery 0.116071, called-only 0.528455,
majority baseline 0.230357. No native class scores/support or AUROC/AUPRC are
invented. Every native projection/call, card and matched baseline replays exactly.
Artifact `65034bf84b7dc87aefee2e21622294cb33f0789b4e9dc41fb81243ee9ec33f07`;
mapping `e1da9a0c47cecdd47087c34026ddc5e5add6d4f950570a81574c481b40654d80`.
**464 relevant checks passed**. Artifact validation now also refuses numeric
outputs as categorical evaluations and categorical outputs as numeric ones.
Initial missing evaluation metadata was refused before any artifact was written;
the executed verification records that diagnostic and corrected scope.

Pending: other missing adapters (`holdout_search`, `multiplex_modules`,
`conformal_calls`), full outer coverage, biological source admission
and calibration. This item does not earn a completion tick from one pilot.
