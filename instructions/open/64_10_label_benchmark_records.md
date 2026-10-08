# 64.10 · Complete the categorical-label benchmark records

Status: OPEN — 35%, 2026-10-08. Native feature-kNN, ortholog-transfer, both native conformal bases and serial random-forest calls verified; remaining adapters/variants and full coverage are pending.

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

Third bounded change: `starplast/conformal_records.py` consumes a verified
train-only native base model, separate calibration/test score populations and
only calibration labels. Native quantiles, per-class thresholds/rare-class overall
fallback, typed set members, native set text and all raw scores are retained.
Unbounded thresholds are explicit status/null, not nonfinite JSON.
`record_scorecards` adds empirical set coverage, mean size, singleton/empty share
and native efficiency at target/class levels, alongside call metrics. Class
score parameters are restricted to each class cohort when deriving its card.

Pilot: `results/conformal_label_pilot_2026_10_08_v2/`, same frozen training state,
2,126 training/569 calibration/560 test genes, alpha 0.1/per-class/kNN fixed
without test-outcome selection. Native function parity is exact on identical
frozen partitions; the released partition chooser is temporarily replaced in
the diagnostic and restored. Every score, set, threshold, card and matched
baseline replays. **908 checks passed**, one optional pdoc module skipped.
Empirical prediction-grade set coverage 0.973214, mean size 24.455357 of 26,
efficiency 0.061786, **zero singleton calls**; 11 rare classes use the overall
fallback. All-training-classes control covers all stored labels with size 26
and efficiency zero. Broad coverage is not a useful-label accuracy claim.
Artifact `a6062f243a4e4fdf4702e9fc9162c438139da54684bf2873948ee1aa04bfa4ad`;
calibration `cc5b0da9c6a5e767376e64cb331aa893bb7bc2d07e52024a5eb018280b7cdc38`.
No independent biology or exchangeability guarantee admitted; no runtime
strategy/data/calibration change.
The canonical v2 also refuses post-load mutations of the base-model payload.
An executed identity-hardening comparison preserves the initial prototype and
verifies identical numerical outcomes with distinct correct code/model lineage.

Pending: missing adapters (`holdout_search`, `multiplex_modules`), full outer coverage, biological source admission
and calibration. This item does not earn a completion tick from one pilot.

Truth-source follow-up **GT-SPATIAL-01**, shared with 65.04:
`results/spatial_truth_source_review_2026_10_08/` verifies seven additional
primary companions and complete spatial marker tables. Toxoplasma's 718 final
markers include all 62 new microscopy outcomes, so those outcomes cannot
independently validate the published classifier trained on them. Plasmodium
Figure 3 provides nine selected IFA outcomes outside both marker fields; all
nine match exact installed gene IDs. Six other attempted targets remain unknown,
not negatives. Coarse microscopy taxonomy, selection bias, context and feature
dependence remain explicit admission gaps. **426 relevant checks passed** and
complete rows/media/mapping replayed. No biological benchmark admitted; the
truth admission gaps and remaining acceptance conditions are unchanged.

## Next frozen partition: logistic conformal, 2026-10-08

Reuse the seed-17 candidate compartment truth, unchanged source SHA and ordered
train/tune/calibration/test cohorts (2,126/572/569/560). Reuse only the verified
training-fitted rank distributions/source exclusions from the kNN pilot; fit the
native balanced logistic model on training labels alone. Fix C=0.5, alpha=0.1
and per-class thresholds before execution. Retain fitted coefficients/intercepts,
all raw calibration/test scores, sets, singleton abstentions and matching cards.
Replay against the native conformal routine on the identical frozen roles and
compare the all-training-class control. No tuning against outer outcomes, source
addition or biological admission. Output: a new immutable
`results/conformal_logistic_pilot_2026_10_08/` packet. Extend meaningful native
parity/boundary tests to both bases; preserve preceding frozen code packets.

Logistic partition completed: `results/conformal_logistic_pilot_2026_10_08/`.
Native training-only logistic coefficients/intercepts, all calibration/test scores,
quantiles, sets and matched cards replay exactly. Same 560 test genes/26 classes;
empirical prediction-grade coverage 0.905357, mean size 6.975, efficiency 0.761;
zero singleton calls/empty sets. Eleven rare classes use overall fallback.
All-training-class control coverage 1, size 26, efficiency 0. No independent
biological truth or exchangeability promise. 1,088 relevant checks passed.
Artifact `3e007f211b59ef86af8ab1ce4c51bc5526f3c79f53ee307bc87091e2e3acb78e`.
Progress 35%; missing adapters/full outer coverage/biology remain open.

Functional complete-Pfam partition, 2026-10-08:
`results/functional_profile_knn_2026_10_08_v4/` reuses the native feature-kNN
adapter with fixed training-only selection/ranks and source/derived-input closure.
All 635 test genes, 333 unsupported profiles, 348 numeric inputs, 1,413 native
training classes and 2,279 reporting classes remain. Exact native/refit/serialized
replay yields 12 correct, eight wrong and 615 abstentions; majority recovers19
profiles. Independent acceptance verifies47 input/47 output/12 artifact receipts.
124 focused adapter/control/exclusion/card/artifact regressions pass. Annotation
recovery remains separate from biological validation. Earlier categorical,
metadata and resource failures are retained; UI packaging and remaining adapters
are pending. Progress stays35%; no full-action completion.

Native random-forest call partition,2026-10-08: forest_records.py invokes the
existing forest with catalogue defaults300trees/leaf2, training-only ranks and
serial n_jobs1 (native4) solely for exact probability addition. Native tree JSON
state, calls/support/scores and unavailable states replay exactly; native statistics
unchanged. Frozen Pf direct EC profile test152genes/four unsupported gives59correct,
93wrong,noabstentions; source grade unresolved. Canonical packet
results/pf_functional_ec_forest_2026_10_08_v2 and independent acceptance retain
source/control/cohort identities. 101 combined focused checks pass; resource and
Qt diagnostics preserved. Adapter excludes importance, tuning, calibration,
deployment and independent biological admission; remaining method/variant coverage
still open. 64.10 stays35%; no whole-action tick.
