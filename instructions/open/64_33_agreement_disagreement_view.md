# 64.33 · Show calibrated agreement and conflicts

Status: OPEN — 50%, 2026-10-08. Paired rates, group intervals and class/profile routes accepted; broader calibrated comparisons remain.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.30, 64.31, 64.32 |
| Estimated remaining engineering time | 3–6 h |
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

### Next bounded partition, 2026-10-08

Expose the accepted Pf paired card in Functional tests, with both calls and all
152 genes visible. Ship the exact pinned report; revalidate its artifact identities,
derived counts and current source context. Preserve individual method views/exports.
Gene clicks and individual paired outcomes must retain the frozen test scope;
missing/corrupt reports or changed contexts clear stale results. No new fit,
calibration, independent truth or ensemble selection. Verify actual desktop routes.

### Accepted paired desktop partition

Discoveries → Functional tests → Comparison: kNN and random forest opens the
externally pinned report from either Pf strategy. All 152 genes/four unsupported
genes retain both calls, named-method outcomes and source profiles. Rates show
fractions/percentages; exact numeric exports are unchanged. A gene click retains
its entity; double-click opens one observed paired outcome without a confidence
probability. Returning to a method restores its original columns/cards/controls.
Missing or corrupt reports and changed source/organism contexts clear stale cards,
rows and exports. Every functional view rechecks current context, so returning to
a method cannot borrow its old accuracy after source changes. The shipped report
exactly matches accepted SHA940da900…8e98b. Canonical desktop V3 packet is under
`results/pf_functional_agreement_ui_2026_10_08_v3/`; independent acceptance V2 verifies
54 output receipts/18 current inputs, including the Tg regression packet. Final
68 focused checks pass. Actual Pf desktop7.43s/810,123,264-byte process peak, serial1900MiB
with workers terminal. Initial deleted-dialog test teardown and import-path
diagnostics are preserved. No original node/claim/benchmark bytes, fits, biological
admission or calibrated probabilities changed. Nightly history records commit/push.
64.33 stays open; whole22/51 unchanged. Next: matched group-aware uncertainty.

### Frozen interval display partition, 2026-10-08

Attach the verified interval JSON SHA0bc3e966…061f9 to the existing Pf paired card.
Preserve exact original rates and all152 genes/four unsupported; show descriptive
95% whole-group percentile bounds and shared-group assumptions. Missing/corrupt or
incompatible interval files must show unavailable uncertainty while retaining a
valid original comparison. Current source/report mismatches still clear the entire
comparison. Cohort intervals are never per-gene confidence. Verify exports and the
actual desktop; no fitting/calibration or biological admission.

### Accepted interval display

The existing Pf comparison now shows five original rates with descriptive95%
whole-group percentile bounds and the forest-minus-kNN all-eligible recovery
difference (12.5percentage points,5.3–20.5). Compact uncertainty text gives500
resamples/145 recorded groups; assumptions and exact values remain in the export.
The pinned interval file exactly matches the numerical packet SHA0bc3e966…061f9.
Missing/corrupt intervals retain a valid original comparison with unavailable
uncertainty; parent/report/source mismatches remain refused. All152 genes/four
unsupported remain, and individual gene outcomes never receive cohort probabilities.
48 focused checks pass, then18 final presentation checks. Actual desktop V2 and
inspected screenshot: `results/functional_group_uncertainty_ui_2026_10_08_v2/`,
6.88s/864,231,424-byte process peak under serial1900MiB, workers terminal.
Independent acceptance verifies48 output receipts/20 current inputs and unchanged
parent/benchmark bundle bytes. Initial display-only empty-reason wording corrected;
no failed runs, fits, source changes or biological/calibration admission.
Nightly history records commit/push; whole22/51 unchanged. Next: bounded class
comparison navigation with exact cohort and annotation-negative semantics.

### Frozen class navigation partition, 2026-10-08

Compare recorded major-class membership or exact complete-profile membership on
the same152 known-source genes. Map a spoken complete-profile call to membership;
abstentions stay unknown. Recorded nonmembership is never biological absence.
Keep all genes for metrics; displayed rows include reference members and either
method's positive calls, retaining false positives. Separate positive and negative
agreements; show native per-method class cards alongside explicitly counted
source precision/recall. Do not reuse pooled profile intervals for this different
question; class intervals/calibration remain unavailable. Preserve current-context
refusals, source identities, exports and gene/individual-outcome navigation.

### Accepted class/profile navigation partition

Both Pf methods expose seven major-class and twelve complete-profile comparisons.
Metrics/exports retain all152 genes/four full-profile-unsupported genes; displayed
rows include reference members or either positive call, preserving false positives.
Recorded nonmembership is not biological absence; abstentions stay unknown.
Class outcomes retain original profiles/outcomes. Native method class cards remain
inspectable; precision/recall show exact denominators. Positive/nonmembership
agreements are separate. Per-gene outcomes have no cohort probability; pooled
profile intervals are not applied to class questions. Paired class uncertainty,
calibration, biology and consensus controls remain unavailable.

19 focused checks and9 final checks pass. Actual desktop packet:
`results/functional_pair_class_ui_2026_10_08/`,9.49s/948,936,704-byte process peak,
serial1900MiB, workers terminal; screenshot inspected. Independent acceptance V2
recounts all19 addresses exactly from original profiles, verifies21 output receipts/
21 current inputs and unchanged source/parent/interval bytes. See
`results/functional_pair_class_acceptance_2026_10_08_v2/`. Preserve eight initial
synthetic-fixture failures (missing group_uncertainty metadata) and the initial
acceptance assertion (compact versus expanded overlap metadata); both corrected
without changing source results. Nightly history records commit/push; whole22/51
unchanged. Other methods, hosts, literature and biological validation remain open.
