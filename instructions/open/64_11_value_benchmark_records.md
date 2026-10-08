# 64.11 · Keep held-out numeric predictions

Status: OPEN — 25%, 2026-10-08. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.05, 64.06, 64.07, 64.08, 64.09 |
| Estimated remaining engineering time | 6–10 h |
| Existing code to inspect | `starplast/strategy_catalog.py`, `starplast/strategy_learning.py`, `starplast/scorecard.py`, `starplast/track_record.py` |

## Deliverable

Row-level numeric evaluation records for regression, imputation, condition shifts, numeric ortholog transfer and conformal values.

## Controlled scope

Numeric-task adapter and its declared pilot target per mechanism; no new measurement ingestion.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Mask real observed values before any dependent preprocessing; group/cell/dataset masks state their intended generalization.
- [ ] Cards expose units, MAE/RMSE, skill against the matched train-only baseline, correlation and coverage; prediction intervals expose width and achieved coverage.
- [ ] Derived shifts preserve both source conditions and their exclusions; out-of-range and unobserved values are not scored as truth.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.

## Frozen prerequisite partition: numeric error and matched-baseline cards

2026-10-08, before numeric pilot fitting. Correct the standard metric path so
constant predictions/truth and small answered cohorts retain identifiable MAE/
RMSE while undefined correlations/normalizations remain unavailable. Add raw-unit
RMSE to the metric glossary. Row cards optionally consume an explicit per-row
training-only baseline prediction, named in parameters; compare MAE/MSE skill
on exactly the strategy's finite answered rows, retain full-cohort baseline
errors separately, and refuse incomplete baselines or contradictory abstention.
No fitted adapter, source addition or biological admission in this partition.
Use analytically known errors, constant baselines, abstentions and zero-error
controls as ground truth; verify current calibration display and legacy records.
Historical frozen calibration values are not rewritten as newly recomputed scores.

Prerequisite validated: `results/numeric_metric_review_2026_10_08/`. Constant
baselines retain MAE/RMSE; matched answered-row baseline skill and explicit
abstention guards are implemented. Analytic known-truth notebook replayed;
544 final checks pass, including calibration/display/legacy record compatibility.
Historical calibration/benchmark snapshots remain untouched. Fitted adapters and
real-data evaluation are next; 10%, no completion tick.

## Frozen first adapter: native ridge / HFF fitness candidate

2026-10-08. Target `fit_invitro_hff` from the existing registered candidate;
7,325 eligible stored observations, seed 17, nested homology groups, inductive
feature access. Preserve direct-experiment grade and every unresolved source/
units/mapping/admission gap; this pilot does not promote independent biology.
Training-only source-family exclusions plus native own-kind family exclusions;
training-only feature observation/variance selection and average-rank/zero
imputation. Native ridge alpha=1 fixed, no tune/calibration/test labels fit
anything. Retain every test prediction/abstention and original score units as
unresolved; matched training mean/median errors/skill. Replay training matrix
against native Context.matrix and test outputs against native _fit_predict.
Output a new immutable `results/ridge_value_pilot_2026_10_08/` packet, or a new
diagnostic version if a guard refuses it. No new source or other numeric
mechanism in this partition.

First ridge execution exposed constant-float std rounding: arbitrary baseline
decile fields were emitted despite constant predictions. Native predictions and
training state replay exactly. Preserve the first packet as diagnostic; use
strict nonzero value range for variation detection, add a 1,101-row regression
control, and write canonical corrected evidence to
`results/ridge_value_pilot_2026_10_08_v2/`. No data/model algorithm change.

Native ridge partition validated: canonical
`results/ridge_value_pilot_2026_10_08_v2/`, 7,325 eligible genes; frozen
train/tune/calibration/test 4,015/1,106/1,103/1,101, 331 features. All training
ranks and all 1,101 native predictions replay exactly. MAE 1.269871, RMSE
1.590402, Spearman 0.701225, R-squared 0.474863, coverage 1; training-mean
matched MAE/MSE skill 0.335439/0.474868. Units remain unresolved and the
direct-experiment-origin candidate is not independently admitted biology.
Constant-float baseline ranking diagnostic retained; exact revision comparison
proves unchanged model/prediction/source values. 1,123 final checks passed.
Artifact `782d9e3745b28552b84f9a97022f270bfcb228e49ce37d1e2da5b1c5af08f1b8`.
25%, no completion tick: other mechanisms/variants, intervals, full coverage
and source/units/biological admission are still pending.
