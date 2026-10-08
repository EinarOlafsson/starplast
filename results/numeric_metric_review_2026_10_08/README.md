# Numeric error/baseline metric prerequisite, 2026-10-08

Constant predictions previously suppressed even identifiable MAE, so training
mean/median controls lacked error cards. The standard metric authority now
computes MAE and raw-unit RMSE whenever at least one finite prediction has truth.
Undefined correlations/rankings/normalizations remain unknown. No answers yields
zero coverage and unknown error. Raw-unit RMSE is registered in the glossary.

Numeric row cards optionally consume explicit per-row baseline predictions and
a named control. Full eligible and model-answered cohort errors remain separate;
MAE/MSE skill compares identical answered rows. Perfect zero-error baselines or
no answers leave relative skill unavailable. Incomplete/infinite baselines,
infinite predictions and contradictory abstention flags are refused. Exact
outcome hashes include the baseline and parameter identities include its name.
Training-only fit/source lineage must still be supplied by the caller/artifact.

The executed notebook checks analytically known fixture values, constant
baselines, one answer and no answers. Example answered-row MAE/MSE skill is 0.5;
full-cohort baseline MAE is 100/3, while matched MAE is 1. These numerical fixtures
are not biological accuracy/admission. No data, fitted model, inference algorithm
or frozen historical calibration/benchmark record was rewritten. Older records
retain their original metric code SHA; new numeric artifacts bind updated code.

Validation: final combined run passes 544 checks, including error/skill edge
cases, tamper identities, baseline, calibration/publication display, strategy
cards/panel, track records and claims. Two existing library warnings. Initial
focused run: 434 passes and one test failure due to a wrong hash-field name;
the corrected test and all final checks pass. Counts overlap, not additive.
Annotated verification, exact controls and code snapshots are frozen here.
64.11 remains open: fitted numeric adapters, real-data pilots, intervals,
full outer coverage and independent biological admission remain pending.
