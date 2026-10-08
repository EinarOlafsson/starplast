# Native logistic conformal pilot, 2026-10-08

Fixed seed-17 Toxoplasma stored-compartment candidate, unchanged ordered source
and exclusion/split lineage. Train/tune/calibration/test: 2,126/572/569/560;
26 native classes. Training-only rank distributions are reused from the verified
kNN packet; the balanced native logistic classifier fits only training labels
with fixed C=0.5/max_iter=500. Alpha=0.1 and per-class thresholds were fixed
before execution. Calibration labels fit only quantiles, never the base model.
The serialized fitted coefficients/intercepts/iteration counts are retained.

All 560 outer-test rows, all calibration/test scores, sets, native set text,
thresholds, class cards and matched train-only baselines are preserved. Native
replay verifies exact scores, fitted coefficients/intercepts, set/call outputs,
thresholds and scorecard aggregates on identical frozen partitions; only the
native partition chooser is temporarily replaced and restored.

Empirical prediction-grade set coverage: 0.905357;
mean set size: 6.975000 of 26;
native efficiency: 0.761000.
No singleton calls, no empty sets; all 560 genes abstain from a singleton answer.
Eleven rare classes use the overall threshold, without a demonstrated class
specific guarantee. All-training-classes control coverage
1.000000, size 26.000000,
efficiency 0.000000. Matched majority surrogate
accuracy 0.230357; it answers a different output
question and is not a prediction-set coverage control.

These labels are stored predictions, not independently admitted biological
truth. Empirical coverage, class discrimination or smaller sets do not establish
biological accuracy, exchangeability or gene-level correctness probabilities.
The logistic pilot was fixed in the action card before execution; comparisons
with the earlier kNN results are descriptive, not outer-test model selection.
No installed data, runtime strategy or calibration algorithm changed. Missing
holdout-search/multiplex adapters, full outer coverage and biological truth
admission remain open under 64.10.

Validation: 466 focused checks and 622 additional disjoint-module checks passed
(1,088 total), including both-base native parity, split/model/calibration guards,
artifacts, scorecards, ground truth/provenance/baselines/capabilities and generated
source scripts. Executed pilot/verification notebooks and original code snapshots
are retained. Packet hashes are recorded in `packet_sha256.json`.

Base model: `6e6c82563b6e51b3eb57256f701ea816272533ee4a1087fbd586aaf63f740247`.
Calibration: `411dd9920add4dcb8f347021ad0a7c8c01c8fe497b52fafffe54799324ef8c1c`.
Held-out artifact: `3e007f211b59ef86af8ab1ce4c51bc5526f3c79f53ee307bc87091e2e3acb78e`.
