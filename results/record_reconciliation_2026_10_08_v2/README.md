# Reviewed scorecard record reconciliation

All 1,038,372 shipped evaluation rows reconcile across 4,112 distinct setting/seed/hold-out-mode/named-set cohorts. There are zero standard-metric or legacy-denominator mismatches. Accuracy among calls and all-hidden accuracy differ in 910 cohorts; 1,852 cohorts have abstentions. Distinct organism/gene addresses total 13,860, not an independent-study count.

`cards.json` retains per-cohort counts, standard metrics, confusion, scope and missing-lineage/uncertainty limits. `cohort_reconciliation.csv`, `scope_coverage.csv`, `summary.json`, manifest and executed notebooks retain the census and exact verification. A separate four-entity/class fixture verifies false-positive denominators and abstention semantics. Existing truth, exclusion and nested-fit/calibration lineage remain unresolved. No algorithm, deployed prediction, numeric data or calibration artifact changed. Future UI and task-specific adapters must use the shared path with verified biological scope.

Superseded by `../record_reconciliation_2026_10_08_v3/`, which adds explicit outcome-row and evaluation-parameter hashes. Original implementation is retained under `code/`.
