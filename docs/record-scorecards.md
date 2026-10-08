# Scorecards from individual evaluation records

`starplast.record_scorecards.aggregate(rows, scope, parameters=...)` is the common
metric path for label calls, numeric estimates, rankings, retrieval, map/module
recovery and replication. `RecordScope` keeps organism, strategy, target,
settings, seed, protocol, partition, benchmark, truth grade, biological unit and
negative semantics explicit. Missing legacy lineage remains a gap.

Every eligible entity appears once within a cohort. Repeated evaluations,
settings, seeds, named sets and hold-out regimes keep separate cards. The
organism overview retains their cards and distinct entity counts; it does not
average unrelated tasks or treat those counts as independent studies.

Label cards distinguish **correct / all eligible hidden entities** from
**correct / answered entities**. They retain eligible, answered, correct, wrong
and abstained counts, coverage, confusion and standard macro metrics. Abstention
is neither a correct nor wrong call. Class precision includes false calls from
other classes, while class recall includes abstentions and wrong calls. A
single gene's recorded outcome is available, but it is not a statistically
estimated probability that the gene's deployment prediction is correct.

Metrics use `starplast.scorecard` for the exact supplied cohort. Unavailable
metrics become null with explicit status; score-free label records cannot
supply AUROC/AUPRC. Positive-only ranking/retrieval can report known-positive
recovery; unknown members cannot become false negatives or verified absence,
and precision/AUROC/AUPRC remain unavailable where not identifiable. Native
numeric units remain declared or unresolved. Cluster placements are distinct
from categorical calls, and replication counts refer to tested findings.

When biological groups are supplied, label cards offer a descriptive group
bootstrap with a fixed seed. Fewer than five groups, missing group identities
and unknown provenance leave intervals unavailable. Grouping does not prove
study independence. Method views link the pipeline card and retain component
accuracy as unavailable until its direct/ablation test is executed.

`legacy_cohorts` separates each stored setting, seed, mode and named set. The
legacy track record's `rate` is accuracy among calls; the standard scorecard's
`accuracy` is all-hidden accuracy. The executed reconciliation at
`results/record_reconciliation_2026_10_08_v3/` checks all **1,038,372 outcomes in
4,112 distinct cohorts**, with zero mismatches under the matching denominator.
The two rates differ in 910 cohorts. There are 13,860 unique organism/gene
addresses across those evaluations, not 1,038,372 independent biological samples.

The legacy file does not pin independent truth, exclusions or nested
fit/calibration lineage. Arithmetic reconciliation does not establish biological
admission or calibrated confidence. The existing legacy display remains a
compatibility path; the shared component under 64.17 will consume these explicit
cards. Task-specific frozen row adapters and biological evidence follow under
64.10–64.16.

Each card also pins its full outcome rows and evaluation parameters. Class-score matrices retain their hash, columns and row count; changed scores cannot reuse a prior metric identity.
