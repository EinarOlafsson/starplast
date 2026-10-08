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

## First frozen categorical adapter

`starplast.label_records.feature_knn` takes an ordered feature population,
the complete training-label cohort, a frozen split and training-derived source
exclusions. It calls the existing native distance-weighted vote and retains
every outer-test row and class-score value, including scores for abstained calls.
It never receives held-out labels or labels support as calibrated confidence.

Training features reproduce `Context.matrix`'s average-rank percentile transform,
centered at 0.5, with missing values filled by zero. Evaluation values use the
frozen training distribution: observed ties have their average rank; unseen
values use the right empirical CDF. All-missing training columns are withheld.
The original whole-context ranks are inappropriate for the declared inductive
protocol because held-out values would change training ranks. The stored model
state includes ordered training distributions and vectors so this distinction
can be inspected and replayed.

`results/label_knn_pilot_2026_10_08_v5/` contains the executed first pilot and
verification on stored compartment calls. The target is prediction-grade;
its recovery measures agreement with those stored predictions, not independent
biological accuracy. Independent biological truth, other label-capable adapters,
full outer-fold coverage and confidence calibration remain open. Earlier pilot
directories retain their original code and explicitly documented serialization
and preprocessing diagnostics.

Artifact JSON and metric identities preserve native floating-point values;
the default rounded pandas JSON representation is not used. A one-ULP score
change invalidates the card identity even if rounded display metrics coincide.

## Frozen ortholog transfer

`starplast.transfer_records.ortholog_transfer` records the native categorical
transfer path. `TransferSource` pins donor species/column/table, the projection
identity, evidence grade and semantic/dependence review. Donor values already
projected through orthogroups are supplied in frozen receiver order. Only
receiver training labels fit the native source-to-label map. Known donor values
derived from receiver truth are refused; unresolved dependence remains a gap.
Receiver column exclusions cannot be applied merely by name to another species.

Every test gene appears, including missing mappings and source categories without
a learned receiver label. The method produces no native class-score matrix,
support or calibrated confidence; AUROC/AUPRC remain unavailable. The spatial
prediction pilot at `results/ortholog_transfer_pilot_2026_10_08/` makes 123 calls
among 560 test genes: 65 agree and 58 disagree with stored receiver predictions,
with 437 abstentions. Called-only agreement of 52.8% and all-hidden agreement of
11.6% describe different denominators; the same-cohort majority baseline is
23.0%. Both endpoints have prediction grade and unresolved context equivalence.
These are surrogate recovery measurements, not independent biological accuracy.

Row-column order is recorded explicitly for replay because JSON object key
ordering is not a dataframe schema. Adaptive numeric and categorical strategies
must use the matching typed output in benchmark artifacts.

## Frozen conformal label sets

`starplast.conformal_records.conformal_calls` requires a verified, complete
train-only native kNN/logistic fitted-model artifact, separate calibration/test
score tables and calibration labels. It keeps native quantiles and rare-class
overall fallback. Empty/multi-label sets abstain from singleton calling. Set
members are typed lists; native display text is retained separately, so a class
name containing a separator does not redefine set membership. Unbounded
thresholds carry explicit status and null value in finite JSON.

The shared card reports empirical set coverage, mean size, singleton/empty
shares and native efficiency, alongside both accuracy denominators. Efficiency
uses the released definition: empty sets count as size one, while set coverage
still counts them as misses. Class cards use the same model class vocabulary
and their own record/score cohort. A singleton is a set of size one, not a
calibrated probability that its label is correct.

The kNN pilot at `results/conformal_label_pilot_2026_10_08_v2/` retains 560 test genes
and all 26 score columns. Set coverage against stored prediction labels is
97.3%, with average set size 24.46/26 and zero singleton calls. Native efficiency
is 0.0618. An all-training-classes control achieves coverage 100%, size 26 and
efficiency zero; high set coverage alone would conceal weak inference capacity.
Eleven rare classes use the overall threshold. Independent biological truth,
exchangeability, the logistic variant and full outer coverage remain open.

## Numeric errors and matched training baselines

Numeric MAE and raw-unit RMSE remain available for constant predictions, constant
truth and even one answered value. Correlations and decile ranking remain
unavailable without adequate variation/sample size; R-squared and normalized
RMSE additionally need nonzero truth spread. No answers means zero coverage
and unknown error, not zero error. Historical frozen calibration values retain
their original code and are not silently recomputed by this metric correction.

A numeric row cohort may include `baseline_prediction` for every eligible entity
with `baseline_name` in its parameters. Cards expose the baseline's full-cohort
errors and errors on exactly the strategy's answered rows. MAE/MSE skill is
`1 - model_error / baseline_error` on those matched rows; a perfect zero-error
baseline or no answers leaves skill unavailable. Missing/infinite baselines,
infinite predictions and contradictory numeric abstention flags are refused.
The outcome hash includes actual baseline values, and the parameter identity
includes its name/units. The caller must still pin the baseline's training-only
fit/source lineage in its artifact; row arithmetic alone cannot establish that.
