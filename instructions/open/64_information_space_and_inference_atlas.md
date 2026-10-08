# 64 · Explore the evidence space and inspect every inference, with ground-truth scorecards

Status: OPEN — 9/40 actions complete (22.5%), 2026-10-08. All legacy outcomes reconciled through shared scorecard semantics; categorical adapters next.

This is the current execution plan under the user's clarified goals. Instructions 53, 60, 62 and 63
retain their completed work, design decisions and evidence; this plan controls their integration order.

## The two goals

1. Explore the information space of published datasets for selected organisms and their hosts.
2. Start from a gene, protein class, localization, label or gene list; inspect every applicable
   strategy's answer; reuse precomputed work; inspect agreement and disagreement; and use tested
   meta-inference with accessible scorecards at every level.

Every inference mechanism needs a clear ground-truth test measuring both accuracy and inference
capacity. Agreement counts alone cannot establish accuracy or calibrated certainty.

The initial end-to-end scope is the existing Tg/Pf spaces and validated human/mouse host gene spaces.
Other organisms under instruction 53 enter through the same admission checklist, one species and
named source/target set per follow-up. Dataset expansion is one verified deposit per bounded task.
The inventory and all views must disclose the admitted scope and the remaining evidence/truth gaps.

## Foundations already implemented

The inspected 0.54.0 runtime has 39 strategies in nine families, six scorecard task types, calibration
and metric definitions, gene/class/target records for 14 label-calling strategies, claims, maps,
the organism registry and the pack framework. Human/mouse protein references exist separately from
future host gene spaces. The completed October 7 literature audit is retained under instruction 63.

Reuse these components. Their existence does not by itself satisfy the new acceptance contracts.
Action rows began at 0% because completion against those contracts had not been verified.

## Common ground-truth and evaluation contract

Every test records organism, target, biological context, evaluation unit, source/version/hash,
eligible population, truth grade, negative/missingness semantics, exclusions, split and settings.

Direct experiment, curation, orthology transfer, prediction, a derived quantity and a synthetic
control remain distinguishable. Recovering a predictor's outputs is not independent biological
validation. Unmeasured genes/pairs are not automatically negatives. Synthetic positives and nulls
test software and method behavior; they complement real-data validation.

Truth is frozen before fitting. Fitted imputation, scaling, feature/representation selection,
setting search, calibration and meta-model training obey the declared split. Gene-group splits
protect the promised homology hold-out. Pair tests state known-node edge completion versus cold-node
generalization. Access to unlabeled test features/graphs is explicitly marked as transductive; it
cannot establish inductive performance. Dataset/context/time hold-outs answer distinct questions
and keep separate test results.

Known-gene profiles use held-out outcomes; unknown-gene deployment predictions have a separate role.
Repeated seeds/settings/folds do not create extra independent biological samples. Uncertainty
respects the relevant biological groups. Small samples, unavailable truth and unsupported questions
remain visible.

### Accuracy and inference capacity by task

| Task | Ground truth and required interpretation |
|---|---|
| Label calls | Hidden known labels; all-eligible accuracy and accuracy-among-calls separately, per-class precision/recall/F1, macro metrics, confusion, coverage and calibration |
| Numeric estimates | Hidden measured values and units; error, skill against train-only baselines, correlation and applicability; intervals add achieved coverage and width |
| Gene ranking / set retrieval | Disjoint query seeds/validation members and a fixed candidate universe; fixed-depth precision/recall, prevalence, AUPRC/lift where identifiable, negative-label assumptions and coverage |
| Relationship inference | Hidden measured/curated relationships, endpoint identity and edge/node regime; matched controls, ranking measures, source exclusions and limits on identifying negatives |
| Map / module recovery | Hidden biological labels/modules after inner selection; recovery, placed/clustered fraction, noise and fragmentation separately from geometry/stability |
| Replication | Findings selected in discovery data and tested on disjoint validation data; tested/replicated counts, effects/direction, uncertainty and null; split-half versus new-dataset replication explicit |
| Set-valued or interval predictions | Disjoint model/calibration/test populations; achieved coverage and set/interval size or singleton efficiency together |

The existing scorecard glossary remains the metric authority. A pooled card must reproduce its
row-level records for the same cohort, target, settings and split. Baselines use the same eligible
population. Thresholds are selected before final testing. Confidence is benchmark-dependent, with
class/context/missingness/understudied strata and outside-range status.

### Cards at every level

| Level | What the user can inspect |
|---|---|
| Dataset / evidence family | Source, mapping quality, context, coverage, missingness, truth grade and lineage; measured performance only where a test supports it |
| Gene / individual prediction | Actual held-out outcome for a known gene; otherwise applicable class/stratum validation, calibrated support, evidence and applicability |
| Class / localization / protein class | Class-specific precision, recall, misses/confusions, prior, coverage and sample counts |
| Label / numeric trait / biological question | Overall and per-class/context metrics, eligible universe, baseline and gaps |
| Method / technique | Directly tested capacity, measured contribution/ablation, or mechanical known-truth/null checks, linked to the pipelines and scope tested |
| Strategy | Task-specific performance per target/organism/settings, failures, coverage and evaluation rows |
| Meta-inference | Paired base-versus-meta results, calibration, dependence, disagreements, coverage and untouched outer-test records |
| Organism | Dataset/target/strategy/test coverage and unresolved gaps; unrelated tasks are not averaged into one accuracy |

The corresponding result opens its scorecard in one action. The card expands to truth source, split,
baseline, method, provenance, failure examples and row-level outcomes. A global rate is not a measured
probability that one unknown gene is correct. A known gene's single right/wrong outcome is not a
statistically reliable per-gene accuracy rate. Composite strategy accuracy cannot be attributed to
each underlying method without testing that method's role.

## Precomputation and meta-analysis contracts

Cache reusable measurements/features, graph operators, representations/maps, permitted fitted models,
held-out outputs, target-specific deployment outputs and scorecards where useful. Artifacts declare
organism/entity/query, dependencies, settings, fit population, role, code/data hashes and version.
Builders use immutable snapshots, resume validated partitions, and invalidate declared descendants.
New user lists, labels/settings and unsupported queries reuse bases plus explicitly tracked on-demand
work; results from a different query cannot be substituted.

First show descriptive agreements/conflicts with measured dependence and evidence lineage. Then
admit a learned combination only through outer-test comparison. Methods sharing source labels,
datasets, derived features, pretrained representations or ensemble bases do not count as independent
votes. Measure matched shared errors/residuals and coverage intersections with uncertainty.

Meta-model training uses inner out-of-fold base predictions; all tuning and calibration exclude the
outer test. Compare its final accuracy, calibration and reach with a base method selected on training,
fixed consensus and matched baseline/null. Keep weak or failed combinations visible. Reuse existing
stacking/claims code where it meets this contract.

## Small, controlled execution and persistent reporting

Each linked card has dependencies, a deliverable, scope, acceptance conditions and a completion
evidence section. Work on one coherent change at a time. Freeze the pilot before expensive runs.
An unexpected source, target, algorithm or prerequisite becomes a bounded follow-up; it cannot
silently expand an item. Additional source/adapter/organism/target partitions are separately validated.

Progress: 0% not started; 25% contract and fixtures fixed; 50% implementation working; 75% evidence
and checks running/reviewed; **✅ only after all acceptance conditions pass, evidence is recorded,
and the change is committed and pushed to nightly**. Move finished cards to instructions/done/
and update parent/index links. Software handling of a missing benchmark may be implemented while
the biological ground-truth acquisition gap remains open.

**User preference:** each time an item is completed, replace its percentage with a green ✅ and
repost the **entire table below**, retaining all unfinished rows and honest time estimates. Keep
this tracker and dated completion evidence current so later sessions follow the same practice.
These estimates are engineering hours, not deadlines. Compute, source access, downloads and hosted
checks add elapsed time; budget them after pilot measurements.

## Action tracker

| Item | Percent done | Time left | Description | Depends on |
|---|---:|---:|---|---|
| [64.01](../done/64_01_entity_query_contract.md) | ✅ | 0 h | Define the entities and questions | — |
| [64.02](../done/64_02_dataset_inventory.md) | ✅ | 0 h | Inventory the available information space | 64.01 |
| [64.03](../done/64_03_measurement_provenance.md) | ✅ | 0 h | Trace each measurement to its source | 64.02 |
| [64.04](../done/64_04_ground_truth_registry.md) | ✅ | 0 h | Register the ground-truth test cases | 64.01, 64.03 |
| [64.05](../done/64_05_heldout_split_protocols.md) | ✅ | 0 h | Freeze the hold-outs and leakage rules | 64.04 |
| [64.06](../done/64_06_benchmark_baselines.md) | ✅ | 0 h | Standardize baselines and negative controls | 64.04, 64.05 |
| [64.07](../done/64_07_strategy_capability_contracts.md) | ✅ | 0 h | Declare what each of the 39 strategies answers | 64.01, 64.04 |
| [64.08](../done/64_08_result_artifact_contract.md) | ✅ | 0 h | Version inference and benchmark artifacts | 64.03, 64.05, 64.07 |
| [64.09](../done/64_09_scorecard_aggregation.md) | ✅ | 0 h | Make all scorecard levels use the same records | 64.04, 64.06, 64.08 |
| [64.10](64_10_label_benchmark_records.md) | 35% | 6–12 h | Complete the categorical-label benchmark records | 64.05, 64.06, 64.07, 64.08, 64.09 |
| [64.11](64_11_value_benchmark_records.md) | 35% | 6–10 h | Keep held-out numeric predictions | 64.05, 64.06, 64.07, 64.08, 64.09 |
| [64.12](64_12_gene_ranking_benchmarks.md) | 0% | 6–10 h | Test gene rankings and set retrieval | 64.05, 64.06, 64.07, 64.08, 64.09 |
| [64.13](64_13_relationship_benchmarks.md) | 0% | 6–12 h | Test relationship and link inference | 64.05, 64.06, 64.07, 64.08, 64.09 |
| [64.14](64_14_cluster_benchmarks.md) | 0% | 6–10 h | Test map and module recovery | 64.05, 64.06, 64.07, 64.08, 64.09 |
| [64.15](64_15_replication_benchmarks.md) | 0% | 4–8 h | Test split-half and conjunction findings | 64.05, 64.06, 64.07, 64.08, 64.09 |
| [64.16](64_16_confidence_applicability.md) | 0% | 6–12 h | Calibrate confidence and measure where it applies | 64.10, 64.11, 64.12, 64.13, 64.14, 64.15 |
| [64.17](64_17_shared_scorecard_component.md) | 0% | 4–6 h | Expose one reusable scorecard component | 64.09, 64.16 |
| [64.18](64_18_dataset_space_browser.md) | 0% | 4–8 h | Browse the selected organism's datasets | 64.02, 64.03, 64.17 |
| [64.19](64_19_gene_evidence_entry.md) | 0% | 4–8 h | Make gene lookup an evidence entry point | 64.01, 64.03, 64.17, 64.18 |
| [64.20](64_20_label_class_entry.md) | 0% | 4–8 h | Make labels and protein classes entry points | 64.01, 64.02, 64.17, 64.19 |
| [64.21](64_21_human_gene_space.md) | 30% | 6–12 h | Build the human gene information space | 64.01, 64.02, 64.03, 64.04, 64.05, 64.08 |
| [64.22](64_22_mouse_gene_space.md) | 0% | 8–16 h | Build the mouse gene information space | 64.01, 64.02, 64.03, 64.04, 64.05, 64.08 |
| [64.23](64_23_space_selection_packs.md) | 0% | 4–8 h | Select and install available organism spaces | 64.18, 64.21, 64.22 |
| [64.24](64_24_host_parasite_evidence.md) | 0% | 4–8 h | Explore typed organism-to-host connections | 64.03, 64.21, 64.22, 64.23 |
| [64.25](64_25_precompute_job_pipeline.md) | 0% | 6–10 h | Build a resumable precomputation pipeline | 64.08, 64.10, 64.11, 64.12, 64.13, 64.14, 64.15, 64.16 |
| [64.26](64_26_toxoplasma_inference_pack.md) | 0% | 4–8 h + compute | Precompute the Toxoplasma inference atlas | 64.07, 64.17, 64.25 |
| [64.27](64_27_plasmodium_inference_pack.md) | 0% | 4–8 h + compute | Precompute the Plasmodium inference atlas | 64.07, 64.17, 64.25, 64.26 |
| [64.28](64_28_human_inference_pack.md) | 0% | 6–12 h + compute | Precompute and validate human inferences | 64.10, 64.11, 64.12, 64.13, 64.14, 64.15, 64.16, 64.21, 64.23, 64.25 |
| [64.29](64_29_mouse_inference_pack.md) | 0% | 6–12 h + compute | Precompute and validate mouse inferences | 64.10, 64.11, 64.12, 64.13, 64.14, 64.15, 64.16, 64.22, 64.23, 64.25 |
| [64.30](64_30_gene_inference_profile.md) | 0% | 6–10 h | Show every strategy's answer for a gene | 64.17, 64.19, 64.26, 64.27 |
| [64.31](64_31_class_label_inference_profile.md) | 0% | 4–8 h | Compare inferences for a label or protein class | 64.17, 64.20, 64.26, 64.27, 64.30 |
| [64.32](64_32_method_dependence_audit.md) | 0% | 6–10 h | Measure dependence between inference mechanisms | 64.03, 64.08, 64.10, 64.11, 64.12, 64.13, 64.14, 64.15, 64.16 |
| [64.33](64_33_agreement_disagreement_view.md) | 0% | 4–8 h | Show calibrated agreement and conflicts | 64.30, 64.31, 64.32 |
| [64.34](64_34_validated_meta_inference.md) | 0% | 8–16 h + compute | Train and test the meta-inference | 64.05, 64.06, 64.08, 64.16, 64.32, 64.33 |
| [64.35](64_35_meta_scorecards.md) | 0% | 4–6 h | Expose meta-inference scorecards at every level | 64.17, 64.30, 64.31, 64.34 |
| [64.36](64_36_navigation_acceptance.md) | 0% | 4–6 h | Connect every evidence and inference level | 64.18, 64.19, 64.20, 64.23, 64.24, 64.30, 64.31, 64.33, 64.35 |
| [64.37](64_37_imported_screen_context.md) | 0% | 4–8 h | Make user data reach the same inference flow | 64.03, 64.04, 64.05, 64.08, 64.23, 64.25, 64.30 |
| [64.38](64_38_prospective_frozen_tests.md) | 0% | 6–10 h | Freeze predictions and score new independent truth | 64.03, 64.04, 64.05, 64.08, 64.16, 64.32, 64.34 |
| [64.39](64_39_reproducible_exports.md) | 0% | 4–8 h | Export the evidence and inference trail | 64.03, 64.08, 64.17, 64.30, 64.31, 64.35 |
| [64.40](64_40_release_goal_acceptance.md) | 0% | 6–10 h + checks | Verify both goals end to end and release | 64.18, 64.19, 64.20, 64.21, 64.22, 64.23, 64.24, 64.26, 64.27, 64.28, 64.29, 64.30, 64.31, 64.33, 64.35, 64.36, 64.37, 64.38, 64.39 |

## Every current strategy's ground-truth contract

Additional authorized actions remain part of **every full completion report**:
**18/51** complete overall after the 68.01 Discoveries browser, with 9/40 of this
instruction's original actions complete. The newest priority is **68.02 function**;
reuse the source/split/record/scorecard contracts above. Do not treat broad annotation
browsing as completion of functional inference or independent biological validation.

| Item | Percent done | Time left | Description | Depends on |
|---|---:|---|---|---|
| [65.01](../done/65_01_dataset_selection_policy.md) | ✅ | 0 h | Implement the dataset-selection rule | — |
| [65.02](../done/65_02_publication_identity_audit.md) | ✅ | 0 h | Audit publication identities and citation rates | 65.01 |
| [65.03](../done/65_03_slot_literature_discovery.md) | ✅ | 0 h | Search literature alternatives for every slot | 65.01, 65.02 |
| [65.04](65_dataset_selection_audit.md) | 4% | 16–30 h + compute | Validate eligibility and justified replacements | 65.01, 65.02, 65.03 |
| [66.01](../done/66_01_source_recovery.md) | ✅ | 0 h | Recover and reconcile missing processed sources | — |
| [66.02](../done/66_02_rbc_candidate_review.md) | ✅ | 0 h | Review the quantitative RBC candidate | 66.01 |
| [66.03](../done/66_03_rhoptry_journal_version.md) | ✅ | 0 h | Verify the host rhoptry journal version | 66.01 |
| [66.04](../done/66_04_host_symbol_projections.md) | ✅ | 0 h | Correct ambiguous host projections | 66.03 |
| [67.01](../done/67_01_work_continuation_watchdog.md) | ✅ | 0 h | Keep continuation active and log actual stops | — |
| [68.01](../done/68_01_discoveries_label_browser.md) | ✅ | 0 h | Browse all available labels and functions in Discoveries | — |
| [68.02](68_discoveries_function_coverage.md) | 25% | 8–14 h + compute | Benchmark and precompute functional inferences | 68.01, 64.04–64.09 |

### Strategy ground-truth matrix

This matrix specifies what must be retained or established. It does not assert that every test
already exists or passes. A strategy can have multiple task contracts for different targets/queries.
64.07 also registers all underlying techniques and their appropriate validation roles.

| # | Strategy | Action | Ground-truth test |
|---|---|---|---|
| 1 | `holdout_search` | [64.14](64_14_cluster_benchmarks.md) | Hidden experimental labels/modules; select map/settings on inner data and score outer recovery. |
| 2 | `geneset_hunt` | [64.12](64_12_gene_ranking_benchmarks.md) | Known set members hidden from the query seeds and map selection; test retrieval on that set's declared population. |
| 3 | `recoverability_atlas` | [64.12](64_12_gene_ranking_benchmarks.md) | Held-out known-label recovery and predictability ranking, retaining target/class and perturbation scope. |
| 4 | `consensus_modules` | [64.14](64_14_cluster_benchmarks.md) | Known-label/module recovery on untouched genes, plus stability and coverage across map settings. |
| 5 | `blind_battery` | [64.15](64_15_replication_benchmarks.md) | Measured findings selected on discovery data and checked on disjoint verification data. |
| 6 | `block_ablation` | [64.10](64_10_label_benchmark_records.md) | Same hidden labels with/without a feature block; paired recovery delta and pipeline contribution. |
| 7 | `feature_knn` | [64.10](64_10_label_benchmark_records.md) | Known labels hidden by gene/homology group; full class calls/scores and abstentions. |
| 8 | `map_neighbours` | [64.10](64_10_label_benchmark_records.md) | Known labels hidden from map construction and neighbor voting; held-out recovery. |
| 9 | `cluster_guilt` | [64.10](64_10_label_benchmark_records.md) | Known labels hidden from cluster enrichment/assignment; per-class calls and coverage. |
| 10 | `label_outliers` | [64.12](64_12_gene_ranking_benchmarks.md) | Synthetic corrupted labels with known originals; separately registered independent adjudication for real label-error claims. |
| 11 | `layer_propagation` | [64.10](64_10_label_benchmark_records.md) | Known labels hidden from graph seeds with target-derived edges excluded. |
| 12 | `layer_vote` | [64.10](64_10_label_benchmark_records.md) | Known labels hidden from every contributing voter and its weight fitting. |
| 13 | `physical_partners` | [64.10](64_10_label_benchmark_records.md) | Known labels hidden from physical-network seeds; explicit network coverage and target exclusions. |
| 14 | `structural_homology` | [64.10](64_10_label_benchmark_records.md) | Measured/curated function or label hidden from fold-neighbor voting; declared truth grade and mapping. |
| 15 | `multiplex_modules` | [64.14](64_14_cluster_benchmarks.md) | Known labels or reviewed modules hidden from module selection and evaluated on untouched groups. |
| 16 | `link_prediction` | [64.13](64_13_relationship_benchmarks.md) | Hidden experimental edges in an explicit known-node/cold-node regime; matched controls and unlabeled-pair limits. |
| 17 | `attention_correction` | [64.13](64_13_relationship_benchmarks.md) | Recovery of separately measured/curated biological relations; source-paper reuse declared and tested. |
| 18 | `unwritten_links` | [64.13](64_13_relationship_benchmarks.md) | Literature recovery is its own test; claims of physical interaction require separately measured physical truth. |
| 19 | `supervised_classifier` | [64.10](64_10_label_benchmark_records.md) | Known labels hidden from all fitted preprocessing, classifier tuning and calibration. |
| 20 | `positive_unlabeled` | [64.12](64_12_gene_ranking_benchmarks.md) | Held-out verified positives; precision/negative-class claims require an independently labelled validation sample. |
| 21 | `trait_regression` | [64.11](64_11_value_benchmark_records.md) | Known measured traits hidden by group/context; train-only baselines and unit-aware residuals. |
| 22 | `masked_imputation` | [64.11](64_11_value_benchmark_records.md) | Measured entries masked before fitting, with entry/group/context mask regimes stated. |
| 23 | `condition_shift` | [64.11](64_11_value_benchmark_records.md) | Hidden condition comparisons derived from measured inputs; both source conditions and training exclusions recorded. |
| 24 | `set_enrichment` | [64.12](64_12_gene_ranking_benchmarks.md) | Disjoint seed/validation set members and known annotation associations; enrichment significance and recovery separated. |
| 25 | `seed_expansion` | [64.12](64_12_gene_ranking_benchmarks.md) | Known members withheld from seeds by group; retrieval on a fixed eligible universe. |
| 26 | `split_clusters` | [64.15](64_15_replication_benchmarks.md) | A discovered biological split repeats with its direction/effect on disjoint validation groups. |
| 27 | `conjunctions` | [64.15](64_15_replication_benchmarks.md) | A discovered joint label/context finding repeats on disjoint validation groups. |
| 28 | `paralog_divergence` | [64.12](64_12_gene_ranking_benchmarks.md) | Ranked paralog differences checked against withheld measured/curated functional differences; both endpoints and proxy limits recorded. |
| 29 | `ortholog_transfer` | [64.10](64_10_label_benchmark_records.md) / [64.11](64_11_value_benchmark_records.md) | Receiver-species measured labels/values hidden; donor-target semantic compatibility and transfer/mapping coverage explicit. |
| 30 | `stratum_focus` | [64.10](64_10_label_benchmark_records.md) | Known labels hidden within declared understudied/lineage/missingness strata, with stratum-specific results. |
| 31 | `triangulation` | [64.10](64_10_label_benchmark_records.md) | Known labels hidden from each base and the combination; shared evidence and abstention visible. |
| 32 | `understudied_first` | [64.10](64_10_label_benchmark_records.md) | Same held-out combination test on the stated understudied population. |
| 33 | `neighbour_space` | [64.13](64_13_relationship_benchmarks.md) | Held-out measured relationship layer/group pairs, excluding its source from fitted features. |
| 34 | `network_training` | [64.13](64_13_relationship_benchmarks.md) | Held-out measured relationship/group pairs, excluding the trained target layer and its sources. |
| 35 | `conformal_calls` | [64.10](64_10_label_benchmark_records.md) | Disjoint model/calibration/test labels; achieved set coverage and set/singleton efficiency. |
| 36 | `graph_convolution` | [64.10](64_10_label_benchmark_records.md) | Known labels hidden from fitted graph model and all target-related feature/edge sources. |
| 37 | `random_forest` | [64.10](64_10_label_benchmark_records.md) | Known labels hidden from train-fitted features, forest selection and calibration. |
| 38 | `stacking` | [64.10](64_10_label_benchmark_records.md) | Outer held-out truth; training stacker sees only inner out-of-fold base predictions. |
| 39 | `conformal_values` | [64.11](64_11_value_benchmark_records.md) | Disjoint model/calibration/test measured values; achieved interval coverage, width and point-error skill. |

## Execution order and gates

| Milestone | Actions | Gate |
|---|---|---|
| Contracts and trustworthy records | 64.01–64.09 | Source, truth, split, capabilities and aggregation contracts agree on small fixtures |
| Ground-truth task coverage | 64.10–64.17 | Mechanisms have appropriate tests; pilot failures/gaps and calibrated scope are accessible |
| Evidence-space exploration | 64.18–64.24 | Parasite and admitted host sources/entities/context are browsable without mixed identifiers |
| Precomputed parasite pilot and profiles | 64.25–64.27, 64.30–64.31 | Gene/class/label routes reach cached outputs and their cards; remaining on-demand work explicit |
| Validated host inference | 64.28–64.29 | Each host uses its own targets, truth, calibration and applicability |
| Agreement and tested meta-inference | 64.32–64.36 | Outer-test evidence shows what the combination adds and where it fails |
| Import, prospective tests and export | 64.37–64.39 | Imported inputs reach the same contexts; frozen predictions are testable; exports reproducible |
| Installed release acceptance | 64.40 | Both product flows pass within their declared supported scope, with required CI/docs/tutorials |

Dependencies control execution; ranges are a browsing aid. The first implementation action is
64.01, followed by its prerequisites for the benchmark and source contracts. Expensive sweeps
require measured pilot parity, runtime/RAM/storage, an immutable snapshot and the shared RAM lease.
Existing environment, GPU-sharing and release instructions remain in force.

## Earlier work incorporated

| Earlier instruction | Where the work is reused |
|---|---|
| 53 — organism spaces | 64.21–64.29: identity, host builders, packs, leakage and bridges; other species enter one at a time |
| 60 — usefulness | 64.03, 64.08, 64.16, 64.19, 64.30, 64.37, 64.39: provenance, digest, target calibration, imports, methods and export |
| 62 — track record | 64.09–64.10, 64.17, 64.19–64.20, 64.30–64.31: existing categorical records, aggregation and drill-down |
| 63 — claims | 64.16, 64.32–64.38: confidence, applicability, dependence, verification and prospective tests |
| 38 / 41 — datasets | 64.02–64.04 define acquisition gaps; each new admission is one verified source with its resolving queries |
| 59 — biological questions | Reviewed question results become entry-point examples and acceptance queries, with matching benchmark evidence |

Prior completed results are preserved. Independent verification, learned combination and prospective
testing have separate meanings and acceptance evidence.

## Completion of the two goals

Goal 1 is met for the admitted release scope when every supported selected organism/host opens its
own information space; the user can browse all admitted datasets/questions, inspect coverage and
context, open measured entities, and trace their sources, including missing/unavailable evidence.

Goal 2 is met for the admitted release scope when a gene, class or label opens a complete account
of applicable strategy outputs/statuses, with reusable precomputation, interpretable agreements and
conflicts, tested meta-inference, and scorecards linked to applicable ground truth at every required
level. A visible ground-truth gap is honest tracking; the missing biological validation stays open.

64.40 records the installed-build evidence matrix and release outcome. Until its acceptance is
verified, this instruction and its implementation actions remain OPEN.

## Planning validation — 2026-10-07

The repository contains 40 individual action cards and 40 matching tracker rows. Their source paths
and new documentation links resolve; dependencies are acyclic. The strategy contract matrix was
checked against the live catalogue and covers all 39 keys. This verifies the action plan's structure,
not implementation or biological performance. All implementation acceptance conditions remain open.
