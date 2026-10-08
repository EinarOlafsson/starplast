# 68 · Broaden Discoveries beyond localization

User priority, 2026-10-08: Discoveries should span all available labels, with
particular emphasis on function. Localization was the first developed inference
path; the existing menu uses only nonempty generated claims and default filters
leave 43 Toxoplasma localization claims visible. Existing cell-cycle claims are
hidden; phenotype recipes are uncalibrated, and domain/EC annotations have not
entered the view. Broaden the actual available information and expose missing
inference/benchmark coverage clearly.

| Item | Percent done | Time left | Description |
|---|---:|---|---|
| [68.01](../done/68_01_discoveries_label_browser.md) | ✅ | 0 h | Browse all available labels and functional classes, memberships and inference/test coverage in Discoveries |
| 68.02 | 35% | 6–12 h + compute | Add functional ground-truth targets, strategy benchmarks and precomputed functional claims |

68.01 is bounded to current categorical annotations/flags plus multi-valued InterPro,
Pfam and EC fields, both existing organisms, and existing claims/recipes. Search labels
and individual terms, retain unclaimed/uncalibrated variables, select class members
and open genes. Known annotations, inferred claims and missing tests must remain
distinct. Existing claim filters/verdicts/calibration are preserved. Source/species
lineage stays explicit; absence of a domain/EC annotation is not adjudicated
biological absence. No new ontology ingestion or algorithm fitting in this item.
Acceptance: exact source membership counts, multi-valued parsing, all existing
recipe/claim targets and annotations visible, functional search, class/gene
navigation, available versus missing tests, offscreen UI/compatibility checks,
executed before/after coverage audit, commit/push and full progress table.

68.02 uses the existing 64.04–64.16 contracts and 64.25 precomputation workflow.
First freeze one named functional target/source/cohort per adapter partition.
Domain/EC source grade, curation/prediction/transfer lineage, source-family closure,
multi-label versus verified negative semantics and homology/source independence
must be reviewed before fitting or declaring biological accuracy. Expand strategy
outputs with accuracy/capacity scorecards; no confidence fabricated from matching
annotations or agreement. Deliver precomputed functional results and explicit
failures/gaps accessible from this view. Whole completion requires actual verified
functional coverage across the declared source/target/mechanism scope.

Next bounded partition: source-grade and ontology/version audit of installed EC,
InterPro and Pfam, then freeze one named functional source/cohort before a native
adapter pilot. Missing functional annotation is unknown; do not train absence as
a verified negative. Retain overlapping EC/domain memberships rather than selecting
the first class. Review obsolete/transferred EC identifiers before merging ontology
classes. Plasmodium domain descriptions are absent in the installed source; any
new descriptions require a pinned authoritative ontology, not guessed transfers.
Link resulting held-out class/label/strategy cards and explicit failures through
Discoveries before counting precomputed results as delivered.

Frozen first partition **FN-EC-01**, 2026-10-08: acquire/pin the current ENZYME
nomenclature metadata and review all installed EC codes (both organisms and
separate orthology-derived fields). Resolve unique transfer chains only; deleted,
unknown, split/ambiguous and preliminary entries keep explicit statuses. Original
gene annotations remain untouched. Build the full set of major enzyme classes for
each completely resolvable annotated gene, including multiple major classes;
freeze this exact profile as a categorical source-recovery target for native kNN.
Unknown/unresolvable gene profiles do not become negatives. Train-only fitting,
source-family closure, protected homology groups, all test rows/native class
scores, profile and member-class cards, matched baselines and exact replay required.
This first recovery test is not independent biological function accuracy, and its
vote shares do not become calibrated probabilities or new verified claims.
Subsequent partitions add independent activity truth, additional functional
targets/adapters and calibration/deployment. First functional result navigation
is now delivered below; full target/mechanism coverage remains open.

These items extend instructions 64/65 and do not replace or shrink them. Repost
the full tracker including both new rows whenever an entire action completes.

68.01 acceptance verified: 41/30 labels, exact replay of every available variable,
21,787/18,879 functional memberships, real class/gene/legacy-scorecard navigation,
413 final checks passed (one optional pdoc skip). Canonical executed evidence:
`results/discoveries_label_coverage_2026_10_08_v4/`. Original claims/recipes/records
unchanged. Functional claims and held-out rows remain zero; **68.02 is next**.
The three earlier audit diagnostics are retained. Current EC ontology versions,
curation/domain prediction/transfer source grades and independent negative semantics
are explicit prerequisites for functional inference, not hidden by this completion.

## 68.02 first functional partition — verified, 25%, 2026-10-08

Current nomenclature is acquired and pinned: ENZYME 02-Sep-2026, InterPro 110.0,
Pfam 38.2. All installed EC fields have explicit unique-transfer/ambiguous/deleted/
missing/preliminary/malformed status; original annotations remain unchanged.
Complete profiles: 1,226 Toxoplasma, 1,050 Plasmodium curated-source, 1,055 separate
orthology-derived genes. Unknown or unresolved annotation never becomes a verified
negative. Current domain names cover 7,926 of 8,043 installed unique identifiers;
112 missing-current InterPro IDs and five officially withdrawn Pfam IDs keep their
original memberships and unresolved/withdrawn statuses. All 34,983 domain memberships
and source descriptions replay exactly. Source assignment releases and gene-level
experimental/curation/prediction grades remain unresolved.

Frozen **FN-EC-01** native numeric kNN: 679 train / 175 tune / 190 calibration /
182 test genes, seed 23, k=15, minimum vote share 0.3. All complete multi-valued
major-class profiles preserved. Training-only ranks/feature selection/imputation;
349 inputs explicitly exclude source/domain/homology and attention counts; no graph.
71/182 exact profiles recovered (39.0%); 171/182 answered (94.0%); 71/171 correct
among answers (41.5%); matched majority baseline 66/182 (36.3%). One unseen profile
and all 11 abstentions retained. Native predictions/scores/model state and all cards
replay exactly. Typed identity:
`349317afd8f24b860d94ae1e11280a6130bf451700b02b9d2145fab85fc99162`.
**Annotation recovery only; zero independent activity benchmarks or new calibrated
unknown-gene deployment claims.** This weak result cannot establish biological accuracy.

Discoveries exposes label, strategy, complete-profile, seven overlapping major-class
and matched-control scorecards, exact gene outcomes and existing gene navigation.
Both organisms support domain-name search; selected classes retain original descriptions,
current-name source/status/release/version caveats and metadata/license URLs. Missing
or corrupt names/results leave browsing usable with explicit unavailable status.
Plasmodium functional benchmarks remain unavailable, rather than inheriting Toxoplasma
performance. Original nodes/claims/recipes/legacy record hashes unchanged.

Canonical executed evidence:

- `results/functional_source_review_2026_10_08_v2/`
- `results/functional_domain_metadata_2026_10_08_v4/`
- `results/functional_ec_knn_pilot_2026_10_08_v4/`
- `results/functional_ec_ui_bundle_2026_10_08_v4/`
- `results/functional_discoveries_audit_2026_10_08/`

Validation: 652 broad checks passed; the remaining new tooltip fixture needed HTML
text normalization, after which all 19 functional reader/UI checks passed. All 653
distinct checks are verified across those runs. Earlier two failures (Qt tuple
lookup and a public-method docstring) are corrected. Prior acquisition/grammar/
context/packaging diagnostics and prototypes remain immutable.

Remaining bounded partitions: independent activity truth and source-release admission;
Plasmodium EC and full InterPro/Pfam target profiles; additional strategy adapters;
source-closure sentinel tests, calibration, applicability and unknown-gene deployment.
The domain-content family and domain-edge family differ, and `search.LAYER_SOURCES`
omits domain: graph adapters must explicitly withhold domain edges and derived
representations. The first pilot has `graph={}`. Attention-count aliases need a
general feature-exclusion follow-up (both counts are already withheld in this pilot).
No full 68.02 completion or new progress-table tick is claimed.

## Source preparation and explicit coverage — 30%, 2026-10-08

The recorded-domain census preserves every original source cell and all
34,983 memberships, including 175 missing-current InterPro and 13 withdrawn Pfam
memberships. The executed domain census exposes complete-profile counts and
singleton profiles without fitting or assigning new functions. Unknown and malformed
profiles cannot become verified negatives; singleton/unsupported test profiles must
remain in future capacity denominators.

The separate immutable prepared-target contract retains normalized complete
profiles, unknown states, universe/order and protected split identities. Its
semantic target identity does not substitute for the original source/census hash.
Capacity reports are reporting-only; test support cannot select settings or features.

The executed coverage census has 880 organism/source/target/strategy/task addresses,
11 installed functional source labels and one verified reference-recovery artifact.
Human/mouse inference adapters remain unsupported/uninventoried; no accuracy transfers
between species or labels. Discoveries now exposes the current organism's 240/320
addresses and links the exact available EC result. Independent biology, calibration
and deployment remain unavailable. Exact installed-node hash and table-value checks
prevent imported or altered contexts from inheriting archived performance.

An opt-in functional exclusion helper closes domain/homology/literature/attention
sources, registered derived columns and explicitly declared representations. Sentinel
checks reject hidden-role selection, cyclic/unresolved lineage and empty encoded
targets. This strengthens future benchmark preparation; it does not retroactively
change old artifacts or prove independence of undeclared inputs.

Shared scorecard, class and host evidence presentation advances 64.17. Source/recovery
snapshots and diagnostics stay immutable. FN-EC-02 stopped with an observed runtime
content-access restriction; it remains unfinished and was not retried. Ordinary
software tasks continued. Independent activity truth and full functional mechanism
coverage remain open; total completion remains 18/51.

Prepared-target evidence is executed in
`results/functional_profile_targets_2026_10_08_v3/pilot.ipynb`. The fixed Toxoplasma
complete-Pfam target retains 8,140 genes, 4,310 eligible profiles, 3,830 unknown
states, 2,279 distinct profiles and 5,952 memberships. Split roles contain
2,381/635/659/635 genes; 333 test profiles lack training support and remain in
the cohort. These are source/split capacity counts, not inference accuracy or
independent biological truth. Exact raw-cell, normalized-target and protected-group
replay is retained with separate source and semantic identities. Earlier failed
serialization attempts remain immutable diagnostics. No model was selected or fit;
68.02 remains 30% and independent truth/source admission remain open.

## Frozen complete-Pfam baseline controls, 2026-10-08

`results/functional_profile_controls_2026_10_08/` consumes the pinned verified
prepared-target V3 with unchanged source, ordered protected-group split and
control seed 20261008. Shared training-majority and seeded training-prevalence
controls are frozen before reading test truth. All 635 test genes remain,
including 333 unsupported profiles; 3,830 unknown source genes stay outside the
eligible cohort. Complete original profile classes are retained, with zero
training prior for unseen profiles. These controls use no features or classifier.

Exact prediction/score/parquet/JSON replay, independently counted class metrics,
full confusion records and shared scorecard snapshots are saved. Majority calls
recover 19/635 profiles (2.9921%); seeded prevalence calls recover 2/635 (0.3150%).
These are reference-annotation control measurements, not biological function
accuracy or calibrated gene confidence. Classes without test truth have explicit
unavailable recall; unsupported test truth still contributes to macro metrics.

16 focused control, split/unknown, test-truth-invariance and input corruption
checks pass. The executed notebook passes under a 400 MiB cap in 6.676 s;
35 input/source/code hashes and 22 output receipts verify. All numeric checks are
exact; no tolerance relaxation. The source release, source/homology independence,
biological admission, actual source-excluded classifier and deployment remain
open. 68.02 stays 30%; whole-action completion remains 22/51 after 64.20.

## Functional input closure prerequisite, 2026-10-08

The opt-in functional-source policy is now v2: publication/paper/citation/fulltext
aliases (including counts, rates and case variants) and declared downstream
representations are withheld. This closes the difference between the numeric
precompute guard and direct categorical adapters before the next Pfam pilot.
64 exclusion/operator checks pass, including 18 alias/derived-lineage sentinels;
ordinary measurement names remain usable. Installed runtime strategies and all
historical source/result packets remain unchanged. New benchmarks pin this policy
code; undeclared biological/source lineage remains unresolved, not certified.

## Frozen complete-Pfam kNN pilot, 2026-10-08

Before the actual fit: reuse only the pinned prepared-target V3 and unchanged
protected train/tune/calibration/test roles (2381/635/659/635). Fixed native
feature-kNN k=15/min_share=0.3, no tune/test selection. Select registered numeric
inputs on training coverage/variation only; close functional/domain, homology,
literature/attention and declared/registered descendants, plus unregistered or
unresolved derivations. Fit native training ranks and frozen query ECDF; use no
graphs/maps. Retain all 635 test genes and 333 unsupported complete profiles.

The frozen pilot records native
raw neighbor votes, method support, class scores/model/exclusions, exact replay,
matched train-only controls, all-eligible/among-call recovery and class/macro
confusions. Unsupported inputs retain explicit abstention/null and unavailable
artifact status. Ground truth remains original source-annotation recovery;
source/homology independence, biological accuracy, calibration and deployment
stay unknown. 15 synthetic source/hidden-role/native/unavailable checks pass.
Do not promote biological claims or complete 68.02 from this one pilot.

Canonical acceptance: `results/functional_profile_knn_2026_10_08_v4/` retains
348 training-selected numeric inputs, 1,413 native training classes and 2,279
full-source reporting classes. Exact native ranks, calls, neighbor weights,
support, scores, refit, parquet/JSON and typed-artifact replay pass. Of 635 test
genes, 20 receive calls: 12 correct, eight wrong, 615 abstentions. All-eligible
recovery is 12/635 (1.89%), coverage 20/635 (3.15%), and correctness among calls
12/20 (60%). Majority/prevalence controls recover 19/2 profiles; the candidate
does not beat majority on all-eligible recovery. These are annotation-profile
metrics, not biological accuracy or confidence in unknown genes.

The accepted executed run finishes in 46.61 seconds with 975.92 MiB process
peak, serial under the existing 1,900 MiB root allowance; all workers terminal.
Independent notebook `results/functional_profile_knn_acceptance_2026_10_08/`
verifies 47 input hashes, 47 output receipts, 12 artifact receipts and exact
denominators. Artifact identity:
`38cfb6bf110cadb559eac41c798298d7d95882efe2524a3591ccf7ad63d0385f`.

Preserved earlier attempts: pandas3 NaN/None categorical representation failure;
two worker signal9 diagnostic stops with unproven cause; root full-matrix
diagnostic proving exact numerical/native results; typed unit/population metadata
refusal; authoritative systemd 800 MiB OOM while saving. Only missing-value
representation, typed scope and retained-memory lifetime changed; settings and
all cohorts remain fixed. Browser packaging, independent truth/source admission,
additional strategies and calibrated deployment remain open. 68.02 stays30%;
whole-action completion remains22/51, with no new tick.

## Complete-Pfam Discoveries integration, 2026-10-08 — partial35%

The offline reader and packager support explicit EC/Pfam namespaces and preserve
all native payload hashes, control names/scopes and protected split roles. Domain
membership cards are independently derived recorded-presence/complement tests,
with their recipe/row identity retained; they do not adjudicate biological absence.
The accepted candidate packet `results/functional_profile_ui_bundle_2026_10_08/`
has 479 native profile cards and 583 domain-membership cards, all635 test genes
and333 unsupported profiles. Packaging0.99s/266.5MiB under400MiB, no fitting.

The executed merge preserves original EC and Pfam entries exactly. Shipped bundle
SHA256 `854d014260330490e2b0f7189d8e6744448ab5d39f7686d75c4aba3cd763a84b`.
Actual desktop audit `results/functional_profile_ui_2026_10_08_v3/` verifies
source→result, strategy, controls, membership/profile exports, abstentions,
gene and coverage navigation, changed-context refusal and no Plasmodium borrowing.
Runtime10.17s, process peak1,067,978,752bytes under serial1900MiB. Screenshot
inspected;118 relevant reader/class/coverage/domain/Discoveries regressions pass.
Independent acceptance verifies30 packet outputs,26 inputs and12 original artifact
receipts, including historical code/bundle snapshots after intentional pin changes.
Preserve first800MiB OOM and second audit's wrong export-envelope key diagnostic.
Independent biological truth, other functional mechanisms/organisms, calibration
and deployment remain open;22/51 whole actions, no completion tick.

## Predeclared Plasmodium EC controls, 2026-10-08

Scope FN-EC-PF-CONTROL-01: direct `ec_number` complete major-class profiles;
verified ENZYME/source review V2 manifest
`5fb85c26e4939c4737a6969c474cd63395ea2faaa1fb111a71ebd241e0aaa81c`.
Keep the whole original Pf gene order, native protected groups and explicit
missing-group fallback/collision checks. Seed20261008 and fractions
(0.55,0.15,0.15,0.15) fixed before execution. Source exclusions and
majority/prevalence controls use training only; unknown and unsupported profiles
remain visible. `ec_number_orthology` is a separate census, never merged into
the target or used as a predictor. No classifier, acquisition, biological
admission, calibration or deployment; direct annotation evidence stays unresolved.

Accepted preparation/control packet:
`results/pf_functional_ec_controls_2026_10_08/`. All5,720 genes retained;
direct EC has1,050 complete profiles/19 classes,170 unresolved and4,500
unannotated. Orthology-derived EC remains a separate census (1,055 complete,
451 unresolved,4,214 unannotated). Roles600/151/147/152 use recorded native
groups with zero missing-group fallbacks or collisions. All four unsupported
test genes across three profiles remain in the denominator.

Majority/prevalence controls recover43/152 (28.29%) and35/152 (23.03%).
Actual control-name scopes, full confusion, profile/seven-major-class metrics,
raw source/normalization/group/split and exact serialization replay retained.
No inference-strategy fit or independent biological accuracy. Runtime4.13s,
261,208KiB process peak under400MiB;128 relevant regressions pass. Independent
acceptance verifies40 output checksums,38 inputs and exact population arithmetic.
The source catalogue now describes direct annotations with unresolved individual
curation/prediction provenance, rather than asserting manual curation. Generated
catalogue refreshed; source/node bytes unchanged. No dataset promoted or replaced.
68.02 stays35%,22/51 whole actions; actual Pf inference and UI result remain next.
