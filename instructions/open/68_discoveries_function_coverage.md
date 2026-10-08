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
| 68.02 | 30% | 6–12 h + compute | Add functional ground-truth targets, strategy benchmarks and precomputed functional claims |

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
