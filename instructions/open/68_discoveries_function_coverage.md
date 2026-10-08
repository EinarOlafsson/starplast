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
| 68.02 | 0% | 8–16 h + compute | Add functional ground-truth targets, strategy benchmarks and precomputed functional claims |

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
