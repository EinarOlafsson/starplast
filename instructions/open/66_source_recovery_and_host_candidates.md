# 66 — Recover source files and verify host candidate upgrades

Status: OPEN — 3/4 controlled items complete, 2026-10-08. Authorized by the user's instruction to finish open
items, download missing datasets and replace sources where warranted.

| Item | Percent done | Time left | Description |
|---|---:|---|---|
| [66.01](../done/66_01_source_recovery.md) | ✅ | 0 h | Reconcile relocated source files and retrieve missing processed sources with explicit failure records |
| [66.02](../done/66_02_rbc_candidate_review.md) | ✅ | 0 h | Review the quantitative RBC candidate against installed fraction/surface evidence |
| [66.03](../done/66_03_rhoptry_journal_version.md) | ✅ | 0 h | Verify host rhoptry journal/preprint data equivalence and version provenance |
| 66.04 | 0% | 4–8 h + validation | Correct ambiguous host symbol projections while preserving original gene-level evidence |

66.01 acceptance: inspect all 162 registered sources against the explicit existing
archive; distinguish raw input from derived and installed files; retain every
missing or ambiguous location and access failure; copy/download only matched
processed sources using recorded URLs, checksums and content validation; never
save an HTML page as assay data. A paper PMID/DOI and its file association remain
separate facts. Do not redownload an existing file merely because it was moved.
All downloads/analyses are kept in executed notebooks. No restricted VEuPathDB
access is bypassed. Do not bulk-download instrument raw files for an unresolved
assay or silently overwrite an existing source.

66.02/66.03 acceptance: primary identity and deposit/sample/table review,
source-specific mapping/QC and quantitative sanity checks, recorded citation-rate,
recency and comprehensiveness factors, and an explicit retain/add/replace decision.
Any actual runtime data change additionally requires the existing non-loss,
leakage, layout, calibration and affected benchmark/tutorial checks. A missing
deposit or an incomparable measurement stays an explicit gap.

Add these four rows to the entire progress report. Completion requires evidence,
relevant checks, a commit and nightly push. This is independent of the cancelled
worker-stagger/retry request; that request remains disregarded.

## Current evidence

The executed initial census covers all 162 sources, identifies stale archive links and
locates GTEx v10 by its original acquisition receipt. The current PMC cloud recovery
retrieves 12 additional named processed inputs, for 57 bindings. Ten attempted sources
still have version/name/identity gaps. See `results/source_recovery_pmc_2026_10_07_v2/`.
The old PMC archive service was retired in August 2026; use public discovered cloud
metadata with exact source identity, not guessed version numbers or filenames.

The RBC candidate MaxQuant ZIP and two publisher supplements have verified deposit
checksums; controlled quantitative review is complete, with retain decision and quantitative-pack
admission remaining under 65.04. PubMed confirms the
rhoptry journal article is an UpdateOf the registered preprint; supplement equivalence
is now verified for all 20,010 genes and three score fields; the canonical citation is
updated, and exact runtime values/v2 provenance are retained.

## 66.04: host symbol mapping correction

The version review found 39 screened gene symbols mapping to multiple reviewed
proteins. `host.uniprot_index()` currently uses the first accession for symbols,
while its Ensembl mapping properly withholds ambiguity. This requires a separate
controlled correction across host sources. Preserve original gene-level values
and identifier/mapping provenance; never silently select one protein, duplicate
an assay value or turn withheld projections into negatives. Audit affected host
deposits, record every ambiguity and distinguish gene perturbation from protein
measurement. Check non-loss against archived gene evidence, leakage, affected
truth/scorecards and any changed inference/calibration artifacts before promotion.
Existing protein projections remain explicitly flagged in the source note meanwhile.

Additional source-review dependency **GT-HOST-TX-02**:
`results/host_expression_source_review_2026_10_08/` verifies GTEx v11 and original
Atlas table, but native source replay differs from installed fields. Preserve
PAR_Y gene identifiers instead of splitting at the first dot. The documented
0.5 CAGE filter per region reproduces all overlapping installed brain values;
mapping/coverage gaps remain. Review original gene/source/mapping lineage with
symbol ambiguity. Original gene rows and replay gaps retained; no numerical
replacement. 534 relevant checks passed.

User steering: unavailable/undisplayable source content is recorded as a gap;
continue the next independent item without waiting for it.

An additional 43 published processed supplement inputs are retrieved and checksummed
for unmatched legacy derivations; two sources remain version/media gaps. See
`results/source_review_input_recovery_2026_10_07/`. They are original review inputs,
not arbitrary renamed legacy outputs. Continue 66.01 source/transform association
and 66.04 host mapping correction while inaccessible content stays recorded.

Final 66.01 evidence: all 162 addresses reconciled, 69 processed bindings, 11 containers and 538 candidate files; hashes verified in `results/source_recovery_final_2026_10_07_v2/`. CSPA original filename recovered with exact installed-field reproduction and protocol-specific non-detection semantics. No numerical runtime promotion. Host mapping correction (66.04) remains open.
