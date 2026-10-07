# 66 — Recover source files and verify host candidate upgrades

Status: OPEN, 2026-10-07. Authorized by the user's instruction to finish open
items, download missing datasets and replace sources where warranted.

| Item | Percent done | Time left | Description |
|---|---:|---|---|
| 66.01 | 60% | 2–4 h + downloads | Reconcile relocated source files and retrieve missing processed sources with explicit failure records |
| 66.02 | 35% | 3–6 h + validation | Review the quantitative RBC candidate against installed fraction/surface evidence |
| 66.03 | 20% | 2–3 h + lookup | Verify host rhoptry journal/preprint data equivalence and version provenance |

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

Add these three rows to the entire progress report. Completion requires evidence,
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
checksums; candidate analysis and admission remain pending. PubMed confirms the
rhoptry journal article is an UpdateOf the registered preprint; supplement equivalence
remains unverified. No installed source was changed.
