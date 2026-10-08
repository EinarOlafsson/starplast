# 66.03 · Verify host rhoptry journal/preprint equivalence

Status: DONE — ✅, 2026-10-07 (local date). Version comparison and citation upgrade
complete; the separate symbol-mapping correction is 66.04.

PubMed PMID 42791346 verifies the journal DOI and UpdateOf PMID 41279910.
The public publisher Dataset EV1 link supplies the verified archive URL.
The archive has one actual screen workbook plus macOS resource forks; the
initial parser failure is retained, and the corrected parser excludes those
non-workbook metadata members and reuses the original checksummed download.

All **20,010 source genes** are identical between journal and preprint. Combined
rhoptry score, beta and Wald FDR match exactly for all 20,010 rows each: no missing
cells and maximum difference zero. The original loader and `%.10g` serialization
reproduce all **18,739 installed protein rows** exactly for each of three fields.
The independent mapping audit finds 18,700 unambiguous gene symbols, **39
ambiguous symbols** and 1,271 unmapped symbols. Existing protein projections
for the 39 ambiguous symbols use the legacy first-accession rule; they are
explicitly flagged in the source note and queued for correction under 66.04.

Decision: retain the exact installed values and original v2 source-file provenance;
upgrade the canonical citation to the peer-reviewed journal, preserving preprint
lineage. The dataset catalogue and generated source script are refreshed. No
measurement, mapping or strategy output is silently replaced.

Primary bioRxiv API confirms v1 first published **2025-10-17**, v2 **2025-10-21**.
The journal publication date is 2026-09-25. Source age must not reset on journal
publication. Journal and preprint are linked versions of one experiment, not two
independent candidate datasets; potentially overlapping citations are not summed.
The recorded primary Europe PMC snapshot has zero citations for each version.
The lineage retrieval crossed UTC midnight (October 8) but remains October 7
in the user's local timezone; all retrieval timestamps are retained without
backdating.

Artifacts: `results/rhoptry_journal_review_2026_10_07_v2/`, the retained initial
diagnostic and primary discovery/linkage responses. **548** relevant host/source
checks and **625** registry/selection/catalogue/review/docstring/organism checks
passed. A stale generated source script caused two checks to fail initially;
regeneration fixed both, and the complete affected check set passed. Exact version
and missingness comparison has regression coverage. The coherent host-review
commit is pushed to `origin/nightly`. Numerical data and calibration are unchanged.
