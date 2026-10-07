# 65.02 — Audit publication identity and citation rates

Status: COMPLETE, 2026-10-07. Supporting publication audit, not dataset admission.

All 162 registry entries are retained in
`results/dataset_selection_2026_10_07/source_publications.csv`. Exact primary
identifier queries resolve 83 recorded PMIDs and 12 recorded DOIs (95 sources).
Dates, counts, annual rates, provider/snapshot and source-link limitations are
recorded. Unknown counts remain unknown, rather than zero.

The remaining 67 source roles/gaps are explicit: 13 local/derived sources without
their own publication; 19 reference resources with unresolved publication lineage;
25 unresolved other publication identities; six accession text matches requiring
origin-paper review; four accession searches without indexed publication hits.
Existence of a recorded PMID is not proof of the registry's source association.

The resumable primary-service audit is `scripts/audit_dataset_selection.py`; its
executed notebook, input/export hashes and 459 used response hashes are frozen
in the result directory. The verification notebook checks every hash and source
identity reconciliation. **576 relevant checks passed, one documentation test
module skipped because optional pdoc is unavailable in the local environment**.
No environment, measurement, strategy, calibration or runtime default changed.

Publication-source lineage review and eligible dataset replacements remain 65.04.
This completion is included in the coherent nightly audit commit/push, under
Einar Olafsson's sole authorship.
