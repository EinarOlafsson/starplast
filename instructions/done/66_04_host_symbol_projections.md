# 66.04 — Correct ambiguous host symbol projections

Status: COMPLETE, 2026-10-08. Acceptance evidence validated; committed and pushed
with the coherent nightly correction.

The host builder withholds symbols mapping to multiple reviewed accessions while
preserving Ensembl behavior and the public mapping schema. The installed human
rhoptry table/cache now withholds only 39 verified ambiguous projections (117
cells); all 20,989 host rows and all other installed values remain exactly intact.
All 20,010 original gene rows are shipped with original scores and mapping
statuses/alternatives, including 39 ambiguous and 1,271 unmapped rows. The
18,700 unambiguous projections are retained. Gene perturbation is explicitly
distinguished from protein measurement in the source catalogue.

The per-cell ledger is checked against immutable source SHA, row/gene identity,
accession alternatives, declared exact decimal encoding and missing installed
cells. Historical non-loss tests reconcile only verified withdrawals; arbitrary
numeric or mapping tampering is refused. Original installed files are backed up
outside the repository, with before/after hashes. Parent frozen snapshots remain
unchanged. New host-dependent artifacts must bind the new source/table hashes.

Evidence: [migration review](../../results/host_symbol_mapping_2026_10_08/migration/README.md),
executed migration and strict verification notebooks, manifest, preserved source,
encoding-row diagnostics and code packet. The review corrects an earlier pandas
tolerance-based “exact” description: pre-existing binary source/cache differences
are explicitly retained; migration preservation uses strict exact comparisons.
GT-HOST-TX-02 retains the distinct expression/mapping/encoding gaps.

Validation: 575 extended follow-up checks passed, two existing skips; initial
795 other passes recorded with two failures subsequently resolved (overlapping
counts are not additive). Host/deposit/provenance/leakage, pipeline, calibration,
track record/claims/scorecards, display and slot checks are covered. Mouse and
parasite feature/graph/track-record/claim SHA values are unchanged. No affected
parasite tutorial input changed; packaging includes both new evidence files.

Host gene spaces, inference packs, independent biological ground-truth admission
and literature-wide dataset replacements remain their separate open actions.
