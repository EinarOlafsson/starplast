# 66.04: installed host projection correction, 2026-10-08

The executed migration withdraws 39 ambiguous human rhoptry protein projections
(117 score/beta/FDR cells), retaining all 20,989 host rows and every other installed
cell exactly. All 18,700 unambiguous projected values per field remain unchanged.
The shipped gene evidence retains all 20,010 primary rows with exact original
scores, mapping status and alternative accessions; ambiguous and unmapped genes
remain accessible rather than becoming negatives.

`migration.ipynb` performs the controlled change; `verification.ipynb` independently
replays preservation with `check_exact=True`, exact retained TSV text and input
hashes. `manifest.json` records before/after hashes, original inputs and external
backups. `starplast.host.projection_withdrawals()` validates the 117-row ledger
against source SHA, gene identity, ambiguity alternatives, declared exact numeric
encoding and actually withheld cells. Historical coverage reconciles through only
these verified withdrawals; its original report remains unchanged. The generated
catalogue and slot views now identify gene perturbation and 18,700 projections.

## Precision correction to the preceding builder audit

The earlier builder audit used pandas' default floating-point assertion tolerance;
its phrase “exact installed scores” must not be read as binary source/legacy
identity. This phase records strict source-versus-legacy binary differences among
unambiguous rows: 7,221 score, one beta and four FDR values (maximum absolute
score difference 2.220446049250313e-16). Of these, 236 score rows are neither
binary identical nor the specific 15-significant-digit decimal encoding tested.
Every source/legacy value is retained in the encoding-row parquet files; no
unambiguous value was normalized or changed. These pre-existing transform gaps
remain GT-HOST-TX-02, separate from ambiguity correction.

Thirteen withdrawn cells differ between source and legacy binary representations;
each exactly matches its explicitly declared 15-significant-digit decimal encoding.
The initial exact-source guard correctly refused these before any mutation. Its
code and diagnostic notebook are retained in the parent preflight packet. A later
verification wrapper completed its assertions but failed to write the notebook
because of wrapper namespace scope; the corrected wrapper reads and compares the
existing precision outputs without overwriting them. Both versions are preserved.

## Validation and scope

The initial broad run passed 795 checks with two failures and two existing skips.
Both failures (a missing local paths import and stale generated slot outputs) were
fixed. The extended follow-up passed 575 checks with two existing skips, including
both affected tests, pipeline, calibration/publishing, track record, claims,
scorecards, display/strategy cards/panels/edges, inventory and slot layouts.
Counts overlap and must not be added as unique tests. Other initial passes covered
host/deposit/provenance, ground-truth/leakage, docstrings and organism invariants.

Mouse tables, parasite feature tables/graphs, track records and claims have
unchanged SHA-256 values. No strategy or calibration algorithm changed; host
strategy adapters and host inference packs remain pending. Tutorial builders use
parasite contexts and do not consume these host rhoptry fields. Historical truth
and benchmark snapshots retain their input versions; any new host-dependent
artifact must bind the changed table SHA. No independent biological accuracy or
GTEx/Atlas replacement is claimed. Packaging includes the new parquet/JSON files.

External exact original backups:
`/media/carruthers/mnt3/claude/toxoplasma_projects/datasets/installed_revision_backups/66_04_2026_10_08/`.
Parent immutable snapshots and hashes remain untouched. `packet_sha256.json`
freezes this migration phase and its diagnostic/preflight inputs.
