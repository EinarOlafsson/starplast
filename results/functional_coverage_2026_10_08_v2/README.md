# Functional source, target and strategy coverage

`census.ipynb` executes a bounded offline inventory and exactly replays every
address in `matrix.json`. Current functional column counts were independently
computed from the installed node tables and exactly matched the previous full
source-membership audit. Historical evidence was scanned in bounded batches:
all **1,038,372** original label-call rows reconcile to **168** scoped summaries.
None describes an installed functional target.

The matrix retains **880** organism/source-label/derived-target/strategy/task
addresses: 240 Toxoplasma, 320 Plasmodium, 160 human and 160 mouse. It inventories
**11** actual functional source labels. Host source slots are explicitly
uninventoried and their gene inference adapters unsupported. No host membership
or other-organism accuracy is inferred from parasite data.

Each row distinguishes declaration-level applicability, known source annotation,
legacy rows, verified reference recovery, independent biological testing,
calibration and unknown-gene deployment. A derived complete EC major-class profile
is separate from the original multi-valued EC field. Unadapted domain/EC sources
explicitly require a multi-valued target adapter. Declaration-only rows never
assert runnable inputs or performance.

Exactly **one** address has a checksum-verified reference-recovery result: the
Toxoplasma EC major-profile feature-kNN label-call pilot. Its metrics retain the
original artifact/scope/settings and never move to Plasmodium, hosts, other source
labels or other tasks. Recovery source populations must fit the current known
annotation inventory and match its gene universe. There are **zero** independently
admitted biological tests, calibrated targets or deployment targets. Pooled
accuracy and pooled annotation-gene counts remain null.

Backend API: `functional_coverage.build(catalogues_by_organism, capabilities=...,
benchmarks=..., legacy_records=..., organisms=...)`. Inputs are current catalogue
rows, declared capabilities, integrity-validated packaged benchmark objects and
exact scoped legacy row counts. Output has `rows`, `source_label_inventory` and
`summary`. Future biological/calibration/deployment admission needs the existing
typed source/benchmark/role contracts; loose flags are deliberately unsupported.

**11 focused checks pass** under 400 MB. Source/code hashes, snapshots and the
executed notebook are retained. The first census is immutable evidence of the
initial schema, before its source-population invariant was strengthened.
