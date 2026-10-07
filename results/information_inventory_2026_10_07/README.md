# Offline information-space inventory — 2026-10-07

All **162 registered sources** appear in **164 source/organism/storage-unit rows**;
two sources declare both gene and pair outputs. **16 retained question-specific
refusals** bring the inventory to **180 rows**. No source acquisition, model fitting
or biological accuracy measurement was performed.

| Availability of declared outputs | Rows |
|---|---:|
| Installed | 157 |
| Source-only, without declared projected columns | 1 |
| Unavailable locally | 4 |
| Bridge records present but source attribution unavailable | 2 |
| Rejected for the named question | 16 |

| Installed reference | Rows | Storage unit |
|---|---:|---|
| Toxoplasma | 8,140 | gene |
| Plasmodium | 5,720 | gene |
| Human host reference | 20,989 | protein |
| Mouse host reference | 15,590 | protein |
| Shared metabolite table | 1,296 | metabolite |
| Toxoplasma host bridges | 733 | pair record |
| Plasmodium host bridges | 10 | pair record |

The human/mouse protein references are not admitted host gene spaces. The shared
metabolite denominator is the installed table universe, not proof that both
organisms were assayed for every metabolite. Each source selects its own columns.

Read [census.ipynb](census.ipynb) for the executed reconciliation and interpretation
limits. [inventory.json](inventory.json) and [inventory.parquet](inventory.parquet)
retain source, organism, storage unit, slot contexts, evidence families, local
availability, declared/missing columns, stored values and missingness counters.
[registry.json](registry.json) freezes the source declarations; [manifest.json](manifest.json)
identifies the installed tables, input/code hashes, software versions and exported
artifact hashes. `base_commit` is the parent baseline; exact worktree code is
identified by its input hashes. Inputs were hashed before and after the census.

Coverage means **any stored nonmissing value in a source's declared columns**, out
of the installed rows in that storage unit. Missing outputs or missing tables
have unknown coverage, rather than zero. Blank cells and nonfinite numeric values
are missing; zero and False are retained. Stored text such as “unknown” remains a
stored annotation, not a validated truth label. Multi-column coverage is a union,
and fractions across different sources or units must not be averaged or summed.

The mouse surfaceome has 1,296 stored rows, including 1,146 False boolean cells,
within the 15,590-protein installed reference. This preserves its stored negatives
without turning the other 14,294 proteins into negatives. `false_cells` is a
storage count; assay interpretation remains tied to the registry/source notes.
Only a supplied measured boolean-negative Observation increments
`explicit_negative_records`.

There are **16 source/unit rows without slot context declarations**. They remain
visible with unspecified context; source notes are retained. Catalog contexts are
question facets, not recovered dose/time/replicate metadata. Missing assay causes
cannot be recovered from merged summaries: no source observations were supplied
for this census, so unmeasured/unmapped/QC/detection counters remain unknown.
The implementation accepts existing typed Observations when those causes are
known. Their counts are records, not independent biological samples.

Declared graph outputs reproduce installed layer record counts (the crosslink
layer has 79 pairs). Multi-layer rows count layer records, not unique biological
relationships. Pair coverage has no assayed-universe denominator. Legacy host
bridges have no source IDs, so all 733/10 records remain unattributed instead of
being assigned to one publication. Remote access was not probed; absent local
data is distinct from a confirmed inaccessible remote source. A refusal applies
to its named question, even if the same publication answers another one.

Checks: **525 passed** for inventory, datasets, slots, queries, observations and
docstrings. Initial synthetic checks passed, but review found an incorrect graph
endpoint suffix in the first census. Correcting it to the shipped `__a`/`__b`
schema and adding a real-graph test exposed a test-helper key error; that was
fixed and the final checks passed. Earlier snapshots are kept outside the tree
under `/tmp/starplast-inventory-first-attempt-20261007` and
`/tmp/starplast-inventory-pre-manifest-version-20261007`.

Reproduce in a fresh directory under a RAM lease and memory cap:

```bash
systemd-run --user --scope -p MemoryMax=4G env QT_QPA_PLATFORM=offscreen \
  OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2 \
  /home/carruthers/anaconda3/envs/starplast/bin/python \
  scripts/build_information_inventory.py --out results/information_inventory_NEW
```

Source-level context, mapping lineage and missingness migration continue in 64.03.
This inventory does not promote any evidence into ground truth or change data,
strategies, shipped claims or calibration.
