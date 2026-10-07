# 64.02 · Inventory the available information space

Status: DONE — ✅, 2026-10-07. Offline inventory and installed-table census verified.

Parent: [64 · Information space and inference atlas](../open/64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.01 |
| Estimated remaining engineering time | 0 h |
| Existing code to inspect | `starplast/datasets.py`, `starplast/slots.py`, `starplast/slot_tree.py`, `starplast/organisms.py` |

## Deliverable

A generated inventory of dataset × organism × unit × biological context × evidence family, with availability and coverage.

## Controlled scope

Inventory existing source records and installed tables; acquisition becomes a separate source-specific task.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] Every registered source appears once per declared species/unit, including unavailable and rejected sources.
- [x] Coverage denominators and missingness distinguish unmeasured, unmapped, inaccessible and explicit negative measurements.
- [x] Inventory row counts reconcile with the registry and installed tables for the existing parasite and host reference data.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Frozen pilot: all 162 registered sources, the existing slot-generator refusals and installed
Tg/Pf gene, human/mouse protein, metabolite, graph and bridge tables. No downloads or strategy
runs. `starplast/inventory.py` uses existing registry and slot declarations;
`scripts/build_information_inventory.py` freezes a fresh snapshot and executed notebook.

Checks: **525 passed** (Python 3.12.13, CPU, offscreen Qt, 4 GB cap per job) for inventory,
datasets, slots, queries, observations and docstrings. Python 3.10 syntax checks passed for
the new modules/builder. The executed notebook and all declared input/output hashes were
verified after generation. Whitespace checks passed.

Artifact: [information_inventory_2026_10_07](../../results/information_inventory_2026_10_07/README.md),
with JSON/Parquet inventory, frozen registry, manifest and executed reconciliation notebook.
All 162 sources reconcile to 164 source/organism/storage-unit rows (two sources have gene and
pair outputs), plus 16 scoped refusals: 180 total. Statuses: 157 installed, one source-only,
four unavailable locally, two unattributed bridge sets and 16 refused questions. Installed
rows reconcile to Tg 8,140, Pf 5,720, Hs proteins 20,989, Mm proteins 15,590, metabolites 1,296
and bridge records 733/10. The mouse surfaceome retains 1,146 stored False cells among 1,296
stored rows without converting the rest of the installed reference into negatives.

Review found an incorrect endpoint suffix in the first census. It was corrected to the actual
`__a`/`__b` schema, a real shipped-graph regression was added, and a test-helper key error was
fixed. Final checks and the regenerated census pass; prior attempts are retained in `/tmp`
as recorded in the artifact README. No incorrect census is published in the tree.

This card's commit, **Inventory registered evidence and reconcile installed storage coverage**,
is pushed to `origin/nightly`. The parent tracker and session handoff identify 64.03 as next.

Coverage is explicitly storage coverage, not an experimental sampling fraction or biological
accuracy. Summary-table null causes remain unknown; optional typed source observations can
separate unmeasured, ambiguous mapping, detection limit, QC and known negatives. No remote
access is claimed or probed: local absence is distinct from confirmed remote inaccessibility.
Legacy bridges lack source IDs, so their records remain unattributed. Sixteen source/unit rows
have no slot context declaration; they remain visible with unspecified context and retained
source notes. These are evidence/lineage gaps for 64.03, not invented measurements. Host
references are protein rows, not yet host gene spaces. Pair counts have no experimental-universe
denominator. No source, strategy, calibration or shipped claim was changed or promoted.
