# 64.18 · Browse the selected organism's datasets

Status: DONE — 2026-10-08. Acceptance verified and committed/pushed on nightly.

Parent: [64 · Information space and inference atlas](../open/64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.02, 64.03, 64.17 |
| Estimated remaining engineering time | 0 h |
| Existing code to inspect | `starplast/slot_tree.py`, `starplast/datasets.py`, `starplast/app.py`, `starplast/organisms.py` |

## Deliverable

A browser over dataset families, biological labels, contexts, coverage and source records for the selected space.

## Controlled scope

Existing datasets/slots in one browser route; no new sources or extra top-level dock by default.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] Filters by organism/host, unit, family, measured context and availability reproduce inventory counts.
- [x] Each source opens its coverage/provenance card and measured entities; unavailable evidence explains its recorded gap.
- [x] Measured, transferred, predicted and user-imported evidence remain identifiable; source links work from a compact summary.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

The Tools → Datasets route shares the existing Slots window. Its qualified inventory
covers all 162 registered sources and 180 scoped rows, including 16 question-specific
refusals. Filters retain organism/host, unit, family, context and availability.
Source cards preserve declared lineage and gaps; immutable exports contain no invented
accuracy. Original gene, host-protein, metabolite and source-attributed pair values
remain accessible in bounded pages. Gene navigation opens the matching organism;
other units retain their own tables. Session imports remain separate sources.
Processed-cache availability is separate from raw-assay provenance, and the slot
viewer no longer confuses stored availability with experimental evidence grade.

Canonical acceptance is the executed notebook in
`results/dataset_browser_2026_10_08_v8/`: 3,900 executed checks, all 180 scoped
rows, 27 current input/code hashes and 24 output receipts verified. Runtime was
98.05 seconds; measured process peak 1,082.8 MiB under a 1,900 MiB cap.
Original numerical values are compared exactly, including zero, False and missingness.
Final screenshots were inspected. All original gene rows remain reachable by paging.

Final regressions: 505 passed, one existing real-GL-context skip. Six observed CI
regressions plus native interval-record checks pass (24 checks). Discoveries now
explains all controls; planted value self-tests retain point checks and explicitly
leave unsupported interval metrics NaN. Package version 0.54.0 and whitespace checks
pass. Full remote CI remains a release gate; no release is claimed.
Earlier successful replays and diagnostic packets are preserved. The first replay
failed while writing a code receipt due to a shadowed variable. A later combined
Qt test process exited 139; its native cause remains unproven. Subsequent checks
use the shared pytest-qt application fixture. Another combined run aborted with
`malloc_consolidate(): unaligned fastbin chunk detected`; explicit native iterator
cleanup before its owning tree destruction was added to the new slot test, and
the same 505-check regression scope then passed. The native crash cause remains
unproven; observed diagnostics are retained. A 400 MiB worker explanation audit
was OOM-killed; serial validation uses the documented 1,900 MiB cap. Invocation
failures are retained separately; the final audit runs as a Python module.

This action admits no new biological truth, sources, calibrated strategy fits,
host gene installation or dataset replacements. Missing provenance, units,
association, mapping and calibration remain explicit gaps.
