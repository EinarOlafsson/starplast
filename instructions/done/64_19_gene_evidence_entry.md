# 64.19 · Make gene lookup an evidence entry point

Status: DONE — 2026-10-08. Acceptance verified and committed/pushed on nightly.

Parent: [64 · Information space and inference atlas](../open/64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.01, 64.03, 64.17, 64.18 |
| Estimated remaining engineering time | 0 h |
| Existing code to inspect | `starplast/identity.py`, `starplast/app.py`, `starplast/organisms.py`, `starplast/slot_tree.py` |

## Deliverable

A gene page with alias lookup, all measured evidence grouped by biological question, source coverage and links to class/label pages.

## Controlled scope

Lookup and evidence navigation for available spaces; inference profiles arrive in item 30.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] Canonical and alias searches select the correct organism/entity and expose ambiguity instead of guessing.
- [x] Every available measured column for the gene is reachable, with context, units and provenance; missing values keep their meaning.
- [x] The initial view is condensed and can expand to all evidence; a measured label and its class both have working scorecard routes.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Exact canonical/table/index aliases use the shared organism-qualified resolver.
Installed identity indices retain their full canonical universe so table subsets
cannot erase collisions. Mapping artifacts and table aliases carry content
identities. Ambiguous text/alias matches require a candidate click; unavailable
index targets retain explicit gaps instead of borrowing another gene’s evidence.

The gene card opens a condensed browser grouped by all matching biological
questions; expansion reaches every original column, including missing values.
Zero, False, exact numerical values and nested source values are preserved.
Contexts, units, declared origins and provenance gaps remain accessible. Scoped
source selection opens the shared immutable source card and the exact dataset
address. Existing label/class routes require the unchanged installed table.

The canonical final packet is `results/gene_evidence_2026_10_08_v5/`: 1,246 executed checks, 22 current input/code
hashes and 20 verified output receipts. Runtime 102.68 seconds, measured process
peak 804.1 MiB under a 1,900 MiB cap. Final screenshots inspected; evidence
columns widened with exact full-value tooltips. The replay covers Tg GRA16 (442 columns)
and Pf PF3D7_0102900 (167 columns), with separate Pf AMA1 alias verification.
Full displayed label/class records equal the original ledger result. Earlier
failed audits are preserved: the first used Tg’s compartment column for Pf;
the second requested an AMA1 class absent from the legacy ledger. Selection
was changed to an explicitly measured class with an existing record, preserving
AMA1’s original annotation and its unscored-class gap; no biology was invented.

473 component/integration/docstring/search/menu checks pass, one existing real-GL
context skip. Four widget checks pass after the final width change. Native alias
ambiguity includes absent candidates, stale selection clears for unavailable
targets, and replacing an evidence window releases its predecessor. The exact
unchanged-table guard prevents importing frozen class accuracy into a changed
context. Version 0.54.0 and whitespace checks pass; replacement full CI remains
a release gate. Earlier successful replays are preserved with their historical
code identities. No source download, new strategy fit, independent
biological validation, host installation or calibrated deployment is admitted.
