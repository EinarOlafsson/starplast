# 64.20 · Make labels and protein classes entry points

Status: DONE — 2026-10-08. Acceptance verified; committed and pushed on nightly.

Parent: [64 · Information space and inference atlas](../open/64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.01, 64.02, 64.17, 64.19 |
| Estimated remaining engineering time | 0 h |
| Existing code to inspect | `starplast/track_record.py`, `starplast/slot_tree.py`, `starplast/app.py`, `starplast/guided.py` |

## Deliverable

Search/browse pages for label variables and their values, including localization, curated protein classes and gene sets.

## Controlled scope

Class/label navigation over declared existing annotations; ontology expansion is a separate registered-source task.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] A class page shows measured membership, definition/context, coverage, per-class precision/recall and confusions for each applicable mechanism.
- [x] A label page shows class balance, overall/macro performance, unmeasured classes and the hierarchy where applicable.
- [x] Class→gene→label→source navigation preserves organism and context; identical text in different labels or species is not merged.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

`LabelSpace` and the function-first label/class browser expose all 41 Tg and
30 Pf native labels, their original values, complete profile memberships,
unknown counts, observed definitions, provenance and source/context gaps.
Overlapping classes use the whole original gene universe; known False values
remain distinct from unannotated genes. Missing parent hierarchies and the
unmeasured class universe are explicitly unavailable rather than invented.

All 39 strategy declarations are visible. Available legacy evaluations retain
separate settings, seeds and partitions. Label cards expose overall/macro metrics;
class records expose precision, recall, F1 and the full cohort confusion.
Precision includes false-positive calls from other classes and recall includes
abstentions. Immutable shared scorecard exports retain original supplied values.
Changed tables, contexts and incompatible raw/derived class scopes cannot borrow
recorded accuracy. Class→gene→label→source routes retain the organism and exact
native target/value. Ambiguous source questions clear the old selection, card
and records; the user must choose the intended address.

Canonical executed acceptance packets:
`results/label_space_2026_10_08_v4_tg/` (418 checks, 41 labels,
112 separate proteasome evaluations; 178.64 s, 1849.4 MiB peak) and
`results/label_space_2026_10_08_v4_pf/` (459 checks, 30 labels,
154 separate proteasome evaluations; 181.66 s, 1373.2 MiB peak).
Each packet has 18 verified current input hashes and 20 verified output receipts.
Independent native TP/FP/FN and exact precision/recall calculations reproduce the
original ledger. Actual navigation and context refusal are exercised; membership
and strategy screenshots reviewed. 494 focused backend, widget, route, source,
docstring, provenance and README checks pass. Version 0.54.0 and whitespace pass.

Initial categorical audit failure, successful V2 arithmetic replay (846 checks)
and combined V3 OOM-kill (279.233 s under 1900 MiB) are preserved with original
code and diagnostics. Final organisms ran serially in separate processes under
the same cap; the lease was not raised. Full nightly CI remains a release gate.
Legacy annotation recovery is not independently admitted biological accuracy.
No source promotion, ontology expansion, new strategy fit, calibrated deployment
or host installation is admitted; broader inference profiles remain 64.31.

## Frozen pilot

Installed Tg/Pf annotations only: all declared native labels/classes, all 39 strategy capability entries, and separately recorded legacy evaluation scopes. Final replay checks exact class arithmetic, changed-context refusal and class→gene→label→source navigation. No new fits, ontology or independent biological truth is admitted. Earlier diagnostic and successful audit packets remain immutable.
