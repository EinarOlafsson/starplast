# 64.39 · Export the evidence and inference trail

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.03, 64.08, 64.17, 64.30, 64.31, 64.35 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `starplast/results.py`, `starplast/report.py`, `starplast/datasets.py`, `starplast/strategies.py` |

## Deliverable

Exportable result tables, scorecards, test records, figures and a concise methods/citation/provenance bundle for the selected query.

## Controlled scope

One common export bundle and adapters for the new profiles; reuse instruction 60's methods/figure work.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] The bundle identifies organism/entity/context, source and artifact hashes, method/settings, benchmark/split and whether outputs are measured/inferred.
- [ ] CSV/table and vector figure exports reconcile with the visible selected data; resolved source citations trace to verified registry records.
- [ ] A reloaded export recovers the same query and scorecard evidence; missing sources stay explicit instead of generating invented citations.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
