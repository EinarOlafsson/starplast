# 64.31 · Compare inferences for a label or protein class

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.17, 64.20, 64.26, 64.27, 64.30 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `starplast/track_record.py`, `starplast/strategy_panel.py`, `starplast/app.py`, `starplast/guided.py` |

## Deliverable

Class/label inference profiles listing measured members, predicted candidates, all applicable methods and per-class/target performance.

## Controlled scope

Class and label profiles over declared annotation schemas; no silent new ontology mappings.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Known membership and inferred membership are distinguished; class-specific precision/recall, priors/lift and missed/confused classes remain accessible.
- [ ] Numeric contexts use their own quantities and value-task cards rather than being collapsed into categorical calls.
- [ ] Query-dependent set results are tied to the exact input members/context; a gene-set scorecard cannot be mistaken for global method accuracy.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
