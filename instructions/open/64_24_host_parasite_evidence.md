# 64.24 · Explore typed organism-to-host connections

Status: OPEN — 0%, 2026-10-07. Existing foundations are reused; completion of this acceptance contract is not yet verified.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.03, 64.21, 64.22, 64.23 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `starplast/host.py`, `starplast/organisms.py`, `starplast/star_edges.py`, `starplast/star_map.py`, `starplast/datasets.py` |

## Deliverable

Explicit cross-space evidence links with organism-qualified endpoints, relation type, context and source quality.

## Controlled scope

Expose existing verified bridges only; new cross-species inference is separately benchmarked.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Measured interaction, orthology, predicted interaction and infection-response links retain distinct labels and provenance.
- [ ] Opening a counterpart selects its own space; tables and measurement denominators remain species-specific.
- [ ] A bridge has mapping/coverage evidence and can be explored without an inference-accuracy claim; inference through it requires a registered benchmark.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Pending: commit, checks, artifact/benchmark identities, measured population, limitations and dated result.
