# 64.08 · Version inference and benchmark artifacts

Status: DONE — ✅, 2026-10-08. Common manifest, data-only writer/reader and dependency invalidation verified.

Parent: [64 · Information space and inference atlas](../open/64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.03, 64.05, 64.07 |
| Estimated remaining engineering time | 0 h |
| Existing code to inspect | `starplast/strategies.py`, `starplast/results.py`, `starplast/searches.py`, `starplast/packs.py`, `starplast/paths.py` |

## Deliverable

A common manifest for outputs, row-level hold-outs, models and scorecards keyed by organism, table/graph/label hashes, code version, settings, evidence exclusions and seed.

## Controlled scope

Schema, writer/reader and invalidation contract for one small artifact of each output task.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] Known-gene held-out outputs and unknown-gene deployment outputs have distinct roles and partition IDs.
- [x] Changed inputs invalidate dependent artifacts; loaders refuse mismatched gene order, organism, target or model/split identity.
- [x] Saved artifacts round-trip with provenance, status, support/calibrated confidence and evaluation scope intact.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

`starplast/artifacts.py` provides schema-1 immutable data-only manifests keyed by typed query/organism/target, strategy/task/output semantics, ordered entities, code version/hashes, table/graph/truth/model/calibration dependencies, settings/seed, exclusion identity, fit population/role, partition/split, evaluation context and truth/negative semantics, status/gaps and confidence scope. Held-out/test and unknown/deployment roles are distinct; known benchmark entities cannot masquerade as unknown deployment. Frozen split validation checks exact ordered test cohort, fit access, benchmark identity and inductive/transductive regime.

Readers require a freshly constructed current expected spec and verify the full manifest and payload hashes/size. Mismatched order, target, settings, split/model/partition or inputs are refused. Changed inputs invalidate declared artifact/model/calibration descendants. Safe finite JSON payloads preserve null metrics; symlinks, unsafe names, nonfinite data and overwriting existing directories are refused. Calibrated confidence requires a calibration identity/applicability scope and a compatible output kind.

Evidence: `results/artifact_contracts_2026_10_08_v2/`, with executed pilot and verification notebooks. **Six scorecard task types and five storage roles**, ten exact round trips, all ten invalidated by changed code and zero by unchanged code; deliberate scope changes refused. **840 relevant artifact/capability/split/truth/provenance/docstring/organism checks passed**. Initial integration incorrectly called a module-level fit guard; the actual split method fixes it, and final checks pass. Python 3.10 syntax and every final input/output/child-artifact identity verified. Original diagnostic implementation retained.

Scope remains manifest/writer/reader/invalidation contracts. Synthetic model metadata is not a fitted estimator; no biological accuracy, new algorithm, numeric dataset or calibration change is claimed. Callers must verify live input hashes; pair/replication endpoint/selection regimes and task-specific row validation remain 64.10–64.15, aggregation is 64.09 and precompute execution is 64.25. Documentation is API-linked and included in the hosted docs build. The coherent validated change is committed and pushed to nightly.
