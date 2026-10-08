# 64.25 · Build a resumable precomputation pipeline

Status: OPEN — 35%, 2026-10-08. Serial resumable framework and one real frozen-result replay verified; production feature/map/model/deployment wiring remains open.

Parent: [64 · Information space and inference atlas](64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.08, 64.10, 64.11, 64.12, 64.13, 64.14, 64.15, 64.16 |
| Estimated remaining engineering time | 4–8 h |
| Existing code to inspect | `scripts/build_track_record.py`, `scripts/build_claims.py`, `starplast/searches.py`, `starplast/packs.py` |

## Deliverable

A dependency-aware job plan for reusable features/operators/maps/models, hold-out records, deployment outputs and scorecard summaries.

## Controlled scope

Pipeline plus tiny fixtures and one bounded pilot; large sweeps begin only with a measured budget and RAM lease.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [ ] Immutable input snapshots and job fingerprints prevent mixed-source runs; interrupted jobs resume by validated completed partitions.
- [ ] A changed dataset/target/settings fingerprint invalidates only declared downstream outputs; cancellation leaves valid earlier partitions.
- [ ] Pilot jobs report runtime/RAM/storage, set a compute budget and package budget, and preserve audit logs and failures.
- [ ] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [ ] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

Verified bounded partition, 2026-10-08: `starplast/precompute_jobs.py` supplies typed
immutable jobs/upstream bindings, acyclic plans, transitive recipe fingerprints,
complete parent artifact identities, explicit byte snapshots, serial execution,
writer locking, an atomic checksummed journal and current-spec/split/payload resume.
Interrupted complete manifests may recover; partial or corrupt outputs remain
preserved and are retried separately. A failed job blocks descendants while
independent jobs continue. Cancellation and explicit attempt/time/current-RSS/
artifact-storage/package limits retain earlier completed partitions. Parent disk
and separately loaded in-memory contents are revalidated at checkpoints.

Executed `results/precompute_ec_views_2026_10_08/pilot.ipynb` consumes only the
verified FN-EC-01 V4 frozen artifact
`349317afd8f24b860d94ae1e11280a6130bf451700b02b9d2145fab85fc99162`.
Three real view jobs retain all182 test rows/scores, recompute matching profile
scorecards, then recompute seven overlapping member-class cards. All six comparisons
are exact, with original split/679 fitting entities/exclusions/unresolved grade.
No fitting, new prediction, source acquisition or biological admission occurs.
The predeclared first-run max_jobs=1 stops after the row partition; resume builds
only two unfinished jobs; an empty-builder replay validates all three cached outputs.
Measured stage times0.350/0.957/0.328s, process peak187.3MiB (not isolated per-job
peak), typed artifacts166,244bytes, complete evidence688,151bytes. External400MiB
cap; per-run120s/3attempt/4MiBstorage/2MiBpackage/currentRSS400MiB budgets retained.

83 focused framework/artifact/split checks passed under400MiB, including34 runner
cases. All34 also passed in the final broad Qt-context run with deterministic
fixture RSS; actual cap enforcement remains external. Cases include invalidation,
failure/blocked descendants, cancellation, corrupt/symlink/partial/interrupted
artifacts, writer lock, snapshot and parent mutation, budget and resume boundaries.
The Linux reader separates executed-process VmHWM from current VmRSS; inherited
launcher getrusage high-water values do not falsely reject a new small worker.
Initial fixture typo, measurement error and missing-method docstring are fixed.

Remaining controlled work: wire production feature/operator/map/model/held-out/
deployment builders to the declared adapter partitions and role-specific input
snapshots, measure those budgets and expose their cache/failure coverage. Scientific
adapters64.10–64.16 remain unfinished; these fixtures and recovery-view jobs do not
complete the broader pipeline or atlas. Builders use cooperative checkpoints and
need an external memory/time cap for preemption. Current-input identities must
be refreshed by the caller or bound through explicit local snapshots.

Next software partition, 2026-10-08: `precompute_features.NumericFeatureBuilder`
provides a bounded intermediate numeric operator with exact ordered source
snapshots, training-only column selection/ranks, frozen held-out ECDF, explicit
missingness and unavailable-state blocking. Split, target, exclusion, recipe,
builder and imported guard-code identities are required dependencies. Complex
values and attention aliases are refused. It produces no labels, confidence or
biological accuracy and checks selected numeric columns, not graph/map lineage.
59 pure feature/runner checks pass under400MB; all738combined integration checks
pass after making Qt-process RSS fixtures deterministic. The first executed
production operator pilot remains pending; overall action stays35%.
