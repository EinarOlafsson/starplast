# 64.05 · Freeze the hold-outs and leakage rules

Status: DONE — ✅, 2026-10-07. Nested partitions, training-only exclusion selection and explicit pair regimes verified.

Parent: [64 · Information space and inference atlas](../open/64_information_space_and_inference_atlas.md).

| Field | Value |
|---|---|
| Depends on | 64.04 |
| Estimated remaining engineering time | 0 h |
| Existing code to inspect | `starplast/strategies.py`, `starplast/track_record.py`, `starplast/search.py`, `starplast/methods.py` |

## Deliverable

Reusable split manifests and exclusion manifests for gene groups, masked values, pairs, datasets and biological contexts.

## Controlled scope

Implement split generation and guard checks; no strategy optimization in this item.

Frozen pilot: reusable training/tuning/calibration/final-test manifests, grouped
entity hold-outs, explicit known-node/cold-node pair regimes, dataset/context/time
group protocols and stage-specific access guards. Freeze two installed parasite
candidate cohorts and exclusion families without fitting or changing runtime tests.
Deliberate group, target-family and derived-layer contamination must fail.
Unknown homology, real pair negatives and external dataset/context ground truth
remain explicit gaps; fixtures test software guards only.

Implement one coherent change at a time. Record the frozen pilot scope before expensive jobs.
If this item needs additional sources, tasks or unplanned algorithm changes, register a bounded
follow-up and its dependency rather than silently widening the item.

## Done when

- [x] New protocol manifests hold out the gene and declared homology group; guards enforce the same split during preprocessing, tuning and calibration. Existing adapters are not certified by this item.
- [x] Pair tests separate edge completion on known nodes from cold-node/generalization tests and state the regime.
- [x] Nested selection has disjoint tuning/calibration/final-test partitions; deliberate target/source/derived-layer contamination is caught.
- [x] Relevant automated checks and any required biological benchmark evidence are recorded, including failures and gaps.
- [x] The result is reviewable, the parent tracker is updated, and the validated change is committed and pushed to nightly.

## Completion evidence

`starplast/splits.py`: immutable content-identified nested partitions, protected
homology/entity/dataset/context/time groups, role-specific fit guards, explicit
transductive joint-feature access, training-only target/source closure, derived-layer
exclusions and separately validated known-node/cold-node pair regimes. Chronological
separation is validated both during generation and manifest reconstruction.

Executed biological candidate pilots in `results/split_protocols_2026_10_07_v2/`:
Toxoplasma HFF growth 7,325 eligible observations in 6,769 declared groups;
train/tune/calibration/test 4,015/1,106/1,103/1,101. Plasmodium insertion index
5,385 observations in 4,619 groups; 2,743/940/917/785. All ten deliberate
leakage probes refused. Source selection cannot inspect held-out target labels.
Both split/exclusion manifests round-trip and retain input/output/cohort identities.
The initial artifact's original code remains archived after a stricter time guard.

532 checks passed, two existing skips, including deliberate leakage, source-family,
endpoint-regime, nested-fit, chronological, organism and docstring checks. Existing
embedding/sklearn warnings are retained, not hidden. The executed verification
checks every input/output hash and Python 3.10 syntax. Commit subject:
**Freeze nested benchmark splits and training-only leakage guards**, on `origin/nightly`.

No model fitted, strategy/data change, benchmark admission or new accuracy claim.
Unknown homology for singleton fallback, assay denominators, real pair negatives,
external dataset/context truth and pretrained representation lineage remain gaps.
The pair examples are synthetic guard fixtures only. Subsequent task-specific
adapters must call these guards on actual inputs; old self-tests and precomputed
maps are not retroactively certified. Continue at 64.06.
