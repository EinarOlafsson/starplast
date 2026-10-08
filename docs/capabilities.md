# Strategy capabilities

`starplast.capabilities.catalog()` declares all 39 current strategies, their
allowed query kinds, deployment outputs, evidence requirements, parameters,
benchmark units/tasks, underlying technique tests and precomputable components.
A new catalogue strategy or technique requires a reviewed declaration.

```python
from starplast.capabilities import resolve

plan = resolve(query, "feature_knn", evidence, settings={"target": "compartment"})
print(plan.applicable, plan.reason, plan.outputs, plan.validation_gaps)
```

`query` is a typed `starplast.query.Query`. `evidence` is an explicit
`EvidenceState` inspected from an admitted gene space **after query-specific
source exclusions**. It records the entity universe, available context,
categorical/numeric columns, exact class values, permitted features and typed
layers. It does not infer availability from a file's name or a protein reference.
The current catalogue adapters cover the two existing parasite gene spaces;
human/mouse gene adapters remain pending.

`resolve` returns a plan or a named reason: unsupported query/output, missing
space, endpoint/context/universe mismatch, missing column/class/evidence,
conflicting targets/seeds, invalid settings or insufficient seeds. It never runs
an inference, guesses a required target, selects an unspecified layer or
manufactures probabilities. Defaults not explicitly supplied remain unresolved
until the context adapter binds them. An applicable plan can still lack enough
training examples, independent biological truth or calibrated support; these
limits remain visible. The later fitting adapters must enforce split guards.

Deployment outputs differ from some legacy self-test task names:

| Example | Actual meaning |
|---|---|
| Recoverability atlas | Class recovery summary |
| Block ablation | Evidence-family prediction loss when removed |
| Link prediction, attention correction, unwritten links | Pair ranking; absent edges are not verified negatives |
| Set enrichment | Feature associations plus gene ranking |
| Masked imputation | Transformed matrix percentiles, not native assay values |
| Condition shift | Rank-scale residual after baseline adjustment |
| Ortholog transfer | Numeric estimate or categorical call according to target type |
| Conformal calls/values | Label sets or intervals; test coverage and size together |

`Capability.validate_output` rejects a changed kind, entity unit or score
meaning. Geometric distances, recovery F1, vote support, model scores and
relationship ranks cannot be relabeled as calibrated gene-label probabilities.

Every technique has a declared direct prediction, pipeline ablation or
known-truth/null-control role. These are **test requirements**, with status
`declared_not_measured`; composite strategy accuracy is not attributed to each
component. Candidate benchmark links must match organism, effective target,
task and biological unit. Matching a candidate does not admit it or establish
context/split applicability. Derived condition residuals require a separately
registered transform and truth target.

The immutable pilot at `results/capability_contracts_2026_10_08_v2/` records
39 declarations, 40 techniques, 78 native plans and 624 query/strategy routes
on synthetic evidence fixtures. It measures software routing, not biological
accuracy or inference capacity. Task-specific execution, truth admission,
precomputation and the interface follow under instruction 64.
