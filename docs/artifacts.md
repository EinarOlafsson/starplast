# Inference artifact identities

`starplast.artifacts` stores immutable, data-only inference snapshots alongside
the existing CSV exports. Every `ArtifactSpec` records a typed query, declared
strategy/task/output meaning, target, ordered population, role, partition,
settings, seed, code version, source/table/graph/model hashes, exclusions,
fitting population, evaluation scope, status and confidence meaning.

```python
from starplast.artifacts import write_artifact, read_artifact

content_id = write_artifact(folder, spec, {"rows.json": rows}, split=split)
restored = read_artifact(folder, expected=current_spec, split=current_split)
```

Construct `current_spec` from freshly verified current inputs. A cache cannot
verify external changes when its caller simply reuses yesterday's expected
hashes. File identities can be obtained with `provenance.SourceFile.inspect`
and verified before building the expected spec. Dependency names retain the
source/version association; checksums alone do not establish biological origin.

The five storage roles are `held_out`, `deployment`, `fitted_model`, `scorecard`
and `reusable_base`. Held-out outputs and scorecards require frozen truth,
benchmark and split identities, plus the complete ordered outer-test cohort.
Unknown-entity deployment has a different evaluation partition and identifies
its model and fitting population. Known benchmark entities cannot be saved as
unknown deployment. Split validation checks fitting access, benchmark identity
and inductive/transductive regime; it does not certify the code that performed
a fit. Pair and replication adapters must separately validate their endpoint
and discovery/validation regimes.

Readers require an exact expected spec, validate the envelope and every payload
hash/size, and refuse mismatched organism, target, entity order, source/settings,
partition, split or model identity. Files remain relative and portable; external
payload symlinks and unsafe names are refused. Payloads are finite JSON data;
missing or unavailable metrics use null. No executable model pickle is loaded.
The writer never overwrites an existing directory and writes its completion
manifest last. An interrupted partial directory stays available for inspection
and cannot be loaded as a completed artifact.

Method support and calibrated confidence are separate. Calibrated probabilities,
sets and intervals require a calibration identity and applicability scope, and
must match the output kind. A ranking or geometric score cannot become a gene
probability. This schema preserves those declarations; calibration quality and
row-level interpretation still require the subsequent benchmark adapters.

`invalidated(artifacts, changed_inputs)` reports changed-input artifacts and
all declared descendants. Named current inputs are `(kind, name): sha256`.
Upstream artifact, model and calibration dependencies pin complete content IDs;
a changed result invalidates descendants even if its settings key is unchanged.
Unlisted inputs are assumed unchanged, so builders must inspect every input
relevant to their declared dependency graph. The later precompute scheduler
uses this contract to decide which partitions can be resumed.

The executed pilot `results/artifact_contracts_2026_10_08_v2/` round-trips all
six scorecard tasks and five storage roles. Its 10 artifacts, model metadata
and predictions are synthetic schema fixtures: no estimator was fitted and no
biological performance was measured. Task-specific row adapters and scorecard
aggregation follow under instruction 64.
