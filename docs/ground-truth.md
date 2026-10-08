# Ground-truth registry

`starplast.ground_truth` records organism-qualified benchmark candidates and
their validation gaps. Each entry declares the target/task, gene or protein unit,
evidence grade and basis, source IDs, frozen values and entity-universe hashes,
eligibility mask, context, quantity units, exclusions and negative semantics.

The current census covers installed parasite targets and human/mouse protein
fields. It does **not** admit independent biological ground truth or change the
shipped self-tests. Source-specific admission and split protocols remain open.
An available installed value is an observation candidate, not certification that
the source assayed the entire table. Measured assay populations remain unknown.

Grades distinguish direct perturbation/abundance assays, curation, orthology
transfer, predictor output, derived quantities and synthetic controls. Unreviewed
fields remain `unresolved`. Native LOPIT classifications are predictor outputs
from measured fractionation profiles; accepting them as independent localization
truth requires another experiment. `compartment_best` mixes native calls and
transfer fallback. In-vivo fitness composites multiply fitness by significance;
AlphaFold confidence is a model output. Predictor agreement cannot substitute for
independent biological accuracy. Synthetic fixtures test these software contracts
and never admit a biological benchmark.

```python
import json
from pathlib import Path
from starplast.ground_truth import read_registry

entries, strategies = read_registry(
    "results/ground_truth_registry_2026_10_07_v2/registry.json"
)
entry = entries[0]
ids = json.loads(Path(entry.universe_file.path).read_text())
eligible = entry.mask(ids)  # rejects reordered or changed entity universes
print(entry.evidence_grade, entry.eligible_population, entry.measured_population)
print(entry.gaps)
```

Missing, unknown and nonfinite values are excluded without imputation. Zero and
False remain stored results; neither automatically means a certified biological
negative. No negative class is currently admitted without assay-specific review.
Known ambiguity and missing assay denominators remain visible. The legacy host
rhoptry protein projections include 39 ambiguous symbol assignments; their
candidate cohort is blocked from admission pending instruction 66.04.

Every one of the 39 strategy tasks has references or an explicit gap in each of
the four requested species. References are **candidate** tests. Relationship
strategies do not inherit gene-label benchmarks: they need pair truth, endpoint
regimes and an assayed negative universe. Replication requires disjoint findings
from discovery and validation; a split-half self-test is not external replication.
Hosts currently have protein references rather than admitted host gene contexts.

Snapshots refuse overwrites. The registry reader validates duplicate addresses,
species/task mismatches, packed-mask counts and padding, checksums and ordered
cohort identity. Files are pinned to absolute paths in this local review artifact;
relocation needs an explicit verified migration. Reproduce the census with:

```bash
python scripts/build_ground_truth_registry.py --out results/ground_truth_NEW
```

The executed notebook, target census, strategy gaps and manifest retain the
evidence. This registry supplies the contracts for subsequent benchmark splits,
records and scorecards; it supplies no new accuracy claim.
