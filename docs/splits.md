# Benchmark hold-outs and source guards

`starplast.splits` freezes four disjoint roles: `train`, `tune`, `calibration`
and final `test`. Explicit protected groups stay together. Gene benchmarks use
homology groups; missing known groups can be declared singleton entities, which
does not establish absence of homology. Dataset/context groups protect those
boundaries and make no additional homology claim. Chronological time groups must
use comparable ISO timestamps; equal times cannot be split across partitions.

```python
from starplast import organisms as O
from starplast.splits import make_split, make_exclusions, write_split

split = make_split(
    eligible_gene_ids, declared_homology_groups,
    organism=O.TOXOPLASMA, benchmark_id=benchmark.benchmark_id, seed=17,
    feature_access="transductive",
)
training = nodes[nodes.gene_id.isin(split.entities("train"))]
excluded = make_exclusions(
    training, benchmark.target, benchmark_id=benchmark.benchmark_id, split=split,
)
excluded.guard_inputs(columns=selected_features, layers=selected_layers)
split.guard_fit("scaling", scaler_fit_gene_ids)
split.guard_fit("setting_selection", tuning_gene_ids)
split.guard_fit("calibration", calibration_gene_ids)
write_split("results/new_split.json", split, excluded)
```

Fitted imputation, scaling, representation, feature selection and initial model
fitting may read training entities only. Setting selection reads training/tuning;
explicit final refitting uses training/tuning, while calibration stays separate.
Calibration reads its own partition. No fitting or selection stage reads final
test labels. An unregistered stage or entity is refused.

Joint feature/graph construction over hold-out entities requires declared
transductive access, checked by `guard_joint_features`. This permission does not
permit fitting on their labels or fitting a scaler on them. Ordinary pointwise
prediction may read held-out features without updating a fitted transform.

Source exclusions reuse Starplast's target-family closure, including source
siblings and layers derived from excluded evidence. Empirical correlations used
to discover exclusions read **training labels only**. The exclusion manifest pins
those entities and the split identity. Passing test rows to exclusion selection
is rejected. Reusing an exclusion manifest selected under another split is also
rejected. Unknown source ancestry remains a provenance gap.

Relationship protocols have separate meanings:

- `known_node_pairs` reserves endpoint-covering training edges, then partitions
  the remaining observed/control pairs. Every held-out endpoint occurs in
  training. Sparse graphs can lack enough edges for nested validation.
- `cold_node_pairs` applies the node/homology partition and retains pairs whose
  endpoints belong to the same role. Both test endpoints are unseen in training.
  Crossing pairs are recorded as excluded, rather than silently discarded.

Neither operation invents negative edges. The caller must provide a reviewed
observed/control universe and distinguish assayed negatives from unobserved pairs.
The current pair examples are synthetic guard fixtures, not biological truth.

Immutable split files carry content identities. Readers validate role separation,
group boundaries and the training cohort used for source selection. Adapters
must call these guards for their **actual** inputs and retain manifests with
their outputs. Merely creating a manifest does not certify an existing self-test,
precomputed map or pretrained representation. No shipped strategy behavior is
changed in this item; adapter integration remains in subsequent action cards.

Executed pilots are in `results/split_protocols_2026_10_07_v2/`. They freeze eligible
cohorts from the ground-truth candidate registry, verify all ten deliberate
leakage refusals and record partition/group counts. Assayed populations and
independent benchmark admission remain unresolved. No model was fitted and no
new accuracy or capacity claim was produced.
