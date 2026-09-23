# 0.43 development benchmark artifacts

See the [benchmark report](../../docs/benchmark-0.43.md) for interpretation and
reproduction commands. These are aggregate results from 29 comparisons, not new
gene annotations.

- `summary.csv`: feature-model metrics for both measured outcomes.
- `per_class.csv`: localization class prevalence, precision, recall and average precision.
- `reliability_bins.csv`: held-out top-class confidence and observed accuracy.
- `coverage_strata.csv`: results using fixed published-feature coverage thirds.
- `paired_group_bootstrap.csv`: 500 paired orthogroup resamples of fixed predictions.
- `network_summary.csv`: fixed-graph transductive comparisons.
- `run_contracts.json.gz`: exact settings, split gene IDs, exclusions and versions.
- `manifest.json`: SHA-256 identities for inputs, implementations and generated artifacts.

The bootstrap intervals omit refitting and model-selection uncertainty. Read the
limitations before comparing inductive feature models with transductive networks.
