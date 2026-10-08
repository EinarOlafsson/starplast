# 64.17 · Expose one reusable scorecard component

Status: COMPLETE, 2026-10-08; acceptance verified, committed/pushed to nightly.
Parent: [64 · Information space and inference atlas](../open/64_information_space_and_inference_atlas.md).

| Item | Percent done | Time left | Description |
|---|---:|---|---|
| 64.17 | ✅ | 0 h | Reusable evidence, performance and individual-outcome scorecards |

The component displays supplied result quality, coverage, uncertainty, baseline,
benchmark population and freshness. Compact cards expand to the original test
record. Definitions, source, split, sample sizes, calibration, controls, failures
and retained outcomes have working detail routes. Supplied calibration metadata
is preserved; missing calibration, dates and source versions are explicit.
Evidence quality and individual outcomes do not acquire performance probabilities.
The snapshot/export retains exact supplied counts, binary64 values and nulls.

Acceptance evidence:

- `results/scorecard_task_views_2026_10_08/audit.ipynb`: 539 software presentation
  checks across six task contracts, ten synthetic aggregate cards, authored
  controls, unavailable states and distinct set/interval definitions.
- `results/desktop_scorecards_2026_10_08_v2/audit.ipynb`: 1,573 checks across
  14 cards in two real Qt hosts. Every metric/detail route, exact snapshot/export,
  expansion, row/outcome anchor dispatch, clearing and source dispatch verified;
  40 output receipts and current code hashes match. Runtime 0.68 seconds and
  measured process peak 150.0 MiB, with an external 1,400 MiB cap. HTTPS browser
  opening was mocked; it is dispatch evidence, not a source-download test.
- `results/functional_coverage_ui_2026_10_08_v4/audit.ipynb`: real Discoveries
  source/card/class/outcome navigation preserves original EC reference-recovery
  outcomes and all original source memberships. The representative human source
  card retains 58,988 canonical genes and all mapping/admission gaps. All
  51 nested output receipts and input hashes match.
- 90 focused shared-view/browser/host/class/Discoveries checks pass; package
  version and whitespace checks pass. Earlier 738 integration checks remain
  recorded in their original scope.

The first desktop attempt is preserved as a diagnostic: sorted JSON mapping keys
changed incidental supplemental-metric list order. The corrected audit compares
complete metrics by identity, with values and definitions exact; it does not
relax numerical comparisons. Screenshots and original code are retained.

This closes the reusable component and representative host presentation scope.
Independent biological accuracy, full strategy adapters, calibrated deployment,
host installation, other entry points and release acceptance remain separate open
items. Missing benchmark data are displayed as unavailable; synthetic fixtures
and annotation recovery do not establish biological validity.
