# 68.01 — Browse functions and every available label in Discoveries

Status: COMPLETE, 2026-10-08; acceptance verified, committed/pushed to nightly.
Parent: [68 — Discoveries function coverage](../open/68_discoveries_function_coverage.md).

| Item | Percent done | Time left | Description |
|---|---:|---|---|
| 68.01 | ✅ | 0 h | Browse all available labels and functional classes, memberships and inference/test coverage in Discoveries |

Acceptance evidence: `results/discoveries_label_coverage_2026_10_08_v4/`.
Actual executed source replay verifies every variable and all 21,787/18,879
functional gene/class memberships for Toxoplasma/Plasmodium. Their menus expand
from 3/2 labels to 41/30, spanning function, phenotype, stage, localization,
structure and other categorical annotations/flags. InterPro/Pfam/EC multi-valued
classes retain overlaps, source descriptions and unknown missingness; replacement
EC identifiers mentioned in notes are not fabricated memberships. All claim,
recipe and legacy-record targets remain visible even with no generated claims.

Qt search/class/member navigation, real application gene opening, label/class
scorecard links, existing claim filters/export/map behaviour and visible verifier
counts verified. Final **413 checks pass; one optional pdoc module skipped**.
Earlier wider application/navigation/help checks passed except one new fixture's
scorecard wording assertion; the final run fixes/verifies it. Earlier executed
source-format diagnostics remain immutable. Source/claims/recipes/track-record
hashes unchanged; no strategy fits or changes to calibration/graphs/layouts.

No functional inferred claims or functional held-out benchmark rows added.
Known annotations and legacy source-recovery tests remain separate from
independent biological accuracy. Missing metrics/tests are explicit and
annotation absence is not a verified negative. Functional source review,
benchmarks, calibration and precomputed inference remain **68.02**, now the
next user-prioritized controlled task. This completes browser coverage only.
