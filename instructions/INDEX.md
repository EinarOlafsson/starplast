# Instructions: what is done and what is left

**New here? Read [START_HERE.md](START_HERE.md) first.**

One file per task. This index is the status table; the files carry the reasoning and the
"done when" condition.

## Done

| # | Task |
|---|---|
| 01 | Nineteen gene codes resolved |
| 02 | Six expression datasets integrated |
| 03 | GUI bugs |
| 04 | `projectionMatrix` compatibility — clicking genes |
| 05 | Visual quality |
| 06 | Studies with no supplementary data |
| 07 | Verify shipped columns against GEO |
| 08 | Standalone data and auto-download |
| 09 | Cell-cycle dataset |
| 10 | Search more targets |
| 11 | README and API reference |
| 12 | Test every aspect of the library |
| 13 | Jobs, stopping, and freeing memory |
| 14 | Tooltips, docstrings, and the scoring explainer |
| 15 | The UMAP gallery — grid, scroll, incremental |
| 20 | Validation: `refit`, candidates, orthogonal evidence, per-category scores in Inference |
| 19 | Logging — opt-in, per-level console control |
| 18 | Navigate / Select modes, 2D and 3D gating |

## Open

| # | Task | % | Blocked by |
|---|---|---|---|
| 21 | Wire the scoring objectives into the interface | 100 (done) | — |
| 16 | Automated walk — UMAP, clustering, per-category scoring | 0 | — (15 landed) |
| 17 | Annotation store | 25 | — (20 landed; it must carry 20's numbers) |
| 22 | Cell diagram under the compartment list, colored to match | 0 | a decision on shared organelles |
| 23 | Left panel becomes "color by": clusterings, runs, binned numerics | 0 | — |
| 24 | Forty black-and-white logo drafts | 0 | — |
| 25 | American spelling — attempted, reverted, read the file first | 0 | — |

## The one ordering constraint that is not negotiable — now satisfied

**17 must not ship before 20**, and 20 landed on 2026-08-12. Annotation is easy to build and easy to
build wrongly. Without validation it gives the project a mechanism for manufacturing unvalidated
claims, which is the one thing it exists to prevent — and a candidate list looks identical whether it
is 90% right or 6%.

The constraint does not disappear now that 20 exists; it becomes a requirement on 17. A saved
annotation carries the validated precision and recall for its category, the cluster's composition,
and the enrichment over prevalence — the numbers the Validation tab now computes. Saving a row
without them would be the same failure by a shorter route.

## Standing constraints on all of it

Measurement, inference, absence and annotation never read as one another.

No candidate is shown without the number that says how much to believe it. On the full proteome one
row out of 7,710 (configuration x compartment) reaches F1 0.5 — `PM - peripheral 1`, at exactly
0.500. A cluster that looks pure is a lead to check, not a result.
