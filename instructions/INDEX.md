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
| 16 | Automated walk — UMAP, clustering, per-category scoring |
| 17 | Annotation store |
| 23 | Left panel becomes "color by" |
| 22 | Cell diagram under the compartment list |
| 25 | American spelling in user-facing text |
| 24 | Forty black-and-white logo drafts |
| 21 | Scoring objectives in the interface |
| 26 | Save and load results, per tab and all at once |
| 27 | Import data, with the preprocessing offered |
| 25b | The identifier rename — job 2 of the spelling task |
| 28 | Every module at 100% coverage |
| 29 | Three more leaks in the circularity guard |
| 30 | Fold in datasets already present and expand missing biological axes |
| 31 | Replace file-shaped feature blocks with biological question slots |
| 32 | GPU paths for k-means, DBSCAN and t-SNE |
| 33 | Crossed-factor findings and fragmentation diagnostic |
| 34 | Repair two metric columns and the CLI finding count |
| 35 | Simplify lighting, make point modes distinct, and add volumetric ray tracing |
| 36 | GPU PBR sphere points and stable GPU density-ray tracing |

## Open

| # | Task |
|---|---|
| 37 | Continuous flashlight, material lab, and ray-rendering comparison |
| 38 | Build a provenance-first pan-Apicomplexan dataset archive |

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
