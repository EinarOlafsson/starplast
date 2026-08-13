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

## Open

**Nothing.** Every task written down here has landed (2026-08-12, v0.16.0).

Two things are known and deliberately not on this list, so that "open" stays a list of work someone
asked for rather than a list of everything imaginable:

- **The identifier rename** — `colour_of`, `colour_mode`, `set_colour_mode`, `categorical_colours`,
  `unknown_colour`, `normalise`, `localisation.py`. Task 25 asked for user-facing text only and says
  this is a separate, deliberate refactor with the tests green on each side. See
  `done/25_american_spelling.md`.
- **The Cryptosporidium VEuPathDB fetch**, which downloaded nothing and reported success. It is
  recorded in `../HANDOFF.md` under "Known outstanding bug", where it belongs: it is a data-fetch
  defect rather than a piece of application work.

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
