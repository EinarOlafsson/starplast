# Structure search, corrected scoring

The same four targets as `../search_2026_08_12/`, re-run after `score_recovery` began excluding labels
that mean *not measured* (`unassigned`, `unknown`, empty, `nan`).

The first run's table was computed under a different rule, so none of its numbers are comparable with
these. It is kept as the evidence that justified the change.

## Why a re-run was needed rather than a recomputation

The exclusion changes which configuration *wins*, not only the winner's score. Once `unassigned` stops
inflating a run's mean, a different combination of blocks and hyperparameters can come out on top — so
the honest comparison required running the whole grid again, 328 runs per target.

## Method

Identical to the first run: 41 dataset combinations × n_neighbors {15, 50} × min_dist {0.0, 0.25} ×
min_cluster_size {25, 60}, on a seeded 3,000-gene subsample (seed 42), median NA policy, rank scaling.
Every scorable run's embedding is stored with its recipe.

## Results

| target | role | before | corrected | change |
|---|---|---|---|---|
| `stage_enriched_derived` | positive control | 0.709 | **0.709** | — |
| `cellcycle_phase` | measured | 0.489 | **0.489** | — |
| `attention_depth` | **negative control** | 0.654 | **0.339** | −0.315 |
| `compartment` | measured | 0.484 | **0.228** | −0.256 |

**The two targets with no absence class did not move at all.** That is the rule behaving exactly as it
should: `cellcycle_phase` is one of G1a/G1b/S/M/C or nothing, and `stage_enriched_derived` is one of
three stages or nothing, so neither ever had a "not measured" class to be rewarded for finding.

## The negative control passes now

Its score fell from 0.654 to **0.339**, below both measured targets, and the per-label table shows why
the first number was never about attention: essentially all of it came from the never-named class at
F1 0.77 over 1,439 genes. With that class excluded, what remains is incidental 0.40, focal 0.27,
substantive 0.15 — weak, which is what a passing negative control looks like.

So the map was **not** substantially organised by study effort among genes that have been studied. It
was separating measured genes from unmeasured ones, and that is a property of missingness rather than of
the literature. The corrected picture is considerably better for the project than the raw 0.654 implied.

## But localisation now sits below the negative control

`compartment` at 0.228 is lower than `attention_depth` at 0.339. Stated plainly: **the map recovers how
much a gene has been studied better than it recovers where the protein is.** That is a real caveat on
any localisation claim and it should travel with one.

Cell cycle at 0.489 is the only target that clears the negative control, and by a comfortable margin.

## The positive control still works

`stage_enriched_derived` reaches 0.709 with oocyst at F1 0.86 — as it must, being a deterministic
function of expression columns the embedding contains. A positive control that stopped scoring high
would mean the pipeline was losing information it had been handed, and every other number would be a
floor rather than an estimate.

## What this does and does not license

    cell-cycle phase   0.489   the one target with a real, checkable signal
    study effort       0.339   the floor any biological claim has to clear
    localisation       0.228   below that floor

And recovery is not prediction: see `../predictions_2026_08_12/`, where the winning cell-cycle structure
turned out to have no cluster pure enough to predict from.

## Files

| file | what |
|---|---|
| `rerun_<target>.csv` | every run's score and settings, with the excluded columns recorded per row |
| `rerunlab_<target>.csv` | per-label precision, recall and F1 for each run |
| `rerun_summary.json` | the winning configuration per target |
