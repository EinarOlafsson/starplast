# Cell-cycle predictions: none, and that is the result

`search.predictions()` was run on the winning cell-cycle structure — the step the whole search exists
for, which turns a recovered structure into named predictions for the genes inside its pure clusters.

**It produced nothing, at the default 80% purity bar.**

## What the winning structure looks like

    mean F1                 0.489
    blocks                  expression_summary + fitness_screens + published_screens
    n_neighbors / min_dist  15 / 0.25
    min_cluster_size        60
    clusters                3, with 50% of genes left as noise
    labelled genes in it    316 of the 3,000 sampled

    cluster 0    44 labelled    purest phase S at 52%
    cluster 1    50 labelled    purest phase C at 32%
    cluster 2    67 labelled    purest phase C at 49%

## It is not one unlucky run

Across the ten best configurations, the highest purity **any** cluster reaches is **52%**. The bar for
emitting a prediction is 80%.

| mean F1 | clusters | best purity |
|---|---|---|
| 0.489 | 3 | 52% |
| 0.452 | 4 | 51% |
| 0.441 | 6 | 43% |
| 0.440 | 4 | 50% |
| 0.437 | 3 | 41% |

## Why F1 0.489 and purity 52% are both true

F1 balances precision and recall, and here **recall is doing the work**. A cluster captures 55% of the
C-phase genes in the sample, which is real structure — while being only 49% C-phase itself.

For scoring "did the map organise anything", recall counts. For *predicting* a phase for an unlabelled
gene, only precision counts, and there is not enough of it. This is exactly why precision and recall are
never blended into one reported number: a single 0.489 hides that half of what a cluster contains is
something else.

## The honest conclusion

The map recovers cell-cycle phase well enough to say it is not random, and cell cycle is the best
biological signal in this project — better than localisation at 0.228. It does **not** recover it well
enough to tell you where an uncharacterised gene sits in the cycle.

Lowering the purity bar would manufacture predictions from clusters that are half something else. The
bar stayed at its default and the answer is no.
