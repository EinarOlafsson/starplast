# Offer the scoring objectives in the Search and Clusters tabs

`starplast/objectives.py` is written, tested and 100% covered, and is currently reachable only from
Python. It should be selectable in the interface, because choosing the objective is the most
consequential decision in a hyperparameter walk and it is presently made for the user.

## What to add

A chooser in the Search tab, and in the automated walk when that exists, offering:

    mean_precision        most clusters are one label
    mean_recall           most labels are in one cluster
    mean_f1               labels and clusters correspond, on average
    best_precision        at least one cluster is mostly one label
    best_f1               at least one label has mostly its own cluster
    n_recovered           how many labels clear a threshold
    v_measure             the whole clustering agrees with the whole labelling
    precision_at_recall   the purest cluster still worth annotating from

Plus: the weighting choice (macro / size), an optional single category to optimise for, and the
`min_recall` floor for `precision_at_recall`.

`OBJECTIVES` maps each name to its one-line description -- use it for the tooltips so they cannot
drift from the implementation.

## What must be shown beside the score

Coverage, cluster count, and chance-corrected agreement (`objectives.agreement`). Every objective
here has a degenerate maximiser, and ARI near zero beside a high score is how you see one. The
explainer is already in Help; the chooser should link to it.

Default to `macro` weighting. `search.score_recovery` weights by label size, so nucleus-chromatin at
769 genes dominates and dense granules at 167 barely registers -- anyone hunting a rare compartment
has been optimising against themselves.

## Done when

A search can be run under any objective, the result records which one was used, and two runs under
different objectives on the same data pick visibly different winners.
