# Three more leaks in the circularity guard — DONE 2026-08-13

Not on anyone's list. Found by running a real search in the real window and reading what won, which
is the same way the display bugs and the settings pollution were found: the suite was green
throughout.

## What was wrong

Recovering `compartment`, the winning configuration was built on the `localization` block —
`lopit_prob_map`, `lopit_prob_mcmc`, `lopit_methods_agree`. Those are hyperLOPIT's own posteriors
and its methods-agreement flag, and the label being recovered is hyperLOPIT's assignment.

**Measured on the shipped cache, those three columns recover `compartment` at mean F1 0.259 — above
interactions (0.192), protein features (0.171), and every expression and fitness block.** Three
columns beating eighteen RNA columns and eight CRISPR screens is not the map finding biology.

Neither existing mechanism could see it:

* **measured association** puts them at 0.29 and 0.43, far under the 0.8 threshold, because a
  posterior does not restate *which* compartment a protein is in;
* **declared derivation** sees nothing, because the label is not computed from the posterior.

They are the same experiment's other outputs. That is a third thing, and asking the same question of
the other targets found two more of it:

* **The negative control was embedding its own inputs.** `attention_depth` is `np.select` over
  `n_papers_focal / substantive / incidental`; none was declared, so none was excluded. They score
  0.25–0.43 against the tiering — a count is not a restatement of a tier while determining it
  completely. A control that sees its own inputs scores too high, and every measured target then
  looks worse than it is by exactly that much.
* **The same quantity, measured another way.** `ortholopit_label` — localization transferred from
  *P. falciparum* and *C. parvum* orthologs — sits at 0.72 against `compartment`; `lit_tier` at 0.71
  against `attention_depth`. Both just under the threshold, which is what a near-copy does.

## What the guard is now

Four tests, each of which was a bug before it was a rule:

    measured association    the undeclared copy: a renamed column, a second inference
    declared derivation     the joint function: a label that is the argmax of three columns
    shared provenance       the same experiment's OTHER outputs, read from the registry
    shared quantity         the same thing estimated another way (`search.SAME_QUANTITY`)

`validate.circularity_error` applies the same tests by name, because the Validation tab was refusing
only the label column and would otherwise have put a validated-looking number on a map built from
that label's own experiment — in the tab whose whole job is to say how much to believe a cluster.

Closing the leak empties whole combinations out of a sweep, so the sweep now reports what it skipped
and why. "8 runs" quietly becoming 7 is the same failure class as a fetch that reports DONE having
downloaded nothing.

## What it costs

The registry now lists all ten columns hyperLOPIT produces rather than five, and declares
`attention_depth` as derived from the three counts. A test asserts every column named in a family
exists in the cache, so a rename cannot silently stop the guard guarding.

## Verified

1,865 tests pass with every module at 100%. Re-checked under pandas 3.0.5 on the real cache: 14
columns excluded for `compartment`, 7 for `attention_depth`, 2 for `cellcycle_phase`; the
`localization` combination is dropped from a sweep and said so; the Validation tab refuses a map
built on `lopit_prob_map` and accepts one built on expression and fitness.

**The published numbers turn out to be unaffected, and that is worth stating precisely.** The
full-proteome battery sweeps a hard-coded six-block base that contains neither `localization` nor
`literature`, so no published row could use the leaking columns. Re-run under the fixed rules, all four
targets reproduce exactly — 0.675, 0.398, 0.322, 0.193 —
`results/full_proteome_2026_08_13_provenance/`.

The leak was reachable from the **interface**, which sweeps every block that has columns: run from
the Search tab on a 600-gene subsample, the localization block won outright, 0.259 against 0.192 for
the best measurement block. A user would have believed it.
