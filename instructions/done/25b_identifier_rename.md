# The identifier rename — job 2 of task 25 — DONE 2026-08-12

Job 1 changed user-facing strings and deliberately left the identifiers, because renaming them is a
wide mechanical refactor that wants its own commit with the tests green on each side. This is that
commit.

## What changed

597 sites across the package, the tests, the dataset scripts and the docs:

    colour   -> color        colour_of, colour_mode, set_colour_mode, categorical_colours,
                             unknown_colour, cluster_colours, and every inflection in prose
    normalis -> normaliz     sources.normalise, embedding.normalise, normalised_by, and the prose
    localis  -> localiz      the module, the feature block, the tests, the generated scripts

`starplast/localisation.py` -> `starplast/localization.py`, `tests/test_localisation.py` ->
`tests/test_localization.py`.

Substring rules rather than whole-word ones, because these words inflect: `coloured`, `recolour`,
`normalising`, `normalisation`, `localised`. Archives are out of scope — `results/` and this
directory record what was done and said at the time, and rewriting them would make the record
disagree with itself.

## The part that could have failed silently

The feature block was called `localisation`, and **a block name is stored inside every saved
recipe**. `columns_for` skips a block it does not recognise — it does not raise — so a recipe
carrying the old spelling would not have failed. It would have rebuilt a *different* embedding, one
whole feature block short, under the name of the run it was supposed to reproduce. There are 1,448
saved embeddings on this machine and every one of them predates the rename.

So the rename ships with a translation, `embedding.RENAMED_BLOCKS`, applied in `EmbeddingSpec`'s
`__post_init__` — the one point every route passes through: a stored recipe, a results row saved
last week and clicked today, a hand-written call. Block weights are keyed by block name and are
translated with it. `EmbeddingStore.list()` translates for display too, so one block does not appear
under two spellings in one column and read as two blocks.

Two tests pin it, one for each route, and both assert the block actually contributed columns rather
than merely that the name changed.

## Verified

1,806 tests pass on pandas 2.3.3. Re-checked under the user's interpreter (pandas 3.0.5, real
17 MB cache, offscreen Qt): the renamed block selects `lopit_prob_map`, `lopit_prob_mcmc` and
`lopit_methods_agree`; a pre-rename recipe builds the same 8,140 x 3 matrix; all 1,448 stored
recipes list `localization` and none lists the old spelling; a pre-rename results row rebuilds its
map; and the window opens with `color_of`, `color_mode` and `color_sources()` under their new names.
