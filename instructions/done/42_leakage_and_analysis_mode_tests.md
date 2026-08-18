# 42 — Test analysis mode for information leakage, aspect by aspect

**Status: open. Requested 2026-08-18 by the user.**

## Why

Analysis mode is where this program makes its claims. Every other part shows data; this part asserts
that a structure is real — and the only thing separating "this cluster is 71% IMC" from a circular
statement is whether the map was built on something that restates IMC.

The guards exist and several are tested. What does not exist is a test suite organised around the
QUESTION rather than around the module: leakage has already got in three separate ways, each blind to
the others, and each was found by accident rather than by a test that was looking for it.

**A guard nobody attacks is a guard nobody has checked.** These tests must try to leak and fail to,
rather than confirm that a function returns what it returned yesterday.

## What to test

### 1. Leakage into the embedding

* A column that RESTATES a held-out target — renamed, rescaled, a second inference over the same
  data — must not survive `excluded_for`. Build such a column deliberately and assert it is excluded.
* `target_family` closure across every slot in the family, and across species for transfers: hold out
  a family, assert no member of it reached the matrix.
* The three mechanisms in `excluded_for` — measured association, declared family, same-experiment —
  each attacked SEPARATELY, since each is blind to what the others catch.
* One-hot categorical inputs: `compartment` fed in as a feature must exclude compartment as a target.

### 2. Leakage into the score

* `battery` marks a used column `used` and a near-copy `derived`; only `held_out` is evidence.
  Feed it a copy of an embedding input under a new name and assert it is not scored as held out.
* The sweep must never hand `battery` a column that is both held out and used — asserted per
  category, not just once.
* A clustering that cannot support a claim must not produce a score (`clustering.degenerate`).

### 3. Leakage through the rows

* A policy that drops genes changes WHICH genes the labels describe. Assert labels and truth are
  aligned on the same index after every na_policy, and that a subset never silently scores against
  the full table.

### 4. Every analysis-mode control does what it says

Each of the seven tabs, driven headlessly: the control changes the spec, the spec changes the matrix,
and the matrix changes the result. A control that is wired to nothing is worse than a missing one.

### 5. The negative controls

* A shuffled target must NOT recover. If it does, the scoring is measuring the clustering's shape
  rather than the label.
* A random column added as a block must not improve any objective beyond noise.
* Holding out EVERYTHING must produce no score rather than a default one.

## Done when

* Each guard has a test that tries to defeat it and fails, named for the attack rather than the
  function.
* The shuffled-target and random-block negative controls are in the suite and pass.
* Every analysis-mode control is driven headlessly at least once.
* 100% coverage of the analysis path, no `pragma`.

## 6. The facets must come from the metadata, not from the label

**Added 2026-08-18, from the user: "all information is not saved in the experiment name necessarily."**

The faceted trees built that day derive stage, host, tissue and condition from the slot's `name` and
`context` — two short strings written by hand. That is the wrong source and it silently bounds how
much can be encoded: a facet nobody happened to type into a label cannot be extracted from it, and the
tree then looks complete while missing exactly the information a hold-out needs.

The authoritative metadata is already in `starplast/datasets.py`: PMID, accession, level, kind,
coverage, and the note. Plus the underlying files themselves — sample names, GEO series titles,
PRIDE deposit descriptions — which carry stage, strain, host cell and timepoint far more reliably
than a slot label does.

To do:

* Derive facets from the REGISTRY entry backing each slot, falling back to the label only where the
  registry is silent, and record which source each facet came from.
* Where a facet is genuinely absent from both, say so in the tree (`stage-unspecified`) rather than
  inferring — and keep the existing `(implied)` marking for the cases where the dataset's stage is
  not in doubt but was never written down.
* Enrich the slot contexts themselves so the label stops being lossy: every slot's context should
  name its stage, host, tissue and condition where those are known.
* A test that every facet a slot's registry entry knows about appears somewhere in its tree address.
  That is the check that stops the tree quietly under-describing the data.
