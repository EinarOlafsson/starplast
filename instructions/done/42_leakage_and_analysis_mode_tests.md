# 42 — Test analysis mode for information leakage, aspect by aspect

**Status: DONE 2026-08-18.** `tests/test_leakage_attacks.py`, plus the panel controls driven in
`tests/test_questions_tab.py`.

## The attacks, and what they caught

Every test is named for the attack rather than the function, and three of them found real defects
rather than confirming a guard.

* **A one-hot categorical cannot also be the target.** The categorical route into the matrix is not a
  block, so the block-dropping guard does not cover it; fed as a one-hot feature, `compartment`
  separates the clusters it is then scored against by construction. Attacked and closed, with the
  mirror test that a categorical which is NOT the target survives — a guard that removed every
  categorical would pass the first test and be useless.
* **Nothing is both held out and used, asserted PER CATEGORY.** The first version of this test
  compared block names against a hierarchy path and "found" a leak that was its own arithmetic; the
  fixed version derives the held blocks from `categories_at`, exactly as the sweep does.
* **Rows stay aligned after every missing-value policy.** A policy that drops genes changes which
  genes the labels describe, and scoring a subset against the full table aligns cluster 3 with the
  wrong genes while every number downstream still looks reasonable.
* **The negative controls** — a shuffled label must not recover, a random block must not improve any
  objective — landed earlier in the day and now have a sibling for every new method: the classifier
  gets the same shuffled-label control, plus one asserting that forty columns of pure noise are not
  "recovered", which is what proves the partition is genuinely out of fold.

## 6, the facet source: built, measured, and REFUSED — with two real bugs found on the way

The instruction asked for facets to come from the dataset registry rather than the hand-written
label. That was implemented at two strengths and both were rejected by what they produced:

* **registry-first** rewrote 142 facets, and `transcription · tachyzoite` — a slot whose name AND
  context both say tachyzoite — became *bradyzoite*. A registry entry describes a DATASET, and a
  stage series names every stage it covers.
* **registry-as-fallback** put `ring`, a PLASMODIUM stage, on the Toxoplasma `shared orthogroup`,
  `shared domain` and `fold confidence` slots. Shared and relational slots declare patterns that
  resolve to registry entries spanning both organisms.

A check that comes back backwards means refuse the source. But writing the invariant that caught it
— **no slot may carry another organism's life-cycle stage** — then found two defects in the SHIPPED
catalogue that had nothing to do with the registry:

1. **`ring` matched inside `conferring`**, so the Toxoplasma slot `resistance-conferring mutation`
   carried a Plasmodium blood stage. Matching is now whole-word.
2. **`sexual` matched inside `asexual`**, so **20 Plasmodium slots that measure the ASEXUAL blood
   stage were filed under the sexual stage** — the biological opposite. Any hold-out addressed at
   the sexual stage was silently taking the asexual slots with it.

Also fixed: the implied-stage vocabulary is Toxoplasma assay language (HFF, BMDM, peritoneum), and
the Plasmodium arm mirrors those slots, so six Plasmodium slots carried `tachyzoite (implied)`.
Stage vocabularies are now per organism.

**One ambiguity is recorded rather than decided by the matcher:** `Pf fitness · in vivo` has the
context "humanised mouse or CHMI", and its host facet moved from `human` to `mouse` when substring
matching went away. Both readings are real — CHMI is infection of actual humans, and a humanised
mouse is a mouse carrying human red cells, which is what the parasite actually inhabits. The right
fix is 42.6's other half, enriching the slot's context so the label states which it means.

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
