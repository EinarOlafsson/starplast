---
name: held-out-recovery
description: Score whether an embedding, clustering or model has recovered something real, without fooling yourself. Use when holding out a label to test a structure, when a cluster-enrichment or feature-importance result looks strong, or when reporting that a map "predicts" a biological property. Covers circularity guards that survive renamed and derived columns, why absence classes inflate every score, and the negative control that catches what the guards miss.
---

# Held-out recovery without fooling yourself

A structure that recovers a held-out label is the standard evidence that an embedding captured
something real. It is also very easy to get a large number for reasons that have nothing to do with the
thing you are claiming. Every failure below produced a confident, plausible, wrong result in a real
project, and none of them raised an error.

## 1. Name the excluded columns and you will still leak

Excluding "the target column" is not enough. In one run, `compartment` fed the embedding and the
enrichment battery reported its exact twin `lopit_map`, plus its derivations `lopit_unified` and
`compartment_best`, as the *top held-out discoveries* at Cramér's V = 0.96.

**Measure association with the inputs instead of trusting names.** Compute it for every column against
every embedding input and exclude anything above a threshold — 0.8 when deciding what an embedding may
see, which is deliberately stricter than the 0.95 used to flag a result as derived after the fact. A
column at 0.85 leaks nearly as much as an identical one.

Use one scale so a single threshold applies across column kinds:

| pair | measure |
|---|---|
| categorical × categorical | Cramér's V |
| continuous × continuous | absolute Pearson correlation |
| categorical × continuous | **correlation ratio (η)** |

That third row is the one people skip, and skipping it is not a small gap: it makes the guard
structurally blind to *every numeric column*. A categorical label computed from continuous features
scored no association with its own sources, because the comparison was never attempted.

## 2. Measured association cannot catch a joint function

A label that is the argmax of three columns associates with each one individually at ~0.6 — under any
workable threshold — while being fully determined by the three together. No pairwise measure can see
this, however good the statistic.

**So declare provenance as well as measuring it.** A derived column should record which columns it was
computed from, and excluding it should exclude those sources *and the rest of their feature block* —
because the sources rarely stand alone either. In the real case, the strongest association to the
derived label (0.72) was a column that was not one of its declared sources but measured the same
biology.

Measurement catches the undeclared leak (a renamed copy). Declaration catches the joint one. Neither
alone is enough.

## 3. An absence class will carry your score

This is the one that changed a published number.

Labels routinely contain a value meaning *not measured*: `unassigned`, `unknown`, empty, `NaN`. If you
score it as a class, you are rewarding the model for separating **measured** entities from **unmeasured**
ones — and that is usually the strongest signal in any real dataset, because an entity absent from most
assays is absent from most of the feature matrix.

Measured, not argued: localization scored mean F1 **0.484**, and its single best-recovered label was
`unassigned` at 0.39. Excluding that pseudo-label, the same run scores **0.207**. The real classes sat
between 0.21 and 0.35 the whole time.

**Exclude absence values from scoring by default**, and keep it overridable — "does the map separate
measured from unmeasured" is a legitimate question, just not a biological one.

## 4. Run a negative control, and believe it

Pick a label that measures *how much attention the entity has received* rather than what it is — number
of publications, depth of literature mention, how many assays include it. Score it exactly like a real
target.

If it scores comparably to your biological targets, the map is substantially organised by study effort.
In the real case the negative control scored **0.654**, above both measured targets, and the per-label
breakdown was decisive: essentially all of it came from one class — the never-mentioned entities at
F1 0.77 — while the genuine attention tiers scored 0.14–0.35.

That single number reinterpreted every other result in the project. Without it, localization's 0.484
would have been reported as a finding.

Run a **positive control** too: a label you know is a deterministic function of the features. It must
score high. If it does not, the pipeline is losing information it was handed, and every other number is
a floor rather than an estimate.

## 5. Never blend precision and recall

Report both, per label. A cluster that is 100% apicoplast while containing 5% of apicoplast proteins
supports no inference, and a single F1 hides which of the two failed. Weight any summary by label size,
so a structure that isolates one small class does not outrank one that organises the whole dataset.

## 6. Record what was excluded, on every row

Reading a results table months later, "was the target actually held out" is the first question, and it
should not require re-running the search. Put the excluded column list in the output alongside the
score.

## The checklist

- [ ] association with inputs **measured**, not just column names listed
- [ ] mixed categorical/continuous pairs actually compared
- [ ] derived columns **declare** their sources; sources excluded with their block
- [ ] absence classes excluded from scoring
- [ ] negative control (study effort) run and reported
- [ ] positive control (known function of the inputs) run and reported
- [ ] precision and recall separate, weighted by label size
- [ ] excluded set recorded on every output row
- [ ] seed, sample and full recipe stored so a hit can be rebuilt
