# 44 — Recipes: a named biological question, its inputs, its holdout, and its control

**Status: open. Proposed by the user 2026-08-18.**

## The change of strategy, and why it is right

Everything so far builds one map from EVERYTHING and holds out one level. That answers *"is this
category redundant given all the others"* — a property of the CATALOGUE. It is a reasonable audit and
it is not a biological question, which is why its output has been hard to read biologically.

A recipe answers *"can THESE measurements predict THIS label"*. That is a hypothesis: it has a wrong
answer, it names what would refute it, and its result is a list of genes rather than a score.

### The user's worked example

> **Which Toxoplasma proteins are on the PV membrane and have a section facing the host?**
>
> UMAP inputs: transcription (all life-cycle stages), translation (all stages), phosphorylation (all
> stages), fitness across all CRISPR screens, gene regulation
> Holdout: LOPIT labels
> Validation holdout: host recognition

## The validation holdout is the important part

A second, independent label that SHOULD map the same way is a positive control inside every run.
Without it a good score says "the map found structure"; with it, "the map found the structure I asked
about" can be told apart from "the map found something else that happens to separate genes".

If PVM-facing is real, host-recognition data should fall on the same clusters as the LOPIT PVM labels.
If the primary holdout recovers and the validation holdout does not, the recovery is probably an
artefact of what was fed in — and that is a result, reported rather than hidden.

## The warning: hand-picking inputs REMOVES a guard

The level sweep closed leakage automatically by holding out a whole class. Choosing inputs by hand
gives that up, and the danger in the worked example is not transcription — it is anything
sequence-derived. Localisation is heavily predictable from signal peptide, TM count and orthoLOPIT
transfer, so a recipe that quietly admits those would "recover" LOPIT by PREDICTING it rather than by
finding biology. The user's input list already avoids them; the builder must enforce it rather than
depend on that.

Therefore, before a recipe runs:

* `excluded_for` runs on the primary holdout AND on the validation holdout, at class scope;
* anything it removes from the chosen inputs is REPORTED, with which mechanism caught it — declared
  family, measured association, or shared experiment;
* a recipe whose inputs still restate a holdout after that is refused, not silently trimmed;
* the validation holdout is checked against the primary holdout too: two labels derived from the same
  experiment are not independent, and a control that is a copy controls nothing.

## What a recipe is

A named, saved, re-runnable object:

    question            free text, and it is the recipe's name
    inputs              category addresses at any level, in any hierarchy
    holdout             the label the map must recover
    validation_holdout  the label that should behave like it, if the answer is real
    excluded            what leakage closure removed, filled in by the run, not by the user
    seed, spec          so a result can be rebuilt exactly

## What running one returns

1. **Mapping quality** — the tuned UMAP and clustering settings, cluster count, coverage, evenness,
   and the degeneracy verdict. A recipe on a degenerate map returns that and stops.
2. **Recovery** — per label and per cluster: precision, recall, F1 for the primary holdout.
3. **The inference** — genes in clusters that are mostly one holdout label but were unlabelled, with
   the cluster's statistics beside each. This is the answer to the question.
4. **The control** — how the validation holdout mapped onto the SAME clusters, and whether it agrees
   with the primary. Reported next to the inference, because it is what says whether to believe it.

## Done when

* A recipe can be written, named, saved, re-run and rebuilt from its seed.
* Leakage closure runs on both holdouts and reports what it removed before anything is built.
* A recipe whose inputs restate its holdout is refused with the reason.
* The worked example above runs end to end and returns genes.
* The validation holdout's agreement is on the same table as the inference it qualifies.
