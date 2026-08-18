# 44 — Recipes: a named biological question, its inputs, its holdout, and its control

**Status: DONE 2026-08-18.** `starplast/recipes.py`, 49 tests, 100% coverage, no `pragma`.

## What shipped

`Recipe` (question, inputs as category addresses, holdout, validation_holdout, seed, spec,
min_precision) -> `close()` returning a `Closure` -> `run()` returning a `RecipeResult` ->
`RecipeStore` for save/reload. `excluded_for` was split into `excluded_detail`, which returns
`{column: mechanism}`, so the closure report is produced by the same code that makes the decision
rather than by a second implementation that could explain a column it did not exclude.

All five acceptance criteria are met, and the run that met the fourth is worth reading twice.

## Six things measurement changed, none of them anticipated when this was written

**1. Class scope, not `target_family`.** The instruction says closure runs "at class scope" and that
turns out to be load-bearing rather than a detail. `target_family` -- the sweep's default -- closes
the family estimating the same quantity and deliberately KEEPS sequence predictors; there is a test
in `test_search.py` asserting exactly that. Asked "which proteins are on the PVM", a recipe scoped
that way accepted `membrane topology` as an input and would have recovered hyperLOPIT by predicting
it from signal peptide and TM count. That is the failure this instruction exists to prevent, and the
first draft shipped it. `scope="biology"` closes the first two levels of what a measurement is ABOUT,
which is what a holdout is a statement about.

**2. A tuned map can be beautiful and useless.** `map_quality` rewards evenness, so tuning chose 121
clusters of ~40 genes: even, non-degenerate, and it named ZERO genes, because no cluster held fifteen
LABELLED ones at 80% purity. A cluster must be big enough that the scoring floor is reachable, and
how big depends on what fraction of genes carry a label at all -- 47% for `compartment`, so 32 genes
minimum. The floor is derived from the holdout's COVERAGE and never from which labels recover:
choosing hyperparameters by the latter is choosing them by the answer.

**3. The confidence threshold belongs to the recipe, not the run.** With the floor in place the map
went to 52 clusters and still named nobody, because the binding constraint was purity, not size: the
purest cluster on this catalogue is 0.74. `min_precision` is now a saved field of the recipe -- a
threshold chosen after seeing the answer is a threshold fitted to it -- and an empty answer reports
how close it came, because "no genes" and "no cluster was 80% pure, the purest was 74%" are different
findings and only the second says whether the question failed or the threshold did.

**4. `unassigned` is a gene waiting to be named, not a class.** hyperLOPIT writes it where it could
not place a protein. Left as a string it is the commonest label in the table, clusters get
"recovered" as unassigned, and -- worse -- the genes this whole module exists to name are counted as
already labelled and never predicted.

**5. A cell-cycle question may not be asked of RNA.** `cellcycle_phase` comes from single-parasite RNA
sequencing and lives under `gene expression > RNA abundance`, so closure empties the entire RNA branch
when it is the holdout. Correct, and it constrains instruction 45: the cell-cycle questions must be
built from protein, fitness or regulation inputs.

**6. `iedb_epitope_count` cannot serve as a validation holdout at all.** 221 labelled genes across
8,140 means no cluster reaches the 15 needed to evaluate one. In the worked example the control was
silent on every cluster the inference came from -- one cluster had 13 labelled genes, the other none.
A control that can never reach the floor is not a control, and the whole `host recognition` axis has
no usable one in this build. That is a finding for the acquisition campaign, and a question dropped
in 45 with a reason rather than shipped with a control that abstains.

## The worked example, run end to end

The user's question, its five input branches, holdout `compartment`, control `iedb_epitope_count`:

* closure removed one column from the inputs the user asked for -- `n_lactylation_sites`, caught by
  measured association with the CONTROL, which no one would have predicted;
* 52 clusters, 55.6% of genes clustered, evenness 0.97, not degenerate;
* 22 labels scored, mean F1 0.222, best 0.505 (60S ribosome). Dense granules 0.73 precision /
  0.27 recall, rhoptries 2 at 0.33/0.62;
* at the default 0.80 confidence: **no genes**, purest scoreable cluster 0.74;
* at 0.70, declared in the recipe: **12 genes** -- 11 mitochondrion-membranes, 1 dense granule --
  each carrying its cluster's precision;
* the control **abstained on all of them**, for the reason in 6 above.

So the machinery answers, and the honest reading of its answer is that this question is not answered
by this data at this confidence: the genes it names are mostly mitochondrial rather than at the host
interface, and nothing corroborates them. That is the result the design is meant to produce -- an
inference with a control beside it saying whether to believe it -- rather than a failure of it.

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
