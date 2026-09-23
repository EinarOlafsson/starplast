# From gene maps to testable gene predictions

Review of Starplast at `ac9019b`, 23 September 2026. This assessment comes from the
current code and bundled tables, with small synthetic checks of scoring behaviour.
It does not establish which model performs best on parasite biology: that requires
a controlled benchmark. The implementation changes below are proposals.

The goal makes sense: collect evidence about genes and proteins, preserve what each
measurement means, and use the combined evidence to propose experiments. I would
keep UMAP as an exploration tool and build the main inference workflow around
specific biological questions, predictions, and evaluation.

## What already exists

The package is substantially more than a UMAP viewer. The current bundled
Toxoplasma table has 8,140 rows and 393 columns, including 360 numeric columns;
Plasmodium has 5,720 rows and 123 columns, including 97 numeric columns. These are
column counts, not counts of independent biological traits. There are 128 dataset
registry entries, which include computed layers as well as source datasets.

| Capability | Current implementation | Next improvement |
|---|---|---|
| Source provenance and gene identifiers | `datasets.py`, `identity.py`, source archive checksums | Connect each observation and derived feature to a versioned source record |
| Biological context | `slots.py`: organism, question, context, entity type, combination policy | Keep replicate and condition detail before producing gene-level features |
| Balanced feature blocks | `embedding.py`: scaling, missingness policies, block variance normalization | Use this machinery consistently, with task-specific feature selection |
| Multiple prediction methods | `methods.py`: logistic regression, boosted trees, single-layer propagation, multiplex communities | Compare them with the same task, data, splits, and prediction metrics |
| Leakage controls | `search.py`, `recipes.py`: target families, derived inputs, excluded edge layers | Perform label-dependent decisions within training folds |
| Validation and candidates | `validate.py`, `holdout_cv.py`, candidate tables and exports | Make one consistent evaluation and prediction contract reach every user-facing result |
| Map quality and rankings | `metrics.py`: trustworthiness, neighbourhood purity, AUROC/AUPRC | Separate projection quality, biological recovery, and prospective prediction quality |

Code references: [embedding](https://github.com/EinarOlafsson/starplast/blob/ac9019b/starplast/embedding.py),
[methods](https://github.com/EinarOlafsson/starplast/blob/ac9019b/starplast/methods.py),
[recipes](https://github.com/EinarOlafsson/starplast/blob/ac9019b/starplast/recipes.py),
[nested validation](https://github.com/EinarOlafsson/starplast/blob/ac9019b/starplast/holdout_cv.py).

## Is the UMAP design rational?

Yes, as a way to explore similarity under a stated selection of features. A
neighbourhood can suggest a useful candidate. It does not by itself identify a
function, an interaction, or a probability that a gene has a trait. An organism
also has no single universal similarity: two genes can be close in expression and
far apart in sequence, localization, or fitness.

The current `recipes.build_map` clusters projected coordinates, and
`EmbeddingSpec` defaults to three dimensions. The headless nested-CV experiment
also constructs three-dimensional UMAPs. This makes a display decision part of
the inference method. I would compare feature-space neighbours, PCA/factor-space
models, and higher-dimensional UMAPs, then display their predictions on a separate
2D or 3D map. Choose the analysis dimension inside validation, rather than assuming
three or prescribing an arbitrary universal number.

UMAP can be useful before clustering. Its own documentation warns that it can
alter densities and split apparent clusters, and recommends exploring embedding
dimensions beyond the display dimension. That argues for a benchmark rather than
removing UMAP. [UMAP clustering documentation](https://umap-learn.readthedocs.io/en/latest/clustering.html).

A specific simplification is needed: there are currently two embedding builders.
`build_graph.embed` uses the legacy feature selection, median fill, z-scoring, and
optional compartment one-hot block weighted by 0.5. For the current Tg table it
selects 16 numeric features. `embedding.build_matrix` has the newer slot resolution,
configurable policies, and block normalization; its default selects 27 numeric
columns before processing. Thus collecting hundreds of columns does not mean the
opening map uses all of them. The graph archives contain `xyz` and edge arrays,
without a stored layout recipe or executed-backend record in the archive itself.

Build the shipped map through the same recipe system as interactive maps, saving
the data snapshot, ordered gene IDs, feature transformations, actual method,
library versions, seed, and fallback outcome. A localization-informed map is fine
for browsing; it cannot independently validate localization prediction.

## Correctness priorities before adding more algorithms

1. **Score classifier outputs as predictions with fixed class meanings.**
   `recipes.build_model` currently passes out-of-fold class codes to
   `search.score_recovery`, which finds the best cluster for each true label.
   That is a descriptive partition score, not classifier accuracy. In a synthetic
   60-gene example with every A called B and every B called A, classification
   accuracy is 0 while partition recovery mean F1 is 1. This does not measure the
   app's biological accuracy; it demonstrates why these metrics cannot share a
   meaning. Use fixed-label out-of-fold precision/recall, confusion matrices and
   rankings. If remapping classes is desired, learn that mapping on training
   data and freeze it before testing.

2. **Abstain when the graph supplies no evidence.** `Propagator.predict` uses
   argmax even when all propagated scores are zero. A synthetic disconnected node
   receives class 0 with scores `[0, 0]`. The downstream enrichment filter may
   reject some groups, but it is not a per-node guarantee. Return an explicit
   unknown/no-support result for these genes.

3. **Separate transductive and inductive evaluation.** Inferring missing labels
   within a fixed known proteome can legitimately use unlabelled genes' feature
   vectors and graph structure. Predicting new genes, strains, or future data is
   a different test. `build_model` builds some transformed matrices before its
   supervised folds; target-dependent exclusion also examines the provided table.
   Fit preprocessing and empirical feature decisions inside each training fold
   for inductive claims. Static exclusion of known target-derived inputs can
   remain a declared rule. [Scikit-learn's leakage guidance](https://scikit-learn.org/stable/common_pitfalls.html).

4. **Do not present cluster purity as a calibrated per-gene probability.**
   `recipes.infer` attaches cluster precision to every candidate in that cluster.
   This is group-level evidence from the labelled subset. It need not apply to a
   poorly measured gene at the edge of the group. Report support, held-out
   performance, disagreement, and missingness; calibrate model probabilities on
   separate validation data where labels permit. Check calibration separately in
   sparse and well-measured genes, since the unlabelled population may differ.

5. **Make fallback behaviour explicit.** Both embedding builders catch UMAP
   failures and can return PCA coordinates. The configurable GPU UMAP call also
   omits the recipe's `metric` argument while the CPU call includes it. Save the
   actual implementation and effective parameters; inference runs should fail
   clearly or record an explicitly accepted fallback.

The existing nested structure-selection helper is a useful foundation, not a
replacement for reviewing how candidate structures, input exclusions, and class
calls were obtained upstream. Random stratified folds alone also do not establish
performance on unrelated protein families or future studies.

## Structure evidence before flattening it

A gene-by-feature matrix is useful as a model input. Keep a richer observation
store from which each task builds that matrix. A minimal observation needs:

`entity ID + entity version + trait + value/unit + organism/strain + stage +`
`host/context + perturbation/dose/time + assay/replicate + uncertainty +`
`evidence status + source accession/version + source location + derivation`

Gene, transcript, protein isoform, residue/site, allele, protein pair, host gene,
and metabolite observations need distinct identities. A phosphosite position or
host-parasite interaction should retain that identity even when a task derives a
per-gene summary. The slot system already recognizes several of these units; extend
it rather than replacing it with an untyped matrix or immediately adopting a
large graph database. Parquet tables plus explicit relations can be sufficient.

Keep measured evidence, computed features, orthology transfers, and model
predictions distinguishable. Preserve disagreements between studies. Distinguish
not assayed, below detection, failed QC, ambiguous mapping, and measured negative;
none of these should silently become an ordinary zero. Avoid treating repeated
representations of one experiment as independent supporting evidence.

The Tg builder includes phosphosite count despite it being missing for 85.6% of
genes in the current cache. That does not prove the map is wrong, but it makes
coverage-stratified evaluation necessary. A model handling NaNs does not by itself
solve the selection bias in which proteins were assayed.

## Additions that could improve inference

| Addition | Why it could help | How to establish whether it helps |
|---|---|---|
| Frozen protein-sequence embeddings | Provide a feature view for proteins with little experimental evidence; test a pretrained encoder rather than training one from scratch | Compare against sequence similarity, domains, and existing features using homology-group holdouts |
| Multi-view latent factors or masked matrix completion | Combine expression, fitness, protein features, and other views while retaining missing-data masks and view contributions | Reconstruct held-out observed measurements, then test independent trait prediction; do not equate reconstruction with function discovery |
| Task-specific network integration | Learn which permitted relation layers support a particular question | Fit weights using training data only; remove target-derived edges and compare against the existing single-layer propagation |
| Multi-label and continuous targets | Proteins can have several roles; fitness, abundance, and drug response need not be forced into categories | Add multi-label evaluation and regression with uncertainty, alongside categorical tasks |
| Structured literature assertions | Turn a mention into a traceable claim about an assay, phenotype, or interaction | Store the exact source passage/table, entity mapping, condition, negation, and extraction status; require review before treating it as a measured label |
| Experiment selection | Suggest the next measurement that would best distinguish competing explanations | Compare predicted information gain with practical assay cost and confirm prospectively |

These are hypotheses to benchmark, not promised accuracy improvements. Pretrained
ESM models provide per-protein representations and extraction tools, making them
a practical candidate input view. [ESM repository](https://github.com/facebookresearch/esm).

MOFA+ is an example of multi-view factor modelling that allows missing values.
Its published applications are not a ready-made parasite gene-function benchmark:
Starplast would need aligned gene rows and appropriate feature views. Missingness
support also does not remove assay-selection bias. [MOFA+ paper](https://link.springer.com/article/10.1186/s13059-020-02015-1).

GeneMANIA illustrates learning network relevance and redundancy before propagating
functional labels. This extends Starplast's current single-layer propagation and
multiplex community consensus. Keeping edge provenance is compatible with an
explicitly learned combination for a particular task.
[GeneMANIA paper](https://link.springer.com/article/10.1186/gb-2008-9-s1-s4).

I would postpone a large graph neural network until the existing linear, boosted,
and propagation baselines have reliable evaluations. A more complex architecture
is useful only if it improves predictions under the same realistic tests.

## Evaluate the question you actually want to answer

Start with three explicit tasks: localization classification, a continuous
fitness phenotype, and a screen-derived phenotype relevant to spaCR. Each needs
its own eligible population, label definition, negative/unknown policy, and
scientifically meaningful baseline.

Use identical splits across methods. Keep related proteins together when measuring
transfer to new families; hold out entire source studies or conditions when
measuring transfer to new experiments. Add a time-based test when dated snapshots
allow it: train with information available at one date and evaluate later
experimental annotations. Time-delayed evaluation is used in CAFA protein-function
assessment. [CAFA evaluation](https://www.nature.com/articles/nmeth.2340).

Report per-trait precision and recall, precision among the top candidates, PR
curves where negatives are justified, and the fraction of genes for which the
method can make a supported call. Treat unknown annotations as unknown, rather
than automatic negatives. Report performance for sparse evidence, rare classes,
and novel families separately. Keep threshold choice and calibration away from
the final test set. Include label-permutation and source-ablation controls.

A useful candidate output is a row naming the gene, proposed trait and condition,
score type, validation support, measured evidence, contrary evidence, missing
measurements, and a plausible discriminating experiment. A cluster label and
purity alone are not enough for deciding what to test next.

## Streamline the product around three tasks

- **Explore a gene:** measurements, evidence timeline, neighbours by evidence type,
  and clear gaps in knowledge.
- **Predict a trait:** choose the biological question; inspect eligible data,
  validated method comparisons, and one ranked candidate table.
- **Compare a screen:** import spaCR results with assay provenance, compare them
  against existing evidence, and export a justified experimental shortlist.

Keep map controls and model tuning available under advanced settings. Reuse the
existing candidate tables, recipe exports and job system. Show the same candidate
selection in the table, evidence panel, and map. Replace several competing notions
of score/confidence with explicitly named quantities shared across methods.

Internally, separate source ingestion, evidence records, feature construction,
evaluation, predictions, and rendering. `app.py` and `analysis_panel.py` together
contain about 6,900 lines; extracting these responsibilities gradually is more
useful than a rewrite or additional visualization modes.

## Suggested order of work

1. Fix prediction scoring and no-support calls; standardize the saved analysis
   recipe and default-map pipeline. Acceptance: swapped labels are penalized,
   disconnected genes abstain, and every result identifies its actual data and method.
2. Establish the three benchmark tasks with group-aware splits and honest
   candidate reporting. Acceptance: methods can be compared on the same held-out genes.
3. Add observation-level provenance, coverage/QC states, and versioned snapshots.
   Acceptance: any candidate feature can be traced to its source measurement.
4. Benchmark sequence embeddings and one multi-view model against existing methods.
   Acceptance: keep additions only when they improve a predefined useful metric.
5. Consolidate the user workflow and spaCR handoff around the validated shortlist.

The strongest product would answer: **what might this gene do, which independent
evidence supports that idea, how reliable is the prediction for a gene like this,
and which experiment would most usefully test it?** The gene map is a valuable
way to explore those answers.
