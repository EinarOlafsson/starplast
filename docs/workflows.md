# Explore, predict and compare

The **Tools** menu has three guided entry points. Advanced map and model controls
remain in the analysis dock. Each workflow uses the organism currently open in
the main window.

## Explore a gene

Search an accession, symbol or product description. Double-click a match to select
it on the map and inspect its measurements, missing values, evidence status and
source links. A gene-table value is often an aggregate: it does not recover the
original replicate, time point or uncertainty when those details were not kept.
Legacy summaries without an explicit evidence classification say **unclassified**.

## Predict a trait

Choose an observed outcome, its type, a method and the grouping used for evaluation.
**Orthogroup** keeps related proteins together; random gene folds answer an easier
question and may overestimate transfer to new families. A study or condition column
can be used when it actually identifies the intended held-out groups. Starplast
does not invent source dates or study membership for a time/study holdout.

The linear and prior models are useful starting points. Boosted trees, feature
neighbours, PCA neighbours, UMAP neighbours and masked factors provide comparisons.
The UMAP predictor fits a separate representation inside each training fold; it
does not classify the opening three-dimensional map. Masked factors are a simple
iterative PCA completion baseline, not an implementation of MOFA+.

For Toxoplasma, **Include frozen protein sequence features** loads the bundled
[ESM-2 representations](protein-sequences.md). It does not download a model.
[Local AF3 features](structures.md) are already present where an exact mapping was
available. Sequence or structure evidence is not automatically an improvement;
evaluate the additions on the same held-out genes.

**Evaluate and predict** runs in the shared Jobs panel. Feature exclusions,
imputation, scaling and dimensionality reduction are fitted on training genes.
The target, its registered relatives and derived inputs are excluded. Labels in
the calibration and evaluation groups are not used for fitting those transforms.
The Stop button requests cancellation between model fits.

Read the held-out metrics before the candidate table. An unlabelled gene is not a
negative example. A gene with insufficient measured features gets no supported
call. Classification can also abstain below the selected model-probability cutoff.
The cutoff is a declared setting, not a threshold optimized on the test labels.

Classification reports fixed-label accuracy, macro F1, per-class precision/recall
and average precision, Brier score and log loss. Unsupported calls count as misses
in accuracy and recall; selective accuracy is separately named. Probabilities use
temperature scaling on reserved groups when enough calibration labels exist. The
table explicitly says when calibration was unavailable. Even a calibrated model
can become unreliable on a less-studied population.

Regression reports RMSE, MAE, R², rank correlation and, when available, intervals
estimated from reserved-group residuals. Their observed coverage is reported;
family or domain shift prevents an unconditional coverage guarantee. Multi-label
tasks are available through the API as separate measured binary columns with
unknown entries preserved.

**Export full run** writes `predictions.csv`, `per_class.csv` and `run.json`.
The JSON records settings, input identity, exact training/calibration/test gene
IDs, excluded features, actual estimators and package versions. Results remain
model hypotheses. Comparing methods here is development validation; evaluating a
selected winner prospectively requires new independent data.

## Compare a screen

Open a gene-level CSV/TSV from spaCR or another measurement pipeline. Choose the
gene identifier and measured-value columns. Identifiers are resolved through the
current table and its identity catalogue. Ambiguous/unresolved identifiers remain
visible; missing or nonnumeric values stay missing, and a measured zero stays zero.

Repeated gene rows, including aliases that converge on one gene, require explicit
aggregation before import. Starplast does not silently average guide-level results.
The exported comparison includes the existing gene evidence and a JSON record of
the source filename, SHA-256 hash and selected columns. Use the advanced importer
when you want to add measurements to an analysis table.

## Evidence and literature records

`starplast.evidence.Observation` stores the biological entity, trait, typed value,
unit, source/version/location, context, replicate, uncertainty, evidence status,
missingness state and derivation. Gene, transcript, protein, residue, allele, pair,
host-gene and metabolite identities stay distinct. Contradictory observations from
different studies or contexts remain separate content-addressed records.

`from_gene_table()` wraps existing summaries with their available provenance; it
does not claim to reconstruct raw experiments. Write source-level observations
directly when loading original tables. Missingness distinguishes unassayed,
below-detection, failed-QC, ambiguous and unknown values.

`propose_assertions()` offers a conservative phrase extractor for source passages.
Its output is always pending review, retains the exact passage and negation, and
omits sentences with ambiguous gene mappings. A named reviewer must accept a claim
before `LiteratureAssertion.observation()` can expose it as a reviewed literature
claim. Acceptance confirms the interpretation of a passage, not experimental truth.
Store and review these records through the API; the desktop does not yet contain
a bulk literature-curation interface.

## Networks and explanations

The `network_prediction.run_network()` API evaluates weighted propagation on a
fixed graph with explicit gene order and edge provenance. Seed labels, layer-weight
validation, probability calibration and outer evaluation use disjoint groups.
The guided method **Weighted networks (fixed graph)** uses the bundled domain,
coexpression and structural-similarity layers with those builders' declared inputs.
Target-derived layers are excluded; disconnected genes abstain. This is explicitly
**transductive** because the fixed graph includes evaluation nodes. It cannot be
reported as performance on newly arriving nodes. Custom network fitting uses the
API so its source declarations remain explicit. Sequence-feature and measured-feature
coverage controls do not apply to network propagation.

`prioritization.explain_candidates()` reports observed source coverage, training
neighbours and, for linear models, feature contributions. It explains unknown-gene
calls using the final fitted model; it does not leak that model into held-out
evaluation explanations. Model contributions describe association, not mechanism.

`prioritization.prioritize()` is a generic budget heuristic for separately
calibrated classification results, using explicit cost, benefit and an optional
uncertainty bonus. It records default assumptions and can limit selection to one
item per supplied group. It is not an optimal budget solver, an information-gain
model or a prospectively validated experimental design system.
