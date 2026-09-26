# Python API

Install the application, analysis modules, and bundled data with `pip install starplast`.
The analysis modules can be used without creating a Qt application.

The [generated module reference](api/starplast.html) includes signatures,
docstrings, and source links. This page covers common entry points. Starplast is
pre-1.0, so pin the package version when preserving a reproducible analysis.

## Read the bundled gene table

```python
import pandas as pd
from starplast import paths

nodes = pd.read_parquet(paths.cache_file("nodes.parquet"))
print(nodes[["gene_id", "product"]].head())
```

Use `pf_nodes.parquet` for *P. falciparum*. `paths.cache_file()` honours
`STARPLAST_CACHE`. This example imports no GUI modules and does not download data.
Keep the table row order when aligning coordinates or cluster labels with genes.

## Build a feature matrix and embedding

```python
from starplast.embedding import EmbeddingSpec, build_matrix, embed

spec = EmbeddingSpec(
    name="fitness-and-expression",
    blocks=("fitness_screens", "expression_summary"),
    na_policy="median",
    scaling="robust",
    method="pca",
    n_components=3,
    random_state=42,
)
X, feature_names, kept_rows = build_matrix(nodes, spec)
coords, feature_names, kept_rows = embed(nodes, spec)
mapped_genes = nodes.loc[kept_rows, "gene_id"].to_numpy()
assert len(mapped_genes) == len(coords)
recipe = spec.to_dict()
```

`kept_rows` is a boolean mask into the input table, including when a missing-value
policy drops genes. `coords` has one row per retained gene. The number of PCA
components is limited by the available matrix dimensions. Use `method="umap"` or
`method="tsne"` for nonlinear projections. `embed(..., return_matrix=True)` appends
the exact feature matrix as a fourth return value.

Use `embed(..., strict=True)` to refuse a failed backend or algorithm. In ordinary
display mode a failed UMAP call can return PCA; the log and saved execution record
identify that fallback. `return_metadata=True` appends a dictionary containing the
requested and executed methods, backend, recipe, ordered gene IDs, matrix and
coordinate hashes, and library versions. Coordinates also carry this record into
`tuning.EmbeddingStore.save()`, which checks gene order and coordinate integrity.
A seed alone does not ensure identical coordinates across software or hardware.

`embedding.default_spec(nodes)` selects the balanced display recipe used by the
packaged map builder. The shipped archives record that execution and the exact
node-file hash. Rebuild just their layouts with `python scripts/rebuild_layouts.py`;
this preserves the existing edge arrays. The display is exploratory: it is not a
held-out prediction of any trait used to construct it.

## Cluster and examine held-out evidence

```python
from starplast.clustering import cluster, battery

labels = cluster(coords, algorithm="hdbscan", min_cluster_size=25)
subset = nodes.loc[kept_rows].reset_index(drop=True)
summary, detail = battery(subset, labels, used_features=feature_names)
```

Noise has label `-1`. `summary` reports feature-level associations; `detail`
contains category or cluster-level results. Supplying `used_features` lets the
battery distinguish input features from held-out evidence. Statistical association
with a cluster does not establish a mechanism.

For label recovery, use `search.excluded_detail()` to inspect excluded inputs and
`search.search()` to run a search with target exclusions. Use
`validate.masked_recovery()` or `validate.validate_all()` to evaluate categories
hidden from cluster selection. Pass `used_columns` and `target_column` so the
circularity checks can run. `holdout_cv.nested_structure_cv()` selects among
precomputed clusterings using inner validation labels and scores untouched outer
folds. A search winner still needs independent confirmation.

## Import a screen table

```python
from starplast.importer import read_any, preprocess, merge_into

screen = read_any("screen.csv")
imported, provenance = preprocess(
    screen,
    gene_column="gene_id",
    columns=["score"],
    quantification="none",
    scaling="none",
    na_policy="median",
    duplicates="mean",
    prefix="spacr_",
)
nodes_with_screen = merge_into(nodes, imported)
```

The current importer recognizes Toxoplasma accessions (`TGME49_`, `TGGT1_`,
and `TGVEG_`). Plasmodium table import is not yet supported. The input must
contain the chosen identifier and measurement columns. Check the
returned provenance and matched genes before interpreting results. Choose a
missing-value policy deliberately; `median` estimates absent values. To include
an imported column in an embedding, add it to `EmbeddingSpec.extra_columns`.

## Run discovery without a window

```bash
starplast-discover --task guilt:compartment_best --budget 10 --name first-pass
starplast-discover --list
starplast-discover --read first-pass
```

Use `--nodes /path/to/pf_nodes.parquet` for a different node table and choose a
layer present in that table. `--out` sets the results directory. Searches can be
long-running; `--help` lists budget, seed, and feature-selection options.

## Evaluate a trait with held-out families

```python
from starplast.prediction import TaskSpec, run

result = run(nodes, TaskSpec("compartment", method="linear", group_column="orthogroup"))
print(result.metrics)
print(result.per_class)
result.save("results/localization")
```

Use `kind="regression"` for continuous measured outcomes, or `run_multilabel()`
for separate observed binary columns. `features` and `exclude` declare inputs;
registered target-derived measurements are still excluded. Calibration and
preprocessing are fitted within training groups. Read the [workflow guide](workflows.md)
for uncertainty, abstention and evaluation limits.

```python
from starplast.evidence import Observation, write_observations

records = [Observation(
    entity_id="example_gene", trait="example_measurement", value=0.0,
    source_id="study_accession", organism="example_organism",
    source_version="v1", source_location="Table 2, row 3",
    unit="relative abundance", context={"condition": "reference"},
    replicate="1", evidence_status="measured",
)]
write_observations(records, "results/observations.parquet")
```

Each record preserves its source and missingness before aggregation. Numerical
zero, measured False and unassayed None have distinct meanings.

## Run, test and trust a strategy

`starplast.strategies` is the Strategies tab without the window: 39 named ways of turning the
tables into a claim, each with a self-test and a scorecard. Nothing here imports Qt.

```python
from starplast import strategies as S

table = S.overview("Tg")            # one row per strategy: name (with its method), task, grade, skill
print(table[["number", "name", "task", "grade", "skill_default"]].head(10))

s = S.get("feature_knn")            # one strategy: its question, explanation and settings
print(s.name, "|", s.method, "|", s.task)
print(s.parameters())               # every setting, its default and why it exists
print(s.techniques_table())         # what the method is built from, each technique explained
print(s.scorecard_table())          # the metrics its self-test reports, each explained
```

Run it, or test it. A run answers the strategy's question on the whole table. A test hides known
information, asks for it back, and judges the answer against the same procedure on shuffled data:

```python
result = S.run("feature_knn", "Tg", target="compartment")
print(result.summary)
calls = result.tables["calls"]      # every output is a DataFrame; result.save(folder) writes them

test = S.test("feature_knn", "Tg", target="compartment")
print(test.summary())               # the verdict and the number it rests on
print(test.card())                  # verdict block, then the task's metrics in their standard order
print(test.scorecard["macro_f1"], test.scorecard["mcc"], test.skill)
```

Every self-test returns the same verdict block (observed, chance, bar, p, skill, hidden) and then
its task's metrics, always in the same order: **label calls** (accuracy, coverage, precision of
calls, macro precision, macro recall, macro F1, weighted F1, Cohen's kappa, MCC, macro AUROC,
macro AUPRC), **ranking** (AUROC, AUPRC, AUPRC lift, prevalence, partial AUROC, R-precision,
precision and enrichment at the top 1%, recall at the top 10%, best F1, nDCG), **set retrieval**,
**cluster recovery** (weighted F1, ARI, NMI, homogeneity, completeness ...), **values**
(Spearman, Pearson, Kendall, R-squared, normalised RMSE, MAE, decile recall) and
**replication**. `S.metrics()` defines each one -- its range, its chance level and how to read
it -- and `S.techniques()` explains every technique. [The scorecards page](scorecards.md) has
the same glossaries and every strategy's measured card.

How far to trust a strategy before running it comes from its calibration: its self-test run over
a grid of settings, several held-out labels and five seeds.

```python
cal = S.calibration("feature_knn", "Tg")
print(cal["grade"], cal["default"]["skill"], cal["default"]["skill_low"], cal["default"]["skill_high"])
print(cal["default"].get("scorecard", {}).get("macro_f1"))  # {mean, low, high, runs}
best = S.tuned("feature_knn", "Tg")                          # the setting calibration chose
result = S.run("feature_knn", "Tg", target="compartment", **best)
```

To run several strategies on one table, or on your own table, keep a context. It caches maps and
matrices, so a second strategy walking the same grid pays nothing:

```python
ctx = S.Context.shipped("Tg").bound(log=print)   # log= shows progress; should_stop= can cancel
for key in ("feature_knn", "supervised_classifier", "stacking"):
    t = S.get(key).test(ctx, target="compartment")
    print(key, t.verdict, round(t.scorecard["macro_f1"], 3))

mine = S.Context(my_nodes, organism="Tg")        # any gene x column table with a gene_id column
```

The scorecard functions work on predictions from anywhere, so a method outside Starplast can be
scored exactly as the strategies are:

```python
from starplast import scorecard as SC

SC.label_calls(predicted_labels, true_labels, positions)   # NaN = no call, counted as wrong
SC.ranking(scores, is_positive)                             # higher score = more likely positive
SC.values(predicted, measured)
```

## Module map

| Task | Modules |
|---|---|
| Resolve data and user state | [paths](api/starplast/paths.html) |
| Inspect datasets and biological evidence groups | [datasets](api/starplast/datasets.html), [slots](api/starplast/slots.html) |
| Import gene-level measurements | [importer](api/starplast/importer.html) |
| Select features and embed | [embedding](api/starplast/embedding.html) |
| Cluster and score associations | [clustering](api/starplast/clustering.html) |
| Search and assess held-out recovery | [search](api/starplast/search.html), [validate](api/starplast/validate.html), [holdout_cv](api/starplast/holdout_cv.html) |
| Rank candidate findings | [discovery](api/starplast/discovery.html), [optimize](api/starplast/optimize.html) |
| Store analyses and recipes | [runs](api/starplast/runs.html), [searches](api/starplast/searches.html), [recipes](api/starplast/recipes.html) |
| Embed the desktop interface | [app](api/starplast/app.html), [analysis_panel](api/starplast/analysis_panel.html) |
| Guided workflows | [workflows](api/starplast/workflows.html) |
| Group-aware classification, regression and multi-label evaluation | [prediction](api/starplast/prediction.html) |
| Weighted transductive networks | [network_prediction](api/starplast/network_prediction.html) |
| Typed observations and reviewed assertions | [evidence](api/starplast/evidence.html) |
| Explain results and compare measurements | [prioritization](api/starplast/prioritization.html) |
| Named inference strategies, their self-tests and calibration | [strategies](api/starplast/strategies.html), [calibration](api/starplast/calibration.html) |
| Standard metrics for every task, and their glossary | [scorecard](api/starplast/scorecard.html) |
| What each strategy's method is built from | [techniques](api/starplast/techniques.html) |
| One neighbour space over every layer; learned edge strengths | [graphspace](api/starplast/graphspace.html) |
| Deposited datasets, derived into columns | [deposits](api/starplast/deposits.html) |

`app.Window(species=...)` requires an existing `PyQt6.QtWidgets.QApplication`.
Call it on the GUI thread. `app.main()` owns the application event loop and is
intended as the console entry point. Scripts that only need tables and analyses
should import the relevant modules directly.
