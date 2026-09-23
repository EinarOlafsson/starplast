# Python API

Install Starplast to use the `starplast` Python package. The distribution carrying
the code is named `starplast-core`; `pip install starplast` installs it for you.
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

The current UMAP implementation falls back to PCA if both accelerated and CPU
UMAP fail. Read the log when checking the method that actually ran. Save the
recipe, input data, package version, backend, and output coordinates with a result;
a seed alone does not ensure identical coordinates across software or hardware.

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

`app.Window(species=...)` requires an existing `PyQt6.QtWidgets.QApplication`.
Call it on the GUI thread. `app.main()` owns the application event loop and is
intended as the console entry point. Scripts that only need tables and analyses
should import the relevant modules directly.
