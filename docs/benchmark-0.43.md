# Prediction benchmark for 0.43

Starplast now evaluates predictions separately from the displayed map. On the
bundled Toxoplasma data, supervised feature models recover localization better
than neighbours in a UMAP. Frozen protein sequence features improve localization
in this development benchmark. The added AF3 summary features do not yet show a
clear improvement on their own.

These are held-out predictions of existing measurements, not newly validated
gene annotations. The results support using the map for exploration and the
**Tools → Predict a trait** workflow for quantitative comparisons.

## Evaluation

The 23 September 2026 run contains 25 feature-model comparisons and four network
comparisons. Localization has 3,827 measured genes in 26 classes; protein abundance
has 748 measured `protein_ibaq_log2` values. Unknown outcomes are not negative
examples. All comparisons use three outer folds grouped by orthogroup, seed 42.
Unassigned genes receive individual groups. This prevents a recorded orthogroup
from crossing folds; it is not a sequence-identity clustering guarantee.

Within each outer training fold, separate groups are reserved for calibration.
Target-family exclusions, association-based exclusions, feature selection,
median imputation, scaling and dimensionality reduction use fitting genes only.
Inputs with less than 5% observed coverage in those genes are removed. Genes need
10% measured coverage for a supported feature-model call. No probability cutoff
was tuned on the evaluation labels; the benchmark cutoff is zero. Each JSON run
contract records the actual columns and fitting/calibration/evaluation gene IDs.

Published-input comparisons exclude AF3 and ESM columns. The next view adds nine
AF3 confidence, coverage and geometry summaries, available for 1,210 genes. The
last adds 320 frozen ESM-2 features for 8,064 proteins. Literature-count predictors
are excluded. ESM is a pretrained model: grouped outcome validation does not imply
that its pretraining corpus lacked these proteins.

The prior, linear model, boosted trees, feature neighbours, PCA neighbours,
UMAP neighbours and masked-factor model all receive the same outer folds for a
given outcome. PCA, UMAP and masked factors use 20 dimensions. The additional
three-dimensional UMAP comparator is trained inside each fold; it does not use
the opening map. Masked factors are iterative PCA completion, not MOFA+.

## Localization

| Model | Inputs | Accuracy | Macro F1 | Macro average precision |
|---|---|---:|---:|---:|
| Class prior | Published | 0.201 | 0.013 | 0.039 |
| Linear | Published | 0.460 | 0.382 | 0.396 |
| Boosted trees | Published | 0.517 | 0.390 | 0.449 |
| Feature neighbours | Published | 0.390 | 0.283 | 0.223 |
| PCA neighbours | Published | 0.359 | 0.254 | 0.194 |
| UMAP neighbours, 20D | Published | 0.314 | 0.228 | 0.169 |
| UMAP neighbours, 3D | Published | 0.305 | 0.215 | 0.166 |
| Masked factors | Published | 0.364 | 0.253 | 0.180 |
| Boosted trees | Published + AF3 | 0.520 | 0.401 | 0.452 |
| Linear | Published + AF3 + ESM | 0.515 | 0.419 | 0.437 |
| Boosted trees | Published + AF3 + ESM | **0.562** | **0.434** | **0.496** |
| Linear, shuffled outcomes | Published | 0.137 | 0.036 | 0.040 |

All feature models supported all 3,827 evaluation genes at the declared coverage
threshold. Accuracy alone hides rare-class failures: even the best model makes
no correct top-class calls for ER 2, PM–peripheral 2 or tubulin cytoskeleton.
The full per-class table retains these zeros. This is not a reliable universal
localization annotator.

Relative to boosted trees on published inputs, adding AF3 changes macro F1 by
**+0.011**, with a paired orthogroup-bootstrap 95% interval of **−0.001 to +0.025**.
Adding AF3 and ESM changes it by **+0.044**, interval **+0.021 to +0.068**.
The latter is evidence for the combined addition, not an estimate of the effect
of ESM alone. These 500-resample intervals condition on the already fitted models;
they omit training and model-selection uncertainty.

Published-feature coverage remains relevant. Boosted-tree accuracy across its
lower/middle/upper coverage thirds is 0.479/0.519/0.554; with AF3 and ESM it is
0.521/0.578/0.588. The thirds use the same published-input coverage for every
model, so dense sequence vectors cannot move a sparse gene into a richer group.

Held-out temperature scaling improves the interpretation of scores but does not
make every probability reliable. For the combined boosted model, the 0.8–0.9
confidence bin has mean confidence 0.850 and accuracy 0.783 (465 genes). The
0.9–1.0 bin has mean confidence 0.937 and accuracy 0.911 (350 genes). Brier score
is 0.588 and log loss 1.433. Reliability may deteriorate further on unstudied genes.

## Protein abundance

| Model | Inputs | RMSE, log2 units ↓ | R² ↑ |
|---|---|---:|---:|
| Mean prior | Published | 3.388 | −0.009 |
| Linear | Published | 3.023 | 0.197 |
| Boosted trees | Published | 2.556 | 0.426 |
| Feature neighbours | Published | 2.598 | 0.407 |
| PCA neighbours | Published | 2.600 | 0.406 |
| UMAP neighbours | Published | 2.650 | 0.383 |
| Masked factors | Published | 2.664 | 0.376 |
| Boosted trees | Published + AF3 | 2.556 | 0.426 |
| Boosted trees | Published + AF3 + ESM | **2.508** | **0.447** |
| Linear, shuffled outcomes | Published | 4.132 | −0.501 |

AF3 gives identical boosted-model predictions here. The combined addition reduces
RMSE by 0.048, but the paired interval for candidate minus baseline is
**−0.125 to +0.022**, which includes no improvement. The best model's nominal
90% residual intervals cover 94.9% of held-out observations, with a mean width of
9.94 log2 units. Those broad intervals limit practical precision; they are not a
guarantee of coverage under a different assay or population.

## Fixed networks

| Layer | Supported fraction | Accuracy, all genes | Macro F1 |
|---|---:|---:|---:|
| Shared domains | 0.159 | 0.058 | 0.111 |
| Coexpression | 0.804 | 0.115 | 0.120 |
| Structural similarity | 0.288 | 0.098 | 0.167 |
| Learned combination | 0.873 | 0.202 | 0.195 |

These are transductive predictions: the fixed graph contains evaluation nodes.
Seed labels, layer-weight fitting, calibration and outer evaluation use separate
groups. Unsupported calls count as misses; accuracy among the combination's
3,341 supported genes is 0.231. Combining layers improves coverage but performs
well below the feature models. It should remain an optional comparison.
Structural edges come from the existing Foldseek layer, not a new all-versus-all
comparison of the locally indexed AF3 models.

## Reproduce and inspect

From a checkout with the application dependencies installed:

```bash
export OPENBLAS_NUM_THREADS=2 OMP_NUM_THREADS=2 NUMBA_NUM_THREADS=2
python scripts/benchmark_inference.py --output /tmp/starplast-feature-runs \
  --sequence starplast/data/esm_features.parquet --controls
python scripts/benchmark_network.py --output /tmp/starplast-network-runs
python scripts/summarize_inference.py --input /tmp/starplast-feature-runs \
  --networks /tmp/starplast-network-runs --output /tmp/starplast-summary
```

The [committed results](https://github.com/EinarOlafsson/starplast/tree/main/results/inference_043)
contain aggregate metrics, per-class results, reliability bins, coverage strata,
paired intervals and compressed full run contracts. The manifest hashes the input
tables, graph, prediction implementations and artifacts. The recorded environment
used NumPy 2.5.3, pandas 3.0.6, scikit-learn 1.9.1, SciPy 1.18.1 and UMAP 0.5.12.
Numeric results can vary with library versions. Unknown-gene candidate tables are
not part of this benchmark release.

There is no prospective test, independent temporal holdout or complete
source-study holdout here. Assay selection, missingness and unrecorded derivations
can still bias evaluation. The models were compared on the same development
folds, so the selected winner needs independent testing. Multi-label prediction,
literature review records, explanations and budget heuristics have software tests;
this report does not establish their biological or decision-making utility.

The practical recommendation is to keep UMAP as a navigable overview, start
prediction with a prior and a linear model, compare boosted trees on the same
groups, and inspect class-level failures and calibration before using a result.
Retain the AF3 data as evidence while distinguishing its availability from proven
predictive value. The shared embedding builder records what actually ran; the
default map remains exploratory because numeric localization confidence can be
among its inputs.
