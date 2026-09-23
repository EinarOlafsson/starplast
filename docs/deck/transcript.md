# Starplast: slide transcript

A guide to the 0.43 release.

## 01. Starplast

Explore what is known.
Evaluate what a model predicts.

A practical guide to the 0.43 release

[Supporting documentation](https://einarolafsson.github.io/starplast/index.html)

## 02. Start with a gene. Keep the evidence.

Bring scattered measurements together without losing what they mean.

Collect: Published measurements, annotations and protein features.

Connect: Match gene identities and retain the source of each value.

Explore: Read evidence, compare relationships and inspect gaps.

Evaluate: Test predictions on held-out genes before interpreting them.

A model suggests a hypothesis. Its evidence and validation determine how useful it is.

[Supporting documentation](https://einarolafsson.github.io/starplast/guide.html)

## 03. Two organisms. Separate maps.

Offline gene tables and networks are included in the installation.

8,140: Toxoplasma gondii genes

5,720: Plasmodium falciparum genes

130: Registered datasets and computed layers

Coverage differs by assay and organism. These counts do not imply complete evidence for every gene.

[Supporting documentation](https://einarolafsson.github.io/starplast/datasets.html)

## 04. Choose the task you came to do.

Explore a gene: Find an accession, symbol or product. Read measurements and their sources.

Predict a trait: Choose an observed outcome. Evaluate a model, then inspect unlabelled-gene hypotheses.

Compare a screen: Join a gene-level results table to the existing evidence. Keep unresolved IDs visible.

Open these guided workflows from the Tools menu.

[Supporting documentation](https://einarolafsson.github.io/starplast/workflows.html)

## 05. The map is one view of the evidence.

One point represents one gene.

Position: Similarity under the chosen feature recipe.

Colour: An annotation or measurement you select.

Lines: One explicitly selected relationship type.

Proximity alone does not establish function or a physical interaction.

[Supporting documentation](https://einarolafsson.github.io/starplast/guide.html)

## 06. Install and open Starplast.

pip install starplast
starplast

Python 3.10+: Use a dedicated environment when possible.

Desktop display: The viewer needs a display and working OpenGL.

CPU included: CUDA is optional; browsing does not require a GPU.

Linux, macOS and Windows. Built tables are available offline.

[Supporting documentation](https://einarolafsson.github.io/starplast/releases.html)

## 07. Your first five minutes.

1  Choose: File > Species: select the organism.

2  Search: Tools > Explore a gene: enter an accession or product.

3  Inspect: Read the evidence and open its source links.

4  Compare: Change a colour or edge layer; keep the question specific.

Start with one familiar gene so you can check that the interpretation makes sense.

[Supporting documentation](https://einarolafsson.github.io/starplast/guide.html)

## 08. Find your way around the workspace.

Map controls on the left. The gene map in the centre. Evidence on the right.

Tools opens guided workflows and analysis docks. View controls appearance and interaction.

[Supporting documentation](https://einarolafsson.github.io/starplast/guide.html)

## 09. Move the map. Keep your bearings.

Rotate and zoom: Left-drag in Navigate mode to rotate. Scroll to move closer or farther away.

Select a gene: Click a point to show its evidence. Search can take you directly to a gene.

Focus a category: Filter checked categories. Double-click a category to centre its genes.

Right-click the map for display options, export and reset.

[Supporting documentation](https://einarolafsson.github.io/starplast/guide.html)

## 10. Search by identity or description.

Tools > Explore a gene

Enter a query: Use a gene accession, symbol or product text.

Check the match: Double-click a result to select the same gene on the map.

Read the table: Inspect available values, missingness and source evidence.

Shown here: an existing alpha-tubulin record, not a model-generated annotation.

[Supporting documentation](https://einarolafsson.github.io/starplast/workflows.html)

## 11. Read the source, not just the value.

What was measured?: Distinguish expression, localization, abundance, sequence and interaction evidence.

In which context?: Check the organism, stage, condition and assay before comparing measurements.

How was it obtained?: Separate measured evidence, computed features, transferred annotations and predictions.

Legacy gene summaries cannot reconstruct replicate-level uncertainty that was never retained.

[Supporting documentation](https://einarolafsson.github.io/starplast/workflows.html)

## 12. The full data catalogue stays one click away.

Evidence family | Examples | Read it as

Expression | Stage profiles and summaries | Condition-dependent measurements

Localization | Compartment annotations | Evidence with assay-specific coverage

Protein features | Sequence, domains, AF3 summaries | Measured or computed descriptors

Relationships | Coexpression and physical assays | Distinct types of connections

Literature | Mentions and reviewed assertions | Traceable source context

Each catalogue record links to its publication, source data or computed inputs.

[Supporting documentation](https://einarolafsson.github.io/starplast/datasets.html)

## 13. Unknown is different from zero.

0: A measured value that really is zero

Unknown: No usable value is available

No call: The model lacks enough support

A missing annotation is not a negative example. Grey marks unknown evidence.

[Supporting documentation](https://einarolafsson.github.io/starplast/workflows.html)

## 14. Different lines answer different questions.

Similarity: Shared domains, orthogroups, expression or structure.

Physical evidence: Crosslink MS and pulldown assays have their own interpretation.

Literature context: Co-mention means appearing together in text; it is not proof of binding.

The Edges controls keep these layers separate.

[Supporting documentation](https://einarolafsson.github.io/starplast/guide.html)

## 15. Inspect one connection layer at a time.

Select: Click a gene you want to inspect.

Choose: Turn on a relevant layer in Edges.

Read: Check neighbouring genes and the evidence behind the relation.

Compare: Switch the layer and see which relationships persist.

Draw all active edges is available, but dense networks are often easier to understand locally.

[Supporting documentation](https://einarolafsson.github.io/starplast/guide.html)

## 16. Select a group without changing the data.

Lasso: screen space: Draw an outline around visible positions. Points behind one another can be included.

Brush: 3D space: Select a spatial neighbourhood around the gene you press.

Choose View > Left mouse button > Select. Check the selected gene list before exporting.

[Supporting documentation](https://einarolafsson.github.io/starplast/guide.html)

## 17. Explore clusters and compare analyses.

Clusters: Group similar coordinates with clustering methods. Noise label -1 means unassigned, not a biological class.

Search and validation: Compare recipes with explicit targets and budgets. Confirm a search winner with independent evidence.

Saved results: Name embeddings and retain their recipes. Right-click result tables to export them; save selected gene lists.

Cluster labels describe the chosen inputs. They do not establish a shared biological function.

[Supporting documentation](https://einarolafsson.github.io/starplast/guide.html)

## 18. Change a map deliberately.

Features: Choose relevant measurement blocks.

Missingness: Decide how absent inputs are handled.

Scaling: Keep feature magnitudes and blocks comparable.

Embedding: Choose a method, dimensions and a seed.

Smaller UMAP neighbourhoods favour local structure. Lower minimum distance packs points more tightly.

[Supporting documentation](https://einarolafsson.github.io/starplast/guide.html)

## 19. Explore in 3D. Evaluate separately.

The displayed map: An exploratory view of selected inputs. Numeric localization confidence can be among them.

A held-out prediction: A separate fitted model with target exclusions, training-only transforms and untouched evaluation groups.

The guided predictor does not classify the opening map.

[Supporting documentation](https://einarolafsson.github.io/starplast/workflows.html)

## 20. Set up a prediction in one place.

Tools > Predict a trait

Outcome and type: Choose the observed target and classification or regression.

Method and groups: Start with a baseline. Keep related genes together with orthogroup.

Run and inspect: Evaluate and predict. Read held-out metrics before the candidate table.

The Jobs panel runs the work. Stop requests cancellation between model fits.

[Supporting documentation](https://einarolafsson.github.io/starplast/workflows.html)

## 21. Keep evaluation genes out of fitting.

Fit: Learn features and model parameters.

Calibrate: Adjust probabilities or residual intervals on separate groups.

Evaluate: Score predictions on groups withheld from both steps.

Outer folds rotate the evaluation groups. Grouped validation is harder and more informative than mixing relatives.

[Supporting documentation](https://einarolafsson.github.io/starplast/benchmark-0.43.html)

## 22. Compare simple methods first.

Method | Useful starting question | Keep in mind

Prior / linear | Is there signal beyond a simple baseline? | Always retain a baseline

Boosted trees | Do nonlinear feature combinations help? | Use the same evaluation groups

Feature / PCA neighbours | Do similar measured genes agree? | Distance depends on preprocessing

UMAP / masked factors | Does a reduced representation help? | A compact map can lose signal

Weighted networks | Which fixed graph layers support labels? | Transductive; disconnected genes abstain

No method is automatically best because it is more complex.

[Supporting documentation](https://einarolafsson.github.io/starplast/workflows.html)

## 23. Read performance and coverage together.

Accuracy: The fraction of evaluation genes called correctly. Unsupported calls count as misses.

Macro F1: Balances precision and recall across classes, so common labels do not dominate the score.

Coverage: The fraction with a supported call. Selective accuracy describes only that supported subset.

Also inspect per-class precision, recall, average precision and calibration.

[Supporting documentation](https://einarolafsson.github.io/starplast/workflows.html)

## 24. What worked in the localization benchmark?

3,827 measured genes · 26 classes · three orthogroup-held-out folds

0.434 macro F1: Boosted trees with published, AF3 and ESM inputs.

56.2% accuracy: Better than the compared UMAP-neighbour models.

Important limit: Several rare classes still have zero top-class recall.

Development cross-validation, not prospective validation of new annotations.

[Supporting documentation](https://einarolafsson.github.io/starplast/benchmark-0.43.html)

## 25. A confident score can still be wrong.

Combined boosted model: observed accuracy in held-out confidence bins

Check calibration: Temperature scaling uses separate groups when enough labels exist.

Read the status: The app explicitly marks runs that could not be calibrated.

Use abstention: A probability cutoff can withhold uncertain calls; it is not a guarantee.

The 0.8-0.9 bin averaged 0.850 confidence but 0.783 accuracy across 465 genes.

[Supporting documentation](https://einarolafsson.github.io/starplast/benchmark-0.43.html)

## 26. Add sequence evidence with one checkbox.

8,064 proteins: Bundled frozen ESM-2 representations.

320 features: A numerical view of protein sequence, separate from experimental measurements.

No download to use: Enable Include frozen protein sequence features in the guided predictor.

Regenerating representations needs the optional sequence extra. Benchmark their value rather than assuming it.

[Supporting documentation](https://einarolafsson.github.io/starplast/protein-sequences.html)

## 27. Your AF3 collection becomes traceable evidence.

1,210 matched genes: Exact sequence or uniquely placed fragment matches in the available local collection.

Nine summary features: Confidence, sequence coverage and descriptive geometry.

Coordinates stay local: The package includes features and provenance, not the private model files.

This is the verified subset in the synced folders, not the total number of proteins ever folded.

[Supporting documentation](https://einarolafsson.github.io/starplast/structures.html)

## 28. Connect local folds to the structure viewer.

pip install 'starplast[structures]'

starplast-index-structures /path/to/models

Index once: Build the local structure catalogue from your model folders.

Select a gene: The structure viewer checks the local index before AlphaFold DB.

Read confidence: A predicted fold and its confidence are not experimental validation.

Recreate the index if the drive moves. Full models take precedence over fragments.

[Supporting documentation](https://einarolafsson.github.io/starplast/structures.html)

## 29. Continuous traits stay continuous.

Protein abundance: 748 measured genes · grouped held-out evaluation

Choose regression: Use Continuous measurement instead of converting an outcome into classes.

Read error units: RMSE here is measured in log2 abundance units; lower is better.

Small added gain: The combined-feature improvement interval includes no improvement.

Nominal 90% intervals covered 94.9%, but their mean width was 9.94 log2 units.

[Supporting documentation](https://einarolafsson.github.io/starplast/benchmark-0.43.html)

## 30. Network coverage is not prediction quality.

87.3% supported: Combining domain, expression and structural layers reaches more evaluation genes.

0.195 macro F1: The combined fixed network remains well below the compared feature models.

The graph includes evaluation nodes: report these results as transductive, not transfer to unseen nodes.

[Supporting documentation](https://einarolafsson.github.io/starplast/benchmark-0.43.html)

## 31. Compare a screen without losing unmatched rows.

Tools > Compare a screen

Open a CSV or TSV: Use a gene-level table and choose its identifier and value columns.

Inspect the join: Mapped and unresolved IDs remain visible. Measured zero stays zero.

Export the comparison: Keep the table and source-file hash together.

Screenshot values are illustrative. They are not biological screen results.

[Supporting documentation](https://einarolafsson.github.io/starplast/workflows.html)

## 32. Move from spaCR results to gene evidence.

spaCR: Generate gene-level measurement or effect estimates.

CSV / TSV: Export identifiers and the measured value you want to compare.

Starplast: Open Compare a screen and inspect the identifier mapping.

Evidence: Read the existing annotations and source context alongside the result.

The spaCR launcher opens Starplast in a separate environment. It does not automatically transfer results.

[Supporting documentation](https://einarolafsson.github.io/starplast/workflows.html)

## 33. Add measurements to a Toxoplasma session.

File > Import data: The advanced importer accepts CSV, TSV, Excel or Parquet for Toxoplasma accessions.

Check the preview: Choose a useful column prefix and explicitly handle duplicates and missing values.

Rebuild deliberately: Imported values do not change map positions until you include them in a new embedding.

Advanced table import is Toxoplasma-only. The guided comparison supports the currently open organism.

[Supporting documentation](https://einarolafsson.github.io/starplast/guide.html)

## 34. Export the evidence behind a result.

predictions.csv: Separate held-out evaluations from unknown-gene hypotheses.

per_class.csv: Inspect performance and support for each class.

run.json: Retain settings, input identity, split genes, excluded features and versions.

Use Export full run. A ranked list alone cannot reproduce an analysis.

[Supporting documentation](https://einarolafsson.github.io/starplast/workflows.html)

## 35. Use the same methods from Python.

import pandas as pd
from starplast import paths
from starplast.prediction import TaskSpec, run

nodes = pd.read_parquet(
    paths.cache_file('nodes.parquet'))
result = run(nodes, TaskSpec('compartment'))
result.save('localization_run')

Headless analysis: Read bundled tables and run methods without opening the desktop.

One contract: The API uses the same evaluation and export rules as the guided workflow.

The generated module reference documents functions, arguments and return values.

[Supporting documentation](https://einarolafsson.github.io/starplast/API.html)

## 36. Keep claims traceable and reviewable.

Observation records: Store entity, trait, value, units, context, source, uncertainty and missingness.

Literature assertions: Retain the source passage and negation. Extracted claims begin pending review.

Explicit limits: The API supports review records; the desktop has no bulk literature-curation interface.

A reviewed interpretation of a passage is not experimental proof. Conflicting observations remain separate.

[Supporting documentation](https://einarolafsson.github.io/starplast/workflows.html)

## 37. Run work in the background. Keep control.

Jobs and Stop: Watch progress in the Jobs panel. Cancellation happens between model fits.

Optional GPU: CPU analysis is the default. NVIDIA CUDA analysis uses the gpu extra on supported Linux systems.

Logs and help: Hover over settings. Enable file logging in Preferences when investigating a problem.

Offline browsing works immediately. Remote structures, publications and optional chat can require network access.

[Supporting documentation](https://einarolafsson.github.io/starplast/guide.html)

## 38. Choose a readable view.

View > Preferences: theme, colours, point rendering, lighting and text size.

Lighting changes presentation only. Use neutral lighting when comparing measurement colours.

[Supporting documentation](https://einarolafsson.github.io/starplast/guide.html)

## 39. Before interpreting a prediction, check four things.

Identity: Are the gene mapping, organism and assay context correct?

Evidence: Which values are measured, computed, missing or transferred?

Validation: Were target-related inputs excluded and related genes held together?

Uncertainty: Do rare classes, calibration and coverage support the intended interpretation?

Model comparisons guide exploration. Independent evidence is needed to establish a new biological claim.

[Supporting documentation](https://einarolafsson.github.io/starplast/benchmark-0.43.html)

## 40. Open a gene. Ask a precise question.

pip install starplast

User guide: https://einarolafsson.github.io/starplast/guide.html

Guided workflows: https://einarolafsson.github.io/starplast/workflows.html

Data sources: https://einarolafsson.github.io/starplast/datasets.html

Python API: https://einarolafsson.github.io/starplast/API.html

Benchmarks: https://einarolafsson.github.io/starplast/benchmark-0.43.html

Issues and source: https://github.com/EinarOlafsson/starplast

Starplast · Einar Olafsson · Open-source software, with original data sources credited.

[Supporting documentation](https://einarolafsson.github.io/starplast/index.html)
