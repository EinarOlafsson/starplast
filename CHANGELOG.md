# Changelog

## 0.44.0

- Add a **Strategies** tab to the right of Evidence and Analysis: 32 named ways of using the combined data for inference, grouped into eight families, each with a tooltip, an explanation, a walkthrough, settings that explain themselves, and Run / Test / Stop.
- Give every strategy a self-test that hides known information -- labels, gene-set members, edges or values -- asks for it back, and compares the answer with the same procedure on shuffled labels, random sets or permuted identities. A strategy passes only above its null's 95th percentile by a stated margin.
- Include the two founding questions as strategies 01 and 02: hold out a category and search maps for the structure that recovers it, and find the map where a gene list forms one cluster with high precision and recall.
- Measure every strategy on the shipped tables and show the verdict beside it: 26 pass, 5 fail and 1 is inconclusive on *T. gondii*; 24, 5 and 3 on *P. falciparum*. The [strategy catalogue](docs/strategies.md) and `results/strategies_2026-09-25/` record every number, failures included.
- Add `starplast.strategies` (context, leakage guard, five holdout test patterns, a planted organism for testing) and `starplast.strategy_catalog`; strategies run headless as well as from the tab.
- Add a `join="louvain"` option to `methods.multiplex_communities`: joining agreed pairs by connected components chains every gene into one community when layers agree about different genes.
- Group evidence into families by the slot catalogue's axis, so knockout screens named for a second background or a genetic interaction count as fitness.

## 0.43.0

- Add a 40-slide introduction and practical guide, presented like spaCR's deck: README cover, linked GitHub slide pages, web viewer, PDF and editable PowerPoint.
- Add guided Explore a gene, Predict a trait and Compare a screen workflows with settings help and full run exports.
- Evaluate feature, linear, boosted, PCA, UMAP, masked-factor and weighted-network predictions with grouped folds, target exclusions, calibration and unsupported-call abstention.
- Add classification, regression and multi-label APIs; preserve unknown outcomes and report fixed-class metrics.
- Include 320 frozen ESM-2 sequence features for 8,064 proteins and nine local AF3 summaries for 1,210 exactly mapped genes. Add model-indexing and sequence-encoding commands.
- Publish 29 reproducible localization and abundance comparisons, source ablations, shuffled controls and uncertainty estimates in the [benchmark report](docs/benchmark-0.43.md).
- Use one balanced embedding builder for packaged and interactive maps; record actual algorithms, input hashes and ordered gene identities.
- Add source-level observation records, reviewed literature assertions, candidate explanations and an explicitly heuristic budget API.
- Restore the Plasmodium evidence ledger and prevent fixture builds from overwriting packaged data.
- Share OpenGL resources between organism windows so switching species retains valid shader programs.
- Expand CI to the complete non-slow test suite, repair documentation and display regressions, and declare numerical thread-control dependencies.
- Repair download badges, retain the linked 130-source catalogue and publish forty additional monochrome logo proposals; keep the approved logo active.

## 0.42.1

- Use the GitHub README as the PyPI project description.
- Use full URLs for README artwork, the rotating gene map, and documentation links.

## 0.42.0

- Simplify the map to individual genes; remove Galaxy and Orthogroup summary modes.
- Adopt the Toxoplasma constellation logo, SVG wordmarks, and window icon.
- Rewrite the README and add a user guide and Python API guide.
- Publish generated API documentation through GitHub Pages.
- Add settings help and document application callbacks.
- Fix importing screen tables whose identifier column is already named `gene_id`.
- Fix Plasmodium gene selection and link its evidence panel to PlasmoDB.
- Publish the application and bundled data as one PyPI project, `starplast`.
- Make CUDA dependencies optional through `starplast[gpu]`.
- Include SVG artwork and compressed sequence tables in wheels; exclude saved embeddings.
- Add synchronized version bumps, build checks, and automatic PyPI publishing.
- Develop on `nightly`; publish version increases on `main` to PyPI and GitHub Releases.
- Link all 128 registered datasets and computed layers to their sources.
- Show a 1440×1080, 30 fps gene-map rotation with selected-gene lighting in the README.

## 0.41.0

Existing development baseline before the packaging and documentation overhaul.
Earlier changes are recorded in [HANDOFF.md](HANDOFF.md) and the Git history.
