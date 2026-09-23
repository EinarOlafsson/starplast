# Changelog

## Unreleased — 0.43.0

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
