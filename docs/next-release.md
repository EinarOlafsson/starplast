# Next release implementation

Target: 0.43.0. Keep the package version at 0.42.1 until the complete release has
passed scientific benchmarks, UI checks, and packaging checks.

- [x] Fixed-class evaluation, unsupported-node abstention, and consistent confidence labels.
- [ ] Shared embedding recipes, strict inference, executed-method provenance, and cache identity.
- [ ] Observation records with source, condition, uncertainty, missingness and review status.
- [x] Local AF3 model inventory, exact gene mapping, fragment handling, structure features and viewer access.
- [ ] Frozen protein sequence embeddings, with reusable caches and model provenance.
- [ ] Group-aware prediction benchmark: feature/PCA/UMAP neighbours, linear/boosted models, multi-view factors and weighted networks.
- [ ] Classification, regression, multi-label support, calibration and explicit abstention.
- [ ] Structured literature assertions with review and source evidence.
- [ ] Candidate explanations and experiment prioritization with stated assumptions/costs.
- [ ] Explore gene / Predict trait / Compare screen workflows and spaCR-compatible imports.
- [x] Forty substantially varied monochrome SVG logo proposals and comparison galleries.
- [ ] Reproducible real-data benchmarks, controls, limitations and release report.
- [ ] Full relevant regression, documentation, wheel and installed-app verification; publish 0.43.0.

The source AF3 directories are inputs, not working directories. Index models and
extract features without altering the shared structure files.
