# Next release implementation

Target: 0.43.0. Keep the package version at 0.42.1 until the complete release has
passed scientific benchmarks, UI checks, and packaging checks.

- [x] Fixed-class evaluation, unsupported-node abstention, and consistent confidence labels.
- [x] Shared embedding recipes, strict inference, executed-method provenance, and cache identity.
- [x] Observation records with source, condition, uncertainty, missingness and review status.
- [x] Local AF3 model inventory, exact gene mapping, fragment handling, structure features and viewer access.
- [x] Frozen protein sequence embeddings, with reusable caches and model provenance.
- [x] Group-aware prediction benchmark: feature/PCA/UMAP neighbours, linear/boosted models, multi-view factors and weighted networks.
- [x] Classification, regression, multi-label support, calibration and explicit abstention.
- [x] Structured literature assertions with review and source evidence (API).
- [x] Candidate explanations and budget prioritization with stated assumptions/costs (API heuristic).
- [x] Explore gene / Predict trait / Compare screen workflows and spaCR-compatible imports.
- [x] Forty substantially varied monochrome SVG logo proposals and comparison galleries.
- [x] Forty-slide introduction and practical guide with spaCR-style GitHub navigation, viewer, PDF and editable PowerPoint.
- [x] Reproducible real-data benchmarks, controls, limitations and [release report](benchmark-0.43.md).
- [ ] Full relevant regression, documentation, wheel and installed-app verification; publish 0.43.0.

The source AF3 directories are inputs, not working directories. Index models and
extract features without altering the shared structure files.

Execution interruptions are recorded in [the task blocker log](task-blockers.md).
An interface restriction was reported on 23 September; its exact triggering task
was not identified. Continue independent software work and retain incomplete
research tasks explicitly in this checklist.
