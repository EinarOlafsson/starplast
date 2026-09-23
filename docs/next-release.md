# 0.43.0 release completion

Release: 0.43.0. Scientific benchmarks, UI checks, documentation, packaging and
the complete regression suite passed before the version bump.

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
- [x] Full relevant regression, documentation, wheel and installed-app verification.
- [x] Publish 0.43.0 and verify the uploaded package and deployed documentation.

The source AF3 directories are inputs, not working directories. Index models and
extract features without altering the shared structure files.

Execution interruptions are recorded in [the task blocker log](task-blockers.md).
An interface restriction was reported on 23 September; its exact triggering task
was not identified. Continue independent software work and retain incomplete
research tasks explicitly in this checklist.

Release acceptance: the complete local suite passed (3,554 tests, nine skips);
subsequent UI and deck checks passed. The final pre-release [GitHub test run](https://github.com/EinarOlafsson/starplast/actions/runs/35912996485)
and [documentation build](https://github.com/EinarOlafsson/starplast/actions/runs/35912996668) passed.
The deck was checked in desktop/mobile browsers and rendered from both PDF and PowerPoint.

Published on 23 September 2026: [PyPI 0.43.0](https://pypi.org/project/starplast/0.43.0/)
and [GitHub release](https://github.com/EinarOlafsson/starplast/releases/tag/v0.43.0),
from commit `1162ebddddd7b95c4dc9f745e767601a2c47b29e`. The public wheel installed
from PyPI, passed dependency checks and rendered both organism maps with working
guided evidence. Its README metadata matches the repository exactly; distribution
checksums match the GitHub release assets. The [deployed 40-slide viewer](https://einarolafsson.github.io/starplast/deck/)
passed navigation, keyboard, mobile swipe, transcript and image-loading checks.
Public PDF and PowerPoint downloads match the committed files byte for byte.

Release automation: [full regression on the release commit](https://github.com/EinarOlafsson/starplast/actions/runs/35915144296),
[build and publication](https://github.com/EinarOlafsson/starplast/actions/runs/35915144281),
and [documentation deployment](https://github.com/EinarOlafsson/starplast/actions/runs/35915144196).
