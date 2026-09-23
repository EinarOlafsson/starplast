# Repository workflow

- Make changes on `nightly`. Keep `main` as the release branch.
- As each coherent change is finished and its relevant checks pass, commit and
  push it to `nightly`. For the current 0.43.0 work, keep development on `nightly`
  until all requested work and release checks are complete. Then bump the version,
  fast-forward `main`, and push `main` to trigger the release and PyPI publication.
- Preserve existing work and resolve branch divergence without force-pushing `main`.
- Attribute commits to Einar Olafsson, `einar.olafsson@gmail.com`. Do not add
  assistant, bot, or other co-author trailers. Einar remains the sole contributor.
- Give each commit a useful subject and a body explaining the problem, the change,
  and the relevant validation.
- Keep logos black and white with thin lines. README badges use standard
  contrasting labels, colours, and appropriate service icons.
- The approved logo is the Toxoplasma silhouette enclosing seven connected stars.
  Regenerate its SVGs with `scripts/generate_constellation_logo.py`; the active
  wordmarks are in `docs/assets/` and the app icon is `starplast/data/icons/starplast.svg`.
- Keep the full linked dataset catalogue in `docs/datasets.md`, generated from
  `starplast.datasets.REGISTRY`; keep a prominent link to it in the README.
- Prepare releases with `python scripts/release.py bump VERSION` and update the
  changelog. Publication runs only when `main` receives a version increase, or
  through an explicit manual release retry. Never publish from `nightly`.
- Publish only the `starplast` PyPI project. Keep optional GPU support in the
  `gpu` extra; do not add separate core, CPU, or GPU distributions.
- PyPI projects and trusted publishers must belong to Einar's own PyPI account.
  Never register the packages under another account or claim publication succeeded
  without verifying the upload. Setup and release commands are in `docs/releases.md`.
