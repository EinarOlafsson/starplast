# Repository workflow

- Make changes on `nightly`. Keep `main` as the release branch.
- As each coherent change is finished and its relevant checks pass, commit and
  push it to `nightly`, merge it into `main`, and push `main`. Return to `nightly`
  for further work. Do not wait for unrelated tasks before publishing finished work.
- Preserve existing work and resolve branch divergence without force-pushing `main`.
- Attribute commits to Einar Olafsson, `einar.olafsson@gmail.com`. Do not add
  assistant, bot, or other co-author trailers. Einar remains the sole contributor.
- Give each commit a useful subject and a body explaining the problem, the change,
  and the relevant validation.
- Keep logos black and white with thin lines. README badges use standard
  contrasting labels, colours, and appropriate service icons.
- Keep the full linked dataset catalogue in `docs/datasets.md`, generated from
  `starplast.datasets.REGISTRY`; keep a prominent link to it in the README.
- Prepare releases with `python scripts/release.py bump VERSION` and update the
  changelog. Publication runs only when `main` receives a version increase, or
  through an explicit manual release retry. Never publish from `nightly`.
- PyPI projects and trusted publishers must belong to Einar's own PyPI account.
  Never register the packages under another account or claim publication succeeded
  without verifying the upload. Setup and release commands are in `docs/releases.md`.
