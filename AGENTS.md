# Repository workflow

- Dataset selection rule (user, 2026-10-07): when multiple biologically suitable datasets
  answer the same information slot, use citation rate per year since first publication as
  the main preference, with a slight recency preference and a comprehensiveness preference.
  Apply the versioned policy in `starplast/dataset_selection.py` and retain its factors,
  publication identity, citation provider/retrieval date and coverage denominator. See
  `instructions/open/65_dataset_selection_audit.md`. Compare matching organism/host, stage,
  assay, quantity and units; check data quality, mapping and access before admission.
  Missing bibliometrics stay unknown; highly cited reviews or unsuitable datasets do not win.
  Additions/replacements require a recorded comparison or explicit unresolved gap, and
  affected data/strategy validation. Keep complementary and context-specific evidence.

- Read `NEXT_SESSION.md` and `instructions/open/64_information_space_and_inference_atlas.md`
  before selecting work. Instruction 64 is the current plan for browsing published organism/host
  evidence, precomputed gene/class/label inferences, tested meta-inference and ground-truth
  scorecards at every level. Reuse completed work from the earlier instructions.
- Keep the instruction-64 action cards and full progress table current. Each time an item is
  completed, replace its percentage with a green tick and show the user the entire table again.
  Completion requires the card's acceptance evidence, relevant checks, commit and nightly push.
- Discoveries priority (user, 2026-10-08): span available labels with function first.
  Follow `instructions/open/68_discoveries_function_coverage.md`: keep annotation
  membership separate from inferred claims; expose unavailable tests and link existing
  scorecards. Prioritize functional source/target benchmarks and precomputed results,
  with source-family exclusions and explicit curation/prediction/transfer lineage.
  Include actions 68.01/68.02 in each full completion progress report (51 rows).
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

- Work continuation/watchdog (user, 2026-10-08): keep the native goal active for
  unfinished instructions 64/65/66. Do not stop at inaccessible content; record
  the gap and advance the next permitted independent task. At each resumed turn,
  read `.starplast-watchdog/state.json` and recent `events.jsonl`, verify actual
  native goal status with the goal tool, then record the current item/progress
  through `python scripts/work_watchdog.py record`. Before a turn boundary or
  explicit stop, record the real reason (and observed goal status); normal final
  reports are continuation boundaries, not full-list completion. Do not infer a
  crash/overload from inactivity, bypass runtime restrictions, or start competing
  agents. See `docs/work-watchdog.md` and completed action 67.01. The systemd
  timer observes actual runtime events independently; preserve its private log.
