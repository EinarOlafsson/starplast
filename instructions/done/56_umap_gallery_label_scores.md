# 56 · A gallery of pregenerated maps, and how well each label maps onto each

**Status: complete (2026-09-29, worktree branch off nightly 0.48.0; not merged, no version bump).**

The user, 2026-09-29:

> the central umap, there should be a panel of pregenerated UMAPs based on all or some of the data,
> so the user can easily pick between them and explore how different labels map onto them. for each
> structure each label should be scored how well it maps onto each UMAP, by several columns: one
> score should reflect how well labels in a category map onto UMAP clusters (precision recall), and
> one for the top mapping label (is there one label that maps really well; this should be weighted by
> how many genes are in the label so I don't get a perfect score for a category with 2 genes with one
> label that maps consistently).

## Plan

1. `starplast/umap_gallery.py` (no Qt): the gallery recipes per organism (all measurements at three
   n_neighbors, all but localization, one per evidence family from the slot catalogue's axes via
   `strategies.Context.family_of`), built WITHOUT any label column, clustered with HDBSCAN; the
   label x map scores; loading and integrity checks.
2. Scores per (label, map): size-weighted mean best-cluster F1 / precision / recall over the label's
   categories; the best single category on the F1 of the Wilson 95% lower bounds of its precision
   and recall; both against a permutation chance level; a circularity flag from the leakage closure.
3. `scripts/build_umap_gallery.py` -> `starplast/data/umap_gallery.npz`, `umap_gallery.json`,
   `umap_gallery_scores.tsv`, recorded in an executed notebook. Heavy steps under a 30 GB scope.
4. `starplast/umap_gallery_panel.py`: a "maps" dock; clicking a map loads it into the central view;
   a sortable score table in both directions; one-click color by label; live scoring of the map on
   screen. Small hook in `app.py`.
5. Tests, `docs/guide.md` + `docs/API.md` sections, CHANGELOG under `## Unreleased`.

## STATE

* **Gallery.** 12 Tg maps, 14 Pf maps: `all_nn10/25/60`, `all_but_localization`, and family maps
  (Tg: fitness, localization, protein abundance, regulation, relation, sequence, transcription,
  translation; Pf adds chemistry and PTM). A family map needs >= 4 source columns and >= 400 genes
  measured in >= 50% of them, and places only those genes. Families refused: Tg PTM (254 genes),
  host effect, metabolism, chemistry, immunity, phenotype; Pf host effect, metabolism, immunity.
  Blocks are `embedding.default_spec`'s slot blocks, not `Context.blocks()`, because for Pf the
  latter falls back to prefix groups (one "family" per column prefix).
* **No label is an input.** Every `ctx.categorical_columns()` column is dropped from the table
  before embedding; Pf had ten 0/1 labels (`is_exported`, `has_domain`, ...) inside slot blocks,
  Tg one (`ortholopit_accepted`). The build raises if one reaches the matrix; a test checks the
  manifest.
* **Clustering: decided NOT to use the Clusters tab's opening setting (25/25, eom).** Measured on the
  shipped maps (`build_umap_gallery.compare_clusterings`, in the notebook): median 4 clusters, ~0%
  unclustered, headline-label skill median 0.01 (Tg compartment) -- the Tg transcription,
  translation, fitness and relation maps each came out as 2 clusters. Leaf selection at 25/5
  (same `clustering.cluster`): median 45-62 clusters, 40-52% unclustered, median skill 0.04-0.05.
* **Scores** (`label_scores`): categories -> clusters = size-weighted mean of each category's best
  cluster F1 (+ precision, recall); unclustered genes count against recall (a stated difference
  from `search.score_recovery`, which drops them). Best category = F1 of Wilson 95% lower bounds of
  precision and recall, chosen by skill over its own shuffled-label chance -- chosen on skill because
  on a one-giant-cluster map the majority value of a yes/no label had a raw bound of ~0.95 and the
  same bound when shuffled (first build: `is_exported` "False" won on the Pf sequence map). 2 of 2
  scores 0.342. Chance = mean of 5 label shuffles over the same genes; skill = (s - c)/(1 - c).
  Circular = map columns intersect `search.excluded_for(label)`.
* **Example (Tg compartment, category skill / best category):** localization-only 0.305 (circular;
  dense granules n=41, 0.60), all_nn25 0.092 (circular; 60S ribosome n=40, 0.54), all but
  localization 0.071 (60S ribosome 0.42), transcription 0.053, fitness 0.046 (dense granules n=180,
  0.33), relation 0.019. Pf lopit_pf_location: every "all" map is circular, including all but
  localization (its closure reaches outside the localization family), as is protein abundance.
* **Files:** `umap_gallery.npz` 1.88 MB, `.json` 0.19 MB, `_scores.tsv` 0.07 MB; wheel 55.3 MB
  (`check_wheel` passes). Build: `scripts/build_umap_gallery.py --workers 8` under
  `systemd-run --user --scope -p MemoryMax=30G`, STARPLAST_GPU=0 (umap-learn on the CPU), ~2 min;
  record `notebooks/umap_gallery_2026_09_29.ipynb`.
* **UI:** `umap_gallery_panel.py`, a "maps" dock tabbed with Analysis (hook: 6 lines in
  `app._gallery`). Thumbnails colored by each map's clusters; click -> `use_embedding` + the map's
  clusters as `cluster_labels` (not kept as a run); table in both directions with plain headers and
  tooltips; Color by label / Color by clusters / Score map on screen.
* **Tests:** `tests/test_umap_gallery.py` (planted scores, the 2-gene case, the giant-cluster case,
  noise as missed, file integrity per organism, no label inputs, every label x map scored, the
  panel both ways, the window loading a map and coloring in one click).
* Not done: the tutorials were not rebuilt; the gallery is not rebuilt automatically when the node
  tables change. The window re-indexes a map by gene id (`Gallery.aligned`), and
  `test_every_map_is_whole_and_over_this_table` fails when the gene table changes, which is the
  signal to re-run the build.
