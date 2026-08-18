# 43 — Sweep every category at every level, tune each map, and infer from what recovers

**Status: DONE 2026-08-18.**

## What closed it

* **Every level of every hierarchy** — `search.sweep_levels` runs `sweep_categories` across all three
  facets and any set of levels, and returns `(table, skipped)` rather than logging the skips away. A
  sweep that covered nine of thirty categories and reported "nine categories" reads exactly like one
  that covered everything, and the commonest reason for a skip here is the circularity guard removing
  every block a category names — which is the fact a reader most needs.
* **Tuning per map** — `tune_umap` by successive halving and `tune_clustering`/`best_clustering` per
  embedding, landed earlier in the day. The finding that came with it: HDBSCAN's default
  excess-of-mass selection collapsed this proteome into 2 clusters on every configuration ever tried;
  with leaf selection the same embedding gives 37.
* **Per label AND per cluster** — `recipes.recovery_by_cluster` returns the full label x cluster
  matrix. `score_recovery` answers "was this label isolated somewhere", which is right for a score and
  wrong for reading a map: a label split cleanly across three clusters is a finding about
  sub-structure and comes back as one mediocre best-cluster F1. Pairs with no overlap are kept as
  zeros, because a missing row reads as "not computed" while a zero says the cluster was checked.
* **Step 4a, the genes** — `search.sweep_inference` returns the sweep table AND a table of named
  genes, each carrying its category, the held-out column, the cluster, the enrichment and the purity
  that produced it. The clustering reaches that function through an `on_clustering` callback rather
  than riding inside the results row: an array hidden in a DataFrame cell would be carried silently
  into every CSV downstream.
* **A biological-relevance score** — `recipes.relevance`, documented as a heuristic, with its four
  terms as COLUMNS beside the ranking so it never overrides the raw numbers it summarises. Strength
  (enrichment, saturating at ten-fold so the ranking does not sort by how rare a label is), reach
  (genes named, saturating at forty so a cluster that swallowed a compartment does not win), novelty
  (imported from `interpret.novelty` rather than reimplemented) and corroboration by the independent
  control — the term no other ranking in this project has. An uncorroborated claim is halved rather
  than zeroed, because several questions have controls too sparse to corroborate anything and
  zeroing would rank answers by whether their control happened to be dense.

**A predicted gene is never written back into the node table.** It comes out as its own table with
its provenance, exactly as this instruction required.

## What the sweep must do

1. Cover **every level** of the hierarchy, not one chosen level.
2. Cover **every information category at each level**.
3. Converge on optimal hyperparameters for **each UMAP**, computationally efficiently.
4. Converge on optimal hyperparameters for **each clustering run on each UMAP**, efficiently.

3 and 4 are the expensive part and the efficiency is the design problem, not an afterthought.

## The procedure

**Step 1.** Build many UMAPs from different category combinations.

**Step 2.** Determine the optimal clustering hyperparameters for each UMAP.

**Step 3.** Map the held-out dataset onto each UMAP and cluster, then ask whether sub-categories or
labels IN the held-out data map onto any of that UMAP's clusters — precision, recall and F1 for each
sub-category/label, for each cluster.

**Step 4.** If held-out labels do map onto a structure, that UMAP scores high, and then:

**Step 4a.** Label the UNLABELLED members of clusters that are mostly one held-out label. Save which
genes were labelled with what, and the precision/recall/F1 for that cluster, for that held-out label,
for all labels in that held-out dataset, and for all clusters. The user reads those inferences and
judges biological relevance. **A scoring system for biological relevance would be ideal.** As a
backup the user inspects cluster/label overlap on the generated UMAPs and infers manually, but the
fully automated route is preferred.

## What this changes about the sweep as built

The current `search.sweep_categories` does one level, one hyperparameter setting, and scores only
whether a held-out category associates with the clusters. It does not tune, does not sweep levels,
does not name genes, and produces no artefact a person can read biologically. It is step 3 with a
fixed configuration and no step 4.

## Constraints that already exist and must not be relaxed

* **A degenerate clustering is not scored** (`clustering.degenerate`). Tuning must optimise toward
  clusterings that pass that guard, not learn to defeat it.
* **The held-out columns must never be among the columns the map was built from**, and `battery`
  marks anything the map saw as `used`/`derived` so only `held_out` rows count.
* **Leakage closure spans the class**, not the leaf (`excluded_for`, scope at the first two levels).
* A gene labelled by inference is a PREDICTION and must never be written back into the node table as
  though it were a measurement. It goes to its own artefact with its provenance.

## Efficiency, since 3 and 4 are the cost

* Successive halving over configurations: score all cheaply on a sample, keep the top fraction, spend
  more only on survivors.
* Cache the built matrix per block set — the same matrix is reused across clustering parameters.
* The clustering sweep is far cheaper than the UMAP; tune clustering fully per UMAP, and be selective
  about which UMAPs get built at all.
* Every run records its recipe and seed, as `search` already does, so a result can be rebuilt.

## Done when

* The sweep covers every level and every category, and says what it skipped and why.
* Each UMAP has tuned clustering parameters chosen by a stated objective.
* Step 3 produces per-label precision/recall/F1 per cluster, not one number per category.
* Step 4a writes a readable table of inferred gene labels with the statistics that produced them.
* A biological-relevance score exists, is documented as a heuristic, and never overrides the raw
  numbers it summarises.
* The UMAPs are inspectable, so manual inference remains possible when the automated score is unsure.
