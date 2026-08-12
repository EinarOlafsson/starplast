# The UMAP gallery, the automated walk, and annotation from clusters

Requested 2026-08-12. This is the largest single piece of work left in the project and it is written
down before any of it is built, because the parts depend on each other and building them in the
wrong order wastes the expensive one.

## The workflow this has to support

The user's own description, kept because it is the acceptance test:

> Run a UMAP hyperparameter search using e.g. all columns except the hyperLOPIT label, then colour
> the graph by hyperLOPIT and see whether any visually-defined cluster is homogeneous. If all the
> GRAs form a cluster, look at the members of that cluster that do *not* carry the GRA label, and
> annotate those as candidate GRA proteins.

So the map is a hypothesis generator, and the held-out label is the check. Nothing here may present
a candidate as a finding: the full-proteome search recovers 1 of 24 compartments, so a homogeneous
cluster is a lead, and the interface has to keep saying so.

## 1. The gallery — see every UMAP, as it is generated

Today "walk hyperparameters" produces a table of numbers and no pictures, and "build this map"
replaces the single central view. Neither lets you look at 288 embeddings.

- **Grid mode**: a wall of thumbnails, one per configuration, each labelled with its parameters.
  Click one to expand it into the central view.
- **Scroll mode**: one large view with a slider stepping through configurations in order.
- **Both modes must populate incrementally.** A UMAP appears the moment it is computed, not when the
  walk finishes. This is the difference between a usable tool and a progress bar: a 288-run walk is
  half an hour, and the user needs to be looking at run 12 while run 13 computes.
- **An expanded UMAP behaves exactly like the main view.** Click a gene and get its evidence, colour
  by anything, draw edges. Not a picture of a map -- a map.

Implementation notes: the walk already stores each embedding through `EmbeddingStore`, so the
gallery reads from the store rather than holding 288 coordinate arrays in memory. Thumbnails can be
rendered offscreen with `grabFramebuffer` at low resolution. The incremental requirement means the
walk must emit a signal per completed configuration, not return a table at the end -- that is a
change to `tuning.walk_umap`'s contract and should be a generator or take a callback.

## 2. Automated mode — UMAP, then clustering, then scoring, per configuration

For every point in the UMAP hyperparameter walk:

1. build the embedding,
2. run a **HDBSCAN hyperparameter search** on it and keep the best clustering,
3. compute precision, recall and F1 for **each category** of the held-out variable against each
   cluster.

The result is a table of "how well does each category map onto a cluster, for each configuration",
which is what lets the user pick a map where the GRAs are clean even if everything else is a mess.

Three optimisation targets, all offered rather than one chosen:

- best single-category mapping (find a map that nails GRAs),
- best average across categories (find a map that is good everywhere),
- both (rank on each and show the frontier).

**Reuse `search.score_recovery`, do not reimplement it.** It already excludes absence labels, and
absence labels were what made the negative control read 0.654 instead of 0.339. It also already
reports `n_labels_recovered`, which is the column that matters -- mean F1 alone is inflated by a
single cluster whenever one class dominates, which is exactly the failure mode a per-category table
invites. See `results/full_proteome_2026_08_12/README.md`.

Clusters must be visualisable in this mode as well, not only categories.

## 3. Annotation — saving a candidate, and the honest label for it

When a cluster is homogeneous for a category, the unlabelled members of that cluster are candidates.
Needs: select a cluster, list its members split by label / no-label, mark the unlabelled ones,
persist.

Persist to a separate annotations file with, for every candidate: the gene, the proposed label, the
configuration and cluster it came from, the cluster's precision and recall for that category, the
date, and free-text reasoning. **It must never be written back into the node table**, and it must
never be indistinguishable from a measured hyperLOPIT call. The project's whole discipline is that
measurement, inference and absence never read as one another; an annotation is a fourth thing and
needs its own colour.

## 4. Validating an annotation — the part that is genuinely unsolved

The user asked for suggestions. Three that would work, in increasing order of strength:

**a. Label-masking recovery (cheapest, do this first).** Hide a known fraction of the labels for one
category -- say 20% of measured GRAs -- rebuild the embedding, and check whether the hidden GRAs
land in the cluster the visible ones form. This is a direct estimate of the precision of exactly the
inference being made, on genes where the answer is known. It costs one extra embedding per fold and
is the closest thing to a ground truth available.

**b. Leave-one-category-out.** Repeat (a) across every category. A method that recovers hidden
dense-granule proteins but not hidden rhoptry proteins is telling you the annotation is only valid
for some compartments -- which is far more useful than one global number, and consistent with the
finding that only `PM - peripheral 1` is recovered at full scale.

**c. Orthogonal evidence.** For a candidate GRA, ask whether it has a signal peptide, whether it
appears in a GRA-baited BioID or IP-MS study, and whether its expression tracks known GRAs. None of
these is proof, but a candidate supported by none of them should be reported differently from one
supported by all three. The IP-MS and BioID corpora are already indexed for this.

Note that (a) and (b) are the ones that produce a *number*, and that number is what any annotation
should be reported with. Without it a candidate list is a list of guesses with no error rate.

## 5. Interaction modes on the map

Left-click currently rotates and selects at once. Split it, as a mode toggle:

- **Navigate**: free rotate, or constrained to x / y / z, plus spin.
- **Select**: draw a gate and take the genes inside it -- 2D (a lasso in screen space, catching
  everything behind it) or 3D (a box or brush in world space). 2D is what you want for "grab that
  visual cluster"; 3D is what you want when the cloud is deep and the 2D lasso would catch the far
  side too.

A gated selection feeds annotation directly: gate a cluster, see its label composition, mark the
unlabelled members.

## Order of work

The gallery (1) is a prerequisite for everything else being usable, and the incremental requirement
is the hard part of it. Selection (5) is independent and small. Automated mode (2) is the biggest
compute change. Annotation (3) is easy to build and easy to build *wrongly* -- it must not ship
before validation (4a), or the project acquires a way to manufacture unvalidated claims, which is
the one thing it has been built to prevent.
