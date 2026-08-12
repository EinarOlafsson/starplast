# A Validation tab, and renaming Meaning

Requested 2026-08-12.

## Rename tab 4: Meaning -> Inference

"Meaning" is a claim the tab does not make. It reports which held-out features distinguish the
clusters, which is an inference and is exactly as strong as the held-out testing behind it. The name
should say so. `Annotation` is the wrong choice here because annotation is what the new tab does.

So the pipeline reads:

    1 Data   2 Map   3 Clusters   4 Inference   5 Search   6 Validation

## Tab 6: Validation

Tests an annotation made from a clustering -- of the Search UMAP or of the Inference UMAP -- by
hiding labels that are already known and seeing whether they come back.

`starplast/validate.py` exists and does the work; this is the interface to it.

### The method

For a category (say dense granules):

1. Hide a fraction -- default 20% -- of the genes that carry it.
2. Pick the cluster holding most of the still-visible members. This is what a user does by eye.
3. Score against the hidden genes, who had no say in that choice.
   - **precision** over the cluster members that do not visibly carry the category: exactly the set
     a user would annotate from.
   - **recall** over the hidden genes.
4. Repeat over folds and average.

Run for every category, ranked. Per-category and never one global number: a method that recovers
hidden dense granules but not hidden rhoptries is not "60% accurate", it works for one compartment
and not the other, and only the per-category table says that.

### What the tab must show and refuse

- The per-category table: n labelled, folds, precision, recall, F1.
- The candidate list for a chosen cluster and category, with `cluster_frac_category` and
  `cluster_frac_contradicting` beside every gene. A candidate must never appear without them.
- **It must refuse to run when the category fed the embedding.** `masked_recovery` already raises;
  the tab must show that as an explanation rather than an error. A cluster matching a feature the
  map was built from is circular and the number would be meaningless.
- The caveat that the current implementation hides labels from the SCORING and not from the
  embedding. Re-embedding per fold is stricter and costs about an hour per fold on this proteome;
  `refit` is the switch, and whether it was used is recorded in the result rather than glossed.

### Saving an annotation

Only from this tab, and only with its numbers attached. Persist to a separate annotations file: the
gene, the proposed label, the configuration and cluster, the cluster composition, the validation
precision and recall for that category, the date, and free-text reasoning.

Never written back into the node table, and never coloured like a measured call. Measurement,
inference and absence never read as one another in this application; an annotation is a fourth thing
and needs its own colour.

### What the first run of this already says

Run against the coarse level-of-detail clustering as a stand-in, precision is 1-6% across
categories -- annotating from those blobs would be almost entirely wrong. That is the expected
result for five huge components rather than real clusters, and it is exactly why the number has to
be shown: the same candidate list, without it, looks like a finding.
