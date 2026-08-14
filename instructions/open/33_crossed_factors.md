# 33 — Crossed factors: is fragmentation a defect or a result?

**Status: open. Not started.** Written 2026-08-14 from an observation made during the first large
search. This one changes what a number means, so read the argument before touching the code.

## The observation

The yield objective keeps winning with very fine clusterings — 60, then 224, 488, 573 clusters over
8,140 genes. The reflex reading is that `discovery.yield_score` sums over findings, so more clusters
means more chances at a claim, and the optimiser has found a way to game it. That reading was
written into HANDOFF.md as a defect.

**It may be wrong, and the reason is biological.** A cluster is not obliged to correspond to one
category of one layer. Golgi proteins in the tachyzoite and Golgi proteins in the bradyzoite are
plausibly two distinct groups — same compartment, different life-cycle stage — and a map with the
resolution to separate them SHOULD return two clusters where a coarser one returns one. Under the
current code both come back as "enriched for golgi", the second looks like a duplicate of the
first, and splitting looks like cheating.

The general form: **the clusters are cells of a cross-tabulation.** Layer A (localisation) has
categories, layer B (stage, cell-cycle phase, fitness regime) has categories, and the map is
resolving the A×B grid rather than either margin. Fragmentation is then not inflation — it is the
map having enough resolution to see an interaction.

## Why the code cannot currently tell the difference

* `discovery.guilt` tests one layer at a time. Two clusters both enriched for golgi are two
  findings, and nothing asks whether they differ in any other layer.
* `discovery.disagreement` is close but asks a different question: it wants a cluster that is
  homogeneous in A and **split within itself** on B. The crossed case is the opposite shape — each
  cluster is homogeneous in BOTH layers, and it is the *set* of clusters that spans B.

So the thing to add is a third claim shape, and a run-level diagnostic that decides which reading of
fragmentation is right for a given run.

## What to build

### 1. `discovery.conjunction(nodes, labels, layer_a, layer_b, ...)`

One row per (cluster, a-category, b-category) where the cluster is enriched for the **combination**
beyond what either category explains alone. The test that matters is not "is this cell enriched" --
with 8,140 genes almost any cell of a large grid is -- but "is the joint enrichment more than the
product of the margins", which is an interaction and should be tested as one. A 2x2 collapse
(this a-category vs rest) x (this b-category vs rest) inside the cluster, against the same table
outside it, is enough; Fisher's exact on the collapsed table, and report the ratio

    lift_joint / (lift_a * lift_b)

as the effect. A value near 1 means the cluster is just the two margins meeting by chance; well
above 1 means the map has resolved the combination.

Guards carry over from `guilt` unchanged: minimum cluster size, minimum category count, `MAX_SHARE`,
BH correction across the whole family — and the family is now large (clusters x |A| x |B|), so the
correction does real work. Circularity marking applies to **both** layers.

### 2. `discovery.explains_fragmentation(nodes, labels, layer_a, layer_b)`

The diagnostic that settles the argument for a given run. Take every cluster enriched for the same
a-category; ask whether those sibling clusters separate on B more than clusters drawn at random
would. Something like: the mean pairwise distance between siblings in B-composition space, against
a permutation null that shuffles B within the a-category. Return one number per a-category and one
for the run.

Read it as: **"of the 573 clusters, N are siblings that differ by stage"**. If that fraction is
high, the fine clustering is resolving structure and the yield score is measuring something real.
If it is near the null, fragmentation is inflation after all and the objective needs the per-gene
normalisation already proposed.

Run it over the searches already saved — `~/.cache/starplast/searches/bigA_*`, `bigB_*` hold the
labels of every configuration, which is exactly what this needs and is why they are stored.

### 3. Wire it through

* `optimize.MODES` gains `conjunction`; `evaluator` calls it with `layers` and `against`.
* `interpret.sentence` gains the wording. Suggested shape:
  *"Cluster 41 is 78% golgi AND 71% bradyzoite-enriched — 3.1x more than either alone would
  predict (q = 2e-9). Cluster 12 is golgi and tachyzoite. These are the same compartment in two
  stages, and the 9 unlabelled genes in cluster 41 are predicted golgi proteins of the bradyzoite."*
* The competing-claims caveat added in 0.27.0 needs revisiting: two claims on one cluster are
  alternatives when they are two categories of the SAME layer, and are not alternatives when they
  are categories of different layers. The current wording says "alternatives" for both.

## What the current data can and cannot support

This matters, and it is the part most likely to be got wrong.

**The exact example does not work today.** hyperLOPIT is a tachyzoite experiment. `compartment` and
`compartment_best` carry one localisation per gene, measured in one stage, so "tachyzoite golgi vs
bradyzoite golgi" cannot be validated against this table — there is no measured bradyzoite
localisation to be right or wrong about. A crossed finding involving localisation x stage is a
**hypothesis generated from expression behaving differently**, not an observation of two
localisations, and the wording must not blur that.

**What is crossable today:** `compartment_best` x `cellcycle_phase` (873 genes measured),
`compartment_best` x `stage_enriched_derived` (1,911 genes, and note it is DERIVED from expression —
`search.SAME_QUANTITY` should already flag the circularity), `cellcycle_phase` x fitness regime,
localisation x fitness.

**What would make the real claim testable:** stage-resolved localisation — a bradyzoite hyperLOPIT.
That is a slot with no dataset in it. See `instructions/open/31_slots.md`; this instruction is the
best argument yet for filling that slot, and it should be recorded there as the reason.

## Acceptance

* A synthetic fixture where the truth is known by construction: four groups built as the 2x2 of
  (compartment A/B) x (stage X/Y), and `conjunction` recovers exactly those four with the
  interaction ratio well above 1, while `guilt` on either layer alone returns two.
* A second fixture where A and B are independent: `conjunction` returns nothing, because a cell
  enriched at exactly the product of its margins is not a finding.
* `explains_fragmentation` run over the four saved searches, with the number reported in the commit
  message. That number decides whether the "yield rewards fragmentation" entry in HANDOFF.md §4
  stays a defect or becomes a feature, and the entry should be rewritten either way.
* 100% coverage, no `pragma`, version bump, no `Co-Authored-By` trailer.
