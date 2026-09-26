# The integrated neighbour space

`starplast/graphspace.py` builds one graph over genes from every permitted evidence layer plus the
measurement table, so that a user can ask *which genes are near this one, and on what evidence*. Two
strategies read it: **33 Put every layer into one space and read a gene's neighbourhood** and **34
Train on the networks and rank the edges they are missing**.

Everything below was measured on the shipped *Toxoplasma gondii* table (8,140 genes, 13 edge layers,
seed 42, default settings) on 2026-09-26. Reproduce it with:

```python
from starplast import graphspace, strategies
sp = graphspace.space(strategies.Context.shipped("Tg"))
print(sp.report.summary())
print(sp.report.frame().to_string(index=False))
print(sp.report.per_layer().to_string(index=False))
sp.gaps(50)
```

## What the space is

The thirteen measured layers are deliberately never merged in this project: each answers a different
question and they disagree with each other. The space does not change that. It is a **separate
object** with its own pairs, and nothing in it is written back into a layer. Every table it produces
carries a `kind` column (`measured` or `inferred`) and a `measured_in` column naming the layers that
record the pair; `measured_in` is empty for exactly the inferred pairs. An integrated edge is a claim
with a probability attached. A crosslink is an experiment.

**Candidate pairs.** For each gene: its partners in every permitted layer (at most 25 per layer,
strongest weight first) plus its 15 nearest genes in measurement space. On the shipped table that is
**276,384 pairs** of the 33 million possible -- 181,054 recorded by at least one layer and 95,330
inferred. The cap has a cost and it is not hidden: a pair in no layer and outside either gene's
nearest neighbours is never ranked, so it can never appear as a gap. `GraphSpace.edge_strength(a, b)`
will still score such a pair on demand; what the candidate set decides is what gets *ranked*, not what
can be *asked*.

**Permitted evidence.** Nine layers on this table (co-expression, co-fitness, shared compartment,
co-translation, shared domain, IP-MS, orthogroup, structure, crosslinks) and 367 measurement columns.
The literature layers are never sources, because co-mention follows attention; the derived layers
(`unwritten_interaction`, `structural_hole`) are never sources, because they are computed from the
others. Holding out a label removes its whole closure first -- the columns through
`Context.banned` and the layers through `Context.banned_layers`, the project's own guard rather than a
copy of it.

## How an edge strength is computed

Per pair, the sources are:

| source | what it is |
|---|---|
| one column per layer | the pair's weight as a **percentile among that layer's own edge weights**, 0 where the layer does not record it. A crosslink count of 5 and a co-expression correlation of 0.95 are not on one scale, and a single weight per source is meaningless until they are |
| one indicator per layer | 1 when the layer does not record the pair at all -- the missing-indicator |
| `measurement_similarity` | the cosine of the two genes' permitted measurement rows. `Context.matrix` is rank-scaled with missing values at the median, so this is a rank correlation across the measured columns in all but name |
| `shared_partners` | Adamic-Adar over the union of the permitted layers: a partner shared through a hub counts for less |
| `same_orthogroup`, `same_pfam` | only where the graph carries no such layer and the closure permits the column. On the shipped table both arrive as layers, so neither appears |

**Degree is not a source.** It is the strongest predictor of a measured edge and the least
interesting one: a model given degree learns which genes are well connected. The degree-matched null
below exists to measure that contamination, and feeding the model degree would be paying for the null
and then defeating it.

**The weights are learned from held-out edges.** One fold per layer worth holding out. In a layer's
fold that layer is removed from its own features *and* the shared-partner count is rebuilt on the
union of the other layers, because a partner shared across the held-out layer is that layer's evidence
arriving by a second route. Edges are held out **by orthogroup** (`Context.groups`), never at random,
so a hidden edge cannot be recovered through a paralog that stayed visible. A source's integrated
weight is then the mean of what it was worth in every fold, **counting as zero the fold where it was
the target** -- so a source that could only predict itself ends at weight zero as arithmetic, not as
good intentions.

What the shipped table's weights say, largest first:

| source | weight | folds it spoke in |
|---|---|---|
| `measurement_similarity` | +1.055 | 5 of 5 |
| `compartment` / `compartment_absent` | +/-0.944 | 5 of 5 |
| `coexpression_absent` | -0.547 | 4 of 5 |
| `shared_partners` | +0.534 | 5 of 5 |
| `cofitness_absent` | -0.524 | 4 of 5 |
| `coexpression` | +0.490 | 4 of 5 |
| `cofitness` | +0.489 | 4 of 5 |
| `domain` / `domain_absent` | +/-0.345 | 5 of 5 |
| `xlms` | +0.256 | 4 of 5 |
| `orthogroup` / `orthogroup_absent` | +/-0.235 | 5 of 5 |
| `struct` | +0.225 | 4 of 5 |
| `cotranslation` | +0.119 | 4 of 5 |
| `ip_ms` | -0.014 | 5 of 5 |

A layer whose edges all weigh the same contributes through two collinear columns (its weight and its
absence), and the fit splits one effect between them -- which is why the tables and
`GraphSpace.evidence()` group contributions per layer when they name the top evidence.

**The probability.** A calibrated logistic: Platt scaling fitted on the pooled held-out pairs, each
scored by the fold that could not see its own layer. It answers one question -- *given the permitted
evidence, is this pair a real relationship of the kind the measured layers record, rather than a
degree-matched non-pair?* -- at one negative per positive. It is therefore **conditional on that
comparison and not an absolute posterior** over all 33 million pairs. Its reliability is measured, not
assumed, and reported in the table below.

## The three nulls, measured

Every number is reported against three negative sets. This is the table to read before quoting
anything from strategy 33 or 34.

| model | null | AUROC | precision@k | Brier | reliability gap |
|---|---|---|---|---|---|
| logistic (**shipped**) | random | 0.678 | 0.954 | 0.213 | 0.025 |
| logistic (**shipped**) | **degree-matched** | **0.670** | **0.931** | **0.216** | **0.022** |
| logistic (**shipped**) | configuration | 0.653 | 0.916 | 0.227 | 0.038 |
| embedding | random | 0.682 | 0.931 | 0.213 | 0.021 |
| embedding | degree-matched | 0.678 | 0.919 | 0.215 | 0.019 |
| embedding | configuration | 0.659 | 0.890 | 0.226 | 0.033 |
| fame (diagnostic, never ranks) | random | 0.682 | 0.954 | 0.212 | 0.022 |
| fame (diagnostic, never ranks) | degree-matched | 0.658 | 0.947 | 0.220 | 0.032 |
| fame (diagnostic, never ranks) | configuration | 0.647 | 0.928 | 0.229 | 0.043 |

9,819 held-out edges over five layers. `random` draws non-pairs uniformly from every gene;
`degree-matched` pairs each positive with a non-pair of similar evidence degree at **both** ends;
`configuration` draws non-pairs from a stub-paired rewiring of the held-out edges, which keeps their
degree sequence exactly and their topology not at all.

**The fame gap** (random minus degree-matched) is printed as its own quantity:

| | fame gap |
|---|---|
| logistic, the shipped model | **+0.008** |
| fame: the same features **plus degree**, trained against random non-pairs | **+0.024** |
| the same recipe on a preferential-attachment graph (`tests/test_graphspace.py`) | **+0.216** |

Read that top row correctly. A small gap is the *result*, not a disappointment: this model is given no
degree and is trained against degree-matched non-pairs, so there was little fame in it to remove.
There are two further reasons it is small on this graph, and both are about the data rather than the
model. The shipped co-expression and co-fitness layers are built with a cap on partners per gene, so
their degree distribution is narrow -- a maximum of 74 partners on 8,140 genes for co-expression -- and
there is less fame available than in a literature-derived network. And the layers that *do* follow
attention -- the two co-mention layers, whose degree reaches 258 -- are excluded from this space
altogether.

The null earns its place by proving that rather than assuming it. On a graph whose edges *are* fame --
preferential attachment, built in the test file for exactly this purpose -- the same machinery reports
the usual recipe at **AUROC 0.922 against random non-pairs and 0.706 against degree-matched ones**,
with the Brier score collapsing from 0.114 to 0.316 because its probabilities were calibrated for a
world in which every non-pair is two obscure genes. The shipped baseline on that same graph scores
0.745 on the hard null: it *loses* to the fame model by 0.38 on the easy null and beats it on the hard
one. Its own gap there is negative, which is what "not exploiting degree" looks like as a number.

**Per held-out layer**, the shipped model against degree-matched non-pairs:

| layer held out | train edges | held-out edges | AUROC |
|---|---|---|---|
| co-translation | 3,418 | 451 | 0.964 |
| crosslinks (`xlms`) | 1,588 | 173 | 0.891 |
| co-expression | 28,091 | 3,023 | 0.779 |
| structure (`struct`) | 7,884 | 563 | 0.706 |
| co-fitness | 49,997 | 5,609 | 0.578 |

Co-fitness at 0.578 is the honest limit of this space: which genes have correlated fitness is nearly
unpredictable from the other evidence, and a gap ranking that rests on co-fitness should not be
trusted. Co-translation at 0.964 is the opposite case and is partly a statement about how that layer
was built. Read the per-layer table, not only the pooled number. The self-test
section below adds the configuration-model bar for each of these layers, and one
of them does not clear it.

For comparison: strategy 16 (`link_prediction`) reports **AUROC 0.97** for held-out crosslinks. It
scores against random non-pairs, keeps the target layer's own visible edges as features and holds out
edges at random rather than by orthogroup. The 0.891 above is the same question asked with those three
things changed.

## What the trained model does and does not add

The learned model is a spectral (truncated-SVD) embedding of the union of the permitted layers,
excluding the held-out one, whose per-pair features are the elementwise product of the two genes'
vectors, added beside the baseline features. Random-walk embeddings were not used because node2vec's
dependency is not installed and this project does not add a dependency to try a model.

**It does not earn its place on the shipped table, and so the baseline is what ships.** Measured gain
on the degree-matched null: **+0.0075 AUROC** pooled, with per-layer gains of +0.002 (co-expression),
+0.010 (co-fitness), -0.005 (co-translation), +0.023 (structure), +0.022 (crosslinks) -- a standard
error of 0.005 across folds. The bar it is held to is 0.01 AUROC *and* twice its own standard error
(0.011), and it clears neither. The reasons for that bar: a gain smaller than its own spread between
folds is a coin landing heads, and a gap table without per-source attribution is worth less than a
slightly worse one that can be read.

This is a negative result, not a failure. It is also not uniform: the embedding gains most on
structure and crosslinks (+0.023 and +0.022), nothing worth the name on co-expression, and loses
slightly on co-translation -- the profile of a model that finds something where the direct evidence is
thin and nothing where it is already thick. On the preferential-attachment graph in the test file it
clears the bar easily (+0.049), so the gate is a real gate and not a way of never shipping it.
Strategy 34 fits it on every run and reports it in
"the three nulls"; a user whose table is mostly crosslinks can set the model to `embedding` and see
what the gate says on their own data. If it clears the bar there, the strategy uses it and says so; if
it does not, the strategy says *that* and ranks with the baseline.

## How to read a gap

A gap is a candidate pair that **no permitted layer records**. Since every layer weight is zero for
such a pair, its probability rests only on `measurement_similarity`, `shared_partners` and any shared
annotation -- which is the right way round: a gap is an inference from everything except a direct
observation. From the shipped table's top of the list:

| gene A | gene B | probability | top evidence |
|---|---|---|---|
| TGME49_204020 (ribosomal protein RPL8) | TGME49_239100 (ribosomal protein RPS7) | 0.854 | `measurement_similarity +1.25`, `compartment -0.56`, `shared_partners +0.28` |

Read it as: the two genes' measurements across 367 columns are nearly identical, they share partners
in the union graph, and **no layer records the pair** -- including the compartment layer, whose absence
indicator is what the negative `compartment -0.56` means. It is a plausible pair of cytosolic
ribosomal subunits that none of these experiments happened to link, which is exactly what a gap should
look like. Three cautions:

1. the probability is against a degree-matched non-pair, so use it as a ranking with a scale, not as a
   chance of being true;
2. check the per-layer table first. If the gap's evidence is carried by a layer the evaluation could
   not predict -- co-fitness on this table -- the ranking is not supported;
3. a gap between two genes with nearly identical measurements may be one experiment's systematic
   effect rather than a relationship. The `top_evidence` column is there so that this is visible
   rather than buried.

## Cost

Building the space on the shipped table takes **30 to 45 seconds** with a peak resident set of
1.8 GB (measured with `/usr/bin/time -v`), most of which is the node table itself. Both strategies
declare `minutes`, since a run builds the space and evaluates it. The self-test of either strategy on
one layer takes 2 to 7 seconds depending on the layer.

## The self-test both strategies use

Pattern 3's shape with the two changes this module exists to make. One layer's edges are hidden by
orthogroup and the layer is removed from its own features. The metric is the AUROC of the hidden edges
against **degree-matched** non-pairs. The null re-run for the bar is **10 configuration-model
rewirings** of the hidden edges, which preserve their degree sequence, so a model reading fame cannot
clear it. A strategy passes above the null's 95th percentile by 0.05. The easy random-null number, the
fame gap, precision@k, the Brier score and the reliability gap are all reported as numbers beside the
verdict, never as the verdict.

On the planted organism at the default setting (hiding co-expression): AUROC 0.902 against
degree-matched non-pairs, against a configuration-model null of 0.555 +/- 0.044 -- PASS by +0.348. On
three null organisms built with every label and every edge dealt out at random: 0.41 to 0.57 against
nulls of the same size -- FAIL on all three, as it must be.

On the shipped table, every layer in turn (each measured with `graphspace.hidden_edge_test`). These
AUROCs differ in the second decimal from the per-layer table above because the self-test draws its own
negatives and fits its own model from its own seed; two honest samples of the same quantity do not
agree exactly, and a document that presented them as one number would be hiding that:

| layer hidden | hidden edges | observed AUROC | configuration null | effect | verdict |
|---|---|---|---|---|---|
| co-expression (**the default**) | 3,023 | 0.788 | 0.492 +/- 0.007 | **+0.296** | PASS |
| co-fitness | 5,609 | 0.581 | 0.500 +/- 0.004 | +0.082 | PASS |
| co-translation | 451 | 0.965 | 0.463 +/- 0.019 | +0.502 | PASS |
| structure | 563 | 0.701 | 0.536 +/- 0.014 | +0.165 | PASS |
| crosslinks | 173 | 0.875 | **0.832 +/- 0.010** | +0.043 | **FAIL** |

**The crosslink failure is the most useful row in this document.** Hidden crosslinks are recovered at
AUROC 0.875 -- the second-best number here -- and a configuration-model rewiring of those same
crosslinks is recovered at 0.832. The reason is structural: a crosslinked complex is a clique, so a
degree-preserving rewiring lands mostly inside the same complexes, and a pair inside a complex that no
layer records looks, to every source except the crosslink layer itself, exactly like a pair that is
recorded. So this space can say which **complexes** the crosslinking missed and cannot say which
**pairs inside them** it missed. A predicted crosslink from strategy 34 should be read as "these two
proteins are in the same neighbourhood", not as "these two proteins touch".

That is also why the default hidden layer is co-expression rather than crosslinks: for crosslinks the
configuration null is not a contrast but nearly the same question. The crosslink layer is one setting
away, and its FAIL is written here rather than avoided.

## Families and numbering

Strategies 33 and 34 join the existing eighth family, *Combine strategies*, rather than opening a
ninth: 31 and 32 combine strategies, 33 and 34 combine the layers, and a new heading would have named
the same idea twice.
