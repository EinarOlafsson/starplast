# Using Starplast

Starplast opens a gene map and an evidence panel. Start with a gene or a table of
screen results, then compare the relevant measurements. The bundled organism maps
are separate: use **File → Species** to open the other organism.

## Search Starplast

The box directly to the right of the **Help** menu searches everything the program
can do, as in spaCR. **Ctrl+Shift+H** (or **Help → Search Starplast…**) puts the
cursor in it. Type a few letters and a list drops down; **Up**/**Down** choose a
result, **Return** or a click goes to it, and **Esc** closes the list.

| Result | Where it takes you |
|---|---|
| command | Carries out the menu command; a submenu opens where it sits in the menu bar |
| panel | Shows the dock on the named tab, or opens Preferences, a guided workflow or the slot tree |
| strategy | Raises the Strategies tab with that strategy selected and its card showing (Details ▸ holds its Guide) |
| setting | Opens the analysis tab or Preferences page holding the control, scrolls to it and outlines it; a display choice such as **Theme: paper** is applied |
| slot | Opens the slot tree on the slot's organism with the slot selected |
| dataset | Shows what the dataset provides, its coverage, publication and source |
| guide | Opens the page of this guide at the heading |
| tutorial | Opens the tutorial in the browser at the section |

Words are matched against names first and descriptions second, and every word
must match, so a second word narrows the list. A method name finds every strategy
that uses it ("hdbscan", "logistic"), and a misspelling close enough still matches.
**Help → Keyboard shortcuts** lists every key the menus bind.

## Map controls

| Control | Action |
|---|---|
| Left drag in Navigate mode | Rotate the map |
| Scroll | Zoom |
| Click a gene | Select it and show its evidence |
| Search box | Find a gene ID or matching product text |
| Category filter | Restrict the map to checked categories |
| Double-click a category | Centre the camera on its genes |
| View → Left mouse button → Select | Draw a selection instead of rotating |
| Lasso (2D) | Select points inside a screen-space outline, including points behind one another |
| Brush (3D) | Select genes within a spatial neighbourhood around the pressed gene |
| Right-click the map | Access display options, export, or reset |

Relationships are selectable in **Edges**. By default, edges are drawn for the
selected gene. **Draw all active edges** can show the broader network, but dense
layers quickly become difficult to read. Co-mention, co-expression, shared
orthogroup, and measured interactions are different kinds of evidence.

## Pregenerated maps

The **maps** tab (beside Analysis) lists maps built ahead of time for the open
organism — 78 for *T. gondii* and 61 for *P. falciparum* — in six groups:

| Group | What the maps are built from |
|---|---|
| All measurements | Every measurement block, at n_neighbors 10, 25 and 60. The three fixed reference maps: their settings are declared, not searched. |
| All but one kind | Everything except one evidence family — except localization, except transcription, except fitness, and so on — so any label can be scored on a map that never saw the kind of measurement that restates it. |
| Evidence families | One family on its own: transcription, translation, protein abundance, fitness, sequence, regulation, relation, and so on. |
| Single experiments | One individual experiment block on its own: each screen, each expression atlas, each proteomics set, the structure and sequence blocks. A block needs at least 3 columns; one or two columns do not make a 3D map, and such a block is covered by its family. |
| Pairs of families | Two families that make a biological question: expression + fitness, fitness + proteomics, sequence + localization, regulation + transcription, and others. Each row says which question. |
| Triples of families | Three families: the whole expression cascade (transcription + translation + protein abundance), cost + amount + place, and a structure-first map with no expression in it at all. |

Maps are built from measurements only. Label columns are removed before
embedding, and the build refuses a map whose matrix contains one. Each map is a
3D UMAP clustered with HDBSCAN with leaf selection. A map built from a subset of
the evidence places only the genes measured in enough of its columns — the
threshold is in each row's tooltip.

### How the settings were chosen

Every map except the three fixed reference maps was **searched** for structure.
Six UMAP settings (n_neighbors 10, 25, 60 × min_dist 0.0, 0.25) and six
clusterings (min cluster size 15, 25, 40 × min samples 5, 10) were tried, and the
combination with the best **structure score** is the one that ships. The search
is successive halving: every UMAP setting is ranked on a 1,500-gene sample, and
only the best two are rebuilt over all the map's genes, because an embedding
costs about a hundred times a reclustering.

The structure score is **label-free** — no label is read while tuning, which is
what lets the labels be scored honestly afterwards. It is the geometric mean of
two things, because either alone can be gamed:

* the project's own `search.map_quality` score — the share of genes clustered
  times the evenness of the cluster sizes, and zero if the clustering merely
  bisected the cloud. This is about the *partition*: it says nothing about
  whether the clusters actually sit apart.
* the mean silhouette of the clustered points on the map's own coordinates,
  rescaled to 0–1 (0.5 is no separation at all). This is about the *geometry*:
  it says nothing about coverage, so two crisp clusters holding 5% of the
  proteome score beautifully on it alone.

Taking the geometric mean means a map has to be good at both, and a zero in
either is a zero overall. Each row of the list shows the score and its parts, and
each row's tooltip gives the full recipe and how many settings were searched.

### Finding a map

With this many maps the list is a tree. The controls above it:

* the **filter** box keeps only maps whose name, group, family or description
  contains every word you type — try `fitness`, `transcription` or `all but`;
* the **sort** box orders them by gallery order, by structure score, by how well
  the label chosen below maps onto them, or by how many genes they place;
* **Group** unticked flattens the tree into one best-first list, which is how to
  find the single best map for a question.

Click a map to show it in the central view. Its clusters become the map's cluster
coloring, and the line under the tree names the map on screen with its recipe and
structure score. **Color by label** colors the map by the label chosen in the
panel, and **Color by clusters** colors it by the clusters.

The table scores every categorical label on every map. Choose **One label, every
map** to see which map a label separates on best — the map on screen is marked
with ▶ and selected. Choose **One map, every label** to see which labels the map
on screen organizes; picking another map replaces the table. Click a header to
sort.

| Column | Meaning |
|---|---|
| Categories → clusters | For each category, the F1 of its best cluster, averaged with weights by category size. Unclustered genes count as missed. |
| Precision, Recall | The size-weighted means behind that F1 |
| Best category, n | The single category that maps onto one cluster furthest above chance, and its size on this map |
| Best category (size-aware) | F1 of the Wilson 95% lower bounds of that category's precision and recall. A 2-gene category alone in a cluster scores 0.34, not 1. |
| Skill, Best skill | The score rescaled against the same score with the label shuffled over the same genes (five shuffles): 0 is chance, 1 is perfect |
| Coverage | Share of the label's labelled genes that are on this map |
| Circular | The map was built from a column in the label's leakage closure, so recovering the label is expected |

**Score map on screen** scores every label on the map in the central view, for
example one built in Analysis. It uses the map's clusters if it has them, and
otherwise clusters it the same way. The result appears as "on screen".

The files are `starplast/data/umap_gallery.npz` (coordinates as float32, cluster
labels as int16), `umap_gallery.json` (every recipe, the structure score and every
setting searched) and `umap_gallery_scores.tsv` (the label × map table). They are
rebuilt by `scripts/build_umap_gallery.py`, and the run is recorded in
`notebooks/umap_gallery_2026_09_30.ipynb`:

```bash
systemd-run --user --scope -p MemoryMax=25G env STARPLAST_GPU=0 \
    python scripts/build_umap_gallery.py --workers 6
```

The build is checkpointed one map at a time into
`results/umap_gallery_checkpoint`, so a run that is stopped resumes where it left
off rather than starting over. To build **more** maps than ship — every block
rather than only those with enough columns, a wider settings grid, or your own
combinations of families — edit the thresholds and `COMBINATIONS` at the top of
`starplast/umap_gallery.py` and write the result somewhere else:

```bash
python scripts/build_umap_gallery.py --no-notebook --out ~/my_gallery
```

`umap_gallery.Gallery("~/my_gallery")` loads any such gallery, and
`MapGalleryPanel(code, gallery=…)` takes one, so a locally built gallery can be
used in place of the shipped one without rebuilding the package.

## Importing results

For *T. gondii*, use **File → Import data** for CSV, TSV, Excel, or Parquet files.
The current importer resolves Toxoplasma accessions; importing Plasmodium tables
is not yet supported. Choose the gene
identifier column and the measurements to import. Inspect the preview before
choosing transformations, missing-value handling, and duplicate aggregation.
Imported columns are joined onto the current organism's gene table. An unmatched
identifier cannot add a new gene to the bundled map.

Use a prefix that identifies the experiment, such as `spacr_egress`. Keep the
original table and the chosen import settings with your analysis. Importing a
column adds data to the session; rebuild the embedding explicitly if you want that
measurement to affect positions.

## Analysis settings

| Setting | What to consider |
|---|---|
| Feature blocks | Choose measurements relevant to the biological question. A block is a group of related input columns. |
| Scaling | Robust scaling reduces outlier influence; z-scores equalize variance; ranks emphasize order. No scaling preserves original units. |
| Missing values | Median imputation estimates missing inputs. Indicators expose missingness to the model. Dropping genes can greatly reduce coverage. |
| Neighbours | Smaller UMAP neighbourhoods emphasize local structure; larger values emphasize broader structure. |
| Minimum distance | Lower UMAP values pack points more tightly. Tighter groups are not necessarily better supported. |
| Random seed | Fix it for a repeatable configuration; compare several seeds to assess stability. CPU and GPU implementations can differ. |
| HDBSCAN minimum cluster size | Sets the smallest retained group. Unassigned points have the noise label `-1`. |
| DBSCAN epsilon | Neighbourhood radius in embedding units. It has no direct biological unit. |
| Held-out target | Exclude the target and related measurements from inputs before evaluating recovery. Inspect the exclusion report. |
| Search budget and sample size | Use small runs for exploration, then confirm candidates using the full gene set and independent evidence. |

The map, clustering, and inference tables can be clicked to revisit the associated
result. Right-click a results table to export it. Save named embeddings and results
to compare configurations later. A high score after trying many configurations
needs confirmation on evidence that was not used to choose the configuration.

## Start here

The **start here** tab, beside Strategies, is for when you do not yet know which
strategy applies. It asks one question at a time, in plain words, and every answer
is a button in the trail along the top that takes you back to change it.

1. **What do you have?** A single gene, a list of genes (a screen's hits, a complex,
   a pathway), your own measurement or screen, a label you want explained or extended,
   or nothing specific.
2. **Which organism?** The installed spaces. Strategies are calibrated per species, so
   this changes what is recommended.
3. **The subject itself.** A gene is searched for by accession, symbol or product; a
   gene list is pasted, loaded from a file, taken from the genes gated on the map, or
   filled with an example set; a label or measurement is chosen from a list that shows
   how many genes carry each one, so an almost-empty column is never picked blind.
4. **What do you want to know?** Find more genes like mine, predict this for the genes
   nobody has measured, explain what defines it, find partners or complex members,
   compare two conditions, or check whether it is learnable at all. Comparing two
   conditions then asks which measurement is the baseline.

It ends at three to six recommended strategies, each with the same card summary the
Strategies tab shows -- name, method, calibration grade and the four headline bars --
one line saying why it is on the list, and the settings your answers decided. **Run it
here** runs it with those settings; **Open in Strategies** opens it there with the same
settings, its guide and its self-test. Where the map gallery or the star map answers
the question better than running anything, it is offered too.

Ranking prefers, in that order, the strategies a goal is for, a scorecard task that
matches the goal, and the calibration grade measured for that strategy *on that
organism*; a strategy this table cannot fill is never offered. When nothing better than
a weak grade exists, the panel says so above the list rather than leaving it to be
noticed.

## Strategies

The **strategies** tab, to the right of Evidence and Analysis, lists named ways of
using the combined data for inference, grouped by how they work: searching the map
space, borrowing from neighbours, walking the measured networks, learning from
examples, starting from a gene list, contrasting two kinds of evidence, crossing
species, and combining strategies. The column beside each name is the verdict its
self-test earned on the shipped data.

1. Select a strategy and read **Guide**: what it infers, why that can work, how it
   fails, a walkthrough, and how it is tested.
2. Set its parameters under **Settings**. Every control has a tooltip. Gene-list
   strategies take pasted accessions, a file, an example set, or the genes gated on
   the map.
3. Press **Test (hold-out)** first. Known information is hidden, the strategy is asked
   for it back, and the answer is compared with the same procedure run on shuffled
   labels, random gene sets or permuted identities. PASS means it beat that null's
   95th percentile by the stated margin *with these settings on this table*.
4. Press **Run**. **Results** shows the summary and tables; click a row to find its
   gene, right-click to save, and use **Show on map** for strategies that build a map.

Every strategy removes the held-out label, anything that restates it, the experiment
that produced it, and any edge layer built from it before it looks at anything else.
The [strategy catalogue](strategies.md) lists all of them with their measured
verdicts, including the ones that fail on this data.

### Discoveries: functions, labels and tested inferences

The **discoveries** tab opens on **Labels and functions**. Browse functional domains
and protein families (InterPro and Pfam), enzyme classifications (EC), phenotypes,
stages, localization, structural labels and other available categorical annotations.
Search a term such as `kinase`, an exact identifier such as `IPR000719`, or a label
name. Select a label and then a class to see its known members; click a gene to open
its evidence. One gene can belong to several functional classes. Descriptions
retain installed source text where available; additional domain names come from
the pinned nomenclature snapshot. Names do not validate gene function or recover
the original annotation release.

Each label reports annotated genes, class count, precomputed claims, and evaluation
availability. **Open held-out scorecard** opens existing label/class results where
they exist. Legacy tests retain their original source-recovery scope. Missing tests
say **not evaluated**; known membership does not create a precision or recall score.
Missing annotation and a false annotation-presence flag do not establish biological
absence. Direct and orthology-derived EC fields remain separate.

For Toxoplasma, **Functional tests** also provides frozen complete-Pfam profile
recovery, domain membership cards, gene outcomes and training-only controls.
Select `pfam_id` and open its functional result, or choose it directly in the
test selector. The fixed candidate recovers 12 of 635 recorded profiles, makes
eight wrong calls and abstains on 615 genes. The majority control recovers 19
profiles. Domain cards measure recorded annotation membership; omissions do not
establish biological absence. Independent biological accuracy and calibration
remain unavailable.

For Plasmodium, select `ec_number` to open the frozen complete EC-major profile
test. It recovers 40 of 152 recorded profiles, makes 103 wrong calls and abstains
on nine genes; all four test genes with profiles unsupported in training remain included. The
majority control recovers 43 profiles and the prevalence control recovers 35.
The test selector also offers the fixed native random forest: 59 profiles
recovered, 93 wrong calls and no abstentions on the same 152 genes. It uses
300 trees and a minimum leaf size of two. Serial tree probability accumulation
keeps the frozen result exactly repeatable. Each strategy retains its own
scorecard; both retain the original training-only controls. Method support is
uncalibrated and is not a probability of biological correctness.
The test uses direct annotations; orthology-derived annotations remain separate.
Individual annotation curation/prediction provenance is unresolved. These scores
measure recorded annotation recovery; independent biological accuracy, calibrated
confidence and predictions for unknown genes remain unavailable.

The old view showed only labels with generated claims. Its default filters left
Toxoplasma localization claims visible, while hiding existing untested cell-cycle
claims. Functional domains and enzymes now remain browsable even where functional
inference and independent validation have not yet been implemented.

The **Inferred claims** view lists what Starplast infers about genes with no label -- a localization, a
phenotype, a stage -- and how each claim was tested. Every claim has:

- a **confidence**: of held-out genes given claims this confident, at least this share were right;
- a **status**: *tested* (evidence measured to be independent of the inference agreed or disagreed),
  *untested* (nothing independent reaches the gene), or *outside tested range* (unlike every gene the
  certainty was measured on, so no number is attached);
- a **lift**: confidence over how common the claimed class is anyway.

The legacy recipe's independence screen admits a check only if it does not repeat the generator's
mistakes much more often than independent evidence would. By default the tab shows tested claims at
80% or more with a lift of at least 2; every filter can be loosened. Click a claim for its reasoning,
**Colour map** to see a label measured and claimed, **Save…** to export. The gene card opens with the
same claims for one gene. These historical calibration and verifier tests are
separate from admission of independent biological benchmarks.

### The track record: would it have known?

Every labelled gene has been hidden once, in five folds that never split an
orthogroup, and each strategy that calls labels was asked what it is. The result
ships with Starplast and reads in four places, each one click from the next:

- **A gene card** ends with *If this gene were unknown*: which strategies would have
  named its class (✓), named another (✗), or said nothing (·). Click the class.
- **The class page** lists, per strategy, how many of that class's genes it got
  right, the rate with its 95% interval, and what it called them instead.
- **The category page**, one click up, lists every class with its best strategy:
  which distinctions the data carries, and which none of the strategies recover.
- **A strategy card** names the classes where it is weakest. Rates are only printed
  from five or more answered genes; below that it says *too few*.

Saying nothing is counted apart from being wrong. Genes are also hidden together, a
whole class at once and in random sets of 1 to 500, so the record shows whether a
strategy still works when a screen's worth of genes is unknown. To hide your own
list (a screen's hits, a complex):

```python
from starplast import strategies as S, track_record as T
ctx = S.Context.shipped("Tg")
rows = T.my_list(ctx, ["TGME49_294550", "TGME49_244470"], "compartment")
T.set_summary(rows)
```

Only the default label of each organism is shipped; `T.evaluate(ctx, key, target)`
builds any other.

## Star map

The **star map** tab, beside Strategies, draws one gene's links to other genes and
says where each link comes from. The chosen gene (the "star") sits in the middle and
the genes linked to it sit on rings around it.

* **Colour** is the source: each measured edge layer of the graph (crosslinks,
  co-fitness, co-expression, shared domain, co-mention and so on), each strategy, and
  one gold colour for your own runs.
* **Solid** lines are measured; **dashed** lines are inferred by a strategy.
* **Width** is strength: the link's score as a percentile within its own run, so
  links from strategies with different score scales share one width scale.
* **Hover** a gene for its id, product and link counts, or a link for its layer or
  strategy, run, setting and score.

Type a gene id and press **Centre**, or select a gene anywhere else (on the 3D map,
by search or in a results table) while **Follow selection** is on. Click a linked
gene to re-centre on it; **Back** (or Backspace) returns. Double-click a gene, or press
**Select on 3D map**, to select it on the 3D map. **1 hop / 2 hops** sets the depth.
**Per source** sets how many of its strongest links each source keeps around each
gene, so a dense layer such as shared compartment cannot hide a sparse one such as
crosslinks. The list on the right is the legend: untick a source to hide it, and read
how many links the centre gene has in each. **Measured** shows the data layers only.

Starplast ships the links from the default-setting runs of the strategies that relate
genes to genes: predicted pairs (16, 17, 18, 33, 34), modules and clusters (01, 04,
09, 15), set expansion from example seed lists (20, 25), and neighbour or partner calls
(07, 13). A module is not drawn as a clique: each member links to its three nearest
co-members, first by how many measured layers join them. Every strategy you **Run**
in the Strategies tab adds its links at once under **your runs**; they are kept in
the user cache (`star_edges/`) and come back in the next session.

### How much of the network to show

The box beside **Centre** chooses how much is drawn.

* **one gene (star)** is the view above: the chosen gene in the middle, its links on
  rings around it, with **1 hop / 2 hops** and **per source**.
* **neighbourhood** grows outwards from that gene for as many **hops** as you ask (up
  to four), taking the strongest links first, and stops at the gene cap (**up to N
  genes**, at most 2,000). Use it to see a pathway or a complex and what sits around it.
* **whole network** lays out every gene the connection rules leave, up to 3,000 genes
  and 30,000 links, and colours the twelve largest clusters so modules are visible at a
  glance. Genes in small clusters are drawn in the muted colour.

Both larger views say on screen what the caps left out, in so many words ("showing the
strongest 30,000 of 112,418 links and the 3,000 best-linked of 7,002 genes — the rest is
left out"). Nothing is silently dropped.

In the larger views: **drag** to pan, **scroll** to zoom around the pointer, **click** a
gene to make it the centre (**Back** returns), **double-click** to select it on the 3D
map, and **Fit** to frame the whole layout again. The small picture in the bottom-right
corner is the whole layout with your view marked on it; click it to jump there.

The layout is computed once for a given centre, mode and set of connection rules, and
then held: hovering changes only the text and the highlight ring, never a position. Pans
and zooms repaint a cached picture rather than re-running the layout.

### Defining the connections

The row under the controls is where you say **what counts as a link**. The resulting
network is counted live beside it — "4,812 links between 3,104 genes from 3 sources" —
together with what each rule removed.

* **The source list on the right** picks which sources may contribute at all: each
  measured layer, each strategy, and your own runs. **All**, **None** and **Measured**
  are shortcuts.
* **strength ≥** keeps a link only when its strength (its score's percentile within its
  own run) is at least that. Raise it to keep each source's best-supported links.
* **≤ N per gene** caps how many links any one gene may contribute, strongest first.
  This is what stops a hub — a gene in hundreds of shared-compartment links — from
  filling the picture by itself.
* **≥ N sources must agree** draws a pair of genes only when that many *different*
  sources say they are linked. Two independent sources agreeing is much better evidence
  than one, and this is the rule that says so.

Type a name in the box beside **Save** to keep the rules now set, and pick that name
again later to load them: "crosslink + co-expression, 2+ sources" comes back exactly as
you left it. **Forget** removes a saved definition. They are kept with your other state,
in `star_edges/star_definitions.json` under the user cache.

## Appearance and performance

**File → Preferences** contains theme, colour maps, point rendering, lighting,
background, text size, window size, and logging. Hover over settings for guidance.
For clearer figures, use neutral lighting and disable distance fading when
comparing colours. Lighting and camera settings change presentation only.

Menus, tooltips, drop-down lists and Starplast's own windows (Preferences, the
guided workflows, the slot tree, explanations) are rounded panes of translucent
black, or translucent white on a light theme, as in spaCR. These windows have no
title bar: drag the background to move one, drag an edge to resize it, and close it
with the **✕** in its corner (**Esc** also closes a dialog). On an X11 desktop without a compositor the
same panes are drawn opaque near-black with their corners cut; set
`STARPLAST_TRANSLUCENT=0` or `1` to override the check.

CPU analysis is available by default. The GPU switch uses installed, available
backends; it cannot install a driver. The GPU comparison tool benchmarks a small
sample and displays the resulting maps. A faster backend can produce different
coordinates for a stochastic embedding.

## Files and network access

Built tables and graphs are read from the package or `STARPLAST_CACHE`. Source data
for rebuilding is resolved through `STARPLAST_DATA`. Saved runs, annotations, and
downloads use the per-user Starplast directory or `STARPLAST_STATE`.

Browsing the bundled map is offline. Fetching protein coordinates, following
publication links, rebuilding from remote sources, and using a configured chat
service may require network access. Chat sends the supplied context to its configured
service; it is optional.

Enable file logging in Preferences when diagnosing a problem. Include the package
version, operating system, reproduction steps, and relevant log excerpt in a
[bug report](https://github.com/EinarOlafsson/starplast/issues).

## Tutorial videos

Six silent screen captures ship with the documentation site, under
`docs/tutorial/video/`, and are linked from the [tutorial index](tutorial/index.html)
and from each written tutorial. They run one to two minutes each, carry no audio, and
burn their captions into the frame.

| Clip | What it shows | Written tutorial |
|---|---|---|
| Find a gene and read its evidence | The find box, the evidence panel, and the rule that a dash is not a zero | 1 Explore |
| Start here, end to end | The guided tab from what you have, through the gene, the goal and the label, to the ranked strategies and Run | 3 Gene list |
| Test before you trust | A strategy card, its four bars, the hold-out test, and the verdict with its chance level | 4 Test and calibrate |
| Pick a map that shows your label | The pregenerated gallery sorted by how well a label maps, and the score table that says so | 2 Hold-out search |
| Walk the network from one gene | The star map: one gene's measured links by source, widened to its neighbourhood, re-centred on a neighbour | 5 Networks and agreement |
| One real biological question, answered | Question 13 of instruction 59, from the question to the named genes | 6 Advanced models |

They are recorded, not acted. `scripts/tutorial_video.py` drives the real `Window`
offscreen, grabs each frame with `QWidget.grab()`, paints the pointer, the highlight
and the caption over it, and encodes the frames with `/usr/bin/ffmpeg`. Every number
and table on screen was computed by the shipped code while the clip was being
recorded, and the two clips that run a strategy wait for the job and show the wait.

Two things are visible in the clips rather than hidden. The central 3D view is a
`QOpenGLWidget`, which the offscreen Qt platform cannot draw, so it records as an
empty strip and says so. And a strategy run takes as long as it takes.

To rebuild them:

```bash
QT_QPA_PLATFORM=offscreen python scripts/tutorial_video.py            # all six
QT_QPA_PLATFORM=offscreen python scripts/tutorial_video.py 4_pick_a_map
python scripts/tutorial_video.py --pages                              # relink only
python scripts/tutorial_video.py --notebook                           # record what they came to
```

`scripts/build_tutorials.py` relinks them at the end of its own run, so rebuilding
the written tutorials does not drop them. The clips are part of the documentation
site, never of the Python package.
