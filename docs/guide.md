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
organism — dozens of them, in six groups:

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
