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
organism. There are three maps of all measurements (n_neighbors 10, 25 and 60),
one of all measurements except localization, and one for each evidence family
that has enough data: transcription, translation, protein abundance, fitness,
sequence, and so on. Maps are built from measurements only. Label columns are
removed before embedding. Each map is a 3D UMAP clustered with HDBSCAN
(min cluster size 25, min samples 5, leaf selection). A family map places only
the genes measured in at least half of its columns.

Click a map to show it in the central view. Its clusters become the map's
cluster coloring. **Color by label** colors the map by the label chosen in the
panel, and **Color by clusters** colors it by the clusters.

The table scores every categorical label on every map. Choose **One label, every
map** to see which map a label separates on best. Choose **One map, every label**
to see which labels a map organizes. Click a header to sort.

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

The files are `starplast/data/umap_gallery.npz`, `umap_gallery.json` and
`umap_gallery_scores.tsv`. They are rebuilt by `scripts/build_umap_gallery.py`,
and the run is recorded in `notebooks/umap_gallery_2026_09_29.ipynb`.

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
