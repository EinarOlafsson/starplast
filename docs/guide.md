# Using Starplast

Starplast opens a gene map and an evidence panel. Start with a gene or a table of
screen results, then compare the relevant measurements. The bundled organism maps
are separate: use **File → Species** to open the other organism.

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

## Appearance and performance

**View → Preferences** contains theme, colour maps, point rendering, lighting,
background, text size, window size, and logging. Hover over settings for guidance.
For clearer figures, use neutral lighting and disable distance fading when
comparing colours. Lighting and camera settings change presentation only.

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
