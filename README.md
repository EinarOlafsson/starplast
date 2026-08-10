# starplast

A 3D browser for the *Toxoplasma gondii* knowledge map. 8,140 genes as points in space, six kinds of
relation as toggleable edges, and one panel per gene showing everything that is actually known about it —
including when the answer is "nothing".

```bash
pip install -e .
python -m starplast.build_graph     # one-off, ~2 min: builds data/graph.npz
starplast
```

Needs a display (PyQt6 + OpenGL). `build_graph` needs the source datasets listed below; the built cache
(`data/graph.npz`, `data/nodes.parquet`) is self-contained, so the app runs without them.

## What makes it different from a network viewer

**Position means something.** Node coordinates are a 3D UMAP embedding of a 43-feature matrix — stage
expression, seven CRISPR fitness screens, hyperLOPIT compartment, paralog number, domain count,
phosphosites, mean AlphaFold pLDDT. Two genes near each other are biologically similar. A force-directed
layout would look similar and mean nothing.

**Six edge types, never merged.**

| edge | n | source |
|---|---|---|
| `comention` | 382 | 33,924 PubMed abstracts (2 shared abstracts minimum) |
| `orthogroup` | 3,452 | OrthoMCL |
| `coexpression` | 49,293 | GSE108740 stage series, top-25 neighbours at r ≥ 0.95 |
| `compartment` | 118,712 | hyperLOPIT |
| `cofitness` | 88,997 | 7 CRISPR screens, top-25 at r ≥ 0.90 |
| `domain` | 10,399 | shared InterPro domain |

They answer different questions and disagree with each other; merging them into one "relatedness" score
would be the single easiest way to make this tool lie.

**Attention correction is on by default.** Raw co-mention edges reproduce the literature's popularity
contest — ROP18 and GRA16 become hubs because they are studied, not because they are central. The default
view shows log2(observed / expected) given each gene's own publication count, so an edge means *more
co-mention than these two genes' individual fame predicts*. Switch it off and the status bar says
`RAW co-mention (attention-biased)`.

**The most important number in the app is 601.** Only 601 of 8,140 genes are named in any of 33,924
abstracts, by symbol or accession. 93% of the proteome has never been the subject of a sentence. Genes below
7 abstracts get an explicit warning in the detail panel: absence of evidence there is absence of attention.

## Interpretation rules the UI enforces

- **"unassigned" is grey, not a 27th compartment.** hyperLOPIT assignment tracks protein abundance, so a
  missing call means unknown, never absent.
- **Fitness scores are competitive growth, not essentiality.** The in vitro screen is predictable from
  protein features (R² = 0.45); the five in vivo screens are not (−0.11 to +0.10). The panel says so.
- **Missing values are rendered grey, never mapped onto the colour scale.** Missingness here is
  informative — large secreted proteins are exactly the ones AlphaFold DB skips.
- **Draw caps are stated.** Above 20,000 drawn edges the status bar reports what was dropped.

## Controls

| | |
|---|---|
| left-drag / scroll | rotate / zoom |
| left-click a node | select; fills the evidence panel |
| search box | gene ID or product substring, then flies to it |
| double-click a compartment | fly to that compartment's centroid |
| level of detail | compartment (galaxy) → orthogroup (system) → gene (planet) |
| edge checkboxes | per-type; edges draw for the selected gene unless "draw all" is on |

## Data sources

Built from `/mnt/firecuda2/Claude/toxoplasma_projects` (home) — see `HANDOFF.md` for the full table and the
work-machine paths. `data/toxodb_gene_names.tsv` is committed and regenerable with
`python -m starplast.fetch_names`.

## Status

v0 + v1 complete: precomputed embedding, GL scatter, picking, evidence panel, edge toggles, attention
correction, level-of-detail. v2 (species switching via cross-species orthology) and v3 (continuous star-map
zoom) are deliberately deferred — see `HANDOFF.md`.
