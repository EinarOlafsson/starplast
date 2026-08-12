# One script per dataset

28 scripts, one for each entry in `starplast.datasets.REGISTRY`. Each fetches its source, reads it,
resolves its identifiers to current ToxoDB ME49 accessions, and reports what that resolution
achieved.

```
python scripts/datasets/xue_singlecell.py
```

```
Single-parasite transcriptional atlas (cell cycle)
  provides   : Measured cell-cycle phase per gene, and pseudotime cluster
  citation   : Xue Y et al. eLife 2020;9:e54129
  note       : Tab-separated despite the .csv extension. RH files use TGGT1_ accessions …

xue_singlecell: fetched 0.1 MB -> …/elife-54129-supp3-v2.csv
read 964 rows x 9 columns from elife-54129-supp3-v2.csv
identifier column 'Unnamed: 0': 879 of 964 rows resolved (91%), 873 distinct genes
```

## What they are for

Inspecting, re-fetching or debugging one source without running a build that touches all 28. Each
script carries that dataset's provenance — citation, PMID, accession, URL, local path — and the
quirks that cost time the first time round, straight from the registry.

## What they are not

They are not a second implementation of the build. The shipped columns are assembled by
`starplast.build_graph.load_nodes()` through the domain modules (`localisation`, `expression`,
`screens`, `cellcycle`), and each script names the one that normalises its columns. A per-dataset
copy of that logic would drift and eventually misrepresent what produced the shipped data.

What these scripts do is the part that genuinely is per-dataset: fetch, read, and standardise the
identifiers. That last step is not a formality — published supplements cite whatever accession was
current when they were written, and `TGGT1_` is more common in them than `TGME49_`. A dataset keyed
on type I accessions joins zero rows without the strain tables, silently, because a join matching
nothing looks exactly like a dataset with no coverage.

## Twelve have no download URL

They report that, name the paper, and say where to put the file. Nothing is silently skipped: a
dataset that cannot be fetched is a gap in a reproduction, and printing nothing would hide it.

## They are generated

From the registry, by `scripts/generate_dataset_scripts.py`. Edit the registry (or `SPECIALS` in the
generator, for things the registry does not carry, like a workbook's sheet name) and regenerate:

```
python scripts/generate_dataset_scripts.py
```

A test fails if regenerating would change anything, so a registry edit cannot leave a stale script
behind.
