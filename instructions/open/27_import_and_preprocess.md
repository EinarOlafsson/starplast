# Import data, with the preprocessing offered rather than assumed

Requested 2026-08-12.

`File ▸ Import` should read a user's own table and put it into the map's feature blocks, offering
the same normalisations this project already applies to published data — and saying what each one
does to the numbers.

## What it has to read

Everything already imported here, which `datasets.py` and `sources.py` between them cover: CSV and
TSV, Excel (`.xlsx`), parquet, and the GEO-style matrices where genes are rows and samples are
columns. `tuning.import_table` already resolves identifiers, including the repair regexes for
`TgME49.208830`-style columns, and routes them through `identity.py` so previous and strain
accessions map forward — that part exists and must be used rather than rewritten.

## The preprocessing to offer

All of it is already in the code, applied to published data. The point is to make the same choices
available and explicit:

- **Quantification type** (`sources.NORMALISATION`): counts, FPKM, TPM, iBAQ, LFQ, log intensity,
  ratio, L2FC — which decides whether a log is taken and whether the column is centred. THE RANGE
  DECIDES, NOT THE FILENAME: `GSE108740_FPKM` reaches 16,520 and is real FPKM, while the columns
  derived from it in the node table top out at 9.7 because they were logged upstream. Show the range
  and let the user confirm.
- **Scaling** (`embedding.SCALINGS`): robust, z-score, rank, none. Rank is the safe default here.
- **Missing values** (`embedding.NA_POLICIES`): indicator, median, drop_columns, drop_genes.
- **Direction**: some screens are sign-inverted relative to others. Offer a flip, with the warning
  that pooling screens without rank-normalising first is how the 64x spread bites.
- **Duplicates**: mean, median, max, or first, for several rows per gene.

## Done when

A user's own CSV can be imported, its columns appear as a feature block that the map can be built
from, and every choice made along the way is visible and recorded with the imported columns.
