# The left panel becomes "color by", and colors by anything — DONE 2026-08-12

The panel said **filter by**, which was half of what it did and the less important half. It says
**color by** now, and it is the one place the map's colouring is chosen — because looking at
structure coloured by many different held-out variables is what this application is for.

## One list, three kinds of thing

Columns, kept clusterings, and quantities cut into bins, in a single chooser. One list rather than
three controls because they answer the same question — what should colour mean right now — and
having to know which of three places to look for the answer was the state this replaced. The value
list under it filters and flies by whatever the chooser names, so the two controls cannot describe
different things while sitting on top of each other.

## Runs are kept

`starplast/runs.py`. A clustering used to be computed, drawn and lost: the second run replaced the
first with no way back, which makes the one comparison this application exists for — does this
structure survive different settings — impossible to make by looking.

A kept run carries three things, each of which was missing:

- **a name**, from the clock to the second so two runs a minute apart are distinguishable without
  anyone typing anything, renameable afterwards (a name already in use is refused, not suffixed);
- **its recipe**, the embedding spec and the clustering parameters, saved beside it — without which
  a name is a label on nothing;
- **which genes it applies to**, not how many. A run over a subsample has fewer labels than the map
  has genes; placed by position, cluster 3's colour lands on whichever gene sits at that index. A run
  whose mask and labels disagree colours nothing at all, because an empty column beats 8,140
  confident mislabels.

Genes outside a run's subsample are **absent**, not unclustered. They were not put in no cluster;
they were not in the map.

## Quantities, binned

Quantile bins, not equal-width: nearly every quantity here is heavy-tailed — the fitness screens span
64x within themselves — and equal-width bins put 95% of the genes in one colour and call that a
colouring. The number of bins is the user's choice.

**Fewer bins than asked for is an answer and is reported, not worked around.** `n_publications` is
zero for most of this proteome, so its quartile edges are all zero and it yields one bin; splitting
that tie by rank or by width would draw four colours over a column with one level. The status line
says what fraction share a value.

That case also turned up a bug worth naming: `qcut` does **not** raise on a constant column, it
returns NaN categories, and casting those to `str` gives the literal `"nan"` — which would have
arrived in the interface as a category called "nan" sitting in the legend beside real ones. Both
paths now end in one honest bin, and a value that lands in no bin is absence.

## Verified

- `runs.py` at 100%; the whole suite green.
- Checked under the pandas 3 interpreter on the real table: colouring by a kept run, by a binned
  quantity at two bin counts, and renaming a run with its recipe intact.
