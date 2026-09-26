# How good each strategy is: calibration over settings, targets and seeds

**Status: complete (v0.45.0, 2026-09-26).**

A strategy's self-test at its defaults says whether it works there. The question the user asked was
broader -- "exhaustively calibrate each strategy so the next time I try the application I know how
good each strategy is, at least on historical data" -- and a single verdict cannot answer it. This
records how the answer was measured and, as much as the numbers, what measuring it uncovered.

## What was run

`scripts/calibrate_strategies.py` runs every strategy's self-test over a grid of that strategy's
own settings (neighbourhood sizes, cluster sizes and selections, thresholds, models, edge layers),
over several held-out labels per organism, and over five seeds -- a seed changes which genes are
hidden and how every map is drawn. **5,940 self-tests, 34 strategies, both organisms.** Results are
in `results/calibration_2026-09-26/`; the shipped summary is `starplast/data/strategy_calibration.json`,
which the Strategies tab, the README table and `docs/calibration.md` are all built from.

## How a number is summarised (`starplast/calibration.py`)

* **Skill** = (observed - chance) / (1 - chance). Zero is the same procedure on shuffled data, one is
  perfect, and a correct-call rate, an AUROC, an F1 and a correlation all read on one scale.
* **Tuned** is chosen on seeds 1-3 and reported on seeds 4-5. Choosing and reporting on the same
  runs reports the luckiest of many settings; a test asserts that a setting lucky only on the
  choosing seeds comes out ordinary.
* **Intervals** are a two-stage bootstrap: held-out targets, then runs within a target. Runs on one
  label share its labels and are not independent, and treating them as such gives intervals several
  times too narrow; a test asserts the interval widens when the variance sits between targets.
* **Grades**: reliable, works when tuned, weak, no skill, untestable.

## The result

On *T. gondii*: 26 reliable, 7 weak, 1 untestable. On *P. falciparum*: 20 reliable, 8 weak, 3 work
only when tuned, 2 untestable, 1 no skill. At their defaults, 28 of 34 pass their self-test on
*T. gondii* and 27 on *P. falciparum*.

The sweep was run twice. The first ran from a code snapshot taken before the second and third data
waves (25 Toxoplasma and 23 Plasmodium columns) and before two screens were consolidated into one
slot each, which changes the blocks maps are built from; since a calibration describes the table it
ran on, it was discarded and the sweep re-run on the final code and tables. The numbers above are
the second run's. The full table is in the README and every number -- per
target, per setting, per parameter value -- is in `docs/calibration.md`.

## What calibrating uncovered

Sweeping settings is also a search, and it found the settings under which a test stops meaning what
it says. Five were not strategies failing; they were tests that could be satisfied without the
inference they claim to test. Each was fixed rather than reported.

1. **A null that picked the biggest cluster.** Strategies 01, 04 and 15 compared their F1 against a
   second search on shuffled labels. A search on random labels picks the largest cluster for every
   label, and a large cluster scores F1 near twice the label's share on size alone, so the null beat
   the real labels on four targets of eight. The null now permutes the hidden genes' labels over the
   SAME chosen clusters, which leaves a large cluster exactly as good under the null and earning
   nothing. Strategy 01 -- the founding question of this project -- passes its *T. gondii* self-test
   under the fair null, 0.132 against 0.039; under the unfair one it failed. Across all targets and
   seeds its calibration still grades it weak: a structure that maps onto a held-out label exists,
   but a cluster-per-label recovery of it is a modest signal on this table.

2. **Tests that ignored their settings.** Strategies 01 and 02 fixed 1,500 genes and the combined
   features and truncated every grid to two values, so every setting the sweep asked about tested
   the same map. They now honour genes per map (to 4,000), feature sets and up to three grid values.

3. **Held-out derived and annotation layers.** Tuning strategy 16 found skill 1.00 on held-out
   `structural_hole` and held-out `domain`. The first is computed from co-expression and co-fitness,
   which stay as evidence; the second joins every gene that shares a label, so it is made of cliques
   and a hidden edge is closed by its neighbours. Both are now refused as targets, by name and reason.

4. **Correlation layers completing themselves.** With those gone, strategy 16 still read AUROC 0.99
   on co-fitness and co-expression -- against random non-pairs, and then, after switching to
   degree-matched non-pairs, still 0.99. A layer drawn as the top correlations among some columns is
   determined by its own visible edges: two genes that correlate at 0.95 with the same partners
   correlate with each other. That is arithmetic, not a missed contact, so strategy 16 now completes
   contact layers only (crosslinks, structure), and strategy 34 answers the question for correlation
   layers from the OTHER evidence. On held-out crosslinks strategy 16 reads 0.78 against
   degree-matched non-pairs, where it used to report 0.97 against random ones.

5. **A layer's own measurements leaking back.** Hiding a layer while its source columns stay in a
   measurement similarity hands the layer back by a second route. `search.layer_measurement_sources`
   now removes the source columns and every column of their experiments whenever a layer is held
   out, in both `graphspace` and strategy 16. Measured honestly, this was not what drove the one
   remaining high number: held-out co-translation still reads 0.958 with all 66 ribosome-profiling
   columns gone and a configuration-model null at 0.53. Co-translation is simply the layer most
   predictable from independent evidence -- genes made together are expressed and needed together.

Two smaller things: five configurations of strategy 05 crashed on a chi-squared table whose row and
column filters could each empty the other's margin (they now run to a fixed point), and 1,000 grid
points asked the graph strategies to hold out layers they rightly refuse, so their grids are now
built from the layers each strategy accepts.

## The machine

The first sweep ran twelve pooled workers that kept every map they drew. They reached 12 GB each,
106 GB in total, and at 21:31 on 2026-09-25 the kernel's OOM killer fired and systemd-oomd killed the
editor with every session inside it. The sweep now schedules against measured memory: each chunk of
work runs in a fresh process that reports its peak, the peak becomes that strategy's estimate, and
a chunk starts only when it fits under a job budget (55 GB) and a system limit (75 GB). It ran as a
capped service (`systemd-run --user -p MemoryMax=62G`) so that even a wrong estimate is stopped by
the kernel inside the job's own cgroup. A machine-wide guard outside the repository warns at 80 GB,
freezes the largest job at 100 GB and kills at 110 GB.

## Where things are

* `scripts/calibrate_strategies.py` -- the sweep (`--workers`, `--job-gb`, `--only`), `--summarize`,
  and `--publish`, which writes the shipped JSON, the README block and `docs/calibration.md`.
* `starplast/calibration.py` -- skill, intervals, the seed split, grades, and reading it back.
* `starplast/strategy_panel.py` -- grades in the list, the calibration paragraph in each Guide, and
  **Use tuned settings**.
* `tests/test_calibration.py` -- 21 tests, including the seed split and the clustered interval.
* `results/calibration_2026-09-26/` -- every run (`runs.jsonl`), per configuration (`summary.csv`),
  the best per strategy (`best.csv`), and per parameter value (`sensitivity.csv`).
