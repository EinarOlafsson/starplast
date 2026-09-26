# The September 2026 data audit: are the information slots filled with the best data there is?

**Status: complete (v0.45.0, 2026-09-26).**

The question was whether every information slot holds the best available data, whether newer
datasets exist, whether the slots are grouped so that holding one out really holds it out, and
whether anything shipped is cited to the wrong paper. The answer to the last one was no, and that
is the most important thing in this record.

## What was searched, and how

Three sweeps, each recording the queries it used so the search can be repeated rather than trusted:

* **ToxoDB / PlasmoDB release 71** — attribute reports re-fetched, which is how the in-vivo CRISPR
  tracks were found to have been REMOVED from the site since the original ingest.
* **NCBI GEO and PRIDE** — index sweeps for both organisms.
* **bioRxiv via Europe PMC** — 2,572 preprints posted 2025-01-01 to 2026-09-25, cross-checked
  against Crossref (member 54368, 835 DOIs; only 3 missing from Europe PMC, all posted in the
  final two days and none of them genome-scale). Worth writing down for the next sweep:
  `api.biorxiv.org/details/...` now returns HTTP 200 with a zero-length body for every interval, so
  the details endpoint is dead and Europe PMC plus Crossref is the working route. bioRxiv also
  changed its DOI prefix from `10.1101` to `10.64898` around 2025-12, and both appear in this
  window.

Every accession was tested for public availability before a dataset was accepted, and every
downloaded file carries its URL and sha256 beside it.

## The citation that was wrong

`crispr_invivo_composite` — the four in-vivo fitness columns, 7,395 genes, in the map since early
on — was cited to PMID 31481656, the 2019 platform paper. It is not that study's data. The columns
are Giuliano et al. 2024 (PMID 38977907, Nat Microbiol 9:2323-2343), and this is not an inference:
against that paper's Supplementary Data 5 the shipped values match at Spearman 0.999 and a maximum
absolute difference of 5e-8, and each column's best match among all 32 numeric columns in the sheet
is its own tissue. The 2019 paper's own supplements are 36- to 620-gene targeted libraries and
cannot be the source.

Two consequences beyond the citation. The same sheet carries **heart and brain**, which nothing in
the map had, and re-reading the supplement rather than the ToxoDB tracks **gained 65 genes** whose
GT1 accessions ToxoDB splits. And a trap: three accessions can resolve to one ME49 gene
(`TGGT1_224540` and `TGGT1_224540A/B`), and averaging all three moved one gene's score from -3.99
to +0.20. The halves are different gene models, not replicates, so the accession that IS the gene
wins and the halves are used only where nothing else reaches that gene.

## What was added

Twenty-one deposits derive into the tables through `starplast/deposits.py`, each with the
computation between the file and the column written down, and each re-derived by
`notebooks/derive_deposits_2026_09.ipynb`. Toxoplasma went from 402 to 438 columns, Plasmodium from
123 to 146, the host table from 34,630 to 36,579 proteins, and slot coverage from 156/204 to 169/215
genes-unit slots.

The additions worth naming:

* **In vivo heart and brain fitness**, and the citation fix above.
* **Serum restriction** (PMID 41407671): fitness in 10% and 1% serum. Fitness in either serum IS
  fibroblast fitness again (rho 0.79-0.85 with the Sidik screen); the DIFFERENCE between them is
  orthogonal to it (rho 0.06) and is the new axis, lipid dependence.
* **Carbon-source withdrawal** (bioRxiv 2025): glucose-only and glutamine-only arms. GDH1 ranks 1
  of 8,155 and PEPCK 2 among genes needed without glucose, which is the paper's result.
* **Translation efficiency, GSE302107**: the most reproducible TE in the organism (replicates rho
  0.97, against 0.87 and 0.67 for the two already shipped).
* **Genome-wide mRNA decay** (GSE329845): 5,997 genes where the shipped column was the 412-gene
  unstable tail of another study. They do not correlate (rho -0.03), which that tail's selection
  explains, so they are kept as separate columns in different units.
* **5' UTR architecture, bradyzoite subtypes, organelle-surface proximity, iron depletion** on the
  Toxoplasma side; **melting temperature, m6A per transcript, chromatin proxiomes, febrile
  phosphorylation, latency, gametocyte proteome, antimalarial target engagement, and male and
  female fertility transferred from P. berghei** on the Plasmodium side.
* **Host**: the fibroblast and hepatocyte responses to infection, a baseline macrophage, and the
  first host-gene CRISPR screen in the map (rhoptry discharge in K562).

## What was refused, and why refusals matter

* **PARFA-Seq (E-MTAB-14557)** — released and well annotated, but no processed data exists: only 36
  FASTQ files. Per-gene values would need a dual-genome alignment build.
* **The 5'UTR reporter assay** — 30,235 scored sequences, but they are variants of twelve endogenous
  UTRs; exact matching recovers 11 genes. It describes sequences, not genes.
* **A translation-efficiency shift under amino-acid starvation (GSE226628)** — replicate agreement
  r -0.36 to +0.27 and zero genes at FDR 0.05, diagnosed to per-library polysome-enrichment
  efficiency varying 0.35-1.50. Shipping it would describe the sucrose gradients. The log fold
  changes from the same deposit, which reproduce the paper at Pearson 0.93, are fine.
* **A thermal melting temperature from the antimalarial screen** — that design holds temperature
  fixed and varies dose, so it cannot give one. Only hit counts ship, and the melting points come
  from a different experiment.
* **The sexual-stage re-ingest** — the journal supplement is byte-identical to the preprint's, and
  the shipped values match the source at Spearman 1.0000, so only the citation changed (now PMID
  42020723). Two further contrasts in it are withheld because their reference group is unstated and
  they anti-correlate at -0.92 with the shipped column: shipping them would mean guessing a sign.

## Numbers a reader should not take on trust

Three published numbers were reproduced and one was found to mean something other than it appears:

* The febrile phosphoproteome's "143 up, 53 down" are counts of ROWS, and a row is a site at one
  phosphorylation multiplicity. Collapsed to unique sites, each taking its strongest change, it is
  **128 and 51**, which is what the shipped columns count.
* The chromatin proxiome's euchromatin counts reproduce exactly (H3K27ac 99, H3K4me3 48). The HP1
  column gives 91, which is NOT the paper's 61: that is its combined heterochromatin group over two
  HP1 baits and the column is one of them.
* The gametocyte translatome is 705 proteins only if "newly made" means found with the label and in
  no control; taking the whole sheet gives 1,179.

## Leakage: the grouping was tested, and two leaks were found

`scripts/leakage_audit.py` measures association between every pair of measured columns, calibrates
what unrelated columns look like, and asks whether any permitted slot predicts a held-out target
like a copy of it. Run on both tables with the new columns
(`results/leakage_audit_2026-09-26/`):

* Between columns from different axes, rank association reaches 0.801 at the 99th percentile and
  0.861 at the 99.9th on Toxoplasma; the closure's 0.8 threshold sits at the 99.0th percentile of
  that distribution, so it is neither arbitrary nor slack.
* **Closure gaps: 0. Residual leaks after the fixes: 0.**

Two real leaks were found and closed:

1. **The berghei transfers were separable.** Holding out the blood-stage phenotype transferred from
   P. berghei and predicting it from the transmission transfer reached 0.58 at nine robust standard
   deviations. They are one knockout library carried across by orthology, so they are now one
   same-quantity family and hold out together.
2. **Withdrawing a nutrient predicted fibroblast fitness at 0.76**, the highest cross-slot
   predictability in the table. Same library, same cells, same passage regime with one nutrient
   removed: an essential gene is essential in both, so the arms now hold out with the fibroblast
   screen. The DIFFERENCE between the arms stays separate, because it is orthogonal to bulk fitness
   and is the question that screen was built to ask.

## Where things are

* `starplast/deposits.py` — the derivations, one function per deposit, with its checks named.
* `scripts/derive_deposits.py` → `notebooks/derive_deposits_2026_09.ipynb` — executed, so the
  notebook in the tree is the run that produced the shipped numbers.
* `scripts/add_deposits.py` — merges derived tables into the caches, and REFUSES any merge that
  would lose a value or write an empty column.
* `scripts/fit_meltome.py` — fits the Plasmodium melting points, which no source publishes.
* `tests/test_deposits.py` — 34 tests: the statistics on planted data, the merge's refusals, and
  every claim above on the shipped tables.
* Raw deposits under `datasets/<level>/<kind>/<accession>/`, each with `URLS.txt` and
  `SHA256SUMS.txt`.
