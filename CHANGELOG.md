# Changelog

## 0.46.0

- **Name each strategy's method.** Every strategy's name now ends with the method it runs in brackets, for example "Hold out a category and search for a map that finds it (UMAP + HDBSCAN)" or "Diffuse a label across one measured network (random walk with restart)". The label is taken from the code each strategy calls, not from its prose. The name appears in the Strategies tab, its Guide, the README calibration table, the strategy and calibration pages and the tutorials. The tab's filter matches it, so typing "HDBSCAN" or "logistic" lists every strategy that uses that method.
- Add `Strategy.method` and `Strategy.name`. A strategy that does not name its method is refused at registration. `strategies.overview()` now has `name` and `method` columns in place of `title`.
- **Give every strategy a scorecard.** Each strategy declares one of six tasks: label calls, ranking, set retrieval, cluster recovery, values or replication. Its self-test reports that task's standard metrics, in a fixed order, on the same hidden genes as its verdict. Label calls report accuracy, coverage, precision of calls, macro precision, macro recall (balanced accuracy), macro F1, weighted F1, Cohen's kappa, MCC, and macro AUROC and AUPRC. Rankings report AUROC, AUPRC and its lift over prevalence, partial AUROC, R-precision, precision and enrichment at the top 1%, recall at the top 10%, best F1 and nDCG. Values report Spearman, Pearson, Kendall, R-squared, normalised RMSE, MAE and top- and bottom-decile recall. Cluster recovery reports weighted F1, precision and recall, ARI, NMI, homogeneity, completeness and the unclustered share. Every metric is checked against scikit-learn and defined once in `starplast.scorecard`, with its range, its chance level and how to read it. The Strategies tab shows the card after every test, with each metric explained on hover. Calibration records every card, and the README and [docs/scorecards.md](docs/scorecards.md) give each metric's mean and 95% interval for every strategy.
- **Explain every technique.** Each strategy lists the techniques its method is built from, and `starplast.techniques` explains each: what it does and why a strategy uses it. The Guide shows them under the strategy's name. `strategies.metrics()`, `strategies.techniques()`, `Strategy.techniques_table()`, `Strategy.scorecard_table()`, `TestResult.card()` and `TestResult.skill` expose the same from Python.
- **Add five strategies (35-39) in a ninth family, Advanced models.**
  - **Conformal label calls** call a gene only where its conformal prediction set holds one label, with a stated error rate that the self-test checks on hidden genes. They also list the genes whose data leaves two labels possible.
  - **Graph convolution** lets a logistic regression learn from each gene's network neighbourhood as well as its own measurements, and reports how much weight it puts on each.
  - **Random forest** ranks measurements by permutation importance on held-out orthogroups.
  - **Stacking** learns out of fold how far to trust measurement neighbours, a linear model and the networks, for each label and class.
  - **Conformal values** wraps a predicted measurement in an interval that holds for a stated share of new genes, and lists the measured genes outside theirs.
- Label-calling strategies now pass their per-class scores to the scorecard (k-nearest neighbours, network vote, random walk, logistic regression, weighted vote, triangulation, map neighbours, cluster guilt), so macro AUROC and AUPRC are measured rather than missing.
- Document the strategies API in [docs/API.md](docs/API.md): overview, run, test, card, calibration, tuned settings, contexts, and scoring predictions made outside Starplast.
- Add `starplast.organisms`, one declaration per species space: code, reference, id pattern, tables, partner, life stages, calibration targets and record links. It declares *T. gondii* and *P. falciparum* today, and tests hold it to every literal it will replace. This is the first step towards separate host and vector spaces ([instruction 53](instructions/open/53_organism_spaces.md)).
- `import starplast` now reaches the main modules as attributes (`starplast.strategies`, `starplast.scorecard`, `starplast.techniques`, `starplast.calibration`, ...), each loaded on first use. Importing the strategies no longer imports Qt: the `Stopped` exception moved to `starplast.stopping`, and `jobs.Stopped` is the same class.
- **Search Starplast from beside the Help menu**, as in spaCR. The box to the right of **Help** (**Ctrl+Shift+H**, or **Help → Search Starplast…**) finds every menu command, panel and tab, strategy (by number, name, method, family or question), setting in the analysis panel and Preferences, information slot of both organisms, registered dataset, heading of the guide and tutorial section, and goes to the exact place: the command runs, the panel is raised on its tab, the strategy is selected in the Strategies tab, the setting is scrolled to and outlined, the slot is selected in the slot tree, the guide opens at the heading. The index is built from the menus, panels and registries themselves (`starplast.help_index`, testable without a display); the field is `starplast.help_search`. **Help → Keyboard shortcuts** lists every key the menus bind, and every menu command now carries a tooltip and a status tip.
- **Format the menus like spaCR's, and draw everything that floats as translucent black glass.** The menu bar is flat and its words light up when pointed at; menus, tooltips and drop-down lists are rounded panes of translucent black (white on a light theme) with rounded item pills, in spaCR's Open Sans. Preferences, the guided workflows, the slot tree, the import and export dialogs and the explanations are rounded glass cards without a title bar: drag the background to move them, an edge to resize, **✕** or **Esc** to close. The explanations and About no longer freeze the map behind a modal box. Without a compositor (X11) the same panes are opaque near-black with their corners cut; `STARPLAST_TRANSLUCENT=0/1` overrides the check (`starplast.glass`). Text typed into fields is now light on the dark field grey in the light themes, where it could not be read.

## 0.45.0

- **Calibrate every strategy.** Each strategy's self-test was run over a grid of its settings, several held-out known labels and five seeds: 5,940 tests in all. The README table gives every strategy's skill at its defaults and at its best setting, with 95% intervals and pass rates. Skill puts every metric on one scale, from 0 for the same procedure on shuffled data to 1 for perfect. The best setting is chosen on seeds 1-3 and reported on seeds 4-5, so the table does not report the luckiest of many settings. Intervals resample held-out targets, then runs within each target, because runs on one label are not independent. On *T. gondii* 26 strategies are reliable, 7 weak and 1 untestable. On *P. falciparum* 20 are reliable, 8 weak, 3 work only when tuned, 2 are untestable and 1 has no skill. [Every number, per target and per setting](docs/calibration.md).
- Show the calibration in the Strategies tab: each strategy's grade in the list, a calibration paragraph in its Guide, and a **Use tuned settings** button. `strategies.overview()`, `strategies.calibration(key)` and `strategies.tuned(key)` expose the same from Python.
- **Correct three self-tests that flattered themselves.** Strategies 01, 04 and 15 were judged against a second search on shuffled labels. That search picks the largest cluster for every label, which scores an F1 near twice the label's share on size alone, and it beat the real labels on four targets of eight. The null now permutes the hidden genes' labels over the same chosen clusters, so a large cluster earns nothing. Strategy 01, the founding question, passes its *T. gondii* self-test under the fair null (0.132 against 0.039), though across all targets and seeds its calibration grades it weak.
- **Stop strategy 16 reporting AUROC 0.99 for arithmetic.** Held-out edges are now scored against degree-matched non-pairs, not random ones. The layer's own source measurements leave its similarity feature. Layers that cannot honestly be held out are refused as targets: derived layers, annotation cliques, literature, and correlation layers whose visible edges determine their hidden ones. On held-out crosslinks it now reads 0.78.
- Make strategies 01 and 02 walk the settings they are given (genes per map, feature sets, grids). Calibrating over settings the test ignored measured nothing.
- Add strategies 33 and 34 and `starplast.graphspace`: one integrated neighbour space from every permitted layer, with per-edge evidence, and a model trained to rank the edges the networks miss. Every number is reported against random, degree-matched and configuration-model non-pairs, with the gap between the first two as its own quantity. A learned embedding was tried and did not beat the interpretable baseline, so the baseline ships. [The measured comparison](docs/graphspace.md).
- **Correct a citation.** The four in vivo CRISPR columns are Giuliano et al. 2024 (PMID 38977907), not the 2019 platform paper: they match that supplement to 5e-8. Re-reading it adds heart and brain fitness and 65 genes.
- Add 21 verified deposits through `starplast.deposits`, each re-derived by an executed notebook (`notebooks/derive_deposits_2026_09.ipynb`). *T. gondii* goes from 402 to 438 columns, *P. falciparum* from 123 to 146, and the host table from 34,630 to 36,579 proteins. Slot coverage rises from 156/204 to 169/215. The audit, with every refusal and its evidence, is [instruction 50](instructions/done/50_data_audit_2026_09.md).
- Measure the slot grouping for leakage instead of asserting it (`scripts/leakage_audit.py`, `starplast.leakage`). The audit found and closed two leaks: the *P. berghei* transfers predicted one another at 0.58, and withdrawing a nutrient predicted fibroblast fitness at 0.76. Afterwards there are 0 closure gaps and 0 residual leaks, and the closure's threshold sits at the 99th percentile of what unrelated columns reach.
- Add five tutorials, each as a GUI walkthrough and a notebook, and a complete guide to every feature ([docs/tutorial](docs/tutorial/index.html)). All are built from the running application by `scripts/build_tutorials.py`.
- Schedule the calibration sweep against measured memory. Each chunk of work runs in a fresh process that reports its peak, and a chunk starts only when it fits under both a job and a system budget. The pooled version grew to 106 GB and took the machine down.
- Fix a chi-squared crash in strategy 05 when filtering a contingency table left a zero margin.

## 0.44.0

- Add a **Strategies** tab to the right of Evidence and Analysis: 32 named ways of using the combined data for inference, grouped into eight families, each with a tooltip, an explanation, a walkthrough, settings that explain themselves, and Run / Test / Stop.
- Give every strategy a self-test that hides known information -- labels, gene-set members, edges or values -- asks for it back, and compares the answer with the same procedure on shuffled labels, random sets or permuted identities. A strategy passes only above its null's 95th percentile by a stated margin.
- Include the two founding questions as strategies 01 and 02: hold out a category and search maps for the structure that recovers it, and find the map where a gene list forms one cluster with high precision and recall.
- Measure every strategy on the shipped tables and show the verdict beside it: 26 pass, 5 fail and 1 is inconclusive on *T. gondii*; 24, 5 and 3 on *P. falciparum*. The [strategy catalogue](docs/strategies.md) and `results/strategies_2026-09-25/` record every number, failures included.
- Add `starplast.strategies` (context, leakage guard, five holdout test patterns, a planted organism for testing) and `starplast.strategy_catalog`; strategies run headless as well as from the tab.
- Add a `join="louvain"` option to `methods.multiplex_communities`: joining agreed pairs by connected components chains every gene into one community when layers agree about different genes.
- Group evidence into families by the slot catalogue's axis, so knockout screens named for a second background or a genetic interaction count as fitness.

## 0.43.0

- Add a 40-slide introduction and practical guide, presented like spaCR's deck: README cover, linked GitHub slide pages, web viewer, PDF and editable PowerPoint.
- Add guided Explore a gene, Predict a trait and Compare a screen workflows with settings help and full run exports.
- Evaluate feature, linear, boosted, PCA, UMAP, masked-factor and weighted-network predictions with grouped folds, target exclusions, calibration and unsupported-call abstention.
- Add classification, regression and multi-label APIs; preserve unknown outcomes and report fixed-class metrics.
- Include 320 frozen ESM-2 sequence features for 8,064 proteins and nine local AF3 summaries for 1,210 exactly mapped genes. Add model-indexing and sequence-encoding commands.
- Publish 29 reproducible localization and abundance comparisons, source ablations, shuffled controls and uncertainty estimates in the [benchmark report](docs/benchmark-0.43.md).
- Use one balanced embedding builder for packaged and interactive maps; record actual algorithms, input hashes and ordered gene identities.
- Add source-level observation records, reviewed literature assertions, candidate explanations and an explicitly heuristic budget API.
- Restore the Plasmodium evidence ledger and prevent fixture builds from overwriting packaged data.
- Share OpenGL resources between organism windows so switching species retains valid shader programs.
- Expand CI to the complete non-slow test suite, repair documentation and display regressions, and declare numerical thread-control dependencies.
- Repair download badges, retain the linked 130-source catalogue and publish forty additional monochrome logo proposals; keep the approved logo active.

## 0.42.1

- Use the GitHub README as the PyPI project description.
- Use full URLs for README artwork, the rotating gene map, and documentation links.

## 0.42.0

- Simplify the map to individual genes; remove Galaxy and Orthogroup summary modes.
- Adopt the Toxoplasma constellation logo, SVG wordmarks, and window icon.
- Rewrite the README and add a user guide and Python API guide.
- Publish generated API documentation through GitHub Pages.
- Add settings help and document application callbacks.
- Fix importing screen tables whose identifier column is already named `gene_id`.
- Fix Plasmodium gene selection and link its evidence panel to PlasmoDB.
- Publish the application and bundled data as one PyPI project, `starplast`.
- Make CUDA dependencies optional through `starplast[gpu]`.
- Include SVG artwork and compressed sequence tables in wheels; exclude saved embeddings.
- Add synchronized version bumps, build checks, and automatic PyPI publishing.
- Develop on `nightly`; publish version increases on `main` to PyPI and GitHub Releases.
- Link all 128 registered datasets and computed layers to their sources.
- Show a 1440×1080, 30 fps gene-map rotation with selected-gene lighting in the README.

## 0.41.0

Existing development baseline before the packaging and documentation overhaul.
Earlier changes are recorded in [HANDOFF.md](HANDOFF.md) and the Git history.
