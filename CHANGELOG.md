# Changelog

## Unreleased

- Add `instructions/open/59_biological_questions.md`: one hundred biological questions about
  *Toxoplasma gondii* and *Plasmodium falciparum* that the shipped data can answer, each with its
  organism, entry point (Start here, Strategies, Analysis, maps, star map), strategy, exact settings,
  what a good answer looks like, and a literature anchor resolved through Europe PMC and NCBI (no
  identifier typed from memory; every query recorded). `scripts/run_biological_questions.py` runs the
  executed subset and writes `results/questions_2026_09_30/` and
  `notebooks/biological_questions_2026_09_30.ipynb`; `docs/questions.md` lists the questions that
  yielded a result, each reproducible as a panel path or a one-line API call. No strategy, calibration
  number, shipped table or UI behaviour changed.

## 0.49.0

- Add a **start here** tab beside Strategies, which is unchanged: a guided path from what a user
  has to the strategies worth running. One question per screen -- what do you have (a gene, a gene
  list, your own measurement, a label, nothing yet), which organism, the gene / gene set / label /
  measurement itself, and the goal (more genes like mine, predict it, explain it, partners,
  compare two conditions, is it learnable) -- with a clickable trail of the answers to go back to
  any of them. Label and measurement choices show their coverage, so an almost-empty column is not
  picked blind. It ends at three to six recommended strategies, each shown with the Strategies
  tab's own card summary (name, method, grade, the four `scorecard.HEADLINE` bars), one line on why
  it is recommended, the settings the answers decided, and **Run it here** / **Open in Strategies**;
  the map gallery and the star map are offered where they answer the goal better. Ranking prefers
  what the goal is for, then a matching scorecard task, then the calibration grade *for that
  organism*, and says plainly when nothing better than weak exists. The question tree and the
  ranking are Qt-free in `starplast/guided.py`, the view is `starplast/guided_panel.py`, and both
  are reachable from the help search.
- Add a **maps** tab beside Analysis with pregenerated 3D UMAPs for each organism. There are three
  maps of all measurements (n_neighbors 10, 25 and 60), one of all measurements except
  localization, and one per evidence family (12 maps for *T. gondii*, 14 for *P. falciparum*).
  Maps are built from measurements only and clustered with HDBSCAN. Clicking a map shows it in
  the central view with its clusters. A sortable table scores every categorical label on every
  map, read by label or by map. It reports categories → clusters (size-weighted best-cluster F1,
  precision and recall) and the best single category, scored on the F1 of Wilson 95% lower
  bounds so a 2-gene category cannot score 1. Both scores come with a shuffled-label chance level
  and skill, plus coverage and a circularity flag from the leakage closure. **Color by label** is
  one click. **Score map on screen** runs the same scoring on a map built in the app.
  `starplast.umap_gallery`, `scripts/build_umap_gallery.py`,
  `notebooks/umap_gallery_2026_09_29.ipynb`; about 2 MB of data.
- Add the **star map** tab (`starplast/star_map.py`): a navigable network centred on one gene,
  with its linked genes on rings around it. Edges are coloured by source (each measured layer,
  each strategy, your own runs), solid when measured and dashed when inferred, and as wide as
  their strength. Hover shows the gene or the link's layer, strategy, run, setting and score.
  Click to re-centre, Back to return, depth 1-2 hops, per-source toggles and link counts, and
  selection follows the rest of the application and can be sent to the 3D map.
- Add `starplast/star_edges.py`, one store of gene-gene links with provenance. Measured links are
  read from each space's graph. Strategy links are bounded: modules link each member to its 3
  nearest co-members, set expansions link each hit to its 3 nearest seeds, neighbour and partner
  calls link to the genes that made them, and no run adds more than 20,000 links. Ships
  `data/star_edges.parquet`, the links of the default-setting runs of 13 strategies on both
  shipped tables, built by `scripts/build_star_edges.py` and recorded in
  `notebooks/star_edges_2026_09_29.ipynb`. Strategy runs made in the application are kept under
  the user cache (`star_edges/`) and appear in the map at once.
- Add the Lourido lab's parasite-density CRISPR screen (Giuliano et al., Cell 2026, PMID 42580337)
  as the Toxoplasma slot "fitness · parasite density". It has five columns: fitness at low (MOI
  0.3) and high (MOI 3) density, the authors' high-versus-low contrast (negative means needed at
  high density), its Bonferroni p, and the paper's 31 density-inhibited mutants. The arms cover
  7,461 genes and the contrast 5,291. Table S1 reproduces the paper's numbers: r = 0.995 between
  the arms, 31 hits, and the 12 high-confidence hits from the stated rule. The arms are held out
  with fibroblast fitness. The contrast is its own leakage family, and the leakage audit finds 0
  gaps and 0 residual leaks (instruction 56).

## 0.48.0

- The Strategies tab opens each strategy on a card instead of prose. The card shows the name,
  method, task and calibration grade, one line saying what the strategy answers, and four painted
  bars in the same places for every strategy: **Better than chance** (skill), **Reach** (coverage,
  or its task's stated analogue) and two task metrics in plain words (for example "Right calls"
  and "Fair across classes" for label calls). Each bar shows its 95% interval from the shipped
  calibration, a marker at chance, the technical name in small type, and a plain sentence on hover.
  **Run**, **Test** and **Details ▸** sit under the bars. The Guide, Settings and Results tabs keep
  all their controls and move behind Details. A test started from the card shows its result on the
  card in the same four bars.
- Add `scorecard.HEADLINE`, `REACH`, `headline()` and `headline_bars()`, which define the four
  headline bars for each task once, including where their chance levels come from.
- Add an "About this test" box for all 39 strategies (`starplast/strategy_explainers.py`). It has
  four short fields: what the strategy does, how it is evaluated, what failure looks like and why,
  and what success looks like and why. Each field was written from the strategy's own explanation
  and test description. The box is collapsed to one line per field and opens on click.
- Add worked examples: one real failure and one real success per strategy and organism
  (`starplast/data/strategy_examples.json`). They are chosen from the 7,640 calibration runs by
  `scripts/build_strategy_examples.py` and recorded in
  `notebooks/strategy_examples_2026_09_28.ipynb`. Each failure explains its failure mode from its
  own numbers. Each success is re-run once and lists its top five calls for genes without a known
  label. Where a strategy never failed on real data, the failure shown is its self-test on the
  noise table, and the card labels it as such.
- Add `docs/strategy_cards.md`, which renders the same cards, explainers and examples for the
  docs. Help search now opens a strategy on its card, and the Guide, Settings and Results entries
  open Details. Tutorial captions now name the moved controls; the tutorials have not been rebuilt.

## 0.47.0

- Add `NEXT_SESSION.md`, the handoff for resuming work: current state, working rules, the prioritised plan for the per-organism spaces and the traps already met. `HANDOFF.md` points to it.

- Rebuild the Pf display from its own feature slots: 118 resolved features now place all 5,720
  genes, instead of the 20-feature legacy fallback. Display recipes are separate from the shipped
  calibration's inference recipe. Seventeen strategy groupings and both compatibility input matrices
  were checked unchanged; edges, node tables, Tg coordinates and calibration values were preserved.

- Add versioned data packs with SHA256 manifests, per-column license declarations, checked HTTPS
  downloads and atomic activation. A shared space builder refuses lost identifiers, columns or
  measured cells and validates graph order and indices. Pack paths resolve from the user's cache;
  organism builders and published download URLs remain pending.

- Publishing a calibration sweep now merges by organism and strategy, preserving unswept results
  and their original provenance. Writes are atomic; malformed existing files are refused. Generated
  calibration pages include every published space and distinguish the latest sweep from older results.

- Split host protein references into separate human and mouse tables, retaining all 36,579
  identifiers, names and measurements. Resolve the 54 identity-only rows from local UniProt and
  bridge records; keep the migration and source checksums in an executed notebook and manifest.
- Require an explicit organism for every dataset and distinguish human/mouse deposits. Provenance
  queries for one species no longer fall back to another. Host loaders, merges, slot coverage and
  the slot audit now read separate species tables; merges refuse lost values even when new rows
  keep the total coverage unchanged.

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
- Add `starplast.organisms`, one declaration per species space: code, reference, id pattern, tables, partner, life stages, calibration targets and record links. It declares *T. gondii* and *P. falciparum* today, and tests hold it to every literal it will replace. The slot tables, the window's species list and record links, the strategy contexts, the calibration targets and the slot generator's life stages now read it. This is the first step towards separate host and vector spaces ([instruction 53](instructions/open/53_organism_spaces.md)).
- `import starplast` now reaches the main modules as attributes (`starplast.strategies`, `starplast.scorecard`, `starplast.techniques`, `starplast.calibration`, ...), each loaded on first use. Importing the strategies no longer imports Qt: the `Stopped` exception moved to `starplast.stopping`, and `jobs.Stopped` is the same class.
- **Search Starplast from beside the Help menu**, as in spaCR. The box to the right of **Help** (**Ctrl+Shift+H**, or **Help → Search Starplast…**) finds every menu command, panel and tab, strategy (by number, name, method, family or question), setting in the analysis panel and Preferences, information slot of both organisms, registered dataset, heading of the guide and tutorial section, and goes to the exact place: the command runs, the panel is raised on its tab, the strategy is selected in the Strategies tab, the setting is scrolled to and outlined, the slot is selected in the slot tree, the guide opens at the heading. The index is built from the menus, panels and registries themselves (`starplast.help_index`, testable without a display); the field is `starplast.help_search`. **Help → Keyboard shortcuts** lists every key the menus bind, and every menu command now carries a tooltip and a status tip.
- Keep the glass style installed once per application, keyed on the Qt application rather than on a Python wrapper that can be collected and remade. The test suite now shares one application for the whole session and closes each module's windows. It had been destroying and rebuilding the application between modules, which left the second one styled with freed memory and segfaulted full runs.
- **Format the menus like spaCR's, and draw everything that floats as translucent black glass.** The menu bar is flat and its words light up when pointed at; menus, tooltips and drop-down lists are rounded panes of translucent black (white on a light theme) with rounded item pills, in spaCR's Open Sans. Preferences, the guided workflows, the slot tree, the import and export dialogs and the explanations are rounded glass cards without a title bar: drag the background to move them, an edge to resize, **✕** or **Esc** to close. The explanations and About no longer freeze the map behind a modal box. Without a compositor (X11) the same panes are opaque near-black with their corners cut; `STARPLAST_TRANSLUCENT=0/1` overrides the check (`starplast.glass`). Text typed into fields is now light on the dark field grey in the light themes, where it could not be read.
- **Add six verified Plasmodium deposits, and the questions they answer that nothing could ask before.** The *P. falciparum* table goes from 146 to 168 columns and slot coverage from 169/215 to 180/221. It now has a subcellular localization for the first time -- 1,646 of 3,000 schizont proteins in 24 niches from hyperLOPIT, with RAP1, MAHRP1, ACP, EXP2 and GAPDH each where they belong -- plus protein abundance in the asexual blood stage and under Hsp90 inhibition, mRNA synthesis and decay rates from hourly 4-thiouracil labelling, RNA dependence of complexes (898 of 3,671 proteins), the sexually committed proteome, field and between-species variation, and an in vitro evolution resistome whose classified top is PfATP4, PfMDR1, the prodrug esterase, cytochrome b, CARL and PI4K. Five empty slots are answered and six new ones exist because the catalogue could not express these measurements. Every headline number is reproduced from the raw deposit in `notebooks/derive_deposits_2026_09.ipynb`; five further candidates were refused, one of them because the paper's own count of 29 drug-sensitive mutants came back as 40. The audit, with every refusal and its evidence, is [instruction 52](instructions/done/52_data_audit_2026_09_b.md).
- Carry LABEL columns through the deposit merge: a compartment has no mean, so where two accessions resolve to one gene the label survives only if they agree. Before this, `parasite_columns` dropped any non-numeric column in silence.
- Declare seven more same-quantity families for the leakage closure and re-run the audit with eight new held-out targets: 0 closure gaps and 0 residual leaks on both tables.

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
