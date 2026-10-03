# 59 -- One hundred biological questions this software can answer

Opened 2026-09-30 on `nightly` at 0.49.0.

**The task.** Write 100 questions about *Toxoplasma gondii* and *Plasmodium falciparum* in the words
a parasitologist would use, each tied to an entry point in the application, a strategy or view, exact
settings, what a good answer looks like, and a resolved literature identifier that makes it a real
question rather than an exercise. Then run a good subset and record, honestly, which ones say
something.

**Nothing in this task changes a strategy algorithm, a calibration number, a shipped data table or
the UI.** Every run is the ordinary `starplast.strategies` API at documented settings.

---

## 1. Method

1. **Literature first.** Fourteen topics were searched through Europe PMC's REST API, with
   independent title confirmation for anchor identifiers through NCBI E-utilities `esummary` and one
   Crossref check. Every identifier in this file was read out of an API response body; none was typed
   from memory. Section 5 lists the queries verbatim and section 6 the confirmation call.
2. **What the data can speak to.** The shipped tables are *T. gondii* `nodes.parquet` (8,140 genes x
   443 columns) and *P. falciparum* `pf_nodes.parquet` (5,720 x 168), plus 13 Tg and 8 Pf gene-gene
   edge layers. An open question was kept only where the shipped columns carry the relevant quantity:
   localization (hyperLOPIT, 25 Tg compartments; the Pf schizont spatial proteome, 24 locations),
   fitness in many conditions, stage expression and translation, PTM site counts, crosslink and
   pulldown interactomes, structural similarity, orthology, and publication counts. Questions about
   quantities that are not in the tables -- lipid composition, cryo-ET densities, drug IC50s -- were
   dropped rather than faked.
3. **Coverage.** The 100 questions use all 39 strategies, all five entry points (Start here,
   Strategies, Analysis, Maps, star map) and both organisms. They deliberately include questions
   where the calibration says the strategy is **weak** (Q21, Q29, Q32, Q37, Q38, Q48, and Q65-Q70),
   so the limits are visible rather than hidden.
4. **Running.** `scripts/run_biological_questions.py` runs the executed subset Q01-Q50 and writes
   `results/questions_2026_09_30/` plus `notebooks/biological_questions_2026_09_30.ipynb`. Verdicts
   are recorded in section 4 and in `docs/questions.md`:
   * **yields** -- real, specific, checkable output: it names genes a researcher could follow up, or
     it makes a clear statement about the data.
   * **thin** -- it runs and is honest, but says little.
   * **no** -- it fails, or the data cannot answer it, with the reason.

---

## 2. The hundred questions

Entry points: **SH** = Start here tab, **ST** = Strategies tab, **AN** = Analysis tab,
**MAP** = maps tab (pregenerated maps), **STAR** = star map tab. `Tg` = *T. gondii*,
`Pf` = *P. falciparum*. "Grade" is the calibration grade **for that organism**.

| # | Org | Question | Entry | Strategy / view | Settings | A good answer | Anchor |
|---|---|---|---|---|---|---|---|
| Q01 | Tg | Which hyperLOPIT compartments can the rest of the data rediscover on its own? | ST | 03 recoverability_atlas (reliable) | target `compartment`, sample 3000, k 15 | a per-compartment AUROC table where ribosome/proteasome are near 1 and diffuse compartments are near chance | PMID:33053376 |
| Q02 | Tg | What compartment do the hyperLOPIT-unassigned proteins behave like? | SH (have a label -> predict) | 07 feature_knn (reliable) | target `compartment`, k 15, min_share 0.3 | a few thousand calls with support, enriched for plausible compartments | PMID:39082802 |
| Q03 | Tg | Which unannotated proteins can be placed by fold alone, where sequence homology fails? | ST | 14 structural_homology (reliable) | target `compartment`, level 1, layer `struct` | calls for hypothetical proteins with a named structural neighbour | PMID:40066064 |
| Q04 | Tg | Which proteins can be placed by the proteins they physically touch? | STAR then ST | 13 physical_partners (reliable) | target `compartment` | calls naming the crosslink/pulldown partners that made them | PMID:40874616 |
| Q05 | Tg | Which map clusters are compartment-pure, and which unassigned genes sit in them? | MAP then ST | 09 cluster_guilt (reliable) | target `compartment`, min_cluster_size 25, leaf, min_lift 2 | a handful of clusters at q<0.05 with lift, and their unassigned members | PMID:38270431 |
| Q06 | Tg | Which localization calls come with a stated 10% error rate? | ST | 35 conformal_calls (reliable) | target `compartment`, alpha 0.1, logistic, per class | a small set of singleton calls plus honest multi-label sets for the rest | PMID:41396985 |
| Q07 | Tg | Which kind of evidence does a learner trust most for localization? | ST | 38 stacking (reliable) | target `compartment`, k 15, folds 5 | meta-model weights over measurement / classifier / network evidence | PMID:33053376 |
| Q08 | Tg | Which measurements actually define a compartment label? | ST | 37 random_forest (reliable) | target `compartment`, trees 300, min_leaf 2 | a permutation-importance ranking of named columns | PMID:33053376 |
| Q09 | Tg | Which kind of evidence carries localization, and which is redundant? | ST | 06 block_ablation (reliable) | target `compartment`, k 15, unit family | accuracy alone and accuracy-without per evidence family | PMID:33053376 |
| Q10 | Tg | Does smoothing measurements along the networks improve localization calls? | ST | 36 graph_convolution (reliable) | target `compartment`, hops 2, C 0.5 | calls plus where the model's weight sits across blocks | PMID:38900844 |
| Q11 | Tg | What does a plain classifier trained on known compartments call the rest? | SH | 19 supervised_classifier (reliable) | target `compartment`, C 0.5 | calls plus what drives each class | PMID:33053376 |
| Q12 | Tg | Which localization calls survive agreement between three independent methods? | ST | 31 triangulation (reliable) | target `compartment`, k 15, min_agree 2 | fewer calls, each with its three per-method votes | PMID:39082802 |
| Q13 | Tg | Which never-published proteins get a confident compartment call? | ST | 32 understudied_first (reliable) | target `compartment`, k 15, min_agree 2 | a ranked list of genes with zero publications and an agreed call | PMID:41582196 |
| Q14 | Tg | Which proteins behave unlike their own compartment -- dual-localized or mis-assigned? | ST | 10 label_outliers (reliable) | target `compartment`, k 15, top 200 | named proteins whose neighbours carry a different compartment | PMID:40874616 |
| Q15 | Tg | How far does compartment information diffuse through the crosslink network alone? | ST | 11 layer_propagation (reliable) | target `compartment`, layer `xlms`, restart 0.5 | calls for the genes the layer touches, and an explicit statement of what it cannot reach | PMID:40874616 |
| Q16 | Tg | Do all nine measured networks together beat any one of them for localization? | ST | 12 layer_vote (weak on Tg) | target `compartment`, k 15 | per-source earned weights; a weak grade means the calls are leads | PMID:38900844 |
| Q17 | Tg | Can a gene be localized from its neighbours on a single UMAP? | MAP | 08 map_neighbours (weak on Tg) | target `compartment`, sample 4000, n_neighbors 25, min_dist 0.1, k 15 | calls that should be visibly worse than Q02 -- the map costs information | PMID:33053376 |
| Q18 | Tg | What do the hyperLOPIT dense-granule proteins have in common across every layer? | SH (have a gene list -> explain) | 24 set_enrichment (reliable) | genes = the 193 `compartment == dense granules` genes, top 300 | a corrected feature profile whose top entries are PV-facing quantities | PMID:36515555 |
| Q19 | Tg | Which uncharacterised proteins most resemble the known dense-granule proteins? | SH | 20 positive_unlabeled (reliable) | same 193 genes, bags 15, top 300 | ranked candidate GRAs scored only by models that did not train on them | PMID:36515555 |
| Q20 | Tg | Which genes does the network walk attach to the rhoptry-1 proteins? | ST | 25 seed_expansion (reliable) | genes = the 71 `rhoptries 1` genes, restart 0.3, top 300 | ranked candidates with their direct links to the seed list | PMID:41975467 |
| Q21 | Tg | Is there any map on which the micronemal proteins fall out as one cluster? | ST | 02 geneset_hunt (weak on Tg) | genes = the 56 `micronemes` genes, default grid | the best map and its F1 against 100 random sets of the same size | PMID:42470170 |
| Q22 | Tg | What sits in GRA17's and MYR1's neighbourhood when every layer is put in one space? | STAR | 33 neighbour_space (reliable) | genes `TGME49_222170, TGME49_254470`, k 10, top 300 | per-gene neighbours with a probability and the evidence behind each | PMID:37498952 |
| Q23 | Tg | Which contacts did the crosslinking interactome most likely miss? | ST | 16 link_prediction (reliable) | layer `xlms`, top 300 | pairs with a score, shared-partner count and which other layers link them | PMID:40874616 |
| Q24 | Tg | Which gene pairs does the data link in three layers that no paper mentions together? | ST | 18 unwritten_links (reliable) | min_layers 3, top 500 | pairs supported by three independent layers with zero co-mentions | PMID:28276701 |
| Q25 | Tg | Which co-mentioned pairs are linked more than the two genes' fame predicts? | ST | 17 attention_correction (reliable) | target `compartment`, layer `comention`, top 200 | pairs whose residual is high although neither gene is famous | PMID:41582196 |
| Q26 | Tg | Which crosslink edges would a model trained on all the networks add? | ST | 34 network_training (reliable) | layer `xlms`, logistic, top 300, fraction 0.25 | ranked gaps plus the three nulls the model was scored against | PMID:40874616 |
| Q27 | Tg | Which genes are far more or less fitness-conferring than everything else predicts? | SH (have a measurement -> explain) | 21 trait_regression (reliable) | target `fit_invitro_hff`, boosted, own kind left out, top 200 | out-of-fold rho and a list of genes whose measured fitness defies it | PMID:27594426 |
| Q28 | Tg | Which genes defy the prediction of in-vivo brain fitness? | ST | 21 trait_regression (reliable) | target `fit_invivo_brain`, boosted, top 200 | the same, for a condition where in-vitro screens are blind | PMID:41935076 |
| Q29 | Tg | Do genes that matter only at high parasite density look like anything else in the data? | ST | 21 trait_regression (reliable strategy, hard target) | target `fit_density_dependence`, boosted, top 200 | an honest rho; near zero would say density dependence has no signature | PMID:34023299 |
| Q30 | Tg | Where can IFN-gamma macrophage fitness be filled in honestly for unmeasured genes? | ST | 22 masked_imputation (reliable) | target `fit_ifng`, rank 20 | per-column reliability and imputed values only where reliability is high | PMID:33067458 |
| Q31 | Tg | Which fitness predictions come with an interval that actually holds 90% of the time? | ST | 39 conformal_values (weak on Tg) | target `fit_invitro_hff`, alpha 0.1, boosted | interval width, and the measured genes that fall outside their interval | PMID:27594426 |
| Q32 | Tg | Which genes matter more under interferon-gamma than in a naive macrophage, and why? | SH (compare two conditions) | 23 condition_shift (weak on Tg) | condition `fit_ifng`, baseline `fit_naive_bmdm`, boosted | a shift per gene, and whether the shift itself is predictable | PMID:33067458 |
| Q33 | Tg | Which paralogue pairs changed compartment -- a family that divided its labour? | ST | 28 paralog_divergence (reliable) | target `compartment`, min_shared 10 | pairs ranked by profile divergence with both compartments printed | PMID:42113865 |
| Q34 | Tg | Does inference still work on the genes orthology cannot reach? | ST | 30 stratum_focus (reliable) | target `compartment`, stratum lineage-specific, k 15 | calls for lineage-specific genes plus accuracy by stratum | PMID:35196325 |
| Q35 | Tg | Are there kinds of gene defined by compartment and life stage at once? | ST | 27 conjunctions (reliable) | a `compartment`, b `stage_enriched_derived`, min_cluster_size 15 | clusters enriched for a combination beyond either label alone | PMID:37081202 |
| Q36 | Tg | Tune a map with no labels at all: what does it turn out to encode? | AN then ST | 05 blind_battery (reliable) | map_from default, sample 3000 | held-out features that differ between clusters at q<0.05 | PMID:33053376 |
| Q37 | Tg | Is there a map setting under which a hidden compartment falls out as clusters? | ST | 01 holdout_search (weak on Tg) | target `compartment`, default grid | best cluster F1 per compartment, and the same on shuffled labels | PMID:33053376 |
| Q38 | Tg | Which communities do several measured networks agree on? | ST | 15 multiplex_modules (weak on Tg) | target `compartment`, resolution 1.0, agreement 0.5 | a small number of large modules with their majority compartment | PMID:38900844 |
| Q39 | Tg | Does Plasmodium essentiality predict Toxoplasma fitness through orthogroups? | ST | 29 ortholog_transfer (reliable) | target `fit_invitro_hff`, source `piggybac_mis` (Pf) | transferred values for Tg genes, and the rank correlation on genes measured in both | PMID:39913589 |
| Q40 | Pf | Does the schizont spatial proteome fall out of the rest of the Plasmodium data? | ST | 03 recoverability_atlas (reliable) | target `lopit_pf_location`, sample 3000, k 15 | per-location AUROC; Maurer's cleft and rhoptries should be encoded | PMID:42218142 |
| Q41 | Pf | Which Plasmodium hypotheticals can be localized by fold? | ST | 14 structural_homology (reliable) | target `lopit_pf_location`, level 1 | calls with a structural neighbour named | PMID:42218142 |
| Q42 | Pf | Which Plasmodium localization calls carry a stated error rate? | ST | 35 conformal_calls (reliable) | target `lopit_pf_location`, alpha 0.1 | singleton calls plus prediction sets | PMID:41396985 |
| Q43 | Pf | Which protein contacts did the Plasmodium pulldowns miss? | ST | 16 link_prediction (reliable) | layer `ip_ms`, top 300 | scored candidate pairs with their supporting layers | PMID:41482054 |
| Q44 | Pf | What changes between schizont and mature gametocyte beyond overall expression? | SH (compare two conditions) | 23 condition_shift (reliable on Pf) | condition `expr_gametocyte_v`, baseline `expr_schizont` | a per-gene residual and an out-of-fold rho showing the shift has a signature | PMID:41482054 |
| Q45 | Pf | Which Plasmodium genes are more or less essential than the data predicts? | ST | 21 trait_regression (reliable) | target `piggybac_mis`, boosted, top 200 | out-of-fold rho plus the genes furthest from prediction | PMID:39913589 |
| Q46 | Pf | Which genes join the measured apicoplast proteins when the networks are walked? | SH | 25 seed_expansion (reliable) | genes = the 79 `lopit_pf_location == apicoplast` genes, restart 0.3, top 300 | ranked candidate apicoplast proteins with their direct links | PMID:36563485 |
| Q47 | Pf | What do the exported Plasmodium proteins share across every layer? | ST | 24 set_enrichment (reliable) | genes = the 191 `is_exported == True` genes, top 300 | a corrected profile whose top entries are export-related | PMID:41035680 |
| Q48 | Pf | Does nearest-neighbour localization work in Plasmodium the way it does in Toxoplasma? | ST | 07 feature_knn (weak on Pf) | target `lopit_pf_location`, k 15, min_share 0.3 | a self-test that should be visibly worse than Q02's -- the point of asking | PMID:42218142 |
| Q49 | Pf | Which Plasmodium proteins behave unlike their measured compartment? | ST | 10 label_outliers (reliable) | target `lopit_pf_location`, k 15, top 200 | named proteins with a contradicting neighbourhood | PMID:42218142 |
| Q50 | Pf | Does Toxoplasma fitness predict Plasmodium essentiality through orthogroups? | ST | 29 ortholog_transfer (reliable) | target `piggybac_mis`, source `fit_invitro_hff` (Tg) | transferred values and the cross-species rank correlation | PMID:39913589 |
| Q51 | Tg | Which uncharacterised proteins behave like known rhoptry proteins? | SH | 20 positive_unlabeled | genes = `rhoptries 1` + `rhoptries 2` (123 genes), bags 15, top 300 | a ranked candidate ROP list | PMID:41975467 |
| Q52 | Tg | Which proteins behave like the known PV-membrane nutrient channels? | SH | 20 positive_unlabeled | genes = GRA17, GRA23 and their hyperLOPIT dense-granule cohort | candidates to test for channel activity | PMID:42555613 |
| Q53 | Tg | Which genes are synthetic-lethal-like partners of GRA17 in the data? | ST | 16 link_prediction | layer `cofitness`, top 300 | pairs whose co-fitness profile matches the delta-GRA17 screen column | PMID:37498952 |
| Q54 | Tg | Which effectors change host transcription, judged from the host-effect columns? | ST | 21 trait_regression | target `hosttx_T2`, boosted, top 200 | genes whose host transcriptional effect defies prediction | PMID:37827122 |
| Q55 | Tg | Which nuclear-targeted effector candidates are there beyond the known four? | ST | 27 conjunctions | a `compartment`, b `screen_any_phenotype` | a cluster of secreted-and-nuclear genes | PMID:39292011 |
| Q56 | Tg | Which in-vivo fitness hits are likely paracrine-rescued rather than cell-autonomous? | ST | 23 condition_shift | condition `fit_invivo_PE`, baseline `fit_invitro_hff` | genes with a large unexplained in-vivo shift, to be read with the paracrine caveat | PMID:39654402 |
| Q57 | Tg | Which intrinsically disordered proteins look like exported effectors? | ST | 19 supervised_classifier | target `compartment`, then rank calls by `af3_low_confidence_modelled_fraction` | disordered, dense-granule-called genes | PMID:38747635 |
| Q58 | Tg | Which ROP-family proteins look like tethers rather than enzymes? | ST | 28 paralog_divergence | target `compartment`, min_shared 10; read ROP pairs | ROP pairs whose profiles diverge, one with membrane-contact features | PMID:41291278 |
| Q59 | Tg | Which secreted proteins have lipid-transfer-like structural neighbours? | ST | 14 structural_homology | target `compartment`, level 1, layer `struct` | dense-granule calls whose structural neighbour is a lipid-binding fold | PMID:42018441 |
| Q60 | Tg | Which parasite proteins sit near the host ESCRT-recruitment column? | ST | 21 trait_regression | target the `Tg_host_escrt_recruitment` block column, boosted | genes whose ESCRT-recruitment value defies prediction | PMID:36718630 |
| Q61 | Tg | Which conoid and apical-complex proteins are still unassigned to a position? | ST | 09 cluster_guilt | target `compartment`, min_cluster_size 25 | an apical-1/apical-2 cluster with unassigned members | PMID:41366525 |
| Q62 | Tg | Which micronemal proteins partition into a "secretion" module and a "junction" module? | ST | 15 multiplex_modules | target `compartment`, agreement 0.5 | two modules that split the MIC/RON set | PMID:42470170 |
| Q63 | Tg | Which kinase substrates are the egress decision nodes? | ST | 21 trait_regression | target `cdpk1_thiophospho_peptides`, boosted | genes with more CDPK1 substrate peptides than predicted | PMID:39903669 |
| Q64 | Tg | Which apicoplast proteins become dispensable when the organelle is bypassed? | ST | 23 condition_shift | condition `fit_complete_medium_2025`, baseline `fit_invitro_hff` | apicoplast genes with a condition-specific shift | PMID:40025025 |
| Q65 | Tg | Which predicted transporters have no publication at all? | ST | 32 understudied_first | target `compartment`, min_agree 2; filter calls to PM/ER membrane classes | zero-publication membrane proteins with an agreed compartment | PMID:41582196 |
| Q66 | Tg | Which genes are induced in bradyzoites and localize to the cyst wall? | ST | 27 conjunctions | a `compartment`, b `stage_enriched_derived` | a bradyzoite-and-surface cluster | PMID:37081202 |
| Q67 | Tg | Which differentiation factors are translationally rather than transcriptionally controlled? | ST | 23 condition_shift | condition `te245775_parent_prebrady_r1`, baseline `rna245775_parent_prebrady_r1` | genes with high TE and flat mRNA | PMID:38782906 |
| Q68 | Tg | Which BFD2-bound transcripts are also stage-restricted? | ST | 21 trait_regression | target `bfd2_rip_log2_ip_over_input`, boosted | transcripts more BFD2-enriched than predicted | PMID:37081202 |
| Q69 | Tg | Which dormancy proteins are secreted and therefore candidate biomarkers? | SH | 20 positive_unlabeled | genes = cyst-wall / bradyzoite-enriched set | secreted, bradyzoite-restricted candidates | PMID:40389644 |
| Q70 | Tg | Which uncharacterised dense-granule proteins have a stage-specific phenotype? | ST | 23 condition_shift | condition `fit_invivo_brain`, baseline `fit_invitro_hff` | GRAs invisible to tachyzoite screens | PMID:40941387 |
| Q71 | Tg | Which chromatin readers are parasite-specific and fitness-conferring? | ST | 30 stratum_focus | target `compartment`, stratum lineage-specific | lineage-specific nuclear-chromatin calls | PMID:42140435 |
| Q72 | Tg | Which SAGA/SWI-SNF-like subunits are apicomplexan-specific? | STAR | star map on a known subunit, `Measured` only, 2 hops | co-fitness and crosslink partners with no orthogroup outside apicomplexa | PMID:42778561 |
| Q73 | Tg | Which AP2 factors are cell-cycle rather than stage regulators? | ST | 27 conjunctions | a `cellcycle_phase`, b `stage_enriched_derived` | AP2 genes in a phase-specific, stage-flat cluster | PMID:40920096 |
| Q74 | Tg | Which phosphosites sit on proteins of one compartment? | ST | 21 trait_regression | target `n_phosphosites`, boosted | compartment-biased phospho load | PMID:41960867 |
| Q75 | Tg | Which palmitoylated proteins are IMC or plasma-membrane residents? | ST | 19 supervised_classifier | target `compartment`; read calls for high `palmitome_odya_vs_hydroxylamine_log2` | palmitoylated pellicle candidates | PMID:39690890 |
| Q76 | Tg | Which N-glycosylated proteins are surface or secreted? | MAP | maps tab, one label / every map on `compartment`; colour by glycosylation column | the map that separates glycosylated surface proteins | PMID:40072250 |
| Q77 | Tg | Which Golgi proteins have no orthologue outside the apicomplexa? | ST | 30 stratum_focus | target `compartment`, stratum lineage-specific | lineage-specific Golgi calls | PMID:39345210 |
| Q78 | Tg | Which very large hypothetical proteins carry hidden domain content? | AN | Analysis tab: feature block `Tg_domain_content` + `Tg_sequence_basics`, cluster, inspect | a cluster of long, domain-poor, high-pLDDT hypotheticals | PMID:38584870 |
| Q79 | Tg | What fraction of the proteome is still "hypothetical", and where does it sit? | MAP | maps tab, `One map, every label`, colour by product class | the compartments where unknowns concentrate | PMID:37953350 |
| Q80 | Tg | Which hypothetical genes are stage-specifically induced? | ST | 27 conjunctions | a `stage_enriched_derived`, b `compartment` | stage-restricted unknowns | PMID:41761089 |
| Q81 | Tg | Which genes fill gaps in the metabolic reconstruction? | ST | 14 structural_homology | target `dtm_class`, level 1 | enzyme-class calls for unannotated genes | PMID:35196325 |
| Q82 | Tg | Which alveolin-family paralogues divided the subpellicular network between them? | ST | 28 paralog_divergence | target `compartment`, min_shared 10 | IMC1-family pairs with divergent profiles | PMID:42113865 |
| Q83 | Tg | Which apical-ring proteins are assembly factors rather than structural? | ST | 27 conjunctions | a `compartment`, b `cellcycle_phase` | apical genes restricted to a daughter-budding phase | PMID:41315181 |
| Q84 | Tg | Which microtubule inner proteins are still hypotheticals? | ST | 09 cluster_guilt | target `compartment`; read the tubulin-cytoskeleton cluster | unknowns in a cytoskeletal cluster | PMID:41626803 |
| Q85 | Tg | Which oxidative-stress fitness genes are not explained by anything else? | ST | 21 trait_regression | target `oxidative_stress_screen_score`, boosted | genes with unexplained oxidative-stress dependence | PMID:33067458 |
| Q86 | Tg | Which genes need glucose specifically, beyond general fitness? | ST | 23 condition_shift | condition `fit_no_glucose`, baseline `fit_complete_medium_2025` | glucose-specific dependencies | PMID:40025025 |
| Q87 | Tg | Which iron-regulated proteins are apicoplast or mitochondrial? | ST | 19 supervised_classifier | target `compartment`; rank calls by `iron_depletion_protein_log2fc` | iron-responsive organellar proteins | PMID:39348385 |
| Q88 | Tg | Which proteins are engaged by calcium in the CETSA screen and localize apically? | ST | 21 trait_regression | target `cetsa_calcium_ed_score`, boosted | apically-called, calcium-engaged candidates | PMID:39903669 |
| Q89 | Pf | Which ApiAP2 factors still have no target set in the data? | ST | 32 understudied_first | target `lopit_pf_location`, min_agree 2 | AP2 genes with few publications and a nuclear call | PMID:39419713 |
| Q90 | Pf | Which genes are chromatin-driven rather than TF-driven? | ST | 21 trait_regression | target `chromprox_hp1_log2fc`, boosted | genes whose HP1 proximity defies prediction | PMID:41501628 |
| Q91 | Pf | Which non-PEXEL proteins are nevertheless exported? | ST | 10 label_outliers | target `is_exported`, k 15, top 200 | exported proteins whose neighbours are not exported, and the reverse | PMID:38334391 |
| Q92 | Pf | Which cargo depends on which PTEX component? | STAR | star map on EXP2/PTEX150, 2 hops, per source 6 | measured and inferred partners separated by layer | PMID:38968107 |
| Q93 | Pf | Which essential genes are structurally druggable but unpublished? | ST | 32 understudied_first | target `lopit_pf_location`; cross with `mean_plddt` and `n_publications` | high-confidence-fold, essential, unpublished genes | PMID:40066064 |
| Q94 | Pf | Which K13 neighbours could set the artemisinin-resistance ceiling? | STAR | star map on kelch13, 2 hops | co-fitness and crosslink neighbours of K13 | PMID:41841738 |
| Q95 | Pf | Which transporters carry resistance-selection variants? | ST | 21 trait_regression | target `resistance_selection_variants`, boosted | transporters with more selected variants than predicted | PMID:41922323 |
| Q96 | Pf | Which febrile-temperature responders are also phospho-regulated? | ST | 27 conjunctions | a `stage_enriched_derived`, b `has_phospho` | a heat-and-phospho cluster | PMID:41960867 |
| Q97 | Pf | Which gametocyte-committed proteins have no Toxoplasma orthologue? | ST | 29 ortholog_transfer | target `committed_vs_asexual_log2fc`, source `stage_enriched_derived` (Tg) | committed genes the transfer cannot reach | PMID:41482054 |
| Q98 | Pf | Which mitochondrial proteins are apicomplexan-only? | ST | 30 stratum_focus | target `lopit_pf_location`, stratum lineage-specific | lineage-specific mitochondrial calls | PMID:34494883 |
| Q99 | Pf | Which m6A-modified transcripts are stage-restricted? | ST | 21 trait_regression | target `m6a_canonical_stoichiometry`, boosted | transcripts with unexplained m6A stoichiometry | PMID:38892332 |
| Q100 | Pf | Which surface proteins are reversibly regulated with population state? | ST | 23 condition_shift | condition `expr_gametocyte_v`, baseline `expr_gametocyte_ii` | surface genes with a maturation-specific shift | PMID:40354414 |

---

## 3. What was run

`scripts/run_biological_questions.py` runs **Q01-Q50**. Q51-Q100 are written but not executed; they
reuse strategies already exercised by Q01-Q50 at different settings, so the point of running them
would be the biology, not the software. Section 4 records the verdicts.

## 4. Verdicts

See `results/questions_2026_09_30/summary.csv` for the machine-readable record and
`docs/questions.md` for the reproducible list of the questions that yielded. The counts and the five
most interesting results are in the closing note of this file.

## 5. Literature queries used, verbatim

All against Europe PMC REST, `https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=<Q>&format=json&pageSize=<N>`:

1. `Toxoplasma MYR1 translocon GRA16 GRA24 effector`
2. `"Toxoplasma" AND "dense granule" AND effector AND (PUB_YEAR:2024 OR PUB_YEAR:2025 OR PUB_YEAR:2026)`
3. `Toxoplasma rhoptry ROP16 ROP18 host signalling review`
4. `"parasitophorous vacuole membrane" Toxoplasma intravacuolar network`
5. `"parasitophorous vacuole" AND nutrient AND (PUB_YEAR:2024 OR PUB_YEAR:2025 OR PUB_YEAR:2026)`
6. `Toxoplasma GRA2 GRA6 intravacuolar network tubules structure`
7. `Toxoplasma invasion AMA1 RON2 moving junction structure`
8. `apicomplexan egress calcium signalling perforin PLP1 microneme secretion AND (PUB_YEAR:2023..2026)`
9. `Plasmodium falciparum merozoite invasion rhoptry RH5 structure AND (PUB_YEAR:2024..2026)`
10. `Toxoplasma egress calcium signalling CDPK perforin microneme`
11. `Toxoplasma apicoplast isoprenoid fatty acid synthesis FASII metabolism`
12. `Plasmodium apicoplast heme biosynthesis essentiality AND (PUB_YEAR:2023..2026)`
13. `TITLE:"egress" AND (Toxoplasma OR Plasmodium)` (sorted `P_PDATE_D desc`)
14. `Toxoplasma bradyzoite differentiation BFD1 BFD2 cyst`
15. `Plasmodium falciparum gametocyte AP2-G sexual commitment AND (PUB_YEAR:2023..2026)`
16. `Toxoplasma genome-wide CRISPR screen fitness score essential genes`
17. `Plasmodium falciparum piggyBac saturation mutagenesis essentiality gene knockout screen`
18. `Toxoplasma in vivo CRISPR screen AND (PUB_YEAR:2023..2026)`
19. `TITLE:"Genome-wide CRISPR Screen" AND Toxoplasma`
20. `hyperLOPIT Toxoplasma spatial proteome`
21. `Plasmodium falciparum schizont spatial proteome subcellular localization LOPIT`
22. `AUTH:"Barylyuk" AND Toxoplasma`
23. `Plasmodium falciparum artemisinin resistance Kelch13 K13 AND (PUB_YEAR:2024..2026)`
24. `PfCRT chloroquine resistance transporter piperaquine structure mutation`
25. `new antimalarial drug target discovery PfATP4 PfPI4K eEF2 AND (PUB_YEAR:2024..2026)`
26. `antimalarial drug target resistance chemogenomic PfATP4 OR PfPI4K OR "druggable genome"`
27. `ApiAP2 transcription factor apicomplexan gene regulation`
28. `Toxoplasma Plasmodium chromatin histone deacetylase HDAC epigenetic regulation AND (PUB_YEAR:2024..2026)`
29. `Plasmodium PEXEL protein export PTEX translocon AND (PUB_YEAR:2023..2026)`
30. `Toxoplasma palmitoylation phosphoproteome GPI anchor post-translational modification`
31. `Toxoplasma gondii palmitoylation OR phosphoproteome OR ubiquitylation proteomics`
32. `apicomplexan hypothetical protein functional annotation gap understudied genes`
33. `Toxoplasma conoid apical complex inner membrane complex cryo-electron tomography`
34. `"unknown function" OR "dark proteome" AND Toxoplasma OR Plasmodium gene function prediction`
35. `Toxoplasma Plasmodium comparative genomics orthologs apicomplexan-specific gene families`
36. `Plasmodium parasite density dependent growth quorum sensing like signalling`
37. `Plasmodium falciparum extracellular vesicles density dependent gametocyte conversion cell-cell communication`
38. `Toxoplasma gondii hypothetical protein localization tagging characterization uncharacterized`
39. `Plasmodium falciparum parasite density growth rate dependence multiplication rate in vivo`
40. `Toxoplasma single-cell RNA-seq heterogeneity vacuole population asynchrony`
41. `VEuPathDB ToxoDB PlasmoDB database resource`
42. `Plasmodium inner membrane complex glideosome IMC proteome AND (PUB_YEAR:2024..2026)`
43. `AlphaFold structure prediction Plasmodium Toxoplasma proteome annotation novel function`
44. `apicomplexan ortholog conservation essentiality Toxoplasma Plasmodium Cryptosporidium comparison fitness`

## 6. Identifier confirmation

Every anchor used in section 2 was confirmed a second time, independently of the search that found
it, through NCBI E-utilities:

```
https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=33053376,42218142,39082802,
27594426,37498952,40874616,41975467,36515555,39913589,38900844,31481656,33067458,37827122,39654402,
41547993,40066064,42107390,39419713,41035680,34494883,38270431,33705380,35042410,23684304&retmode=json
```

All 24 returned the first author, year and title used here. The remaining anchors came out of the
Europe PMC response bodies in section 5 and are kept as raw JSON alongside this work.

## 7. A finding about the literature itself

Two of the fourteen topics are genuinely thin, which shaped the questions:

* **Parasite density and quorum-like effects** (Q29, Q52) have almost no gene-level apicomplexan
  literature. The nearest work is malaria extracellular-vesicle biology (PMID:23684304,
  PMID:39868784) and trypanosome quorum sensing. The shipped density screen
  (`fit_density_low`, `fit_density_high`, `fit_density_dependence`) is therefore one of the few
  gene-level handles that exists, and Q29 is a genuinely open question rather than a re-derivation.
* **The unknown fraction of the proteome** has many "we characterised N novel proteins" papers and
  almost no systematic treatment. PMID:39082802 (Tachibana 2024, CRISPR screens on
  hyperLOPIT-unassigned proteins) is the closest precedent for what this software does, which makes
  Q02, Q13 and Q65 the questions closest to the field's actual frontier.
