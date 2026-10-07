# Instructions: what is done and what is left

**New here? Read [START_HERE.md](START_HERE.md) first.**

One file per task. This index is the status table; the files carry the reasoning and the
"done when" condition.

## Done

| # | Task |
|---|---|
| 01 | Nineteen gene codes resolved |
| 02 | Six expression datasets integrated |
| 03 | GUI bugs |
| 04 | `projectionMatrix` compatibility — clicking genes |
| 05 | Visual quality |
| 06 | Studies with no supplementary data |
| 07 | Verify shipped columns against GEO |
| 08 | Standalone data and auto-download |
| 09 | Cell-cycle dataset |
| 10 | Search more targets |
| 11 | README and API reference |
| 12 | Test every aspect of the library |
| 13 | Jobs, stopping, and freeing memory |
| 14 | Tooltips, docstrings, and the scoring explainer |
| 15 | The UMAP gallery — grid, scroll, incremental |
| 20 | Validation: `refit`, candidates, orthogonal evidence, per-category scores in Inference |
| 19 | Logging — opt-in, per-level console control |
| 18 | Navigate / Select modes, 2D and 3D gating |
| 16 | Automated walk — UMAP, clustering, per-category scoring |
| 17 | Annotation store |
| 23 | Left panel becomes "color by" |
| 22 | Cell diagram under the compartment list |
| 25 | American spelling in user-facing text |
| 24 | Forty black-and-white logo drafts |
| 21 | Scoring objectives in the interface |
| 26 | Save and load results, per tab and all at once |
| 27 | Import data, with the preprocessing offered |
| 25b | The identifier rename — job 2 of the spelling task |
| 28 | Every module at 100% coverage |
| 29 | Three more leaks in the circularity guard |
| 30 | Fold in datasets already present and expand missing biological axes |
| 31 | Replace file-shaped feature blocks with biological question slots |
| 32 | GPU paths for k-means, DBSCAN and t-SNE |
| 33 | Crossed-factor findings and fragmentation diagnostic |
| 34 | Repair two metric columns and the CLI finding count |
| 35 | Simplify lighting, make point modes distinct, and add volumetric ray tracing |
| 36 | GPU PBR sphere points and stable GPU density-ray tracing |
| 37 | Continuous flashlight, material lab, and ray-rendering comparison |
| 48 | The orphan alarm never looks at the host table |
| 49 | A Strategies tab: 32 ways to infer something, each testing itself |
| 50 | The September 2026 data audit: a wrong citation, 21 deposits, two leaks closed |
| 51 | How good each strategy is: calibration over settings, targets and seeds |
| 52 | The second September 2026 data audit: six Plasmodium deposits, five refusals |
| 54 | [Scorecards, techniques and five advanced strategies](done/54_scorecards_and_advanced_models.md) |
| 55 | [Strategy cards: four standard bars, the test explained, a real failure and success](done/55_strategy_cards.md) |
| 56 | [Pregenerated map gallery and label x map scores](done/56_umap_gallery_label_scores.md) |
| 56 | [Star map: gene-gene links by source, with provenance, including your own runs](done/56_star_map.md) |
| 56 | [The parasite-density screen (Giuliano et al., Cell 2026): fitness at low and high density, and what crowding needs](done/56_lourido_density_screen.md) |
| 57 | ["Start here": a guided path from what a user has to the strategies worth running](done/57_start_here_tab.md) |
| 58 | [The blank panel: a tab switched away from comes back as bare background](done/58_blank_panel_bug.md) |
| 61 | [Tutorial videos: six screen captures of the real application](done/61_tutorial_videos.md) |
| 54 | [Many more pregenerated maps, searched for structure, and a navigable maps panel](done/54_many_umap_maps_and_the_score_table.md) — shipped in 0.50.0 |
| 58 | [Star map: stop the hover twitch, show the larger network, define the connections](done/58_star_map_steady_larger_defined.md) — merged and shipped in 0.50.0 |
| 64.01 | [Define organism-qualified entities, questions and provenance-preserving alias resolution](done/64_01_entity_query_contract.md) — 472 focused checks passed |
| 64.02 | [Inventory registered evidence and reconcile installed storage coverage](done/64_02_dataset_inventory.md) — 162 sources, 180 rows, 16 scoped refusals; 525 focused checks passed |
| 65.01 | [Versioned citation-rate, recency and comprehensiveness preferences](done/65_01_dataset_selection_policy.md) — 389 focused checks passed; publication audit follows |

## Open

The current execution order is [instruction 64](open/64_information_space_and_inference_atlas.md):
the user's October 7 goals, 40 controlled action cards, a ground-truth contract for all 39 current
strategies, and the full progress table. Existing instructions below retain their evidence and
remaining work; instruction 64 integrates their priorities. Repost the complete progress table
whenever an action is completed, replacing its percentage with a green tick.

| # | Task |
|---|---|
| 65 | [Audit every dataset choice and retain the citation-rate/recency/comprehensiveness rule](open/65_dataset_selection_audit.md) — four bounded policy, publication, challenger-search and validated-promotion items |
| 64 | [Explore published organism/host evidence and inspect precomputed inferences with ground-truth scorecards](open/64_information_space_and_inference_atlas.md) — 2/40 complete (5%); continue with 64.03 |
| 63 | [Claims: generate knowledge with a measured certainty, then test it independently](open/63_claims_generate_and_verify.md) — steps 1+2 shipped in 0.54.0; October 7 literature audit complete (no Tg candidate passes, Pf candidate adds only two claims; no promotion). Next: orthology/experimental verifier coverage, frozen prospective tests, genes outside the tested range |
| 62 | [The hold-out track record](open/62_holdout_track_record.md) — shipped 0.52.0-0.53.0 for 12 labels; left: pooled record vs scorecard check, the questions page |
| 38 | Build a provenance-first pan-Apicomplexan dataset archive |
| 41 | Fill as many slots as possible: Toxoplasma, Plasmodium, and host |
| 60 | [What would make Starplast most useful](open/60_usefulness_review.md): ranked plan; weeks 1-2 = fix the dock bug, result provenance, methods paragraph + citations, figure export |
| 59 | [One hundred biological questions this software can answer](open/59_biological_questions.md): 100 questions with entry point, strategy, settings and a resolved literature anchor; Q01-Q50 executed and graded yields/thin/no |
| 53 | [One space per organism](open/53_organism_spaces.md), then links between pathogen, vector and host. R0-R3, Pf display, pack framework and lossless calibration publishing implemented. Remaining: organism builders and published packs · Hs/Mm spaces · Space menu · Cp/Pb, Ag/As, Rn/Fc spaces · per-space leakage and calibration · cross-space bridges (version numbers to be assigned; 0.47-0.54 went to other work). Small follow-ups: retry active learning with conformal uncertainty; next Tg/host data audit; `spaces/As/abundance/PXD001647/` is junk to delete. See WORK_LOG.md for active progress and blockers. |

## The one ordering constraint that is not negotiable — now satisfied

**17 must not ship before 20**, and 20 landed on 2026-08-12. Annotation is easy to build and easy to
build wrongly. Without validation it gives the project a mechanism for manufacturing unvalidated
claims, which is the one thing it exists to prevent — and a candidate list looks identical whether it
is 90% right or 6%.

The constraint does not disappear now that 20 exists; it becomes a requirement on 17. A saved
annotation carries the validated precision and recall for its category, the cluster's composition,
and the enrichment over prevalence — the numbers the Validation tab now computes. Saving a row
without them would be the same failure by a shorter route.

## Standing constraints on all of it

Measurement, inference, absence and annotation never read as one another.

No candidate is shown without the number that says how much to believe it. On the full proteome one
row out of 7,710 (configuration x compartment) reaches F1 0.5 — `PM - peripheral 1`, at exactly
0.500. A cluster that looks pure is a lead to check, not a result.
