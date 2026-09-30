# 56. The parasite-density screen (Lourido lab, Cell 2026)

**Status: done (2026-09-29).** The user asked for "the new Cell paper from the Lourido lab ... low
vs high abundance of toxo, which genes are important at high density".

## The paper, resolved (never typed from memory)

| | |
|---|---|
| Title | Convergent evolution of metabolic regulation governs redox adaptation in Toxoplasma |
| Authors | Giuliano CJ, Kalluraya CA, Kloehn J, Sloan MA, Bunkofske ME, Hunter CA, Soldati-Favre D, Harding CR, Lourido S |
| Journal | Cell, 2026 Aug (first published 2026-08-11) |
| PMID | 42580337 |
| PMCID | PMC13531851 (NIHMS2204255, author manuscript) |
| DOI | 10.1016/j.cell.2026.07.029 (PII S0092867426008275) |

Queries:

1. Europe PMC `AUTH:"Lourido S" AND PUB_YEAR:[2025 TO 2026]`: 19 hits, one Cell paper (above).
2. Europe PMC `EXT_ID:42580337 AND SRC:MED`, `resultType=core`: authors, PMCID, abstract ("We
   screened ... to identify genes that support parasite fitness during crowding. NAD(P)+
   biosynthesis was required at high parasite density ...").
3. Crossref `works/10.1016/j.cell.2026.07.029`: title, author list, and the PII in the full-text link.
4. NCBI efetch `db=pmc&id=13531851`: the methods, the Figure 1 legend and the supplementary-material
   labels (Table S1 = "CRISPR screening data, related to Figure 1").
5. Europe PMC `AUTH:"Giuliano CJ" AND Toxoplasma`: there is no preprint of this paper.

This is the density screen: one barcoded genome-wide library, selected for four passages at MOI 1,
then split for four passages at MOI 0.3 (low density) and MOI 3 (high density). No other candidate
matched the description.

## The data

The article is not open access. Europe PMC `supplementaryFiles` refuses it, and PMC `/bin/` links
return a proof-of-work page to scripts. The publisher's CDN serves the table directly:
`https://ars.els-cdn.com/content/image/1-s2.0-S0092867426008275-mmc2.xlsx` (mmc2 = Table S1; mmc3 =
Table S2, RIP-Seq, not taken). It is stored at
`<STARPLAST_DATA>/DNA/CRISPR_screen/42580337/mmc2.xlsx` with `URLS.txt`, `SHA256SUMS.txt`
(sha256 `03e0d9a1...a1ef`) and an entry in the folder's `SOURCES.md`. The screen's reads are not in
GEO. GSE327264 and GSE329845 hold the RIP-Seq and RNA-Seq, and PXD077563 holds the proteomics.

Sheets used: `GWS mean L2FC to input` (the arms), `GWS UMI-filtered scores` (the contrast) and
`GWS DIM Hits` (the call).

The same paper was already in the map. Its actinomycin-D RNA-Seq (GSE329845) is
`mrna_log2_remaining_4h_actinomycin` (instruction 50), which is already cited to PMID 42580337.
The screen is the new part.

## Verification: the paper's numbers, recomputed from Table S1

| Paper states | Recomputed | |
|---|---|---|
| Arms agree at r = 0.995 (Fig. 1B) | Pearson r = 0.9951 (8,331 genes, passage 8) | reproduced |
| 31 density-inhibited mutants (DIMs) | 31 in the hit sheet. All 31 are among the 32 genes the contrast puts at Bonferroni p < 0.05 on the depleted side. The 32nd, TGGT1_264610, sits at adj. p 0.048 on 0.5 clones at high MOI | reproduced (1 borderline gene the authors dropped) |
| 13 of them hypothetical proteins | 13 | reproduced |
| 12 top-scoring DIMs: at least four-fold at adj. p < 1e-4 | 12, identical to the sheet's "high confidence" set | reproduced exactly |
| 4 of 5 NAD(P)+ biosynthesis genes are DIMs; NAPRT trends but misses the cut-off | NMNAT, NAD synthetase, NAD kinase and nicotinamidase are DIMs. NAPRT has a contrast of -2.24 on a single clone, so it has no p and no call | reproduced |
| NMNAT and NAD synthetase lead | They are the two most negative contrasts (-4.90, -4.08) | reproduced |

Sign check: cytosolic ribosomal proteins have a median of -9.8 in both arms, against -4.1 for other
genes. Negative means needed, as in every other fitness column.

## Columns (`starplast/deposits.py::density_screen`, deposit `crispr_parasite_density`)

| column | what | Tg genes |
|---|---|---|
| `fit_density_low` | gene score at MOI 0.3, passage 8 (mean gRNA log2 FC to input) | 7,461 (91.7%) |
| `fit_density_high` | the same at MOI 3 | 7,461 (91.7%) |
| `fit_density_dependence` | authors' gRNA-UMI contrast, log2(high/low). **Negative = needed at high density** | 5,291 (65.0%) |
| `fit_density_dependence_log10padj` | -log10 Bonferroni-adjusted two-sided t-test against non-targeting clones | 4,173 (51.3%) |
| `fit_density_dim` | 1 = one of the paper's 31 DIMs, 0 = scored and not called, empty = not scored | 5,379 (66.1%) |

GT1 accessions resolve to ME49 through the identity layer (`add_deposits.toxo_resolver`). 8,331
rows reach 7,461 genes; 802 GT1 rows have no ME49 gene. 29 of the 31 DIMs land. The two that do not
are TGGT1_314770 (tRNA m1G methyltransferase, -3.48) and TGGT1_254255 (hypothetical, -1.39): the
GT1 strain table maps both to N/A. TGME49_500365 carries TGME49_314770 as a previous id and the
same product, so it is probably the same gene. It is not mapped here, because mapping by accession
number is not a rule the identity layer uses. This is a candidate for an identity-layer fix.

The 266-gene targeted follow-up (with a nicotinamide arm) is not shipped because it is not genome-wide.

## Slot and leakage

* Slot: new **`fitness · parasite density`** (`Tg_fitness_parasite_density`), context "HFF, MOI
  0.3 and MOI 3", policy `separate`, one slot per screen as for serum and carbon source. Slots go
  from 296 to 297, and gene-unit coverage from 180/221 to 181/222.
* `search.SAME_QUANTITY`:
  * The two arms join **`fitness in vitro`**. They correlate with `fit_invitro_hff` at rho 0.70.
    In a counterfactual audit that left them out of that family, the density slot was the best
    permitted predictor of held-out fibroblast fitness (0.716, robust z 3.2). With them grouped,
    the best permitted predictor is bradyzoite transcription at 0.491.
  * A new family, **`parasite density`**, holds all five columns. The contrast is computed from
    the arms' clones, and the call from the contrast and its p.
  * The contrast is orthogonal to bulk fitness (rho -0.08) and stays out of `fitness in vitro`.
* `scripts/leakage_audit.py` (target `fit_density_dependence` added):
  `results/leakage_audit_2026-09-29/`. Tg: **0 gaps, 0 residual leaks** (400 columns). Pf: 0 and
  0. The strongest permitted predictor of the new target is in vivo spleen fitness at 0.150.

## Rebuilt

`scripts/derive_deposits.py` (notebook section 20, re-executed), `scripts/add_deposits.py` (Tg 438
→ 443 columns; nothing lost; the host tables were only rewritten byte-for-byte, so they were
restored), `scripts/generate_slot_table.py`, `scripts/update_datasets.py`,
`scripts/generate_dataset_scripts.py` (`scripts/datasets/crispr_parasite_density.py`) and
`scripts/rebuild_layouts.py --organism Tg` (314 → 319 layout features; edges unchanged).
The calibration sweep (recipe step 9) was **not** run. It needs a multi-hour sweep and the user's
go-ahead.

## Other new genome-scale Toxoplasma screens seen (listed only, not added)

Found with Europe PMC queries `TITLE_ABS:"Toxoplasma" AND (genome-wide | CRISPR screen | genetic
screen | pooled screen) AND FIRST_PDATE:[2025-09-01 TO 2026-09-29]` and a second query on screen +
CRISPR/knockout/knockdown:

* **Host N-glycan / rhoptry discharge** (K562 host screen), now published: Valleau et al., EMBO J
  2026, PMID 42791346, DOI 10.1038/s44318-026-00911-z. It is already shipped as
  `host_k562_rhoptry_screen` and cited to the bioRxiv preprint, so the citation should be updated.
* **ELFN2 host autophagy screen**: Autophagy 2026, PMID 42329082, DOI 10.1080/15548627.2026.2693781.
  A host-side genome-wide screen.
* **EAF1 pooled image-based CRISPR screen** (ESCRT subversion): bioRxiv DOI
  10.64898/2026.07.08.737057.
* **MIC11 in vivo fitness gene**: Nat Commun 2026, PMID 41935076, DOI 10.1038/s41467-026-71423-x.
  Scope (targeted or genome-wide) not checked.
* Serum restriction (PMID 41407671) is already shipped.
