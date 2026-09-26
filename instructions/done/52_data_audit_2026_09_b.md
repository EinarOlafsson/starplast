# The second September 2026 data audit: what the Plasmodium table was still missing

**Status: complete (2026-09-26).**

Instruction 50 finished earlier the same day, and it left a question behind: it filled Toxoplasma
slots well and Plasmodium slots where it found data, but the Plasmodium arm still had no
localization, no protein abundance for the stage the parasite lives in, no RNA half-life, no field
variation and no resistance data — the questions that organism's literature has been answering for
twenty years. This pass went looking for those, found six datasets, refused five more, and wrote
down why each time.

Everything here is Plasmodium. That is the finding, not the plan: the Toxoplasma side was swept
again with the same queries and produced nothing new that could be verified. Four Toxoplasma
candidates and one host candidate were read and refused, and they are listed below with their
reasons.

## What was searched, and how

Every query is written out so the sweep can be repeated rather than trusted.

**NCBI GEO** (E-utilities `esearch` against `db=gds`, then `esummary`, then the series' own FTP
supplementary directory):

    Toxoplasma gondii[Organism] AND ("2025/01/01"[PDAT] : "2026/12/31"[PDAT])          640 hits
    (Toxoplasma[All Fields]) AND (Homo sapiens[Organism] OR Mus musculus[Organism])
        AND ("2023/01/01"[PDAT] : "2026/12/31"[PDAT]) AND (gse[Entry Type])             50 hits

**Europe PMC** (`/search?query=...&resultType=core`), title-restricted because a free-text term for
an assay returns reviews. The queries that produced the accepted datasets:

    (TITLE:"Plasmodium falciparum" AND TITLE:"proteome") AND (FIRST_PDATE:[2023-01-01 TO 2026-12-31])
    (TITLE:"Plasmodium" AND (TITLE:"mRNA decay" OR TITLE:"mRNA stability" OR TITLE:"half-lives"
        OR TITLE:"half-life" OR TITLE:"transcriptional dynamics"))
    ("Plasmodium falciparum") AND (TITLE:"chemogenomic" OR TITLE:"resistome"
        OR TITLE:"in vitro evolution" OR TITLE:"drug resistance selection")
    (TITLE:"Plasmodium falciparum" AND (TITLE:"ATAC" OR TITLE:"chromatin accessibility"
        OR TITLE:"regulatory landscape"))
    (TITLE:"Toxoplasma" AND (TITLE:"proteome" OR TITLE:"proteomic" OR TITLE:"spatial proteome"))
        AND (FIRST_PDATE:[2024-01-01 TO 2026-12-31])
    (TITLE:"Toxoplasma" AND (TITLE:"turnover" OR TITLE:"half-life" OR TITLE:"proteome dynamics"
        OR TITLE:"protein degradation"))
    ("Toxoplasma gondii") AND (TITLE:"THP-1" OR TITLE:"monocyte" OR TITLE:"iPSC"
        OR TITLE:"induced pluripotent") AND (FIRST_PDATE:[2022-01-01 TO 2026-12-31])
    ("Toxoplasma") AND (TITLE:"brain" AND (TITLE:"proteome" OR TITLE:"proteomic"
        OR TITLE:"transcriptome")) AND (FIRST_PDATE:[2022-01-01 TO 2026-12-31])

Every identifier was then resolved rather than typed: `EXT_ID:<pmid> AND SRC:MED` with
`resultType=core` for each PMID, `api.crossref.org/works/<doi>` for each preprint, and
`ebi.ac.uk/pride/ws/archive/v3/projects/<accession>` for each PRIDE deposit. Supplementary files
came from `europepmc.org/…/<PMCID>/supplementaryFiles`, and bioRxiv's own
`/content/biorxiv/early/<date>/<doi>/DCn/embed/media-n.xlsx` for the preprints — the route
instruction 50 established, still working, including the changed `10.64898` DOI prefix.

Three things about the search route worth keeping:

* **Europe PMC's supplementary endpoint serves open-access articles only.** `PMC13531851` — the Cell
  paper whose GEO series 50 took for mRNA decay — answers `Article with id … is not open access
  one`, which is why the genome-wide crowding screen in that same paper could not be taken.
* **A PRIDE record can disagree with the paper's own table.** PXD079493's description says 133
  Hsp90-dependent hits where the supplementary table lists 131. The reproducible number ships.
* **Two of the strongest candidates were already on disk.** `datasets/incoming_2026_09_25/` held the
  Hsp90 chemoproteomics and the M2K1 seasonal atlas, downloaded during instruction 50 and never
  processed. One is accepted here and one refused; a folder in `incoming/` is not a decision.

## What was added

Six datasets, 22 columns, all Plasmodium -- seven deposit entries, because the spatial proteome's
variation table is registered separately from its localization. `starplast/deposits.py` derives
each one and `notebooks/derive_deposits_2026_09.ipynb` (sections 13-19) re-derives every number
below from the raw file. The table goes from **146 to 168 columns**, and slot coverage from **169/215 to 180/221**
genes-unit slots: five slots that were empty are answered, and six new ones exist because the
catalogue could not express what these datasets measure.

| dataset | what it adds | verified by |
|---|---|---|
| **Spatial proteome of the schizont** (PMID 42218142, Nat Commun 2026;17:6192) | `lopit_pf_location`, `lopit_pf_svm_score` — the FIRST subcellular localization in the Plasmodium table | 3,000 proteins mapped, 1,646 classified, 24 niches, all three the paper's own; RAP1 in the rhoptries, MAHRP1 in the Maurer's clefts, ACP in the apicoplast, EXP2 and HSP101 at the PVM, GAPDH in the cytosol |
| **Field and between-species variation**, same paper's Supplementary Data 3 | `field_pnps_adj`, `field_variant_fraction` (the empty field-variation slot) and `dnds_laverania`, `dnds_plasmodium` (a new slot) | 5,232 genes; field pN/pS agrees with PlasmoDB's laboratory-strain SNP ratio at rho 0.29, so it is a different population and not the column the table already had |
| **mRNA synthesis and decay rates** (PMID 29985403, Nat Commun 2018;9:2656) | `mrna_decay_rate_4tu` (the empty RNA-stability slot), `transcription_rate_4tu` (a new slot) | the paper's claim of transcription in every stage: 616–962 genes peak in each of the six windows; 4,373 of 4,373 synthesis rates positive, 4,401 of 4,420 decay rates negative; of the genes whose transcription peaks in a ring window, 87–92% also peak in the shipped ring expression column |
| **Blood-stage proteome and Hsp90 dependence** (bioRxiv 10.64898/2026.08.28.747854, PXD079493) | `proteome_blood_log2` and the two inhibitor responses (two empty slots), `hsp90_dependent` (a new slot) | the paper's 131 dependent proteins reproduce EXACTLY from its stated rule, and all 131 are the same proteins; 124 survive gene mapping |
| **RNA dependence** (PMID 38355719, Nat Commun 2024;15:1365) | `rna_dependent`, `rna_dependence_qvalue` — a new slot | 898 of 3,671, the paper's own number, from the deposit's q-values; RNA helicases enriched (28 of 61, odds 2.7, p 2e-4) and the proteasome 0 of 14 |
| **Sexually committed proteome** (PMID 41482054, Mol Cell Proteomics 2026;25:101505) | `committed_vs_asexual_log2fc`, `committed_vs_asexual_fdr` — a new slot | MSP1 +0.34 at FDR 0 and MSP2 +0.67, the paper's merozoite-surface finding; MSRP1 at +1.68 is only the sorting marker |
| **In vitro evolution resistome** (PMID 39607932, Science 2024;386:eadk9893) | `resistance_target_compounds` and three selection counts (the empty resistance slot), `pf6_field_dnds`, `pf6_field_nonsyn_snvs` | 724 clones and 118 compounds, both exact; the classification's top is PfATP4 and PfMDR1 at 5 compounds, the prodrug-activating esterase 4, then cytochrome b, CARL, PI4K beta and NPC1 at 3 |

Six new slots, each with data: `localization · schizont`, `chaperone dependence`,
`RNA dependence of complexes`, `transcription rate · asexual blood stage`,
`protein abundance · sexually committed`, `between-species selection (dN/dS)`.

### Three decisions inside those datasets that changed a number

* **The label had nowhere to go.** Every deposit before this one produced numeric columns, and
  `parasite_columns` merged them with `mean(numeric_only=True)` — which drops a text column in
  silence. `lopit_pf_location` is a compartment name, so the merge now carries label columns too,
  and where two accessions resolve to one gene the label survives **only if they agree**: a
  compartment has no mean, and voting for the first row would attach a location to a gene nobody
  measured there.
* **A count of mutations is not a resistance gene.** Ranked by how many compound selections mutated
  it, the resistome's top gene is AP2-G with 13, followed by PfEMP1 — genes that mutate under
  prolonged culture whatever the drug, because losing gametocyte production is cheap in a flask.
  Shipping only the counts would have made AP2-G the parasite's leading resistance gene. The
  paper's own classification ships beside them and is what answers the slot.
* **The atovaquone gene nearly vanished.** The resistome names the three mitochondrial-genome genes
  by pre-2010 identifiers (`mal_mito_3`), which match no accession pattern and were being dropped —
  one of them is apocytochrome b, with 32 selected clones. Each is mapped on the deposit's own gene
  description matching exactly one product in the shipped table, written out in
  `PF_RESISTOME_LEGACY` rather than pattern-matched.

## What was refused, and why

Five datasets were read in full and not shipped. None of these is a gap in the search; each is a
decision about what the source can support.

* **The M2K1 seasonal single-cell atlas** (bioRxiv 10.1101/2025.04.14.648697, downloaded, parsed).
  889 genes differentially expressed between asymptomatic dry-season infections and clinical
  malaria, tabulated per 10 hours-post-invasion bin. Two things stop it: only values passing the
  paper's threshold appear, and **45% of the genes change sign between bins**, so a per-gene mean
  would average opposite directions. The paper assigns an overall direction to 42 of the 889, and on
  those 42 the mean sign computed here agrees with it in every case — which is the evidence that the
  parse is right and the summary is still not one the source makes.
* **The chemogenomic piggyBac screen** (PMID 37067430): 553 isogenic mutants against
  dihydroartemisinin and bortezomib, which would have filled the empty drug-sensitivity slot. Its
  headline could not be reproduced: the paper's Figure 3 describes 29 drug-sensitive mutants
  (log2FC < 0, P < 0.05 in at least one screen, with growth intact in the no-drug control), and the
  deposit gives **40** under that rule on the raw p-value and 12 on the adjusted one. A dataset
  whose stated result cannot be recomputed does not ship.
* **The temperature ATAC-seq** (PMID 41721429), which would have filled the empty chromatin-
  accessibility slot. The accessible regions are genomic intervals: 1,089 peaks with coordinates and
  no gene, and the only gene-keyed sheet is 186 filtered peak–gene pairs. Assigning peaks to genes
  needs a P. falciparum annotation this tree does not carry offline, and the VEuPathDB service
  answers 401 by design — which is not something to work around.
* **The TgPRO RIP-seq** (GSE327264, from PMID 42580337, whose GEO series title states that TgPRO
  binds transcripts encoding metabolic proteins). The deposit is raw input and IP counts, so the
  enrichment has to be computed here, and two things came out wrong. The RNA-binding-deficient point
  mutant shows **more** apparent enrichment than wild type (median +0.64 against −0.01, and the two
  agree at rho 0.61), and the series' own claim does not reproduce: among the top decile of
  enrichment the share of genes carrying an EC number is 0.198 against 0.190 elsewhere (odds 1.05,
  p 0.66). The paper's own processing is behind a paywall, so there is nothing to check the
  normalization against.
* **The cross-lineage m5C methylome** (PMID 40830525), which would have thickened the Toxoplasma RNA
  modification slot. Its supplements are differential methylation between the three clonal lineages;
  there is no per-gene methylation level for a reference strain, which is what the slot asks for.

Three more were read and set aside without a full attempt: the **culture-adaptation proteome** (PMID
41036225) publishes per-strain MaxQuant dumps of up and down lists rather than a per-gene abundance;
the **BKI-1708 affino-proteomics** (PMID 42654955) is a drug-affinity pulldown of 453 parasite and
831 human proteins, so it cannot answer the host-proteome slot its title suggests; and the **early
monocyte response** (PMID 41246354), which would have answered two empty human-monocyte host slots,
ships one PDF data sheet and no machine-readable table. A search for a mouse-brain host proteome or
transcriptome during chronic infection returned nothing genome-scale with a per-gene table.

## Leakage: the new columns were tested, and nothing had to be closed

`scripts/leakage_audit.py` was re-run on both tables with eight new Plasmodium held-out targets
(`results/leakage_audit_2026-09-26/`), and seven same-quantity families were declared in
`starplast/search.py` before it ran:

    Hsp90 inhibition response        the hit call is a function of the two responses
    4-thiouracil mRNA dynamics       synthesis and decay come out of one model fit
    field variation                  two papers reading overlapping MalariaGEN releases
    between-species selection        one estimator on two nested ortholog sets
    in vitro evolution selection     one resistome: counts and the call made from them
    RNA dependence                   a flag and the q-value it is thresholded from
    sexual commitment proteome       a contrast and its FDR
    localization (extended)          the spatial proteome's label and its SVM confidence

**Closure gaps: 0. Residual leaks: 0**, on 161 measured Plasmodium columns over 12,171 pairs and on
395 Toxoplasma columns over 76,218. Between columns from different axes, Plasmodium rank association
reaches 0.586 at the 99th percentile, so the closure's 0.8 threshold is well outside what unrelated
columns do in this table.

Two numbers from that run are worth reading rather than filing:

* **Ribosome footprints predict the new proteome at 0.686**, the highest permitted predictability in
  the table — and it is not a leak (robust z 2.2). Footprint density predicting protein abundance is
  the biology that ribosome profiling exists for, measured on different molecules in different
  laboratories.
* **Field pN/pS is predicted by the laboratory-strain SNP columns at 0.407.** That is the decision
  to keep them in separate families, measured: related, not a restatement.

## Two things a reader should not take on trust

* **`mrna_decay_rate_4tu` is not a half-life.** It is transcripts lost per minute at the gene's own
  peak, so it tracks abundance (rho 0.42 against blood-stage expression) and an abundant stable
  transcript can out-decay a scarce unstable one. It is in a different unit from the Toxoplasma
  actinomycin columns that answer the mirrored slot, which is why that slot's policy is `separate`
  and not `average`.
* **`rna_dependent` is not 'binds RNA'.** RNA helicases are enriched and the proteasome is absent,
  as they should be — but ribosomal proteins are *depleted* (23 of 137), because a ribosome is
  mostly RNA and does not come apart in this assay the way an mRNA-bridged complex does.

## One thing this audit found and did not change

**The Plasmodium layout is built from 20 of 168 columns.** `scripts/rebuild_layouts.py` was re-run
and produced coordinates bit-identical to the shipped ones, with only the recorded node-file hash
changing. The reason is `embedding.BLOCKS`: its patterns and its slot blocks are built from the
*Toxoplasma* catalogue, so the Plasmodium table contributes whatever happens to match — the
expression series and a few protein features. Twenty-two new columns therefore cannot move that map,
and neither could the twenty-three instruction 50 added. Fixing it means deciding what the
Plasmodium embedding should be built from, which changes the shipped map and is not an audit's call.

## Where things are

* `starplast/deposits.py` — the six derivations, each with its check named, and the label-carrying
  merge rule.
* `notebooks/derive_deposits_2026_09.ipynb` sections 13–19 — executed, so the notebook in the tree
  is the run that produced the numbers above.
* `starplast/datasets.py` — seven registry entries (the spatial proteome and its variation table are
  separate entries, so holding out a localization does not withdraw a dN/dS).
* `scripts/generate_slot_table.py` — six new Plasmodium slots, five new pattern lists for slots that
  were empty, and five references resolved through Europe PMC.
* `tests/test_deposits.py` — 48 tests now, 14 of them new: the planted-label rule, every headline
  number above, and the families.
* Raw deposits under `datasets/<level>/<kind>/<PMID or accession>/`, each with `URLS.txt` and
  `SHA256SUMS.txt`.
