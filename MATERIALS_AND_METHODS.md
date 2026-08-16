# Materials and Methods

Draft for the PLOS ONE manuscript. Prose is written to be adapted directly; **bracketed items in bold
need a citation or a value confirmed before submission** — `python -c "from starplast import datasets;
print([d.key for d in datasets.unresolved()])"` lists them programmatically.

Figures were read from the built cache at commit `29b0bd9`.

---

## Implementation and availability

starplast is a desktop application written in Python 3.10 using PyQt6 (6.7.1) and pyqtgraph (0.13.7) for
OpenGL rendering, with pandas, NumPy, scikit-learn and umap-learn for data processing. It is distributed
with a precomputed cache (30 MB) containing all 373 columns for 8,140 genes and all 12 relation
types (296,412 edges), so the application requires neither a network connection nor any source dataset at
runtime. The cache is installed inside the package, so it is present in a wheel and resolves without
configuration. Source
code, the cache, and a notebook that reproduces every download are available at
**[repository URL]** under **[license]**.

Genes are rendered either as flat point discs or as glossy/metallic OpenGL sphere impostors. For
each fragment, the shader reconstructs a hemisphere normal from point coordinates, evaluates a GGX
microfacet response against up to eight active lights and a procedural studio environment, and
writes the curved sphere surface to the depth buffer. Scene illumination can optionally upload a
48³ point-density volume as an `R32F` texture and march 24 shadow samples per gene in the vertex
shader. This provides soft cluster occlusion rather than triangle-surface reflection or hardware
path tracing. Rendering comparisons, static-frame hashes, and full-map timings are retained under
`results/pbr_lighting_2026_08_14/` with the active OpenGL renderer recorded alongside each benchmark.

A machine-readable registry of every input dataset — its provenance, accession, the node-table columns it
produces, its coverage, and its known limitations — is included as `starplast/datasets.py` and is the
single source from which this section, the download notebook, and the application's provenance display
are generated.

## Gene universe and identifier resolution

The gene universe comprises 8,140 *Toxoplasma gondii* ME49 genes, obtained by de-duplicating an 8,227-gene
annotation table on gene identifier. Gene products, orthogroup assignments and InterPro domain annotations
were taken from ToxoDB and OrthoMCL **[cite ToxoDB release]**.

Because the literature does not use a single identifier for a gene, we implemented an explicit identifier
resolution layer. Each gene is associated with its current ME49 accession, its pre-2012 accessions, the
corresponding GT1 and VEG strain accessions, its ToxoDB symbol, and Tg-prefixed variants of that symbol.
Symbols and previous identifiers were retrieved from the ToxoDB REST API (accessed 2026-08-11). Cross-strain
accessions were mapped by numeric suffix; we validated this correspondence rather than assuming it, finding
that 99.53% of ME49/GT1 gene pairs sharing a numeric suffix and an OrthoMCL release also share an
orthogroup (VEG, 99.46%).

Ambiguity was recorded rather than resolved. A string claimed by two genes at the same confidence tier was
withdrawn from the index and retained for reporting (153 strings); across tiers, the more specific
identifier takes precedence, so a gene's own accession is never withdrawn because another gene formerly
carried it. Three precision constraints were applied to symbol matching: symbols containing no digit
require an upper-case surface form, because several valid gene symbols (HOOK, CLAMP, CLIP, SPARK, REMIND)
are also ordinary English words; tokenisation is Unicode-aware, since an ASCII-only character class
truncates accented words (French *Santé* → `Sant`, itself a gene symbol); and an accession prefix is not
re-interpreted as a symbol. *Toxoplasma* strain designations were excluded from symbol matching. The
resulting index resolves 19,368 strings plus 15,192 strain accessions.

This layer proved necessary rather than precautionary: one published in vivo CRISPR screen cites pre-2012
accessions for every gene it reports and contributed no rows at all until its identifiers were resolved.

## Subcellular localization

Compartment assignments derive from hyperLOPIT **[cite Barylyuk et al.]**, which assigns 3,827 of 8,140
genes (47.0%) to one of 26 compartments. Both the MAP and MCMC assignments were retained as separate
variables together with their posterior probabilities, because the two inference procedures disagree for
980 of the 3,827 assigned genes (26%); reporting a single assignment would conceal this.

To extend coverage we implemented orthoLOPIT, a conservative transfer procedure. For each gene lacking a
native assignment, the measured hyperLOPIT labels of its *Plasmodium falciparum* and *Cryptosporidium
parvum* orthologs were consulted **[cite Pf LOPIT; Guérin et al. 2023]**. Because compartment vocabularies
differ between species, all terms were first mapped through a curated dictionary of 87 species-specific
terms onto 12 unified categories. A label was transferred only when every available donor agreed;
disagreement was left unresolved rather than settled by majority. *Trypanosoma brucei* was excluded as a
donor: it is not an apicomplexan and its compartment set lacks the apical secretory organelles. This added
126 genes, for a total of 3,953 of 8,140 (48.6%).

Measured and transferred labels are stored in separate variables, with a provenance field recording which
applies to each gene, so that an inference cannot be mistaken for an observation. Genes with neither are
rendered as unknown rather than assigned a default compartment, because hyperLOPIT assignment tracks
protein abundance and the unassigned set is therefore biased toward low-abundance proteins.

## Transcriptomic, genetic and protein-level data

Stage-resolved transcriptomes were taken from GEO accessions GSE108740 (tachyzoite, day 3/5/7 and in vivo
tissue cyst; 12 columns; 7,739 genes) and GSE206344 (oocyst sporulation series; 6 columns; 7,974 genes).
All 18 raw FPKM columns are distributed; summary variables are log2(mean FPKM + 1). Three additional
transcriptional series already present in the shipped cache were retained as their own biological
slots rather than pooled: acute and chronic mouse-brain infection with purified bradyzoites
(PMID 31726967; 14 columns; 7,663 genes), alkaline-stress differentiation (GSE132248; 10 columns;
7,880 genes), and MORC/BFD1 perturbation (PXD058095 supplementary RNA-seq workbook; 18 columns;
7,841 genes). The first of these is FPKM transcript abundance despite residing in a mixed
transcriptome/proteome source directory; it is never treated as a fitness screen.

Eight CRISPR fitness screens were incorporated: a genome-wide in vitro screen in human foreskin
fibroblasts (Sidik et al. 2016; PMID 27594426), four in vivo composite scores (peritoneum, lung, liver,
spleen) obtained from ToxoDB, naive bone-marrow-derived macrophage and IFN-γ screens (PMID 33067458),
and a further in vivo screen covering 115 genes.

Two of those eight carry no confirmed citation and are marked as such in the registry rather than
attributed by inference. The ToxoDB in vivo composite scores are recorded against PMID 31481656, which
is also the source of the targeted in vivo platform below; whether the composite scores derive from that
same study or were computed by ToxoDB from another has not been established, and the two must not be
cited as one without checking. The 115-gene in vivo screen has no recorded publication at all.
`starplast.datasets.unresolved()` returns exactly these entries, so an unverified attribution cannot
reach a manuscript unnoticed. Four additional published screens
were incorporated from publisher supplementary material: a genome-wide screen for genes synthetically
lethal with GRA17 (7,553 genes; PMID 37498952), two targeted screens across *T. gondii* strains and mouse
subspecies (236 and 232 genes; PMID 40240328), a targeted in vivo platform (168 genes; PMID 31481656), and
a pooled single-cell screen scoring effectors by their perturbation of host transcription (252 genes;
PMID 37827122).

We emphasise that these screens share neither a sign convention nor a scale — the macrophage screens are
inverted relative to the others and standard deviations differ approximately 15-fold — and that the two
GRA12 screens are not replicates (in vivo log fold-changes correlate at *r* = 0.41). They are therefore
retained as separate variables and are not pooled. Coverage is highly uneven, and genes untested by a
targeted library are represented as missing rather than as null effect.

Protein abundance is represented by the median log2 iBAQ across replicates from two immunoprecipitation
experiments (PRIDE PXD043808 and PXD065585), covering 748 genes (9.2%). **This is enrichment, not a deep
proteome**, and should not be described as proteome-wide. A distinct total-proteome series from
AP2XII-1/AP2XI-2 perturbation (PXD039400/PXD042658) contributes 12 replicate-abundance and three
log2-fold-change variables for 3,005 genes; abundance and fold-change were normalized separately.
Oocyst developmental-stage protein abundance is represented by eight iTRAQ ratios from PXD003765
(2,079 genes), retained as ratios rather than re-centered.

Phosphosite counts, without positions, were available for 1,175 genes (14.4%); this column is inherited
from the upstream node table and its
originating publication is not recorded, so it is listed by `datasets.unresolved()` and must be
confirmed before citing. (A separate phosphosite dataset with positions and ratios, PXD017032, is
represented in the gene table by five variables: measured-site counts and median up/down ratios over
1,603 genes. Residue positions remain in the source workbooks because they are not a per-gene
quantity.) Per-gene AlphaFold model confidence (mean pLDDT) was
available for 6,480 genes (79.6%); coordinates are not distributed but are retrieved on demand from the
AlphaFold Database and cached locally.

## Relations between genes

Twelve relation types are stored separately and never merged, since they answer different questions and
disagree with one another. Ten are derived directly from data: co-mention in abstracts (435 edges,
minimum two shared abstracts) and in open-access full texts (7,733; two shared paragraphs); crosslinking
mass spectrometry (2,842 pairs over 1,630 genes) **[cite StarPath]**; immunoprecipitation of tagged baits
(64 pairs); structural similarity by Foldseek TM-align at TM ≥ 0.7 over 6,900 *T. gondii* AlphaFold models
(11,684 pairs); shared orthogroup (3,452), shared InterPro domain (10,399), shared hyperLOPIT compartment
(118,712), and correlation-based co-expression (49,293; top-25 neighbours at *r* ≥ 0.95) and co-fitness
(88,997; *r* ≥ 0.90).

Crosslink identifiers required care: the source resource uses RH88 accessions whose numbering does not
correspond to ME49 (`TGRH88_016370` is `TGME49_210408`, not `TGME49_216370`), so mapping used the alias
column distributed with the export rather than the numeric-suffix rule valid for GT1 and VEG.

## Attention correction

Raw co-mention counts reproduce the field's publication history rather than biology. We therefore report
an attention-corrected residual by default. For two genes appearing in *n*₁ and *n*₂ of *N* units, the
expected co-occurrence under independence is *n*₁*n*₂/*N*, and the residual is
log2((observed + 0.5)/(expected + 0.5)). The unit is the record for abstracts and the paragraph for full
texts, since two genes named in one paragraph plausibly stand in a relation whereas two genes named
anywhere in a long article often do not. Units naming more than 12 genes were excluded as lists or hit
tables, and per-gene counts were computed over the same population, so that observed and expected derive
from one denominator.

We additionally classify each gene by *where* it is named, a distinction read from document structure
rather than assigned as a weight: focal (named in a title; 286 genes), substantive (named in an abstract;
464), incidental (named only in a body paragraph or caption; 1,816), or absent (5,574). This distinction
is material. Adding full texts and the identifier resolution layer raised the number of genes named
anywhere from 601 to 2,566 (31.5% of the proteome), but the number named in a title or abstract rose only
from 601 to 750 (9.2%): almost all of the apparent gain consists of genes listed once in a screen's
supplementary table and never discussed.

## Derived relations

Two relation types are derived rather than observed, and are labelled as such.

**Structural holes** (255 pairs over 457 genes) are gene pairs that co-express across the stage series and
co-behave across the seven CRISPR screens, yet appear together in no abstract and no open-access
paragraph. Two constraints were necessary to prevent the measure from detecting paralogy. First, homology
(shared orthogroup or shared domain) constitutes a single evidence family and cannot serve as either of
the two required independent lines, because paralogs almost always share domains; counting them separately
made 53 of the first 66 candidates pure paralogy, and permitting homology to pair with co-expression
admitted a further 291 pairs of which 76% were same-orthogroup. Second, shared compartment was excluded
entirely, being both unspecific (118,712 edges) and abundance-dependent, which would preferentially link
well-expressed and therefore well-studied genes.

**Unwritten interactions** (2,673 pairs over 1,602 genes) are pairs measured to interact physically, by
crosslinking mass spectrometry or immunoprecipitation, that appear together in no abstract and no
open-access paragraph. This describes 2,673 of 2,906 measured pairs (92%); in 147 of these, both genes are
individually well studied.

## Visualisation

Node coordinates are a three-dimensional UMAP embedding (n_neighbors = 25, min_dist = 0.25, Euclidean
metric, random_state = 0) of a 43-column matrix: 16 numeric features (three expression summaries, seven
fitness screens, mean pLDDT, paralog number, InterPro count, phosphosite count, domain presence and
lineage specificity), z-scored following median imputation, concatenated with 27 hyperLOPIT one-hot columns
scaled by 0.5. Points are rendered with translucent blending and depth testing, and edge opacity scales
with edge weight.

## Structure search and held-out label recovery

Positions in the map are an embedding of measured features; whether that embedding organises anything
biological is a separate question, and it is answered by searching for structures that recover a label
the embedding was never given.

For a target label, every node-table column that substantially restates it is excluded from the feature
matrix before embedding. Exclusion is by measurement rather than by name: association between each
column and the target is computed (Cramér's V for categorical pairs, Pearson correlation for continuous
pairs, and the correlation ratio for mixed pairs, all on 0–1 so one threshold applies), and any column
at or above 0.8 is removed along with the target itself. Naming the columns is not sufficient — in an
earlier version `compartment` fed an embedding and its exact twin `lopit_map` together with its
derivations `lopit_unified` and `compartment_best` were duly reported as the top held-out discoveries at
V = 0.96.

Measured association is necessary but not sufficient. A label computed as a function of several columns
is a *joint* function of them: `stage_enriched_derived` is the argmax of three expression columns and
associates with each individually at only 0.56–0.66, below any workable threshold, while being fully
determined by the three together. Derived columns therefore additionally declare their sources in the
dataset registry, and a declared source is excluded together with the rest of its feature block.

Each combination of dataset blocks, missing-value policy, scaling and UMAP and HDBSCAN hyperparameters
is scored by how well the resulting clusters recover the held-out label. For each label the single best
cluster is taken, and precision and recall are reported separately and never blended, because a cluster
that is purely apicoplast while holding 5% of apicoplast proteins supports no inference. The summary
statistic is the label-size-weighted mean F1, so a structure that isolates one small class does not
outrank one that organises the proteome. Every run records its full recipe, seed, sample and the exact
set of columns excluded, so a result can be re-derived rather than trusted.

**Labels meaning "not measured" are excluded from scoring.** This is not a detail. Scored with them
included, localization reached mean F1 0.484 and its single best-recovered label was `unassigned`
(F1 0.39); excluded, the same run scores 0.207, with real compartments between 0.21 and 0.35. Because
hyperLOPIT assignment tracks protein abundance, `unassigned` is largely "too scarce to call", so a
structure separating it is separating measured genes from unmeasured ones.

A negative control makes the point unarguable. Scoring depth of literature attention — how much a gene
has been studied — gives mean F1 0.654, higher than either measured target, and essentially all of it
comes from the never-named class at F1 0.77 over 1,439 genes, while the genuine attention tiers score
0.14–0.35. The embedding encodes which genes have been measured, because a gene absent from most assays
is absent from most of the feature matrix. Any target carrying an absence class inherits that signal.

Under the corrected rule, cell-cycle phase is the target with a clean result: it has no absence class,
and every one of its labels is a real phase. A positive control behaves as required, with the derived
stage label — a deterministic function of expression columns the embedding contains — recovered well
above the measured targets.

## Limitations

Several limitations follow from the above and should be read alongside any use of the map.

Missing values are imputed to the column median before embedding, so the pattern of *which genes were
measured* is partially encoded in node position: genes lacking fitness data lie 0.38 map-radii from those
possessing it (mean pLDDT 0.21; phosphosite count 0.22). The structure search quantifies the consequence:
a clustering recovers "named in no paper at all" at F1 0.77, which is higher than it recovers any
biological label, so this is the dominant organising signal in the map rather than a marginal one. Since phosphosite counts are missing for 85.6% of
genes, that variable functions largely as an indicator of inclusion in a phosphoproteomics experiment.
Relatedly, the hyperLOPIT block contributes only 1.1% of the feature matrix's variance and therefore has
little influence on position despite being included. UMAP preserves local neighbourhoods rather than global
distances, so proximity is interpretable while inter-cluster separation is not; both derived relation types
are computed from edges rather than from embedding distance and are unaffected.

Predicted complexes for crosslinked pairs are frequently inconsistent with the measurement they are meant
to explain: of 2,397 scored pairs, 60% place no crosslink within reach, the median interface ipTM is 0.16,
and only 162 pairs are both confident and crosslink-consistent. Many failures involve dense granule
proteins, which are largely disordered. The crosslink is the measurement; the model is a hypothesis about
the pose.

Finally, the two literature sources are different populations. Abstracts cover the field, whereas full
texts comprise only those articles deposited open access, so full-text coverage cannot be quoted as
coverage of *Toxoplasma* research generally. The full-text corpus grew during development (5,493 to 6,667
articles); the build log records the count actually used.

## Data and code availability

All input datasets, their accessions and their download locations are enumerated in
`starplast/datasets.py`. The notebook `notebooks/download_datasets.ipynb` reproduces every download from
that registry.

**[Add repository DOI / Zenodo archive.]** This is the one placeholder that cannot be filled from the
project: minting a DOI requires an account and an archived release, and is the author's to do.
