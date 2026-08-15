# Toxoplasma dataset gap audit — 2026-08-14

## Outcome

The rebuilt cache contains 8,140 genes and 295 node columns. Every cache column belongs to exactly
one leaf slot; there are no unassigned or multiply owned columns. The 112-slot Toxoplasma catalogue
now has 47 well-covered, 8 partial, 12 thin and 45 still-empty slots. Empty relationship slots
often require pair-shaped data and therefore cannot be filled honestly by adding a gene-level summary.

## Newly downloaded and integrated

| accessions | biological slot(s) | cache coverage |
|---|---|---:|
| [GSE99395](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99395) | intracellular/extracellular RNA, ribosome occupancy, relative translation efficiency | 7,437 genes |
| [GSE129869](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129869) | host-context RNA, ribosome occupancy, relative translation efficiency | 7,880 genes |
| [GSE19092](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19092) | synchronized tachyzoite cell-cycle transcription | 6,504 cache genes |
| [GSE51780](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE51780) | merozoite and tachyzoite transcription | 6,504 cache genes |
| [GSE168155](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168155) | transcription under CPSF4/RNA-processing perturbation | 7,869 genes |
| [GSE200962](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200962) | transcription under restriction-checkpoint/cyclin perturbation | 7,880 genes |
| [GSE253884](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253884) + [GSE253885](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253885) | in-vivo fitness of hyperLOPIT-unassigned proteins | 551 genes |

## Evidence audit

The accession records and primary papers were checked as assay evidence, not merely as keyword
matches. The integration preserves parasite stage, host context, assay type and perturbation in the
column names and slot facets. It does not infer strain, host-cell, time-point or causal mechanism where
the deposited processed table does not establish one. In particular:

- RPF, RNA abundance and their within-sample relative translation-efficiency contrast are three
  related measurements, not interchangeable replicates.
- Extracellular tachyzoites in GSE99395 represent that assay's lytic-stress context; they do not stand
  in for a developmental bradyzoite measurement.
- GSE168155 measures the transcript response to CPSF4 depletion. It supports an RNA-processing
  perturbation slot, not a direct per-gene RNA-modification slot.
- GSE253884/5 are targeted libraries of proteins lacking a hyperLOPIT assignment. Their untested genes
  remain missing, and copied comparator-screen columns are not imported a second time.
- The legacy microarrays use probe-to-old-model mappings followed by the explicit ToxoDB previous-ID
  resolver. Ambiguous or unresolved probes are excluded rather than guessed.

These are observational feature integrations. They do not by themselves establish gene function or
subcellular localization; those proposed outcomes require the independent holdout validation reported
by the GPU walk.

The downloaded-file manifest and checksums are in
`/mnt/firecuda2/Claude/toxoplasma_projects/datasets/toxoplasma_acquisition_2026_08_14/MANIFEST.md`.

## Public datasets that can fill remaining slots, but need a second integration layer

| slot gap | public source | why it was not forced into the gene table in this pass |
|---|---|---|
| histone marks / chromatin state | GSE87834 and GSE104347 | processed peaks/tracks are available, but they are coordinate-shaped and use older assemblies; gene overlap requires an explicit assembly/feature mapping |
| chromatin accessibility | GSE222832 (ATAC-seq) | broadPeak and bigWig outputs are available; promoter/gene assignment and genome-version provenance must be defined first |
| direct RNA modification | GSE294543 (m5C) | GEO sample metadata describes methylated-RNA sites, but no series-level processed supplement is exposed; source retrieval and gene/site aggregation remain |
| bradyzoite subtypes | GSE311669 | a raw single-cell matrix is public; a cluster-aware pseudobulk representation is needed rather than treating cells as gene columns |
| TF binding | GSE106864 and related ApiAP2 ChIP series | binding peaks must be mapped to promoters/genes without collapsing factor and condition |
| newer high-resolution translation | GSE302107 | raw sequencing is public but no processed series supplement is supplied; reproducing the published gene-level processing is a separate workflow |
| PTM classes (acetylation, lactylation, nitrosylation, ubiquitination, glycosylation) | PRIDE accessions listed in the slot table | PRIDE commonly holds instrument/search outputs; a verified processed, site-to-gene supplement is required for each modification before integration |

## Interpretation rules retained

- GSE168155 is a response to depletion of an m6A-reading/polyadenylation component. It does not fill
  the direct per-gene m6A slot.
- GSE99395 extracellular parasites are a lytic-stress context, not bradyzoites.
- GSE253884/5 tested a selected library. Missing genes are **not** assigned a zero fitness value.
- Microarray probes from GPL7186 are resolved through explicitly registered legacy ToxoDB gene models;
  unresolved or ambiguous probes are dropped.
- Coordinate tracks and pairwise relations stay out of the node table until their mapping is explicit.
