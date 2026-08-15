# Build a provenance-first pan-Apicomplexan dataset archive — open

## State

Starplast currently integrates Toxoplasma-shaped functional, expression, localization, interaction,
phenotype, and sequence/structure evidence. Future versions should be able to extend the same
leakage-aware slot hierarchy to other Apicomplexan parasites without repeating dataset discovery.

## What to do

1. Search primary repositories and literature for downloadable gene-level datasets from any
   Apicomplexan parasite, prioritizing Plasmodium, Cryptosporidium, Eimeria, Babesia, Theileria,
   Neospora, Sarcocystis, Cyclospora, Gregarina, and related experimentally studied genera.
2. Cover as many slot-compatible evidence families as practical: transcript abundance and response,
   proteomics and PTMs, localization and organellar targeting, essentiality/fitness, interactions,
   host response, life-stage timing, comparative genomics, and structural/sequence annotations.
3. Download public source files into
   `/mnt/firecuda2/Claude/toxoplasma_projects/datasets/` without overwriting existing data.
4. Keep every source in a species/study directory with a machine-readable manifest recording DOI or
   accession, repository URL, download URL, retrieval date, checksums, file sizes, license/access
   status, biological scope, and likely Starplast slot family.
5. Verify archives and table readability. Record inaccessible, controlled-access, missing, or
   supplement-less studies in the manifest rather than representing them as downloaded.
6. Do not integrate these pan-Apicomplexan files into the current Toxoplasma graph; this task builds
   a future-ready source archive, as explicitly requested.

## Done when

- Discovery is documented and reproducible from a checked-in catalog/downloader.
- Every successful download has a checksum and provenance record, and existing files are preserved.
- The archive spans multiple genera and multiple non-leaking biological evidence families.
- Failed or gated sources are explicit, so a later session does not rediscover the same dead ends.
