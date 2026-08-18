# Build a provenance-first pan-Apicomplexan dataset archive — MACHINERY DONE 2026-08-18

**Status: the archive and its discovery are built, tested and run. The catalogue is 248 sources and
grows by re-running rather than by editing.** `starplast/archive.py`, `scripts/build_archive.py`,
30 tests, 100% coverage.

## Discovery is resolved, not remembered

The only thing checked in is which hostname serves which clade — the one fact about the WEB rather
than about the data. The release, the organisms and the files are all read from each site's own index
at run time, so re-running next year finds next year's release without an edit. A release number
typed into a URL is correct until the next release, and then it is a 404 that looks like a missing
dataset.

Resolved on 2026-08-18: **211 organism directories across four sites**, yielding **248 gene-level
sources** — Plasmodium and other haemosporidia 118, Toxoplasma/Neospora/Eimeria/Sarcocystis/Besnoitia
64, Cryptosporidium/Cyclospora/Gregarina 36, Babesia/Theileria 30. That covers every genus this
instruction names except those with no VEuPathDB presence.

## Three properties, each because its absence has cost time

* **Nothing is overwritten.** A file already on disk is recorded as `present` and left exactly as it
  is. The archive is a store of what was retrieved and when, and a re-run that silently replaced a
  2019 download with a 2026 one would destroy the provenance it exists to keep.
* **A failure is recorded, not omitted.** `unreachable`, `restricted` and `absent` are three
  different facts and each goes in the manifest. An archive listing only what worked sends the next
  session to rediscover the same dead ends — the same service the empty slots' `blocked_by` verdicts
  perform in the slot atlas. A source declared controlled-access is not even attempted, because
  retrying a licence does not change it and reporting a network failure sends somebody to debug the
  wrong thing.
* **Every file is opened, not sized.** sha256, byte count, and a shape read by actually parsing:
  sequences for FASTA, feature rows for GFF, columns for tables, members for archives. A 404 page
  saved as `.csv` is a few kilobytes of valid HTML, and every archive that checked only sizes has one.

## Not integrated into the graph, deliberately

As this instruction requires. Integrating it now would quietly make Toxoplasma claims out of
Cryptosporidium files, and the whole point is a future-ready store for extending the leakage-aware
slot hierarchy to another parasite without repeating discovery.

## What remains, and it is acquisition rather than engineering

The 248 sources are sequence and annotation — the layer everything else is keyed against, and
correctly the first thing archived. The other evidence families this instruction lists (proteomics
and PTMs, localisation, essentiality, interactions, host response, life-stage timing) live in GEO,
PRIDE, ArrayExpress and supplementary tables, and each needs a resolved accession per study rather
than a directory walk. The machinery takes them unchanged — a `Source` names its repository, licence
and scope — so that is a curation campaign against a working archive, not more building.

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
