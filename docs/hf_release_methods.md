# Methods — Toxoplasma gondii tagged-protein interaction studies

The methods for the derived dataset published at
`huggingface.co/datasets/einarolafsson/toxoplasma-interaction-studies`. This file ships inside the
release as `METHODS.md` and is kept here in the repository, so the two cannot drift apart.

Reproduced end to end by [`notebooks/build_hf_release.ipynb`](../notebooks/build_hf_release.ipynb).

## What the dataset is

Which genes appear in the supplementary tables of 97 published proximity-labelling (BioID, TurboID,
APEX) and pulldown (IP-MS, co-IP) studies that used a tagged *Toxoplasma gondii* protein.

68,085 rows, one per (study, gene). 97 studies, of which 51 carry a confirmed CC-BY licence.

**It is membership, not interaction.** This is the single most important thing about the table and
the reason it is published in this form rather than as an interaction network. See below.

## How studies were identified

1. Screen 33,924 *Toxoplasma* PubMed abstracts (the corpus described in the main README) for
   proximity-labelling and pulldown vocabulary — BioID, TurboID, APEX, miniTurbo, IP-MS, co-IP,
   immunoprecipitation followed by mass spectrometry.
2. Keep those where the tagged protein is a *Toxoplasma* protein rather than a host protein, since a
   host-baited pulldown answers a different question.
3. Record the PMID, PMCID, year and method for each. The result is
   [`notebooks/catalogue_bioid_ipms.tsv`](../notebooks/catalogue_bioid_ipms.tsv), which is committed
   so the screen is inspectable rather than a number in a sentence.

## How supplements were retrieved and parsed

Supplementary files were fetched from the publisher for each study with a PMCID. Formats that could
not be read are recorded per study in `unreadable_formats` rather than dropped silently — a study
with zero genes and no recorded reason is indistinguishable from a study with genuinely no hits.

Tables were parsed from `.xlsx`, `.xls`, `.docx` and `.csv`. Every sheet of a workbook is read, not
only the first: the relevant table is frequently not sheet 1.

## How gene identifiers were resolved

Identifiers pass through the same identity layer the map uses:

- `TGGT1_` and `TGVEG_` strain accessions map forward to current ToxoDB **ME49** accessions.
  `TGGT1_` turned out to be *more* common in these supplements than `TGME49_`, so a pipeline that
  only recognised ME49 accessions would silently lose most of the data.
- Pre-2012 accessions map forward through the identity table.
- Where an identifier is ambiguous, the ambiguity is recorded rather than resolved by guessing.

## The membership caveat, stated properly

A study's supplement is usually its **complete quantification table**, not its hit list. The median
parsed study contributes 599 genes; eleven contribute more than 2,000, and the largest contributes
7,866 — essentially the entire proteome.

Treating a row in this table as an interaction would manufacture tens of thousands of false edges.
The rows say *this gene appears somewhere in the supplementary data of this study*, which is a fact
about a document, not about a protein.

Converting membership into interactions requires per-paper curation of which sheet and which column
mark enrichment, and what threshold that paper applied. That work is not done here, and this dataset
should not be used as though it were.

## Licensing

Only derived facts are redistributed. Of the 97 source articles, 51 carry a confirmed CC-BY licence;
the rest are non-commercial, ambiguous, or have no locally readable licence block. **No raw
supplementary file is mirrored**, whatever its licence — the dataset carries the extracted membership
facts, and every row carries its source PMID so the original is one click away.

Cite the original studies. This table is a finding aid, not a source.

## Files in the release

| file | rows | what |
|---|---|---|
| `study_gene_membership.parquet` | 68,085 | one row per (study, gene) found in a supplement |
| `studies.parquet` | 97 | per-study metadata, licence, whether it parsed, and why not |
| `METHODS.md` | — | this document |
| `build_hf_release.ipynb` | — | the notebook that reproduces both tables |
| `README.md` | — | the dataset card |
