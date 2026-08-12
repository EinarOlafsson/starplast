---
license: cc-by-4.0
task_categories: [tabular-classification]
tags: [toxoplasma, proteomics, interactome, bioid, ip-ms]
---

# Toxoplasma gondii tagged-protein interaction studies

Which genes appear in the supplementary tables of 97 published proximity-labelling
(BioID / TurboID / APEX) and pulldown (IP-MS / co-IP) studies with a tagged *Toxoplasma gondii* protein.
Studies were identified by screening 33,924 PubMed abstracts; supplementary files were retrieved from
the publishers.

## This is membership, not interaction

A study's supplement is usually its **complete quantification table**, not its hit list. Median genes per
parsed study is 754; eleven list more than 2,000 and the largest lists 7,866 — essentially the whole
proteome. Treating a row here as an interaction would manufacture tens of thousands of false edges.
Converting membership to interactions needs per-paper curation of which sheet and column mark enrichment.

## Contents

| file | rows | what |
|---|---|---|
| `study_gene_membership.parquet` | 67,099 | one row per (study, gene) found in a supplement |
| `studies.parquet` | 97 | per-study metadata, license, and whether it parsed |

Gene identifiers are resolved to current ToxoDB ME49 accessions. `TGGT1_` accessions turned out more
common in these files than `TGME49_`, so identifiers pass through an identity layer that maps strain and
pre-2012 accessions forward.

## Licensing

Only derived facts are redistributed here. Of the source articles, 51 carry a confirmed CC-BY
license; the remainder are either non-commercial, ambiguous, or have no locally readable license block,
and their **raw files are not mirrored**. Every row carries its source PMID so the original is one click
away. Cite the original studies, not this table.

Produced by [starplast](https://github.com/EinarOlafsson/starplast).
