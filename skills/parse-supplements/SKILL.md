---
name: parse-supplements
description: Extract text and gene identifiers from published supplementary files — PDF, DOCX, XLS/XLSX, TXT — when a paper's hit list is not in a machine-readable table. Use when supplementary data must be mined from documents, when accession extraction from spreadsheets has returned nothing, or when a study lists proteins by name rather than by ID. Covers Toxoplasma, Plasmodium and Cryptosporidium accession formats.
---

# Parsing published supplementary files

Supplementary data arrives in whatever the publisher felt like. In one real corpus of 97 tagged-protein
interaction studies, only 32 shipped accessions in a spreadsheet; the rest published their hit lists as
PDF or DOCX, or by gene name alone. This skill reads all of those uniformly.

## Run it

```bash
python scripts/extract.py <dir-or-file> --organism toxo --out found.tsv

# also match proteins named by symbol, and keep the extracted text for inspection
python scripts/extract.py studies/ --symbols toxodb_identity.tsv \
       --text-dir /tmp/extracted --out found.tsv
```

Output is a TSV of `(file, kind, identifier)` where `kind` is `accession` or `symbol`.

No pip installs. It uses `pdftotext` (poppler), Python's `zipfile` for DOCX, `libreoffice --headless`
for legacy `.doc`/`.xls`, and pandas for modern spreadsheets.

## What matters when doing this

**Use `pdftotext -layout`, never the default.** Without it a two-column table is reflowed into prose and
the columns interleave, so every extracted row is a blend of two different rows. The corruption is
invisible unless you read the output — the text looks fine, the data is wrong.

**DOCX is a zip of XML; convert cell and row boundaries before stripping tags.** Replace `</w:tc>` with a
tab and `</w:tr>` with a newline first. Strip tags without doing that and every cell in a row runs
together, which turns a gene/description table into one unusable line.

**Accessions get mangled by typesetting.** `TGME49_208830` appears as `TGME49-208830`, `TGME49.208830`,
`TGME49 208830` and `TGME49208830` depending on the journal's line-breaking. Normalize separators before
matching, and re-insert the underscore when it has been dropped entirely.

**A digit-free symbol must appear in upper case to count.** `HOOK`, `CLAMP`, `CLIP`, `SPARK` and `REMIND`
are all real *Toxoplasma* gene symbols and all ordinary English words. Matching them case-insensitively
fills the output with prose. Symbols containing a digit (`GRA16`) are safe either way.

**Tokenise Unicode-aware.** An ASCII-only character class truncates accented words mid-token — French
*Santé* becomes `Sant`, which is a real gene symbol — and manufactures matches from nothing.

## After extraction

Route every identifier through a gene-identity layer before use. Papers cite whatever accession was
current when they were written: in the same corpus, `TGGT1_` strain accessions were **more common than
`TGME49_`**, and one 2019 screen used pre-2012 `TGME49_0xxxxx` ids for every gene it reported, which
resolve to nothing without the mapping.

## What this does not do

It extracts *identifiers mentioned in a document*, which is not the same as a hit list. A supplement is
usually the paper's complete quantification table — in that corpus the median parsed study named 754
genes and the largest named 7,866, essentially the whole proteome. Deciding which rows are hits needs the
paper's own threshold column, and that is named differently in every paper. Treat the output as
membership and curate from there; treating it as interactions would invent tens of thousands of them.
