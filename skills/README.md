# skills

Reusable techniques this project had to work out, packaged so the next one does
not rebuild them. Each folder is a Claude Code skill: a `SKILL.md` describing
when to reach for it and what goes wrong, plus any scripts.

These are also installed under `../.claude/skills/` so they are available to a
session automatically. The copies here travel with the repository.

| skill | for |
|---|---|
| `parse-supplements` | pulling text and gene identifiers out of published supplementary files — PDF, DOCX, XLS/XLSX — when the hit list is not in a machine-readable table |

## What makes something worth saving here

Not "code that worked". A skill earns a folder when the *knowledge* is the
expensive part and would otherwise be rediscovered by trial and error:

- `pdftotext` needs `-layout`, or a two-column table reflows into prose and
  every extracted row silently blends two different rows
- DOCX is a zip of XML; convert `</w:tc>` and `</w:tr>` to tab and newline
  before stripping tags or a whole table collapses into one line
- a digit-free gene symbol must appear in upper case, because HOOK, CLAMP,
  CLIP, SPARK and REMIND are all real symbols and all real English words

None of those is obvious, each cost time, and each produces output that looks
correct while being wrong.

## Candidates not yet packaged

Techniques used here that would generalise, if a second project needs them:

- **identity resolution for published gene lists** — papers cite whatever
  accession was current when written; strain accessions may outnumber the
  reference ones; ambiguity must be recorded rather than guessed
  (`starplast/identity.py`)
- **the circularity guard for cluster enrichment** — naming the columns an
  embedding used is not enough, because a renamed or derived column defeats it;
  association with the inputs has to be measured (`starplast/clustering.py`)
- **normalising by quantification type** — the range decides, not the filename;
  the same series can be raw FPKM in one file and already-logged in another
  (`starplast/sources.py`)
