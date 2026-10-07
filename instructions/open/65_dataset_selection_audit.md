# 65 · Audit dataset choices and retain a reproducible selection rule

Status: OPEN, 2026-10-07. User-requested priority alongside instruction 64.

User rule: when multiple datasets answer the same information slot, prefer the
publication with the highest citation rate per year, with a slight preference
for newer datasets and a preference for more comprehensive datasets.

Biological admission precedes ranking: organism, host/cell type, stage, assay,
quantity, units, measured versus inferred evidence, mapping, quality checks and
data access must fit the question. A highly cited review is not a dataset.
Incompatible cohorts remain separate. Popularity cannot establish ground-truth
accuracy. Unknown dates/counts/coverage remain unknown, not zero or guessed.

## Controlled action items

| Item | Percent done | Time left | Description |
|---|---:|---|---|
| [65.01](../done/65_01_dataset_selection_policy.md) | ✅ | 0 h | Implement and test a versioned selection policy |
| 65.02 | 0% | 4–8 h + lookup | Audit publication identity and citation rates for all 162 registered sources |
| 65.03 | 0% | 6–12 h + review | Search every current slot for literature challengers and expose comparison gaps |
| 65.04 | 0% | 8–16 h + compute | Validate superior alternatives, promote eligible replacements and rerun affected checks |

65.02 depends on 65.01; 65.03 depends on 65.01/65.02; 65.04 depends on all three.
Each completion requires recorded evidence, relevant checks, a commit and nightly
publication. Repost the full progress table (instruction 64 plus these four rows)
when an item completes. Keep uncertain comparisons and unsuccessful replacements
visible rather than declaring every incumbent the best in the literature.

## Frozen audit scope

All 162 registry entries, all 297 slot declarations, current installed storage
coverage and candidate/refusal records. Host status is measured separately from
the pending host gene-space builders: seven human and three mouse dataset sources,
20,989/15,590 protein reference rows, 8/26 Tg and 4/24 Pf host slot views populated.
These are storage/slot counts, not host inference or benchmark completion.

Publication metadata comes from recorded primary-service queries and frozen
responses. Only identifiers verified by those responses may enter the audit's
resolved publication fields. Dataset accession links remain distinct from papers.
No wholesale data replacement follows from citation ranking alone; each proposed
replacement must pass its biological admission and affected validation checks.

## Completion evidence

65.01: `starplast/dataset_selection.py` implements the versioned, admission-gated
policy; `docs/dataset_selection.md` documents factors and unknowns. **389 focused
checks passed** for policy, docstrings and organism invariants. AGENTS.md and the
session handoff persist the user rule. Publication metadata and literature search
are next; no dataset replacements have been made by this policy-only step.
