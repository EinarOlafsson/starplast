# Repository review

Review of the local 0.41.0 development baseline, with the packaging and interface
changes prepared for 0.42.0. This is a code and product review, not a reanalysis of
the underlying studies.

## What Starplast is trying to do

Starplast is a gene-evidence workbench for *Toxoplasma gondii* and
*Plasmodium falciparum*. It combines heterogeneous measurements and annotations,
projects selected inputs into a map, and lets a researcher investigate genes or
groups using the underlying evidence. Its search tools try feature combinations
and clusterings to recover withheld labels or identify candidate findings.

That is a useful companion to spaCR. The practical workflow is microscopy screen,
gene-level result table, contextual evidence, then follow-up experiment. Keeping
the two applications in separate environments is sensible: the graph and analysis
dependencies need not constrain spaCR's imaging installation. The current launcher
opens Starplast but does not automatically transfer a screen or its provenance.
The importer currently recognizes Toxoplasma identifiers only; the Plasmodium
map can be browsed, but importing Plasmodium screens needs additional work.

## What works well

- Dataset and biological-slot registries record where measurements came from and
  what question they address. This is more useful than an undifferentiated matrix.
- Different relationship types stay separate. A shared publication is not treated
  as the same observation as a measured interaction.
- Missing measurements, inferred annotations, and known labels have distinct roles.
- Feature selection, scaling, missing-value policy, and random seed can be saved
  in an embedding recipe.
- Search and validation code contains explicit target exclusions, circularity
  checks, and leakage tests. These are necessary for the stated use case.
- CLI discovery and analysis modules can run without opening a window.
- Built caches support offline browsing, and there is substantial regression coverage.

## Where the design needs care

**The 3D map should support an evidence workflow.** It is useful for exploring
neighbourhoods and selecting candidates. It does not directly represent cellular
geometry, interaction probability, or biological distance. Filters, evidence
inspection, and export are the parts that turn a view into a research tool.
A sortable candidate table with source links would be a valuable next addition.

**Search scores can reward the search procedure.** Holding a label out of the
input is necessary, but repeatedly selecting configurations using that same label
can still overfit the evaluation. The existing
`holdout_cv.nested_structure_cv` helper separates inner model selection from outer
label evaluation; that is a good foundation for reporting generalization. Make
that distinction clear wherever search scores are displayed, and confirm candidates
with independent evidence. Existing leakage tests are useful checks, not proof
that all biological proxies have been excluded.

**The executed method needs to travel with the result.** In `embedding.embed`,
UMAP exceptions can lead to PCA fallback. The log reports it, but the supplied
recipe still requests UMAP. Persist the actual method and backend, and offer a
strict mode that fails instead of changing methods. This remains a follow-up;
the release polish does not change search semantics.

**The GUI is carrying too many responsibilities.** `app.py` handles loading,
rendering, preferences, import, export, evidence formatting, and job controls.
`analysis_panel.py` coordinates several scientific workflows. Keep their existing
interfaces while extracting data loading/export, preferences, and evidence
presentation into smaller modules. A wholesale rewrite would add risk without
improving the biological task.

**Data snapshots need explicit identity.** Bundling a modest cache is convenient.
For future growth, separate versioned data snapshots from code and attach checksums,
source versions, build parameters, and coverage reports. The current source
registry is a good foundation. A wheel should never include a developer's saved
embeddings or analysis output.

**The README overstated coverage.** A gene absent from the indexed corpus is not
necessarily unstudied, and the app cannot show everything known about a gene.
The new README describes the indexed evidence and links to a generated catalogue
instead of making those broader claims.

## Why the aggregation views were removed

The Galaxy mode found connected components in a coarse occupancy grid over the
embedding, then drew their centroids. It summarized the shape of a particular
projection; it did not identify a biological compartment or a validated module.

The Orthogroup mode averaged the positions of visible genes in families with at
least four members. A dispersed family's mean can lie between unrelated regions.
The colour was taken from its first member, so it did not summarize the family's
annotation distribution. The summary dots had no membership drill-down; picking
still operated on gene coordinates.

These modes could become useful with explicit membership tables, spread measures,
and navigation back to the constituent genes. In their existing form, they added
ambiguity without supporting a distinct task. They have been removed from the app.
The individual-gene map, orthogroup edges, filters, and selection tools remain.
The standalone `lod` helpers remain available for older analysis code.

## Changes in this pass

The release adds a logo and application icon, a shorter README, user and API guides,
a generated documentation site, missing application docstrings, settings help,
optional GPU installation, package-content checks, and automated version publishing.
The `starplast-core` distribution name is retained for compatibility with spaCR.

The next work should prioritize recording the executed analysis method, transferring
screen tables and provenance from spaCR, and giving candidate evidence a table-based
review workflow. More rendering modes would be a lower priority.
