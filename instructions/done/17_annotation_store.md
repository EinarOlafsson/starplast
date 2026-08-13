# Saving an annotation, with the numbers that justify it — DONE 2026-08-12

Shipped after task 20, which is the ordering constraint the index calls non-negotiable, and it is
the reason this task is mostly a set of refusals rather than a feature.

## `starplast/annotations.py`

An `Annotation` carries: the gene, the label proposed, the target column, the cluster and the full
configuration behind it, the cluster's composition (`cluster_frac_category`,
`cluster_frac_contradicting`, `enrichment`), the **validated precision and recall for that
category**, how many folds and whether they re-embedded, the date, and free-text reasoning. Every
field is there because leaving it out makes the row unfalsifiable.

`AnnotationStore` is a CSV of its own — readable by a person, because the reasoning field is the part
a collaborator will want to argue with, and re-read before every write, because a store that cached
a shared file would overwrite whatever someone else added while the window was open.

## Three refusals

- **No validated precision, no save.** An error, not a warning. The message names the tab to run.
- **Precision below 10%, no save.** Measured on this map's coarse components, annotating from them
  would be right 1-6% of the time. A store that accepted those rows would be a machine for producing
  them.
- **One bad row refuses the whole save.** A partial save that dropped the unjustified rows silently
  would leave the user believing they had saved what they were looking at.

Re-annotating the same gene for the same target replaces rather than duplicates, and is logged; the
same gene under a different target is a different proposal, not a changed one. Withdrawing is as easy
as proposing.

## Never confused with a measurement

Never written into the node table — a test loads the table, saves an annotation for one of its genes,
and asserts the table is byte-identical afterwards. On the map, `annotations` is its own colour mode
in a colour used for nothing else: a test measures the distance from every categorical palette in
every theme, from the attention tiers and from both greys, and fails if the annotation colour comes
near any of them. Genes with no annotation are grey — not a category, not zero: nobody has proposed
anything for them.

## The interface

Saving is possible **only from the Validation tab**, only from a validated run, and only for the
candidate list on screen. `annotations.from_candidates` joins the candidates to their error rate in
one place, so a precision from a different category or a different clustering cannot travel with a
row. A refusal appears as an explanation in the same place the circularity refusal does.

## Verified

- `annotations.py` at 100%, `analysis_panel.py` back to 100%.
- The refusals are tested by their effects: a refused save leaves no file at all, and an unvalidated
  candidate list leaves the store empty and the reason on screen.
