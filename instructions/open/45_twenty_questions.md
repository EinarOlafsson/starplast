# 45 — A hundred questions, the best twenty, shipped with the program

**Status: open. Requested by the user 2026-08-18. Execute after 44 builds the recipe object.**

## What to do

1. Think carefully about what questions a recipe can actually answer, and write **100** of them.
2. Take the best **20**.
3. For each of the 20, work out its UMAP inputs, its inference holdout, and its control holdout.
4. Ship those 20 in the program, so a user picks a question rather than assembling one.

## What makes a question worth shipping

Write these down against each candidate; the ones that fail are the ones to drop.

* **It has a wrong answer.** "Which proteins are on the PVM facing the host" can come back empty or
  contradicted. "What is the structure of the proteome" cannot. Only the first is a question.
* **The holdout has enough labelled members to score.** A label with eight genes is recovered or not
  by luck; the scoring floor is already set at 15 elsewhere in this codebase and that is the number to
  respect here.
* **The inputs do not restate the holdout.** This is where these go wrong, and the trap is
  predictable per question type: for anything about LOCALISATION the danger is sequence-derived
  columns -- signal peptide, TM count, orthoLOPIT transfer -- because localisation is heavily
  predictable from them, and a recipe admitting them recovers the label by PREDICTING it rather than
  by finding biology. For anything about FITNESS the danger is the other screens. For anything about
  STAGE the danger is the stage-derived labels.
* **The control is independent of the primary holdout.** Two labels from one experiment or one
  curation source are not two observations. A control that is a copy controls nothing.
* **The answer would change what someone does next.** A gene list somebody would actually pick from
  for a knockout, an antibody, or a localisation experiment.
* **It is answerable with the data that is IN the application.** 159 of 223 slots are filled; a
  question needing one of the empty ones is a question for the acquisition campaign, not for a recipe.

## Spread, not variations

Twenty recipes that are all "where does this protein family live" is one recipe run twenty times. Aim
across the axes the catalogue already separates: localisation and export, life-cycle stage and
conversion, fitness and essentiality, host interaction, metabolism, regulation and chromatin,
immunity and antigenicity, and the relationships between genes.

Say for each of the 20 which of these it belongs to, so the spread is visible rather than assumed.

## Deliverables

* The list of 100, with a one-line reason each was kept or dropped. The dropped ones are the useful
  half of that document: they record which questions this data cannot answer and why, which is the
  same service the empty slots' `blocked_by` verdicts perform.
* The 20, each as a runnable recipe: question, inputs, holdout, control, and the expected shape of a
  believable answer.
* Those 20 reachable in the program -- pick a question, run it, read the table.
* At least one recipe expected to FAIL, kept deliberately and labelled as such. A library where
  everything works is a library that has been fitted to its answers.

## Done when

* The 100 exist with their keep/drop reasons.
* The 20 run end to end and return genes, with their control's agreement beside the inference.
* A user can choose one in the application without writing a recipe by hand.
* The deliberate failure fails, and its failure is legible.
