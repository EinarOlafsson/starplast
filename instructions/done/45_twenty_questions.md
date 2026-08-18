# 45 — A hundred questions, the best twenty, shipped with the program

**Status: DONE 2026-08-18.** `starplast/data/questions.json` (100), `starplast/questions.py`,
`scripts/generate_question_table.py`, `scripts/run_question_catalogue.py`, and tab **8 · Questions**
in the analysis panel.

## What shipped

**100 candidates, 20 shipped, 1 kept because it fails.** The hundred were written against a
generated inventory of what this table can actually be asked -- 50 categorical columns with two or
more classes of at least fifteen members, and 290 numeric columns that reach that once binned --
rather than against anyone's memory of the catalogue.

**Every verdict is measured rather than asserted.** Each of the 100 is run through the real
`recipes.close` against the shipped cache by `scripts/generate_question_table.py`, so a question its
author kept but closure refuses is dropped with the closure's own words, and the published document
regenerates from the catalogue instead of being maintained beside it.

## The selection rule, which the data chose

Instruction 45 warns that twenty variations of one question is one recipe run twenty times. The rule
that prevents it is **one question per distinct holdout**: twenty recipes predicting the SAME label
from twenty input sets differ only in their inputs, while twenty distinct holdouts are twenty
questions. Applied to the 26 candidates that survived closure it yields exactly 20, and the six it
sets aside stay in the file marked `covered` with their reasoning intact.

Spread of the twenty, by axis: fitness and essentiality 5, relationships between genes 4, life-cycle
stage and conversion 3, localisation and export 2, metabolism 2, regulation and chromatin 2, host
interaction 2.

## Immunity and antigenicity ships NOTHING, and that is the most useful thing here

Every candidate on that axis that survived closure turned out to hold out a **fitness screen** rather
than an immunity measurement -- "which genes are required to survive in IFN-gamma-activated
macrophages" is a fitness question wearing an immunity label. They were re-axed to what their holdout
actually measures rather than shipped under a label that would have made the spread look complete.

The reason the axis is empty is measured: `n_bcell_epitopes` has 34 labelled genes and
`iedb_epitope_count` has 221, out of 8,140. A cluster needs fifteen labelled genes before it can be
evaluated at all, and in the instruction-44 worked example no cluster reached that -- the control
abstained on every one. **A label that cannot reach the floor is not a holdout and not a control.**
That is an acquisition target stated precisely enough to act on, and it belongs beside the empty
slots' `blocked_by` verdicts rather than hidden inside a question that quietly returns nothing.

## The two questions their authors got wrong, and how we know

Both were kept by the agent that wrote them and both were refused by closure, for the same reason and
neither by hand:

* a PVM question controlled by `compartment_best` against a `pvm_proximity_positive` holdout --
  measured association, so the control is a copy;
* a cell-cycle question controlled by `cellcycle19092_7h_r1` against `cellcycle_phase` -- declared
  family, the same experiment's own replicate.

The first is now **shipped deliberately as the recipe expected to fail**, labelled as such in the
program and in the document. A library where every question works is a library that has been fitted
to its answers.

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
