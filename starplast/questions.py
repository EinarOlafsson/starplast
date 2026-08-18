#!/usr/bin/env python3
"""The shipped question catalogue: what this data can be asked, and what it cannot.

Instruction 45 asked for a hundred questions and the best twenty. The hundred are kept rather than
discarded because the DROPPED ones are the useful half: each records something this data cannot
answer and why, which is the same service the empty slots' `blocked_by` verdicts perform in the slot
atlas. A question dropped for "the holdout has 221 labelled genes and a cluster needs 15 to be
evaluated" is a target for the acquisition campaign, stated precisely enough to act on.

**A verdict written by hand is a claim; a verdict from `recipes.close` is a measurement.** Every
question in this catalogue is run through the real leakage closure against the real table, and the
closure's refusal -- not the author's reasoning -- is what marks a question undeliverable. The two
disagree often enough to matter: a question can look circular and survive closure, and a question can
look clean and be refused because a control shares an experiment with its holdout three steps away.

The catalogue lives in `starplast/data/questions.json` and the published tables regenerate from it
with `scripts/generate_question_table.py`, so the document and the program cannot drift apart.
"""
from __future__ import annotations

import json
import os

import pandas as pd

from .recipes import Recipe, close

#: Where the shipped catalogue lives.
CATALOG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "questions.json")

#: The axes instruction 45 asked the twenty to spread across. Named here so "the spread is visible
#: rather than assumed" is checkable by a test rather than by reading the list.
AXES = ("localisation and export", "life-cycle stage and conversion", "fitness and essentiality",
        "host interaction", "metabolism", "regulation and chromatin", "immunity and antigenicity",
        "relationships between genes")


def load(path: str = CATALOG) -> list:
    """Every candidate question, kept and dropped, as plain dicts."""
    if not os.path.exists(path):
        return []
    with open(path) as fh:
        return json.load(fh)


def as_recipe(question: dict) -> Recipe:
    """One catalogue entry as a runnable recipe.

    The catalogue carries fields a recipe does not (the verdict, the reason it was kept or dropped),
    so this is a projection rather than a cast: `Recipe.from_dict` already ignores what it does not
    know, and doing the same here would let a typo'd key ride along unnoticed in both directions.
    """
    fields = ("question", "inputs", "holdout", "validation_holdout", "holdout_bins", "control_bins",
              "min_precision", "axis", "expectation", "seed", "scope", "organism", "expect_refusal")
    return Recipe.from_dict({k: v for k, v in question.items() if k in fields})


def verify(nodes: pd.DataFrame, questions: list | None = None) -> pd.DataFrame:
    """Run every question through the real closure and report what it actually does.

    The columns worth reading are `declared` against `closure_ok`: a question its author kept that
    closure refuses is a question the author was wrong about, and it is dropped with the closure's
    own words rather than argued with. The reverse -- dropped by hand, accepted by closure -- is
    worth a look, because it is usually a criterion the closure cannot see (whether the answer would
    change what anyone does next) rather than a mistake.
    """
    rows = []
    for q in (questions if questions is not None else load()):
        recipe = as_recipe(q)
        result = close(nodes, recipe)
        holdout = result.holdout_column
        labels = 0
        if holdout and holdout in nodes.columns:
            from .recipes import label_series
            labels = int(label_series(nodes, holdout, recipe.holdout_bins).notna().sum())
        rows.append({
            "question": q.get("question", ""),
            "axis": q.get("axis", ""),
            "declared": q.get("verdict", ""),
            "shipped": bool(q.get("shipped")),
            "expect_refusal": bool(q.get("expect_refusal")),
            "holdout": recipe.holdout,
            "holdout_labelled": labels,
            "control": recipe.validation_holdout,
            "closure_ok": result.ok,
            "blocks": len(result.blocks),
            "columns": len(result.columns),
            "removed_from_inputs": len(result.removed),
            "refusal": result.refusal,
            "reason": q.get("reason", ""),
        })
    return pd.DataFrame(rows)


def shipped(nodes: pd.DataFrame | None = None, questions: list | None = None) -> list:
    """The questions that may be offered in the program: kept by their author AND passed by closure.

    `shipped` rather than `verdict == "keep"`, because more questions are sound than are worth
    offering: twenty recipes predicting the SAME label from twenty input sets is one recipe run
    twenty times, so the catalogue ships one question per distinct holdout and marks the rest as
    kept-but-covered. The unshipped ones stay in the file with their reasoning intact.

    `expect_refusal` questions are the exception and they are deliberate. Instruction 45 asked for at
    least one recipe expected to FAIL, kept and labelled as such, because a library where everything
    works is a library that has been fitted to its answers.
    """
    questions = questions if questions is not None else load()
    kept = [q for q in questions if q.get("shipped") or q.get("expect_refusal")]
    if nodes is None:
        return kept
    verified = verify(nodes, kept)
    ok = set(verified.loc[verified.closure_ok, "question"])
    return [q for q in kept if q["question"] in ok or q.get("expect_refusal")]
