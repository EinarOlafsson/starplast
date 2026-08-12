# American spelling, done safely

Requested 2026-08-12. User-facing text should use American spelling throughout.

## This was attempted and reverted, so read this before trying again

A regex over every `.py` file replacing British forms with American ones broke the suite badly:
13 failures and 121 errors. It renamed **function names**, not only prose -- `normalise()` became
`normalize()` while every caller and test still said `normalise` -- and it rewrote a state
constant's VALUE, `CANCELLED = "cancelled"` to `"canceled"`, which desynced it from the tests that
compared against the literal.

The lesson is specific: **substitute inside string literals only, never across whole lines.** Parse
the file and rewrite only `ast.Str`/`ast.Constant` nodes, or at minimum only the text between
quotes, and never touch a `def`, a class name, an import, or the right-hand side of a constant that
other code compares against.

## Two separate jobs, and only the first was asked for

**1. User-facing text.** Labels, tooltips, status messages, dialog text, docstrings. This is the
request and it is safe if done through the AST.

**2. Identifiers.** `colour_of`, `colour_mode`, `colour_by`, `set_colour_mode`, `categorical_colours`,
`unknown_colour`, `localisation.py`, `normalise()`. Renaming these is a wide mechanical refactor
touching the app, the panel, the theme module and their tests. It is worth doing once, deliberately,
with the tests run between each rename -- not folded into another change.

Do 1 first and confirm the suite is green before starting 2.

## Words in play

optimise/optimize, visualise/visualize, normalise/normalize, summarise/summarize, colour/color,
behaviour/behavior, centre/center, labelled/labeled, cancelled/canceled, neighbourhood/neighborhood,
favour/favor, modelling/modeling, recognise/recognize.

`localisation` appears as a module name, a feature-block name and a scientific term. Decide once
whether the block name changes, since it is stored in saved embedding recipes -- renaming it silently
invalidates every stored recipe that names it.

## Done when

The suite is green, no user-facing string uses a British form, and any identifier rename was a
separate commit with the tests green on each side of it.
