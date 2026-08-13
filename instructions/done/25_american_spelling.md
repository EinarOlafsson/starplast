# American spelling, done safely — DONE 2026-08-12

Job 1 of the two the task file separates: **user-facing text only**. Done through the AST, as the
file instructed after the first attempt was reverted.

## How

A script parses each module and rewrites only string literals, never source lines. Three rules keep
it from repeating the failure that reverted it:

- **Only strings containing a space.** That one rule puts every state value (`"cancelled"`), dict
  key, column name and feature-block name out of reach — the previous attempt rewrote
  `CANCELLED = "cancelled"` and desynced it from the tests.
- **Identifiers inside prose are protected.** Backticked names and any snake_case token are masked
  before substitution, so `` `colour_of` `` in a docstring still names the attribute that exists.
- **f-strings are handled whole, with the braces masked.** Before Python 3.12 the Constant parts of
  a JoinedStr report the position of the entire f-string; rewriting them rewrote the interpolations
  too. The first run of the corrected script did exactly that — `{colour}` became `{color}` and the
  variable no longer existed — which is why the braces are masked and the suite was run before
  anything was committed.

131 strings across 28 modules. The knock-on changes: five tests asserting on the old spellings, and
the README's dataset table plus the generated dataset scripts, both regenerated from the registry.

## Deliberately not changed

`localisation` — module name, feature-block name **and** scientific term. The block name is stored
inside every saved embedding recipe, so renaming it invalidates them; and changing the prose while
the control still reads `localisation` would make the documentation disagree with the interface. It
stays until job 2 decides that deliberately.

## Job 2, not done here, as the task file requires

Renaming the identifiers — `colour_of`, `colour_mode`, `set_colour_mode`, `categorical_colours`,
`unknown_colour`, `normalise`, `localisation.py` — is a wide mechanical refactor across the app, the
panel, the theme module and their tests. The task says it is worth doing once, deliberately, in its
own commit with the tests green on each side. It is left for that.
