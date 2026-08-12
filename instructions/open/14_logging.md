# A real logger, opt-in, with per-level console control

Requested 2026-08-12. Modelled on the logger in spaCR: levels, and separate control over what
reaches the console versus what is kept.

## Why this is not just `print`

The console pane already tees stdout, so a walk's progress is visible. What it cannot do is
distinguish a progress line from a warning from a traceback, which means it also cannot filter them,
and it keeps nothing after the session ends. Two concrete failures this project has already had that
a log would have caught immediately:

- every AlphaFold fetch 404ing for months, reported as the ordinary "no model available",
- 40 of 288 search configurations completing with no record of which ones or how long each took.

## What to build

`starplast/logging_util.py`:

- Standard levels: DEBUG, INFO, WARNING, ERROR.
- **Opt-in.** Off by default, enabled from Preferences. A tool that starts writing files to a user's
  disk without being asked is a tool people stop trusting.
- **Per-level console control.** Independent choice of which levels print to the console pane and
  which only go to the file, because DEBUG is what you want in the file during a long walk and never
  what you want scrolling past on screen.
- File output under the cache directory, rotating so a long session cannot fill a disk.
- The console pane gains a level filter alongside its existing text filter, and colours by level.

## Preference

Preferences gains: enable logging (off by default), file level, console level, and where the file is
written. The setting persists.

## What must be logged

Every job start, finish, failure and cancellation, with elapsed time. Every dataset fetch, with the
URL and the outcome. Every embedding and clustering with its full recipe. Warnings for the cases
that currently pass silently: a dataset that could not be fetched, a join that matched zero rows, a
version-pinned URL returning 404.

That last category is the point. A silent failure is the expensive kind, and this project has shipped
two.
