# A real logger, opt-in, with per-level console control — DONE 2026-08-12

`starplast/logging_util.py`, wired into the console pane, Preferences, jobs, fetches and the
embedding.

## Two switches, because they answer different questions

**The file is opt-in and off by default**, as asked. **What reaches the console is a separate
choice, and it is never off** — it defaults to WARNING. That is a deliberate departure from reading
"opt-in" as "silent until enabled": a warning that is only written to a file nobody enabled is a
silent failure with extra steps, which is the exact failure mode this task exists to end. Nothing is
written to disk unless asked; warnings are still audible.

## What it does

- Levels DEBUG / INFO / WARNING / ERROR, one logger root (`starplast`), never propagating to the
  root logger so another library's `basicConfig` cannot print everything twice.
- Rotating file (2 MB, 3 backups) under the platform cache directory. A directory that cannot be
  written costs the file and reports it, not the console.
- `configure()` is idempotent — it replaces handlers rather than adding them, which is how a logging
  system starts being called noisy and then gets turned off.
- The console handler resolves `sys.stdout` **at emit time**, because the console pane installs its
  tee afterwards; a handler holding the stream it was built with would write to the terminal and
  never appear in the pane.
- `Timed` logs start, duration and failure-with-traceback for a block.

## The console pane reads the level

Console lines are `[LEVEL] message`, a prefix so a wrapped line still starts with its level. The
pane colours by level and filters by minimum level, alongside the existing text filter. **Levelless
output is never hidden**: most of what this application writes is a `print`, and folding it into
INFO would make an INFO filter useless.

Making that affordable required fixing a real wart: the pane re-rendered the entire document on
every line, which is quadratic on the one screen whose whole purpose is watching something long run.
Lines are appended as they arrive now, and only a filter change re-renders.

## What is logged

Every job start, finish, cancellation and failure with elapsed time (failures with the traceback, so
it outlives the session). Every fetch with its URL and outcome — `sources._get` at INFO on success
and WARNING on failure, and `structures.fetch_alphafold` warning with every URL it tried when none
worked, which is precisely the 404-for-months failure named in the task. Every embedding logs its
full `EmbeddingSpec`; every clustering logs its parameters.

## Preferences, persisted

Enable, file level, console level, and the path being written, in Preferences and stored in
`QSettings`. Applied at startup before anything can want to log. Persisted only when a setting is
actually changed, so opening the application does not rewrite the user's settings file.

Tests now run against an isolated settings directory. Without that, a developer with logging enabled
would have had every `Window` in the suite writing to their real log file, and the suite would have
been rewriting their real preferences.

## Verified

- 1,423 tests pass headless. `logging_util.py` at 100%, `console.py` at 99% across the two files
  that exercise it.
- Run under the pandas 3 interpreter: log enabled at DEBUG, a job run and a job failed, the file
  holding both timings and the `ZeroDivisionError` traceback, the pane receiving level-tagged lines.
