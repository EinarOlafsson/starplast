# 57 · "Start here": a guided path from what a user has to the strategies worth running

**Status: complete (nightly, 2026-09-30).** The tab is beside Strategies, which is untouched.
`tests/test_guided.py` (29 tests) passes, and so does the full suite.

The user, 2026-09-30:

> "instead of presenting users with strategies, ask them about what data they have and go from there
> and end in the strategies. so in addition to the strategies, start from the data or gene list or
> whatever the user has, keep the strategies tab but add another tab. this tab should start with one
> question and lead the user through and end at the strategies. during this process the user will
> choose one gene or a set of genes they care about, choose a label they care about, etcetera."

## What was built

* `starplast/guided.py` — **Qt-free**: the question tree, its options, the branching, and
  `recommend(answers, ctx)`. Testable headless and usable from a notebook.
* `starplast/guided_panel.py` — the view: one question per screen, a clickable trail of answers, and
  the recommendation cards. It holds no window; it emits `run_strategy`, `open_strategy` and
  `open_view`, and `install(window)` wires those. `Window._guided()` is the whole hook in `app.py`.
* `starplast/strategy_panel.py` — `headline()` and `shipped_card()` became module functions that the
  methods now call, so the Start-here cards show the Strategies tab's numbers rather than a second
  computation of them. Nothing else in that tab changed.
* `starplast/help_index.py` — `guided_entries()`, derived from `guided.STEPS` and `guided.GOALS`, so
  the tab, every question and every goal are reachable from the search beside Help.

## The tree

| step | asks | options |
|---|---|---|
| `have` | What do you have? | one gene · a gene list · my own measurement · a label · nothing yet |
| `space` | Which organism? | the installed spaces (`organisms.codes(available=True)`) |
| `gene` / `genes` / `label` / `measure` | the subject itself | search by id, symbol or product; paste / file / map gate / example set; a column list showing each column's coverage |
| `label` (again) | which label to give one gene | only for a single gene whose goal is `predict` |
| `goal` | What do you want to know? | more like mine · predict it · explain it · partners · compare · learnable — the ones that make sense for what they have |
| `baseline` | Compared with what? | only when the goal is `compare` |

"Nothing yet" skips the goal question and implies the tour. Every first answer reaches an end state,
and every option leads somewhere: `_walk` in the test walks the whole tree and checks both, plus
that every end state produces at least one recommendation.

## How the ranking works

`3.0 − 0.25 × rank` for a strategy the (goal, subject) pair is FOR, else 0.6; `+0.8` when its
scorecard task is one the goal asks for; the calibration grade **for that organism**
(reliable +1.5, works when tuned +0.8, weak −0.4, no skill −3.0, untestable −1.0, unknown 0); and
cost as a tie-break only. Dropped before ranking: anything this table cannot fill (`usable`, asked of
the strategy's own defaults, not of a list kept here) and anything that cannot take the subject
chosen (`fits_subject` — a label-calling method is never offered for someone's numeric screen).
`caveat()` says plainly when nothing better than weak exists. `prefill()` fills only parameters the
strategy has, from an answer that means the same thing, and every prefill in the whole tree is
handed to `Strategy.settings` in a test.

## Verified

* 29 tests in `tests/test_guided.py`: the tree from every first answer, every option, the ranking
  rules, per-organism grades, prefilled settings accepted by `Strategy.settings`, the weak-only
  caveat, the panel headless with a tooltip on every control, the window hook, and the help search.
* Screenshots of every step and of the recommendations for four branches, in
  `results/guided_start_here_2026_09_30/` (rebuilt by walking the panel offscreen). Looked at and
  iterated on: the cards collapsed on top of one another in the scroll area until height-for-width
  was declared on them, and a strategy question ending in "?" read as "?." until the full stop was
  made conditional.
* No new "Tg"/"Pf" literals; no strategy algorithm, calibration number or data table was touched.
