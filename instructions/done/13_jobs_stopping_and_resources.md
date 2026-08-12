# Jobs, stopping, and freeing memory — DONE 2026-08-12

Requested: a way to stop an analysis or all analyses, per-job quit/force-quit, and a way to clear
RAM / VRAM / CPUs as in spaCR.

## What the bug turned out to be

`AnalysisPanel` ran its work on a private QThread, separate from the JobRunner. That single fact
produced three complaints at once: its jobs never appeared in the Jobs panel, none could be stopped,
and `_run` refused to start while its thread was busy — so for the several minutes a search takes,
every button on Data, Map, Clusters and Inference did nothing but write "a job is already running"
into the status bar. Four tabs that look broken, and the explanation in the one place nobody is
looking while watching a map.

## Done

- The panel submits through the window's runner. Jobs are named for what they are, and several run
  at once.
- Stop selected, stop all, right-click to stop or copy a failed job's traceback.
- Cooperative cancellation via the progress callback, which every long analysis calls once per
  configuration. Measured: a job cancelled during a 10,000-step loop stopped at step 194.
- A stop records as `cancelled` with no error and no traceback. A stop reported in red teaches people
  to ignore the failure list.
- A resource line — process RSS, system memory, GPU memory where visible, CPU count — polled every
  three seconds, and a "free memory" button that drops the level-of-detail cache, collects garbage
  and empties the GPU cache, then reports what it released.
- Free memory deliberately does **not** touch running jobs, and a test asserts a job in flight
  survives it. A "clear RAM" button that silently killed a search would be a data-loss button wearing
  a housekeeping label.
