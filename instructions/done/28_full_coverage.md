# Every module at 100% — DONE 2026-08-13

The standing goal, reached: 7,389 statements, 0 uncovered, 1,849 tests.

The starting point was 99% — 76 lines across five modules. What they were is the interesting part,
because "99%" hid three different things.

## Two of them were dead code

`celldiagram.recolor` carried a nested `paint()` that nothing called: the loop below it does the
substitution inline. `celldiagram._greyscale` was left from the version of the diagram that greyed
the artwork, and the diagram has been hollow outlines since 0.15.0. Both deleted rather than
covered — a test for a function nobody calls pins behaviour nobody depends on, and makes the module
look larger than it is.

## One was a real defect, found by writing the test

`objectives.adjusted` reported `n_labels` as every distinct string in the label column. `score`
excludes the absence labels (`unassigned`, `unknown`, `nan`, empty) and anything below `min_label`,
so on the real compartment column it said "24 classes" for a score computed over eleven — and
`n_labels` exists precisely to say how many classes the null is correcting for. It now counts the
classes the score was actually computed over.

## The rest was the branch nobody drives

The uncovered lines were almost entirely of two kinds:

* **The file dialogs.** Every save and load has a "no path given, so ask" branch and a "the user
  pressed cancel" branch, and the tests all passed a path. Cancel is now tested for all five —
  results, bundle, CSV, import — because a cancel that writes a file under a default name is a file
  nobody will find again.
* **The failures that must not take the session down.** A fetch that raises (logged at WARNING and
  re-raised, because a quiet failed download is indistinguishable from a source with nothing to
  give); `/proc` that is not there; a torch that will not answer; a framebuffer grab with no GL
  context; a group the artwork does not contain; a column a later cache no longer has. Each is now
  driven for real — a fake clipboard, a raising `open`, a fake `torch.cuda`, a supplied `QImage` —
  rather than marked as unreachable.

No `pragma: no cover` was added anywhere. The rule stands: drive the branch, or delete it.

## Verified

1,849 pass on pandas 2.3.3 with every module at 100%. Re-checked under the user's interpreter
(pandas 3.0.5, real cache, offscreen Qt): `adjusted` reports 2 classes for a column of a, b and
unassigned, and the window still opens.
