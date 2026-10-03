# 61 · Tutorial videos: six screen captures of the real application

Opened and finished 2026-10-03 on a worktree of `nightly` at 0.50.0.

**The task.** spaCR ships a silent mp4 and a poster image per tutorial. Do the same here, but
generate the clips from the running program rather than filming a desktop: drive the real `Window`
and panels offscreen the way `scripts/build_tutorials.py` already drives them, grab frames with Qt,
and encode with the system `ffmpeg`. Six clips, each a different entry point and a different
strategy, each a real workflow end to end. Publish them with the tutorials.

**Nothing here changes a strategy algorithm, a calibration number, a shipped data table or the UI.**
This is a recording job. The only new code is `scripts/tutorial_video.py` and four lines each in
`scripts/build_tutorials.py` and `scripts/build_docs.py`.

---

## 1. The harness

`scripts/tutorial_video.py`.

* **`Studio`** — one offscreen `QApplication` and one real `Window`, at 1680x920, with `QSettings`
  and `STARPLAST_STATE` redirected to a scratch directory, exactly as `tests/conftest.py` does and
  for the same reason. One application films all six clips: building six windows would cost six
  times as long, and Qt keeps style state across applications.
* **`Clip`** — the recorder. A frame is the cached window grab, plus an overlay this module paints:
  a synthetic arrow pointer, a rounded highlight around the control about to be used (flaring on
  the press), a label across the empty central strip, and a 130-pixel caption band under the window
  carrying the beat number, a heading and a sentence. The frame is 1680x1050, both even, as yuv420p
  requires.
* **Beats** — `say`, `note`, `move_to`, `press`, `type_into`, `scroll`, `hold`, `wait`,
  `mark_poster`. The window is re-grabbed only when something changed, so a four-second pause costs
  one grab rather than fifty and a whole clip encodes in about half a minute.
* **`wait`** — a strategy started from the window goes to the job runner, so the panel returns at
  once and the answer arrives later. The wait is pumped and **recorded** (the first seven seconds of
  it) rather than cut out: the jobs panel is doing something, and a tutorial that hides the wait
  teaches the wrong thing.
* **Encoding** — PNG frames into a temporary directory (`STARPLAST_VIDEO_SCRATCH` to put them off a
  small `/tmp`), then `/usr/bin/ffmpeg` to H.264, yuv420p, 12 fps, CRF 27, `+faststart`; a marked
  frame becomes `<slug>_poster.jpg` at 1120 px.
* **Determinism** — fixed window size, fixed beats, isolated settings and state, shipped tables and
  calibration. `reset_window` puts the window back between clips (docks hidden, filters cleared, the
  opening map restored, the guided tab restarted), and the two clips that run a strategy are
  recorded *after* the star-map clip, because a run leaves its predicted links in the star map's
  store and the star map would otherwise open showing links the viewer never made.

## 2. The six clips

Measured by `ffprobe` on the committed files (`notebooks/tutorial_videos_2026_10_03.ipynb`).

| Clip | Shows | s | kB | Tutorial |
|---|---|---|---|---|
| `1_find_a_gene` | Find box → evidence panel for GRA16, scrolled through the measurements; a dash is not a zero; a product-text search; colour-by | 60.3 | 759 | 1 Explore |
| `2_start_here` | The guided tab: what you have → organism → gene → goal → label → six ranked cards → Run, the job, the result | 63.4 | 1301 | 3 Gene list |
| `3_test_before_you_trust` | Filter by method → strategy 07's card and its four bars → the worked failure beside the worked success → Test → the verdict and the scorecard | 61.6 | 821 | 4 Test and calibrate |
| `4_pick_a_map` | 78 pregenerated maps grouped → sorted by how well the label maps → the top row's own warning that it restates the label → the best map that does not → colour by label → the score table | 62.2 | 689 | 2 Hold-out search |
| `5_walk_the_network` | Star map centred on GRA16 → measured sources only, then all → the neighbourhood view → re-centred on GRA24 → Back | 64.6 | 1026 | 5 Networks and agreement |
| `6_one_question` | Question 13 of instruction 59 run live: strategy 32, target `compartment`, `min_agree` 2 → 1,723 candidates, the top three named | 62.0 | 826 | 6 Advanced models |

Six posters, 83–128 kB each. **6,068,075 bytes (6.1 MB) in all**, against a 40 MB budget.

Every number on screen is the one the program produced while the clip was being recorded. Two are
worth naming because they are checkable:

* Clip 3's verdict: `PASS` — correct calls per hidden gene 0.334 against 0.052 ± 0.007 under 10 runs
  on shuffled labels, effect +0.283, 951 hidden.
* Clip 6's top three candidates: `TGME49_209000`, `TGME49_225745`, `TGME49_266010`, which are the
  first three rows of `results/questions_2026_09_30/Q13_understudied_first_candidates.csv`. The
  notebook asserts the match.

## 3. Two limitations, shown rather than hidden

* The central 3D view is a `QOpenGLWidget` and the offscreen Qt platform does not support one, so it
  records as an empty strip. The strip carries the sentence "the 3D map needs OpenGL, which a
  headless recording cannot draw" rather than being left as an unexplained black column. Compositing
  a committed screenshot into it was considered and rejected: it would be staging.
* A strategy run takes as long as it takes, and the clip holds on the panel before and after.

## 4. Publishing

* The clips go to `docs/tutorial/video/`.
* `publish_pages()` writes a player into each written tutorial's settings panel and a grid of
  posters into `docs/tutorial/index.html`, between `<!--starplast-video-->` markers so that running
  it twice leaves one copy. Both pages are generated, so `scripts/build_tutorials.py` calls
  `publish_pages()` at the end of its own run and a rebuild cannot drop them.
* `scripts/build_docs.py` now copies `docs/tutorial/` into `docs/site/` and carries a **Tutorials**
  link in the site navigation. The written tutorials were built and committed but never published;
  they are now.
* They are documentation, not package data. `pyproject.toml` names its packages explicitly and sets
  `include-package-data = false`, so nothing under `docs/` can reach the wheel; no exclusion was
  needed. Verified: the 0.50.0 wheel is **65.0 MB** (limit 75), `scripts/check_wheel.py` passes, and
  the wheel contains no `docs/`, `tutorial` or `.mp4` entry.

## 5. What was verified

* Frames extracted from every clip and looked at, one every 90 frames (7.5 s), through several
  rebuilds. Four things were found by looking and fixed: raw HTML markup burned into a caption
  (labels are rich text — `first_line` strips tags); a row number typed into the star map's gene box
  instead of an accession (`last_nodes` carries row numbers, not identifiers); a gallery map and a
  jobs dock left on screen by the previous clip (`reset_window`); and the maps table scrolled
  sideways far enough to push its column headings out of frame.
* The caption on clip 4 originally claimed no gallery map is built from the label being scored. The
  panel's own sentence says otherwise for the top-scoring map, so the clip was rewritten to read
  that warning aloud and then pick the best map that does **not** restate the label. A tutorial that
  contradicts the program it is teaching is worse than no tutorial.
* Full suite: see the commit message.

## 6. Rebuilding

```bash
QT_QPA_PLATFORM=offscreen python scripts/tutorial_video.py             # all six, about 4 minutes
QT_QPA_PLATFORM=offscreen python scripts/tutorial_video.py 4_pick_a_map
python scripts/tutorial_video.py --pages                               # relink without recording
python scripts/tutorial_video.py --notebook                            # record what they came to
```
