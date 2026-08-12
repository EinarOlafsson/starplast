================================================================================
README (NARROW SCOPE), API REFERENCE, AND TOOLTIPS
================================================================================

STATE
    README.md is 298 lines and is mostly an essay about the science: "What
    makes it different from a network viewer", "Interpretation rules the UI
    enforces", "Finding structures that predict something they were never
    told". Installation and use are a small part of it.

    Tooltips exist only in the Preferences dialog. Every module carries a full
    docstring holding the reasoning, not just the signature, so the raw
    material for an API page is already written; nothing generates it.

THE README'S SCOPE IS NOW FIXED (user instruction, 2026-08-11)
    The README contains THREE things and nothing else:

      1. what the program does          -- short. What it is, what you see.
      2. how to run it                  -- install, launch, troubleshoot.
      3. which datasets it includes     -- AS A TABLE, with for each dataset:
                                           its reference (citation / PMID /
                                           accession) and the TYPE OF DATA.

    Everything else moves out. It does not get deleted -- the reasoning is the
    most valuable writing in the repo -- it goes where it belongs:

      the science, the design decisions, why the map is a UMAP  -> HANDOFF.md
      methods prose for the manuscript                          -> MATERIALS_AND_METHODS.md
      interpretation rules, circularity, attention correction   -> HANDOFF.md
      per-function detail                                       -> the API page

    Rationale: a README is read by someone deciding whether to run the thing
    and then running it. Argument belongs in the paper, not in the front door.

WHAT TO DO
    - Rewrite README.md to those three sections. Keep it short enough to read
      in full before deciding to install.

    - GENERATE the dataset table from `starplast.datasets.REGISTRY`, do not
      hand-write it. The registry already carries name, level, kind, provides,
      coverage, pmid, accession, citation and url for all 26 datasets. A
      hand-written table drifts the moment a dataset is added, and a README
      that misstates the data is worse than one that omits it. Add
      `datasets.readme_table()` returning markdown, and a test asserting the
      committed README table matches what the registry currently produces --
      so adding a dataset without updating the README fails the suite.

    - Group the table by level (DNA / transcription / translation /
      post-translation / reference), which is the same taxonomy datasets/ uses
      on disk. Columns: dataset, type of data, coverage, reference.

    - Tooltips on every control in the main window, not just Preferences. The
      useful ones say WHY, not what: "occlude: nearer points hide farther ones
      -- additive blending saturates dense regions to white and destroys the
      colour encoding" beats "blending mode".

    - An API page generated from the docstrings (pdoc is one file of config;
      Sphinx is heavier and buys little here), published to GitHub Pages by an
      Actions job.

HOW TO KNOW IT WORKED
    A person who has never seen the repo can read the README in three minutes,
    install it, launch it, and see in a table exactly which published datasets
    are inside and where each came from -- without meeting a single argument
    about method.

    And: add a 27th dataset to the registry without touching README.md, and
    the test suite fails.

TRAPS
    Do this AFTER task 02. It documents a UI that 02 changes, and screenshots
    taken now are screenshots to retake.

    Do not let the dataset table become a list of file paths. The reference is
    the PUBLICATION -- that is what a reader needs to judge the data. The path
    is an implementation detail and belongs in the registry, not the README.


================================================================================
DONE -- 2026-08-12
================================================================================

README: 298 lines to 125, holding the three things the instruction fixed --
what the program does, how to run it, and which datasets are inside with their
references and data types.

    The dataset table is GENERATED from starplast.datasets.REGISTRY. A test
    asserts the committed README matches what the registry currently produces,
    so adding a 27th dataset without updating the README fails the suite. A
    hand-written table drifts the moment a dataset is added, and a README that
    misstates which data is inside is worse than one that omits it.

    The reference column names the PUBLICATION, never the file path -- a path is
    an implementation detail that changes with the layout, a PMID is not. Where
    a citation is not confirmed the cell says so rather than being blank,
    because a blank hides it and unresolved() exists precisely so those cannot
    reach a manuscript unchecked.

    A test also fails if a heading about circularity, interpretation or
    structural holes reappears, which is how the scope crept the first time.

API REFERENCE: pdoc over all 24 modules, published to Pages by an Actions job on
every push to main. Not committed -- a checked-in copy drifts from the
docstrings it came from.

    Two tests guard it. Every module must have a docstring of at least forty
    words, since a one-line summary is a signature restated and leaves a hole in
    the published page. And pdoc is actually run, which imports every module --
    catching anything that cannot be imported headlessly, or a docstring whose
    markup breaks the renderer. Both otherwise surface only in CI.

TOOLTIPS: on every main-window control, not only in Preferences. Search, level
of detail, colour mode, spin, the attention correction, draw-all-edges and the
compartment list. Tested for substance rather than presence -- a tooltip under
fifteen words is failing, because "colour mode" restates the label while "grey
always means unknown, never a category and never zero" is worth reading.

NOT DONE
    Screenshots. The instruction to take them came before task 02 changed the
    UI, and they would want retaking after the visual work; they are also the
    one thing here that cannot be verified by a test.
