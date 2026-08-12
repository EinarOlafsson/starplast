================================================================================
API REFERENCE, README, AND TOOLTIPS THROUGHOUT
================================================================================

STATE
    Asked for: "tooltips and an API and detailed readme on github".

    Tooltips exist only in the Preferences dialog. Every module has a full
    docstring already -- they carry the reasoning, not just the signature -- so
    the raw material for an API page is written; nothing generates it.

    README.md is thorough on the science and thin on use: no screenshots, no
    install troubleshooting, no worked example.

WHY IT MATTERS
    The repo is going into a PLOS ONE submission. A reviewer who cannot see
    what the tool looks like or how to call it will judge the manuscript on the
    manuscript alone.

WHAT TO DO
    - tooltips on every control in the main window, not just Preferences.
      The useful ones say WHY, not what: "occlude: nearer points hide farther
      ones -- additive saturates dense regions to white and destroys the colour
      encoding" beats "blending mode".
    - an API page generated from the docstrings (pdoc is one file of config;
      Sphinx is heavier and buys little here)
    - README: screenshots of the four themes and three LOD tiers, a worked
      example from `search()` to a prediction, the install troubleshooting
      table, and a pointer to MATERIALS_AND_METHODS.md
    - a GitHub Actions job that publishes the API page to Pages

HOW TO KNOW IT WORKED
    A person who has never seen the repo can install it, open it, and find one
    gene's evidence without asking a question.

TRAPS
    Do this AFTER task 02. It documents a UI that 02 changes, and screenshots
    taken now are screenshots to retake.
