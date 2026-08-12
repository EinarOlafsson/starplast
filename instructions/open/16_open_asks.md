# Everything still asked for, in one place

The running list, so nothing depends on being remembered from a conversation. Items that are done
are recorded in `../done/`; this file is only what is outstanding, with what exists already so the
next person does not start from nothing.

## From 2026-08-12, still open

### A. Selection and navigation modes on the map  (not started)
Left-click currently rotates and picks at once. Split into a mode toggle:
- **Navigate** — free rotate, or constrained to x / y / z, plus spin.
- **Select** — draw a gate: 2D lasso in screen space, or a 3D box/brush in world space. 2D is for
  "grab that visual cluster"; 3D is for when the cloud is deep and a lasso would catch the far side.
A gated selection is the natural input to annotation: gate a cluster, see its label composition,
mark the unlabelled members. See `13_umap_gallery_and_annotation.md` §5.

### B. The UMAP gallery  (not started — the big one)
Grid and scroll modes, populating **incrementally**, expanded views behaving like the main map.
Full specification in `13_umap_gallery_and_annotation.md` §1.

### C. Automated walk: UMAP → HDBSCAN search → per-category F1  (not started)
With three optimisation targets: one category, the average, or both.
Full specification in `13_umap_gallery_and_annotation.md` §2.

### D. Annotation persistence  (not started)
Save a candidate with its cluster composition and validation numbers attached, to a separate file,
never into the node table. `validate.candidates()` already returns the rows with the composition
columns; what is missing is the store and the UI. `13_umap_gallery_and_annotation.md` §3.

### E. Logging  (not started)
Opt-in, per-level console control, rotating file. `14_logging.md`.

### F. Validation tab  (partly done)
`starplast/validate.py` is written and tested, and tab 6 exists and runs it. Outstanding:
- `refit` — re-embed per fold rather than hiding labels only from the scoring. This is the stricter
  test and the one whose number is worth quoting; it costs roughly an hour per fold on this proteome.
- Candidate list in the tab, with `cluster_frac_category` and `cluster_frac_contradicting` shown
  beside every gene.
- Orthogonal-evidence check for a candidate: signal peptide, appearance in a bait study of that
  compartment, expression tracking known members. `15_validation_tab_and_rename.md`.

### G. Retire or keep the HuggingFace dataset  (decision, not work)
The derived release is 228 KB and has been committed at `starplast/data/hf_release/` since
`0d6c13d`, so the HuggingFace copy is redundant for distribution. Keep it only if a citable DOI-able
mirror is wanted. Nothing depends on the answer.

## Standing constraints these must all respect

Measurement, inference, absence and annotation never read as one another. An annotation is a fourth
kind of thing and needs its own colour and its own file.

No candidate is shown without the number that says how much to believe it. On the full proteome,
one row out of 7,710 (configuration × compartment) reaches F1 0.5 — `PM - peripheral 1`, at exactly
0.500. A cluster that looks pure is a lead to check, not a result.
