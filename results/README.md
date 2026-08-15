# Results

Analysis output, kept as evidence for the claims made in `instructions/done/` and in the manuscript.

Deliberately **outside** `starplast/data/`. That directory is the application's cache and ships inside
the wheel, so everything in it is carried by every install; these are results a reader may want and the
program never reads.

| directory | what |
|---|---|
| `search_2026_08_12/` | the structure search that found the negative-control failure — the tables that justify excluding absence labels from `score_recovery` |
| `search_2026_08_12_corrected/` | the same battery with absence labels excluded — the numbers 3f quotes for the 3,000-gene subsample |
| `full_proteome_2026_08_12/` | all four targets over all 8,140 genes; the subsample turned out to be flattering |
| `predictions_2026_08_12/` | what the winning cell-cycle structure predicts: nothing clears the 80% purity bar |
| `full_proteome_2026_08_13_provenance/` | the full-proteome battery re-run after three leaks were found in the circularity guard — every number identical, which is the result |
| `lighting_2026_08_14/` | fixed-camera OpenGL comparisons for point modes, light moods and volumetric ray-traced shadows, plus renderer-labelled frame timings |
| `pbr_lighting_2026_08_14/` | GPU PBR sphere-impostor comparisons, GPU-vs-CPU density-ray timings, visual pixel differences, and repeated-frame stability hashes |
