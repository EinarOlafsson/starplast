# Initial all-source location census

All 162 registered sources were inspected against the explicit existing archive
`/media/carruthers/mnt3/claude/toxoplasma_projects/datasets`. Existing files are
hashed; installed caches, source directories, missing paths and stale links are
kept separate. GTEx v10 is bound by matching registry URL, organism, filename
and the original acquisition receipt hash; it was not downloaded again.

Status totals: 44 present processed files, one verified moved file, 40 installed
caches, ten source directories requiring member manifests, 15 without declared
paths, 23 outside-root paths/stale links, 20 requiring source/repository review,
nine failed supplement requests. There are 45 processed-source bindings.
The nine Europe PMC requests timed out. All failures are in the executed notebook
and source-location JSON; no HTML error page was installed as data.

The follow-up `../source_recovery_pmc_2026_10_07_v2/` uses the current public PMC
Cloud Service and adds 12 verified inputs. Missing sources and stale directory
links still need source-specific review under 66.01. This census is not a claim
that every raw dataset is recovered or every registry association correct.
