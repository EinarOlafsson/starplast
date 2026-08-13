# Offer the scoring objectives in the interface — DONE 2026-08-12

Closed out. The Search tab has carried the chooser since the objectives landed: all eight objectives
from `objectives.OBJECTIVES` with their explanations, `macro`/`size` weighting, the category list
(select none to score every one), `min_recall`, and the two scoring floors — `min_cluster` and
`min_label` — which are distinct from HDBSCAN's `min_cluster_size` and were only half exposed
before. Every run records the objective it was scored under, along with chance-corrected agreement,
because a score whose objective is not beside it cannot be compared with another.

The one part that was outstanding is the clause "and in the automated walk when that exists". It
exists now (task 16), and it uses the objective: the walk keeps, per embedding, whichever clustering
scores best **under the objective in force** rather than under a fixed metric, and the frontier
column marks the configurations nothing beats on both mean and best F1.

The index had this at "100 (done)" while the task file sat in `open/`. It is in `done/` now.
