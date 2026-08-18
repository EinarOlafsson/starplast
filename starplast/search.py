#!/usr/bin/env python3
"""Searching for structures that recover a held-out label.

This is the point of the whole apparatus. Everything else — identity, missing-value policy, block
weighting, clustering walks — exists so that this question can be asked without fooling ourselves:

    Is there a combination of datasets and hyperparameters whose structure recovers a label that was
    never given to it?

If a map built from expression and fitness alone puts 97% of apicoplast proteins in one cluster, then
compartment is *predictable from that structure*, and the unlabeled genes in that cluster acquire a
testable prediction. The same question for cell-cycle or stage tells you when a hypothetical gene is
likely active. That is inference, not description — and it is only valid because the target was held out.

Three rules make the search trustworthy, and each of them was a bug before it was a rule:

1. **The target and everything that restates it are excluded from the embedding.** Naming the target
   column is not enough: `compartment` fed one embedding and its twin `lopit_map` was duly reported as
   the top discovery at V = 0.96. `excluded_for` resolves the full set by measured association.
2. **Precision and recall are scored separately and both are kept.** A cluster that is 100% apicoplast
   but holds 5% of apicoplast proteins is useless for inference; one holding 97% at 90% precision is the
   goal. A single blended number hides exactly that difference, so the objective is per-label F1 with
   both parts stored beside it.
3. **Every run records what it would take to reproduce it** — the full `EmbeddingSpec`, the seed, the
   clustering parameters, the excluded columns, and the resulting scores. A hit nobody can rebuild is
   not a result.

The search is deliberately exhaustive rather than clever. This space has not been mapped, so the useful
first move is to cover it and look, not to optimize into one corner of it.
"""
from __future__ import annotations

import itertools
import json
import os
import time
from dataclasses import asdict, dataclass, replace

import numpy as np
import pandas as pd

from .clustering import NOISE, _association_with_inputs, cluster
from .embedding import BLOCKS, SLOT_BLOCKS, EmbeddingSpec, build_matrix, columns_for

DEFAULT_SEED = 42

# Targets worth asking about, and the columns that carry them.
TARGETS = {
    "compartment": "compartment",           # hyperLOPIT, 3,827 labelled
    "compartment_best": "compartment_best",  # + orthoLOPIT transfer
    "lopit_unified": "lopit_unified",       # 12 cross-species categories
    "attention_depth": "attention_depth",

    # Cell cycle. `cellcycle_phase` is the second MEASURED target this project has (873 genes from
    # single-parasite sequencing) and the only one that is not about location -- which is the point:
    # a structure that recovers compartment and a structure that recovers cell-cycle phase are
    # different claims, and until this column existed only the first could be made.
    "cellcycle_phase": "cellcycle_phase",
    "cellcycle_pseudotime": "cellcycle_pseudotime",

    # DERIVED. Kept as a target on purpose, as a control rather than a discovery: it is computed from
    # expression columns, so any embedding containing expression should recover it easily, and if one
    # does NOT, that says the embedding lost information it was handed. Reading a good score here as a
    # finding would be the circularity mistake this project has already made once.
    "stage_enriched_derived": "stage_enriched_derived",
}


@dataclass
class RunStep:
    """One scored configuration of a search: its numbers, its map, and its clustering.

    The automated walk's output is a table of (configuration x category) scores, and a table of
    scores is not a result -- the structure it scored is. A step carries the coordinates and the
    labels so the configuration can be shown as it finishes, which is the difference between a walk
    that can be watched and one that reports at the end. Same contract as `tuning.WalkStep`, plus
    the clustering, because here the clustering is half of what was scored.
    """
    index: int                   # 1-based, so it reads as "7 of 20"
    total: int                   # configurations this search expects to run, at most
    row: dict                    # the best-scoring row for this embedding
    per: pd.DataFrame            # its per-category precision / recall / F1
    coords: np.ndarray           # the embedding, at the same scale as any other map
    genes: np.ndarray            # boolean mask over the node table: which genes have a position
    labels: np.ndarray           # the clustering that was scored
    spec: EmbeddingSpec          # the full recipe, with this configuration's hyperparameters

    @property
    def label(self) -> str:
        """The configuration in one line, for a thumbnail caption or a status message."""
        return (f"{self.row.get('blocks', '?')}  nn={self.row.get('n_neighbors', '?')} "
                f"md={self.row.get('min_dist', '?')} mcs={self.row.get('min_cluster_size', '?')}")


def frontier(R: pd.DataFrame, columns=("mean_f1", "best_f1")) -> pd.Series:
    """Which rows are on the Pareto frontier of two objectives: nothing beats them on both.

    The task asked for the two targets "offered rather than chosen", and ranked on each. Ranking on
    each is two sorts; the frontier is the thing neither sort shows -- a configuration that is second
    on both is often the one to use, and it appears at the top of neither list. Rows not on the
    frontier are beaten outright by some other row, which is a fact about them worth having.
    """
    cols = [c for c in columns if c in R.columns]
    if R.empty or len(cols) < 2:
        return pd.Series([True] * len(R), index=R.index)
    v = R[cols].to_numpy(dtype=float)
    on = np.ones(len(R), dtype=bool)
    for i in range(len(R)):
        if not np.isfinite(v[i]).all():
            on[i] = False
            continue
        # Dominated: some other row is at least as good on both and strictly better on one.
        better = (v >= v[i]).all(axis=1) & (v > v[i]).any(axis=1)
        on[i] = not better.any()
    return pd.Series(on, index=R.index)


#: Columns that estimate the same quantity, grouped by the quantity rather than by the experiment.
#: Recovering a label is only a claim about the map if the map was not shown that label by some other
#: route, and this project measures several things twice: localization is measured by hyperLOPIT and
#: transferred from two other species; attention is tiered twice. Membership is written out rather
#: than pattern-matched, because a rule that guessed at family from a column NAME would eventually
#: throw away a real measurement for looking like the target.
SAME_QUANTITY = {
    "localization": ("compartment", "compartment_best", "compartment_source", "lopit_map",
                     "lopit_mcmc", "lopit_unified", "lopit_prob_map", "lopit_prob_mcmc",
                     "lopit_methods_agree", "lopit_confident", "ortholopit_label",
                     "ortholopit_donors", "ortholopit_accuracy", "ortholopit_accepted"),
    "attention": ("attention_depth", "lit_tier", "n_publications", "n_fulltext",
                  "n_papers_focal", "n_papers_substantive", "n_papers_incidental"),
    "cell cycle": ("cellcycle_phase", "cellcycle_pseudotime"),
}


def excluded_group(nodes: pd.DataFrame, hierarchy: str, path, organism="Tg") -> set:
    """Every declared column below a hierarchy branch.

    This is the class-level operation the slot catalogue previously promised but did not implement.
    ``hierarchy='evidence'`` omits an assay/data class, ``'biology'`` omits everything about a
    biological subject, and ``'context'`` omits a system or stage.  Dependencies are closed in both
    directions by :func:`excluded_for` when a target is supplied.
    """
    from . import slots
    out = set()
    for slot in slots.slots_in_group(path, organism=organism, hierarchy=hierarchy):
        out.update(slots.declared_columns(nodes, slot))
    return out


def excluded_for(nodes: pd.DataFrame, target: str, threshold=0.8,
                 scope="target_family") -> set:
    """The target, everything that restates it, and everything the same experiment produced.

    Three mechanisms, because each is blind to what the others catch and this guard has now leaked
    three separate ways:

    * **measured association**, for the undeclared copy -- a renamed column, a second inference over
      the same data. The threshold is stricter than the 0.95 used to flag derived columns after the
      fact, because here we are choosing what an embedding may *see*, and a 0.85-associated column
      leaks nearly as much as an identical one: `lopit_mcmc` is a second inference over the same
      experiment at 0.74;
    * **declared derivation**, for the joint function -- a label that is the argmax of three columns
      shows 0.56-0.66 against each of them separately, so no pairwise statistic can see it;
    * **shared provenance**, for the same experiment's OTHER outputs, which is neither of the above.
    """
    if scope not in ("direct", "target_family", "evidence", "biology", "context"):
        raise ValueError("scope must be direct, target_family, evidence, biology or context")
    assoc = _association_with_inputs(nodes, {target})
    out = {target} | {c for c, v in assoc.items() if v >= threshold}

    # The generated hierarchy is the primary declaration.  A target family is the narrow,
    # leakage-safe default (all direct estimates of the same quantity).  Broader scopes are explicit
    # sensitivity analyses: omit the whole assay class, biological subject or context branch.
    from . import slots
    matched = slots.target_slots(target)
    if scope != "direct" and matched:
        if scope == "target_family":
            selected = {s.key: s for seed in matched
                        for s in slots.family_slots(seed.target_family)}.values()
        else:
            selected = {}
            for seed in matched:
                path = getattr(seed, f"{scope}_path")
                # The final component is assay-specific.  Its parent is the useful class-level
                # holdout (e.g. all spatial-localisation evidence, not only hyperLOPIT posterior 1).
                selected.update({s.key: s for s in
                                 slots.slots_in_group(path[:-1] or path, hierarchy=scope)})
            selected = selected.values()
        for slot in selected:
            out.update(slots.declared_columns(nodes, slot))

    # Shared provenance. `lopit_prob_map` is the posterior of hyperLOPIT's own assignment: not a
    # restatement of the compartment (association 0.29, far under any workable threshold), and not a
    # declared source of it either, since the label is not computed from the posterior. It is the
    # same experiment's other output, and a map built on it is being scored against a label that
    # experiment also produced.
    #
    # Not a small effect. On the shipped cache the three-column `localization` block recovered
    # `compartment` at mean F1 0.259, above interactions (0.192), protein features (0.171) and every
    # expression and fitness block. Three columns beating eighteen RNA columns and eight CRISPR
    # screens is not the map finding biology.
    from . import datasets
    same = datasets.provenance(target)
    if same is not None:
        out |= {c for c in same.columns if c in nodes.columns}

    # The same QUANTITY, however it was arrived at. Provenance is about which experiment produced a
    # column; this is about what the column is an estimate OF, and the two come apart wherever the
    # project estimates something twice. `ortholopit_label` is a localization label transferred from
    # P. falciparum and C. parvum orthologs -- a different experiment, in a different species,
    # estimating the same thing, at association 0.72 with `compartment`; `lit_tier` is a second
    # tiering of how much a gene has been written about, at 0.71 with `attention_depth`. Both sit
    # just under the exclusion threshold, which is what a near-copy does.
    for family in SAME_QUANTITY.values():
        if target in family:
            out |= {c for c in family if c in nodes.columns}

    # Anything the target was DECLARED to be computed from, plus the rest of that column's block.
    # Measured association is pairwise and cannot see a label that is a joint function of several
    # columns: stage_enriched_derived scores 0.56-0.66 against each of its three sources, under any
    # workable threshold, while being a deterministic function of all three together. The block is
    # taken as a whole because the sources do not stand alone either -- the strongest association to
    # that derived label, 0.72, is a tissue-cyst FPKM column that is not one of its declared sources
    # but measures the same biology.
    from . import datasets
    declared = set(datasets.derived_sources(target))
    if declared:
        # A declared source may itself be a summary of a registered experiment.  Now that the raw
        # GSE columns are selectable, excluding expr_cyst while leaving the twelve measurements it
        # summarizes would put the same quantity straight back into the map by a longer route.
        # Provenance closes that route without guessing from column-name prefixes.
        for source in declared:
            origin = datasets.provenance(source)
            if origin is not None:
                out |= {c for c in origin.columns if c in nodes.columns}
        for block in BLOCKS:
            cols = set(columns_for(nodes, EmbeddingSpec(blocks=(block,))).get(block, []))
            if cols & declared:
                out |= cols
        out |= declared

    # And the reverse direction: a summary or label computed from anything already banned is also
    # banned.  Repeat to a fixed point because a derived output can itself feed another derivation.
    changed = True
    while changed:
        before = len(out)
        out |= {c for c in datasets.derived_dependents(out) if c in nodes.columns}
        changed = len(out) != before
    return out


def _spec_without(spec: EmbeddingSpec, nodes: pd.DataFrame, banned: set) -> EmbeddingSpec:
    """Drop any block or categorical that would feed a banned column into the embedding."""
    keep_blocks = []
    for b in spec.blocks:
        cols = columns_for(nodes, EmbeddingSpec(blocks=(b,))).get(b, [])
        if cols and not (set(cols) & banned):
            keep_blocks.append(b)
    keep_cat = tuple(c for c in spec.categorical if c not in banned)
    return EmbeddingSpec(**{**asdict(spec), "blocks": tuple(keep_blocks), "categorical": keep_cat})


# Label values that mean "not measured" rather than naming a class. A structure that separates these
# has recovered which genes were MEASURED, not what they are -- and measurement tracks study effort, so
# the score is about the literature rather than the biology.
#
# This is not a hypothetical. Scored with them included, `compartment` reaches mean F1 0.484 and its
# single best-recovered label is `unassigned` at 0.39, while the real compartments sit at 0.21-0.35;
# excluding it, the same run scores 0.207. The negative control makes the point unarguable:
# `attention_depth` scores 0.654, of which essentially all comes from recovering the never-named class
# at F1 0.77, its real tiers scoring 0.14-0.35.
ABSENCE_LABELS = frozenset({"unassigned", "unknown", "", "nan", "none", "unlabelled", "unlabeled"})


def score_recovery(labels: np.ndarray, truth: pd.Series, min_label=15,
                   exclude_labels=ABSENCE_LABELS) -> tuple:
    """How well does the clustering recover `truth`? Per-label precision, recall and F1.

    For each label the best single cluster is taken -- the question is whether the structure isolates
    the label somewhere, not whether the clustering happens to use the same number of groups.
    """
    ok = (labels != NOISE) & truth.notna().to_numpy()
    if ok.sum() < 50:
        return {}, pd.DataFrame()
    lab, tru = labels[ok], truth[ok].astype(str).to_numpy()
    rows = []
    for t in pd.unique(tru):
        if str(t).strip().lower() in exclude_labels:
            continue
        n_t = int((tru == t).sum())
        if n_t < min_label:
            continue
        best = None
        for cl in np.unique(lab):
            inter = int(((lab == cl) & (tru == t)).sum())
            if not inter:
                continue
            prec = inter / int((lab == cl).sum())
            rec = inter / n_t
            f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
            if best is None or f1 > best["f1"]:
                best = {"label": t, "cluster": int(cl), "n_label": n_t,
                        "n_in_cluster": inter, "precision": prec, "recall": rec, "f1": f1}
        if best:
            rows.append(best)
    per = pd.DataFrame(rows)
    if per.empty:
        return {}, per
    # Weighted by label size, so a structure that nails one tiny class does not outrank one that
    # organises the whole proteome.
    w = per.n_label / per.n_label.sum()
    from . import gpu
    summary = {"backend": gpu.backend_id(),
               "mean_f1": float((per.f1 * w).sum()),
               "best_f1": float(per.f1.max()),
               "best_label": str(per.loc[per.f1.idxmax(), "label"]),
               "n_labels_scored": int(len(per)),
               "n_labels_recovered": int((per.f1 >= 0.5).sum())}
    return summary, per


def _flushing(log):
    """A long search redirected to a file shows nothing until it ends unless output flushes."""
    def out(m):
        log(m)
        try:
            import sys; sys.stdout.flush()
        except Exception:
            pass
    return out


def search(nodes: pd.DataFrame, target: str = "compartment",
           block_sets=None, na_policies=("median",), scalings=("rank",),
           n_neighbors_values=(15, 50), min_dist_values=(0.0, 0.25),
           min_cluster_sizes=(25, 60), sample_size: int | None = None,
           seed: int = DEFAULT_SEED, store=None, save_above: float | None = None,
           objective: dict | None = None, on_run=None, should_stop=None, log=print) -> tuple:
    """Walk combinations of datasets and hyperparameters, scoring recovery of a held-out target.

    `objective` selects what "good structure" means, as the keyword arguments `objectives.score`
    takes -- which of precision, recall or both, averaged or best-case, and over every category or
    one. Omitted, the per-label F1 summary is used, which is what this always did.

    Whichever is chosen, the per-label table is computed regardless and returned alongside, because
    no single number survives contact with a real result: the winning configuration is chosen by the
    objective, and then read per label to see whether it earned it.

    `should_stop` is asked before every configuration, and a search that is stopped RETURNS what it
    has -- ranked, with its per-category table -- rather than raising. A stop is a decision that
    enough has been seen, not an error, and throwing away twenty minutes of finished configurations
    because the user pressed stop on the twenty-first would make the button unusable. Asked per
    configuration rather than per log line, because the log is throttled to every fortieth run and a
    stop that takes four minutes to land reads as a button that does not work.

    `on_run` is called with a `RunStep` as each EMBEDDING finishes -- carrying its best clustering
    among the `min_cluster_sizes` tried, its per-category scores, its coordinates and its labels. Per
    embedding rather than per row, because the rows for one embedding differ only in the clustering
    and a gallery of the same map five times is not a gallery. Every row still reaches the returned
    table: what is streamed is what can be looked at, not a subset of what was scored.
    """
    log = _flushing(log)
    truth_col = TARGETS.get(target, target)
    if truth_col not in nodes.columns:
        raise ValueError(f"target {truth_col!r} not in the table")
    banned = excluded_for(nodes, truth_col)
    log(f"target {truth_col!r}: {int(nodes[truth_col].notna().sum()):,} labeled genes")
    log(f"  excluded from every embedding ({len(banned)}): {', '.join(sorted(banned))}")

    if block_sets is None:
        # Six biological questions, not six file families. Localization remains held out by default
        # because it is commonly the target; literature attention is never a feature.
        base = ["Tg_transcription_tachyzoite", "Tg_transcription_bradyzoite_tissue_cyst",
                "Tg_transcription_oocyst_sporozoite", "Tg_fitness_hff_in_vitro",
                "Tg_protein_abundance_tachyzoite", "Tg_fold_confidence_disorder"]
        base = [b for b in base if columns_for(nodes, EmbeddingSpec(blocks=(b,))).get(b)]
        # A minimal/imported table may predate the slot catalogue and carry only a legacy family
        # such as ``fit_*``. Keep that table searchable; saved recipes retain their old block name.
        if not base:
            legacy = ["expression_summary", "expression_raw", "fitness_screens",
                      "published_screens", "protein_features", "interactions"]
            base = [b for b in legacy if columns_for(
                nodes, EmbeddingSpec(blocks=(b,))).get(b)]
        block_sets = [tuple(c) for r in (1, 2, 3) for c in itertools.combinations(base, r)]
        catalog = SLOT_BLOCKS if any(b in SLOT_BLOCKS for b in base) else BLOCKS
        left_out = [b for b in catalog if b not in base and columns_for(
            nodes, EmbeddingSpec(blocks=(b,))).get(b)]
        log(f"  {len(block_sets)} dataset combinations from {len(base)} blocks: "
            + ", ".join(base)
            + (f"  (not swept: {', '.join(left_out)})" if left_out else ""))

    rows, per_label, runs, last_reported, emitted = [], [], 0, 0, 0
    # Combinations the sweep never ran, and why. A skip that is not counted turns "8 runs" into "7
    # runs" with nothing to explain the difference, and the commonest reason for one is now the
    # circularity guard removing every block a combination names -- which is exactly the thing a
    # reader of the table needs told.
    skipped = {}
    total = (len(block_sets) * len(na_policies) * len(scalings)
             * len(n_neighbors_values) * len(min_dist_values) * len(min_cluster_sizes))
    # How many EMBEDDINGS the walk expects, which is what a step counts against: the clusterings of
    # one embedding are the same map scored several ways.
    embeddings = max(total // max(len(min_cluster_sizes), 1), 1)
    log(f"  {total} runs")
    t0 = time.time()
    sub = None
    if sample_size and sample_size < len(nodes):
        sub = np.random.default_rng(seed).choice(len(nodes), sample_size, replace=False)
        sub.sort()

    stopped = False
    for blocks, pol, sc in itertools.product(block_sets, na_policies, scalings):
        if should_stop is not None and should_stop():
            stopped = True
            break
        spec0 = EmbeddingSpec(blocks=blocks, na_policy=pol, scaling=sc, random_state=seed)
        spec0 = _spec_without(spec0, nodes, banned)
        if not spec0.blocks:
            skipped.setdefault("every block they name feeds a column the guard excluded",
                               []).append("+".join(blocks))
            continue
        try:
            X, names, keep = build_matrix(nodes, spec0, log=lambda *a: None)
        except ValueError:
            skipped.setdefault("no usable feature matrix", []).append("+".join(blocks))
            continue
        idx = np.arange(len(nodes))[keep]
        if sub is not None:
            sel = np.isin(idx, sub)
            X, idx = X[sel], idx[sel]
        if len(X) < 200:
            skipped.setdefault("fewer than 200 genes survive the missing-value policy",
                               []).append("+".join(blocks))
            continue
        try:
            import umap
        except ImportError:
            log("umap-learn not installed"); return pd.DataFrame(), pd.DataFrame()

        for nn, md in itertools.product(n_neighbors_values, min_dist_values):
            if should_stop is not None and should_stop():
                stopped = True
                break
            if nn >= len(X):
                continue
            Y = np.array(umap.UMAP(n_components=3, n_neighbors=nn, min_dist=md,
                                     metric="euclidean", random_state=seed).fit_transform(X), copy=True)
            made = []
            for mcs in min_cluster_sizes:
                if mcs >= len(Y):
                    # HDBSCAN raises rather than returning all-noise when min_cluster_size exceeds the
                    # sample. One unusable grid point should cost that point, not the whole walk.
                    continue
                lab = cluster(Y, algorithm="hdbscan", min_cluster_size=mcs)
                truth = nodes[truth_col].iloc[idx]
                summary, per = score_recovery(lab, truth)
                runs += 1
                if not summary:
                    continue
                if objective:
                    # Scored under the chosen objective, and the choice is recorded on the row. A
                    # score whose objective is not beside it cannot be compared with another.
                    from .objectives import agreement, score as score_objective
                    obj = score_objective(lab, truth, **objective)
                    summary = dict(summary)
                    summary["objective_score"] = obj["score"]
                    summary["objective"] = obj["objective"]
                    summary["objective_detail"] = obj["detail"]
                    summary["weighting"] = obj["weighting"]
                    # Chance-corrected agreement travels with it, because every objective here has a
                    # degenerate maximiser and an ARI near zero beside a high score is how you see one.
                    summary.update(agreement(lab, truth))
                row = {"target": truth_col, "blocks": "+".join(spec0.blocks),
                       "na_policy": pol, "scaling": sc, "n_neighbors": nn, "min_dist": md,
                       "min_cluster_size": mcs, "seed": seed,
                       # The subsample this run drew, recorded so the run can be rebuilt exactly.
                       # `n_genes` is what survived the missing-value policy afterwards and is not
                       # the same number, so reconstructing the draw from it silently produced a
                       # different set of genes -- and therefore a different map from the row.
                       "sample_size": int(sample_size or len(nodes)),
                       "n_genes": int(len(X)),
                       "n_features": X.shape[1],
                       "n_clusters": int(len(set(lab[lab != NOISE]))),
                       "noise_frac": float((lab == NOISE).mean()),
                       # Recorded on every row so a saved table can be audited on its own. Reading a
                       # results CSV months later, "was the target held out" is the first question,
                       # and it should not require re-running the search to answer.
                       "n_excluded": len(banned),
                       "excluded": ";".join(sorted(banned)),
                       **summary}
                rows.append(row)
                per = per.assign(**{k: row[k] for k in
                                    ("blocks", "na_policy", "scaling", "n_neighbors",
                                     "min_dist", "min_cluster_size", "target")})
                per_label.append(per)
                made.append((row, per, lab))
                # Default to saving the best runs relative to what this search actually found,
                # not an absolute bar. A fixed 0.6 threshold saved nothing at all on the first
                # real search, whose best was 0.592 -- the gate silently discarded the answer.
                keep_at = save_above if save_above is not None else 0.0
                if store is not None and summary["best_f1"] >= keep_at:
                    name = (f"{truth_col}_{'+'.join(spec0.blocks)}_{pol}_{sc}"
                            f"_nn{nn}_md{md}_mcs{mcs}")
                    # The stored spec has to carry the hyperparameters this run ACTUALLY used. spec0
                    # holds the block, policy and scaling, but n_neighbors and min_dist are loop
                    # variables -- saved unchanged, every stored recipe claimed the EmbeddingSpec
                    # defaults (25 and 0.25) whatever the run did. The name encoded the truth and the
                    # machine-readable field did not, which is the worse half to get wrong: reopening
                    # a saved embedding by its recipe would have rebuilt a different map.
                    used = EmbeddingSpec(**{**asdict(spec0), "n_neighbors": nn, "min_dist": md})
                    store.save(name, Y, used, gene_ids=nodes.gene_id.iloc[idx], features=names,
                               extra={"target": truth_col, "excluded": sorted(banned),
                                      "clustering": {"algorithm": "hdbscan",
                                                     "min_cluster_size": mcs},
                                      "scores": summary})
            if made and on_run is not None:
                # The best clustering of this embedding, by whichever objective is in force -- the
                # HDBSCAN search is a search, and what it is for is keeping the winner.
                key = "objective_score" if objective else "mean_f1"
                best_row, best_per, best_lab = max(made, key=lambda m: m[0].get(key, float("-inf")))
                emitted += 1
                genes = np.zeros(len(nodes), dtype=bool)
                genes[idx] = True
                used = EmbeddingSpec(**{**asdict(spec0), "n_neighbors": nn, "min_dist": md})
                from .embedding import normalize
                on_run(RunStep(index=emitted, total=embeddings, row=best_row, per=best_per,
                               coords=normalize(Y), genes=genes, labels=best_lab, spec=used))
        # Report on crossing each multiple of 40 rather than on exact equality. The check sits at the
        # end of a dataset combination, so `runs` jumps by however many hyperparameter points that
        # combination had: equality only ever fired when that stride happened to divide 40, and with a
        # different grid shape a long walk printed nothing at all and looked like a hang.
        if runs - last_reported >= 40:
            last_reported = runs
            log(f"    {runs}/{total} runs, {time.time() - t0:.0f}s")

    R = pd.DataFrame(rows).sort_values("mean_f1", ascending=False) if rows else pd.DataFrame()
    if not R.empty:
        # Marked rather than filtered: a dominated configuration is still a result, and the column
        # says which ones nothing beats on both objectives at once.
        R["on_frontier"] = frontier(R)
    if store is not None and not R.empty and save_above is None:
        log(f"  saved embeddings for every run; the top result is the first row of the table")
    P = pd.concat(per_label, ignore_index=True) if per_label else pd.DataFrame()
    if stopped:
        log(f"  STOPPED after {runs} of {total} runs -- what finished is kept, ranked and saved")
    for why, which in skipped.items():
        names = sorted(set(which))
        log(f"  skipped {len(names)} combination(s) -- {why}: "
            + ", ".join(names[:6]) + (f" and {len(names) - 6} more" if len(names) > 6 else ""))
    log(f"  done: {runs} runs in {time.time() - t0:.0f}s")
    if not R.empty:
        b = R.iloc[0]
        log(f"  best: mean F1 {b.mean_f1:.3f} (best label {b.best_label} at F1 {b.best_f1:.3f}) "
            f"from {b.blocks} nn={b.n_neighbors} md={b.min_dist} mcs={b.min_cluster_size}")
    return R, P


def rebuild(nodes: pd.DataFrame, row, log=print) -> tuple:
    """Rebuild one row of a search result: the same map, on the same genes, with its clustering.

    The inverse of what `search` records, and it lives here so it cannot drift from the loop that
    wrote the row. Every step is taken the way the search took it, because a rebuild that differs in
    any of them puts a different map on screen from the one the row's numbers describe:

    * the feature matrix is built over the WHOLE table and then subset, not built over the
      subsample -- rank scaling over 3,000 genes is not rank scaling over 8,140 restricted to them;
    * the subsample is redrawn from the recorded seed and `sample_size`, then sorted;
    * the columns the run excluded are excluded again, from the row rather than recomputed, since
      the exclusion is what makes the score mean anything;
    * the clustering uses the row's own `min_cluster_size`.

    Returns `(coords, genes, labels, features)`, where `genes` is a boolean mask over the node table
    saying which genes have a position -- the rest are not in this map and must not be drawn in it.
    """
    get = (lambda k, d=None: row.get(k, d))
    blocks = tuple(b for b in str(get("blocks", "")).split("+") if b)
    if not blocks:
        raise ValueError("that row names no feature blocks")
    seed = int(float(get("seed", DEFAULT_SEED)))
    spec = EmbeddingSpec(blocks=blocks, na_policy=str(get("na_policy", "median")),
                         scaling=str(get("scaling", "rank")),
                         n_neighbors=int(float(get("n_neighbors", 15))),
                         min_dist=float(get("min_dist", 0.1)), random_state=seed)
    banned = {c for c in str(get("excluded", "")).split(";") if c}
    spec = _spec_without(spec, nodes, banned)
    if not spec.blocks:
        raise ValueError("every block in that row feeds a column the run excluded")
    X, names, keep = build_matrix(nodes, spec, log=lambda *a: None)
    idx = np.arange(len(nodes))[keep]
    size = int(float(get("sample_size", 0) or 0))
    if size and size < len(nodes):
        sub = np.random.default_rng(seed).choice(len(nodes), size, replace=False)
        sub.sort()
        sel = np.isin(idx, sub)
        X, idx = X[sel], idx[sel]
    log(f"rebuilding {'+'.join(spec.blocks)} over {len(X):,} genes, seed {seed}")
    import umap
    from .embedding import normalize
    Y = np.asarray(umap.UMAP(n_components=3, n_neighbors=spec.n_neighbors,
                             min_dist=spec.min_dist, metric="euclidean",
                             random_state=seed).fit_transform(X))
    mcs = int(float(get("min_cluster_size", 25)))
    log(f"clustering at min_cluster_size={mcs}")
    # Clustered on the raw coordinates, displayed normalized: normalizing is a uniform scaling, so
    # it cannot change the clustering, and doing it in this order keeps that guarantee obvious.
    labels = cluster(Y, algorithm="hdbscan", min_cluster_size=mcs)
    genes = np.zeros(len(nodes), dtype=bool)
    genes[idx] = True
    return normalize(Y), genes, labels, names


def predictions(nodes: pd.DataFrame, labels: np.ndarray, truth: pd.Series, gene_index,
                min_precision=0.8, min_cluster_labelled=10) -> pd.DataFrame:
    """Turn a recovered structure into predictions for the unlabeled genes in each cluster.

    Only clusters that are already purely one label are used, and the precision achieved on the
    *labeled* members is carried through as the prediction's confidence -- it is the only honest
    estimate of how often the prediction will be right.
    """
    out = []
    tru = truth.to_numpy()
    for cl in np.unique(labels[labels != NOISE]):
        m = labels == cl
        known = m & pd.notna(tru)
        if known.sum() < min_cluster_labelled:
            continue
        vals, counts = np.unique(tru[known].astype(str), return_counts=True)
        top, n_top = vals[counts.argmax()], counts.max()
        prec = n_top / known.sum()
        if prec < min_precision:
            continue
        for gi in np.asarray(gene_index)[m & pd.isna(tru)]:
            out.append({"gene_id": gi, "predicted": top, "cluster": int(cl),
                        "cluster_precision": float(prec),
                        "n_labelled_in_cluster": int(known.sum())})
    return pd.DataFrame(out)


# --------------------------------------------------------------------------- the category sweep
def categories_at(hierarchy: str = "evidence", level: int = 1, organism: str = "Tg",
                  blocks=None) -> list:
    """Every category address at one LEVEL of a hierarchy, with the blocks beneath it.

    Returns ``[(path, [block, ...]), ...]``. A level rather than the whole tree, because the question
    "does this class of evidence recover on a map built without it" is asked of a class, and which
    level counts as a class is the reader's choice: level 1 is "molecular measurements", level 3 is
    "transcript abundance". Both are real questions and they have different answers.
    """
    from . import slots
    keep = set(blocks) if blocks is not None else None
    out = {}
    for slot in slots.all_slots(organism):
        path = slots.hierarchy_path(slot, hierarchy)
        block = path[-1]
        if keep is not None and block not in keep:
            continue
        if len(path) <= level:                       # a slot shallower than the level asked for
            continue
        out.setdefault(tuple(path[:level]), []).append(block)
    return [(path, sorted(blocks)) for path, blocks in sorted(out.items())]


def sweep_categories(nodes: pd.DataFrame, spec, hierarchy: str = "evidence", level: int = 1,
                     organism: str = "Tg", only=None, algorithm: str = "hdbscan",
                     min_cluster_size: int = 25, should_stop=None, on_result=None,
                     log=print) -> pd.DataFrame:
    """Hold out each category in turn, rebuild the map without it, and ask whether it comes back.

    This is the project's whole argument, run as a procedure instead of by hand. For each category:
    build the embedding from every OTHER block, cluster it, then test the held-out category's columns
    against those clusters with `clustering.battery` -- which marks anything the map actually saw as
    `used` or `derived`, so only genuinely held-out evidence is scored.

    A category that recovers strongly is structure the rest of the evidence already implies. One that
    does not is either independent information or noise, and the two are told apart by looking, not by
    this function -- which is why the per-feature table comes back rather than only a score.

    `should_stop` is asked between categories and a stopped sweep RETURNS what it has, for the same
    reason `search` does: a stop is a decision that enough has been seen, not an error.
    """
    from .clustering import battery, cluster
    from .embedding import EmbeddingSpec, columns_for, embed
    groups = categories_at(hierarchy, level, organism, blocks=spec.blocks)
    if only:
        wanted = {tuple(p) if isinstance(p, (list, tuple)) else (p,) for p in only}
        groups = [(path, blocks) for path, blocks in groups if path in wanted]
    rows = []
    for path, held in groups:
        if should_stop is not None and should_stop():
            log(f"sweep stopped after {len(rows)} of {len(groups)} categories")
            break
        remaining = tuple(b for b in spec.blocks if b not in set(held))
        if not remaining:
            log(f"{' > '.join(path)}: skipped, it is everything that was ticked")
            continue
        without = replace(spec, blocks=remaining)
        held_columns = sorted({c for cols in
                               columns_for(nodes, replace(spec, blocks=tuple(held))).values()
                               for c in cols})
        if not held_columns:
            log(f"{' > '.join(path)}: skipped, no columns in this cache")
            continue
        used = [c for cols in columns_for(nodes, without).values() for c in cols]
        log(f"{' > '.join(path)}: holding out {len(held_columns)} columns, "
            f"building from {len(used)}")
        # `embed` returns (coords, feature names, kept rows). The kept rows matter: a policy that
        # drops genes with missing values means the labels describe a SUBSET, and scoring them against
        # the whole table would align cluster 3 with the wrong genes.
        coords, _names, kept = embed(nodes, without, log=lambda *a: None)
        labels = cluster(np.asarray(coords), algorithm=algorithm, min_cluster_size=min_cluster_size)
        sub = nodes.loc[kept] if kept is not None else nodes
        scored, _detail = battery(sub, labels, used_features=used, features=held_columns,
                                  log=lambda *a: None)
        held_rows = scored[scored.evidence == "held_out"] if len(scored) else scored
        best = float(held_rows["score"].max()) if len(held_rows) else float("nan")
        row = {"category": " > ".join(path), "hierarchy": hierarchy, "level": level,
               "blocks_held_out": len(held), "columns_held_out": len(held_columns),
               "columns_used": len(used), "clusters": int(len({int(x) for x in labels if x >= 0})),
               "tested": int(len(held_rows)), "best_score": best,
               "recovered": bool(len(held_rows) and best >= 0.30)}
        rows.append(row)
        if on_result is not None:
            on_result(row, held_rows)
    return pd.DataFrame(rows)
