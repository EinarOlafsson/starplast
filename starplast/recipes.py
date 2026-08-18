#!/usr/bin/env python3
"""A recipe: one biological question, the measurements allowed to answer it, and two holdouts.

The level sweep in `search` builds one map from EVERYTHING and holds out one category at a time. That
asks *"is this category redundant given all the others"* -- a property of the CATALOGUE. It is a fair
audit and it is not a biological question, which is why its output has been hard to read biologically.

A recipe asks *"can THESE measurements predict THIS label"*. That is a hypothesis: it has a wrong
answer, it names what would refute it, and what comes back is a list of genes rather than a score.

**The validation holdout is the part that makes it trustworthy.** A second, independent label that
SHOULD map the same way is a positive control inside every run. Without it a good score says "the map
found structure"; with it, "the map found the structure I asked about" can be told apart from "the map
found something else that happens to separate genes". If the primary recovers and the control does
not, the recovery is probably an artefact of what was fed in -- and that is a result, reported beside
the inference rather than hidden.

**Hand-picking inputs REMOVES a guard, so this module adds it back.** The level sweep closed leakage
automatically by holding out a whole class; choosing inputs by hand gives that up. The danger is
predictable per question type -- for LOCALISATION it is anything sequence-derived, because signal
peptide, TM count and orthoLOPIT transfer predict localisation well enough that a recipe admitting
them "recovers" LOPIT by PREDICTING it rather than by finding biology. So before a recipe builds
anything:

* `excluded_detail` runs on the primary holdout AND on the validation holdout, at class scope;
* everything it removes from the chosen inputs is REPORTED, with the mechanism that caught it;
* a recipe whose inputs still restate a holdout after that is REFUSED, not silently trimmed;
* the two holdouts are checked against each other, because two labels from one experiment are not two
  observations and a control that is a copy controls nothing.

That last check is not hypothetical on this table: `compartment`, `compartment_best`,
`lopit_unified`, `lopit_map`, `lopit_mcmc` and `ortholopit_label` are six views of one experiment, and
any pair of them would look like a question with a control.
"""
from __future__ import annotations

import json
import os
from dataclasses import asdict, dataclass, field

import numpy as np
import pandas as pd

from .embedding import EmbeddingSpec, columns_for, embed
from .search import (ABSENCE_LABELS, CLUSTER_GRID, DEFAULT_SEED, TARGETS, best_clustering,
                     excluded_detail, map_quality, predictions, score_recovery, tune_umap)

#: How much a control label must be enriched in a cluster, over its own rate across the whole map,
#: before it counts as agreeing with the primary. Two-fold is a low bar deliberately: the control is
#: asked to CORROBORATE a cluster the primary already picked out, not to recover the map by itself,
#: and a control label that is merely at its base rate in the cluster says nothing either way.
CONTROL_ENRICHMENT = 2.0

#: The floor for scoring anything, matched to `score_recovery`'s own. A label with eight members is
#: recovered or not by luck.
MIN_LABEL = 15


# --------------------------------------------------------------------------- the object
@dataclass
class Recipe:
    """A named biological question, saved so that its answer can be rebuilt rather than remembered.

    `inputs` are category ADDRESSES rather than column names: ``("evidence", ("molecular
    measurements", "transcription"))`` names a branch of one hierarchy and takes every slot beneath
    it, so a recipe keeps meaning something after a slot is added. A bare string is taken as a block
    name, which is what an imported table or a legacy recipe uses.

    `excluded` is deliberately NOT a field a user fills in: it is what closure removed, recorded by
    the run. A recipe that carried its own exclusions could disagree with the guard that produced
    them.
    """
    question: str
    inputs: tuple = ()
    holdout: str = ""
    validation_holdout: str = ""
    organism: str = "Tg"
    #: Class scope, and NOT the `target_family` default the sweep uses -- this was measured, not
    #: assumed. `target_family` closes the family that estimates the same quantity, and it
    #: deliberately KEEPS sequence predictors: there is a test asserting exactly that. Asked "which
    #: proteins are on the PVM", a recipe scoped that way accepted `membrane topology` as an input --
    #: signal peptide, TM count -- and would have "recovered" hyperLOPIT by predicting it from
    #: sequence, which is the failure instruction 44 exists to prevent. `biology` closes the first two
    #: levels of what a measurement is ABOUT, which is what a holdout is a statement about.
    scope: str = "biology"
    seed: int = DEFAULT_SEED
    spec: dict = field(default_factory=dict)
    #: The purity a cluster must reach before its unlabelled genes are named. Part of the RECIPE and
    #: not of the run, because it is a claim about how confident an answer has to be before it counts
    #: -- and a threshold chosen after seeing the answer is a threshold fitted to it.
    min_precision: float = 0.8
    holdout_bins: int = 0          # >0 bins a QUANTITY into that many classes before scoring
    control_bins: int = 0
    axis: str = ""                 # which of the catalogue's axes this question belongs to
    expectation: str = ""          # what a believable answer looks like, written before the run
    expect_refusal: bool = False   # a recipe kept BECAUSE it is refused; see instruction 45

    def __post_init__(self):
        self.inputs = tuple(_normalize_address(a) for a in self.inputs)
        self.spec = dict(self.spec)

    @property
    def name(self) -> str:
        """The question is the name. A recipe with a second name would soon have two meanings."""
        return self.question

    def to_dict(self) -> dict:
        """The recipe as JSON-safe data. Addresses become lists, since JSON has no tuples."""
        d = asdict(self)
        d["inputs"] = [list(a) if isinstance(a, tuple) else a for a in self.inputs]
        return d

    @classmethod
    def from_dict(cls, d) -> "Recipe":
        """Rebuild a recipe from stored data, ignoring fields a later version added.

        Ignoring rather than raising: a recipe saved by a newer build must still load, or upgrading
        the program destroys the questions somebody already asked.
        """
        d = dict(d)
        d["inputs"] = tuple(_normalize_address(a) for a in d.get("inputs", ()))
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


def _normalize_address(address):
    """One address shape, whatever it arrived as.

    JSON has no tuples, so a saved recipe comes back as lists; a hand-written one is a tuple; and a
    bare block name is a string. All three become either ``str`` or ``(hierarchy, (part, ...))``.
    """
    if isinstance(address, str):
        return address
    hierarchy, path = address
    if isinstance(path, str):
        path = (path,)
    return (hierarchy, tuple(path))


def address_label(address) -> str:
    """An address as one readable string, for a report a person reads."""
    if isinstance(address, str):
        return address
    hierarchy, path = address
    return f"{hierarchy}: {' > '.join(path)}"


def address_blocks(address, organism: str = "Tg") -> tuple:
    """The blocks one address names. A slot key IS a block, which is what lets an address be input."""
    if isinstance(address, str):
        return (address,)
    from . import slots
    hierarchy, path = address
    return tuple(s.key for s in slots.slots_in_group(path, organism=organism, hierarchy=hierarchy))


def label_series(frame: pd.DataFrame, column: str, bins: int = 0) -> pd.Series:
    """One holdout column as labels, with absence expressed as absence.

    Two conversions, and both of them decide what the answer means rather than how it looks.

    **A quantity is binned or it cannot be a holdout at all.** Most of the targets worth asking about
    are numeric -- every fitness screen, every abundance, every epitope count -- and scored raw they
    produce thousands of one-member classes, none of which clears the floor, so the recipe returns
    "nothing recovered" for a question the data could have answered. Quantile bins via
    `runs.bin_column`, which reports fewer bins than asked for rather than inventing distinctions in a
    mostly-constant column.

    **"unassigned" is not a class, it is a gene waiting to be named.** hyperLOPIT writes `unassigned`
    where it could not place a protein, and left as a string it becomes the commonest label in the
    table: clusters get "recovered" as unassigned, and -- worse -- the genes this whole module exists
    to name are counted as already labelled and never predicted. Mapped to absence, they become
    exactly what they are: the candidates.
    """
    s = frame[column]
    if bins:
        from .runs import bin_column
        s = bin_column(s, bins=bins)
    text = s.astype("object").where(s.notna(), None)
    return pd.Series([None if (v is None or str(v).strip().lower() in ABSENCE_LABELS) else v
                      for v in text], index=frame.index, dtype="object")


# --------------------------------------------------------------------------- closure
@dataclass
class Closure:
    """What a recipe is allowed to build from, and everything taken away from it, with the reason.

    Reported BEFORE anything is built, because the interesting failure is a recipe that would have
    scored well on inputs it should never have been given.
    """
    blocks: tuple = ()
    columns: tuple = ()
    excluded: dict = field(default_factory=dict)      # every banned column -> why
    removed: dict = field(default_factory=dict)       # of those, the ones the inputs asked for
    dropped_blocks: dict = field(default_factory=dict)
    emptied: tuple = ()
    holdout_column: str = ""
    control_column: str = ""
    refusal: str = ""

    @property
    def ok(self) -> bool:
        return not self.refusal

    def report(self) -> str:
        """The closure as text, which is what goes in a log and in a methods section."""
        if self.refusal:
            return f"REFUSED: {self.refusal}"
        lines = [f"{len(self.blocks)} block(s), {len(self.columns)} column(s) may build this map;"
                 f" {len(self.excluded)} column(s) excluded"]
        for block, why in sorted(self.dropped_blocks.items()):
            lines.append(f"  dropped {block}: {why}")
        return "\n".join(lines)


def close(nodes: pd.DataFrame, recipe: Recipe, threshold: float = 0.8) -> Closure:
    """Resolve a recipe's inputs, remove everything that restates either holdout, and say what went.

    Four refusals, and each of them was a way this could have produced a confident wrong answer:

    * a holdout that is not a column here -- the question cannot be scored at all;
    * a control that is not independent of the primary, which is the same experiment twice;
    * an input address emptied entirely by closure, which means the user asked to predict a label
      from itself. Trimming that silently would leave a smaller recipe answering a question nobody
      asked;
    * nothing left to build from.
    """
    out = Closure()
    holdout = TARGETS.get(recipe.holdout, recipe.holdout)
    out.holdout_column = holdout
    if not holdout or holdout not in nodes.columns:
        out.refusal = f"the holdout {recipe.holdout!r} is not a column in this table"
        return out

    excluded = {c: f"primary: {why}" for c, why in
                excluded_detail(nodes, holdout, threshold=threshold, scope=recipe.scope).items()}

    control = TARGETS.get(recipe.validation_holdout, recipe.validation_holdout)
    out.control_column = control
    if control:
        if control not in nodes.columns:
            out.refusal = (f"the validation holdout {recipe.validation_holdout!r} is not a column "
                           f"in this table")
            return out
        # The control is checked against the primary's own closure. Two labels one experiment
        # produced are not two observations, and this is the mechanism that already knows it.
        if control in excluded:
            out.refusal = (f"the validation holdout {control!r} is not independent of the primary "
                           f"{holdout!r} -- {excluded[control].split(': ', 1)[1]}. A control that is "
                           f"a copy controls nothing")
            return out
        for c, why in excluded_detail(nodes, control, threshold=threshold,
                                      scope=recipe.scope).items():
            excluded.setdefault(c, f"control: {why}")
    out.excluded = excluded

    if not recipe.inputs:
        out.refusal = "the recipe names no inputs"
        return out

    kept, emptied, absent, unknown = [], [], [], []
    for address in recipe.inputs:
        blocks = address_blocks(address, recipe.organism)
        if not blocks:
            # A misspelt address must be loud. Resolved to nothing and passed over quietly, it would
            # produce a map built from the OTHER inputs and an answer to a question nobody asked.
            unknown.append(address_label(address))
            continue
        survivors, leaked = [], False
        for block in blocks:
            cols = columns_for(nodes, EmbeddingSpec(blocks=(block,))).get(block, [])
            if not cols:
                out.dropped_blocks.setdefault(block, "no columns in this cache")
                continue
            leaking = {c: excluded[c] for c in cols if c in excluded}
            if leaking:
                first = sorted(leaking)[0]
                out.dropped_blocks[block] = f"{first} -- {leaking[first]}"
                out.removed.update(leaking)
                leaked = True
                continue
            survivors.append(block)
        if not survivors:
            # Emptied by the guard and emptied by an empty cache are different facts about the
            # recipe, and reporting the second as the first would tell a user their question is
            # circular when the truth is that this build has no such data.
            (emptied if leaked else absent).append(address_label(address))
        kept.extend(survivors)

    out.emptied = tuple(emptied)
    # Order-preserving deduplication: two addresses may legitimately name the same slot, and a block
    # fed twice would be weighted twice.
    out.blocks = tuple(dict.fromkeys(kept))
    if unknown:
        out.refusal = "these inputs name nothing in this catalogue: " + "; ".join(unknown)
        return out
    if emptied:
        out.refusal = ("these inputs restate a holdout and were removed entirely rather than "
                       "trimmed: " + "; ".join(emptied))
        return out
    if absent:
        out.refusal = "these inputs have no columns in this build: " + "; ".join(absent)
        return out
    out.columns = tuple(c for cols in columns_for(
        nodes, EmbeddingSpec(blocks=out.blocks)).values() for c in cols)
    # Belt and braces, and it costs nothing: assert the survivors really are clean. A block is
    # dropped whole, so a banned column reaching this point would mean the two column lookups
    # disagree -- which is exactly the kind of silent leak this module exists to stop.
    still = sorted(set(out.columns) & set(excluded))
    if still:
        out.refusal = f"columns that restate a holdout survived closure: {', '.join(still)}"
    return out


# --------------------------------------------------------------------------- the control
def dominant_by_cluster(labels: np.ndarray, values: pd.Series, min_label: int = MIN_LABEL,
                        enrichment: float = CONTROL_ENRICHMENT) -> pd.DataFrame:
    """Which label dominates each cluster, how much of it that label holds, and whether that is more
    than the label holds everywhere.

    Written for the validation holdout and used for the primary as well, because "how pure is this
    cluster" is the same computation as "does the control corroborate it". Sharing it is what makes
    the near-miss report honest: the purity quoted when nothing clears the confidence threshold is
    computed by the same code that decides whether the control agrees, so the two can never drift.

    Per cluster: the commonest control label, how much of the cluster it holds, how much of the whole
    map that label holds, and the ratio. Enrichment rather than raw share, because a label carried by
    60% of all genes is not evidence when it is carried by 60% of a cluster.

    `agrees` is False for two quite different reasons and they are told apart in `why`: too few
    control-labelled genes in the cluster to say anything, or enough of them sitting at their base
    rate. Reporting both as "no" would let a cluster nobody could evaluate read like a cluster the
    control refuted.
    """
    labels = np.asarray(labels)
    values = values.astype("object").where(values.notna(), None).to_numpy()
    known = np.array([v is not None for v in values])
    rows = []
    base = pd.Series([str(v) for v in values[known]]).value_counts(normalize=True) if known.any() \
        else pd.Series(dtype=float)
    for cl in sorted({int(v) for v in labels if v >= 0}):
        m = (labels == cl) & known
        n = int(m.sum())
        if not n:
            rows.append({"cluster": cl, "label": "", "n_labelled": 0,
                         "share_in_cluster": float("nan"), "base_rate": float("nan"),
                         "enrichment": float("nan"), "agrees": False,
                         "why": "no gene in this cluster carries a label"})
            continue
        counts = pd.Series([str(v) for v in values[m]]).value_counts()
        top, n_top = str(counts.index[0]), int(counts.iloc[0])
        # The base rate cannot be zero here: `top` was counted among this cluster's genes, and this
        # cluster's genes are part of the base. No guard, because there is nothing to guard against.
        share, rate = n_top / n, float(base[top])
        ratio = share / rate
        if n < min_label:
            agrees, why = False, f"only {n} labelled gene(s); the floor is {min_label}"
        elif ratio >= enrichment:
            agrees, why = True, f"{top} is {ratio:.1f}x its rate across the map"
        else:
            agrees, why = False, f"{top} sits at {ratio:.1f}x its rate across the map"
        rows.append({"cluster": cl, "label": top, "n_labelled": n,
                     "share_in_cluster": share, "base_rate": rate, "enrichment": ratio,
                     "agrees": agrees, "why": why})
    return pd.DataFrame(rows)


def recovery_by_cluster(labels: np.ndarray, truth: pd.Series,
                        min_label: int = MIN_LABEL) -> pd.DataFrame:
    """Every scoreable label against every cluster: precision, recall and F1 for each pair.

    `score_recovery` answers "was this label isolated SOMEWHERE", which is the right question for a
    score and the wrong one for reading a map. A label can be split cleanly across three clusters --
    a real finding about sub-structure -- and come back with one mediocre best-cluster F1 that hides
    it. This is the table a person actually looks at, and the one instruction 43 asks for.

    Pairs with no overlap are KEPT as zeros rather than omitted. A missing row reads as "not
    computed"; a zero says the cluster was checked against that label and holds none of it, and the
    difference matters when the question is where a label is absent.
    """
    labels = np.asarray(labels)
    text = truth.astype("object").where(truth.notna(), None).to_numpy()
    known = np.array([v is not None for v in text])
    counts = pd.Series([str(v) for v in text[known]]).value_counts() if known.any() \
        else pd.Series(dtype=int)
    keep = [lab for lab, n in counts.items() if n >= min_label]
    clusters = sorted({int(v) for v in labels if v >= 0})
    rows = []
    for lab in keep:
        is_label = np.array([str(v) == lab for v in text]) & known
        for cl in clusters:
            in_cluster = labels == cl
            inter = int((is_label & in_cluster).sum())
            n_cluster_labelled = int((in_cluster & known).sum())
            prec = inter / n_cluster_labelled if n_cluster_labelled else 0.0
            rec = inter / int(is_label.sum())
            rows.append({"label": lab, "cluster": cl, "n_label": int(is_label.sum()),
                         "n_labelled_in_cluster": n_cluster_labelled, "n_in_cluster": inter,
                         "precision": prec, "recall": rec,
                         "f1": 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- running one
@dataclass
class RecipeResult:
    """Everything one recipe produced, including the reasons it produced nothing.

    A refused recipe and a degenerate map both come back as a result rather than an exception: both
    are answers to the question that was asked, and a caller that has to catch an exception to learn
    "this map cannot support a claim" will sooner or later report the exception as a failure of the
    program instead of a finding about the data.
    """
    recipe: Recipe
    closure: Closure
    quality: dict = field(default_factory=dict)
    settings: dict = field(default_factory=dict)
    summary: dict = field(default_factory=dict)
    inference_note: str = ""
    recovery: pd.DataFrame = field(default_factory=pd.DataFrame)
    per_cluster: pd.DataFrame = field(default_factory=pd.DataFrame)
    inference: pd.DataFrame = field(default_factory=pd.DataFrame)
    control: pd.DataFrame = field(default_factory=pd.DataFrame)
    control_summary: dict = field(default_factory=dict)
    genes: np.ndarray = field(default_factory=lambda: np.array([], dtype=bool))
    labels: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
    stopped_because: str = ""

    @property
    def ok(self) -> bool:
        return not self.stopped_because

    def tables(self) -> dict:
        """The result as saveable tables, for `results.save_bundle`."""
        return {"recovery": self.recovery, "per_cluster": self.per_cluster,
                "inference": self.inference, "control": self.control}


def run(nodes: pd.DataFrame, recipe: Recipe, tune: bool = True, sample: int = 1500,
        log=print) -> RecipeResult:
    """Answer one question: close leakage, tune a map, cluster it, score it, and name genes.

    The order matters and it is the argument of the whole module. Closure runs FIRST and reports
    before anything is built, so a refusal costs no compute and is legible. Tuning runs on what
    survived. The clustering is checked for degeneracy before it is scored, because a recovery number
    computed on a bisection reads exactly like a real one. Only then is the primary holdout scored,
    the unlabelled genes in pure clusters named, and the control put on the SAME clusters beside them.
    """
    closure = close(nodes, recipe)
    log(closure.report())
    if not closure.ok:
        return RecipeResult(recipe=recipe, closure=closure, stopped_because=closure.refusal)

    spec = EmbeddingSpec(blocks=closure.blocks, random_state=recipe.seed, **recipe.spec)
    settings = {}
    if tune:
        walk = tune_umap(nodes, spec, sample=sample, seed=recipe.seed, log=lambda *a: None)
        usable = walk[walk["usable"]] if len(walk) else walk
        if not len(usable):
            return RecipeResult(recipe=recipe, closure=closure,
                                stopped_because="no UMAP setting gave a map worth clustering")
        best = usable.sort_values(["stage", "score"], ascending=[False, False]).iloc[0]
        settings = {"n_neighbors": int(best["n_neighbors"]), "min_dist": float(best["min_dist"])}
        spec = EmbeddingSpec(**{**asdict(spec), **settings})
    coords, _names, kept = embed(nodes, spec, log=lambda *a: None)
    coords = np.asarray(coords)

    # A cluster that cannot hold MIN_LABEL labelled genes cannot support a claim about them, and
    # the tuner does not know that: `map_quality` rewards evenness, so on this catalogue it chose 121
    # clusters of ~40 genes -- a beautifully even map that named ZERO genes, because no cluster held
    # fifteen labelled ones at 80% purity. The floor is derived from the holdout's COVERAGE (what
    # fraction of genes carry any label at all), never from which labels recover: tuning on the
    # latter would be choosing hyperparameters by the answer.
    # Measured on the genes that are actually IN the map, not on the whole table. A missing-value
    # policy that drops rows changes which genes the labels describe, so the labelled fraction of the
    # table is the wrong denominator for a floor about the map's clusters.
    sub = nodes.loc[kept] if kept is not None else nodes
    labelled = float(label_series(sub, closure.holdout_column, recipe.holdout_bins).notna().mean())
    floor = int(np.ceil(MIN_LABEL / labelled)) if labelled else MIN_LABEL
    sizes = tuple(v for v in CLUSTER_GRID["min_cluster_size"] if v >= floor) or (floor,)
    log(f"clusters must hold >= {floor} genes for {MIN_LABEL} of them to carry a label "
        f"({100 * labelled:.0f}% of genes are labelled); trying {sizes}")
    chosen = best_clustering(coords, grid={**CLUSTER_GRID, "min_cluster_size": sizes},
                             log=lambda *a: None)
    from .clustering import cluster
    from .search import TUNE_CLUSTERING
    params = {k: chosen[k] for k in ("min_cluster_size", "min_samples")
              if k in chosen and chosen[k] != "auto"}
    params = {k: int(v) for k, v in params.items()}
    labels = cluster(coords, **{**TUNE_CLUSTERING, **params})
    quality = map_quality(labels)
    settings.update(params)

    genes = np.zeros(len(nodes), dtype=bool)
    genes[np.arange(len(nodes))[kept] if kept is not None else np.arange(len(nodes))] = True
    result = RecipeResult(recipe=recipe, closure=closure, quality=quality, settings=settings,
                          genes=genes, labels=labels)
    if not quality["usable"]:
        # A recipe on a degenerate map returns that and stops. There is no score to report and a
        # number here would be read as one.
        result.stopped_because = f"the map cannot support a claim: {quality['why_not']}"
        log(result.stopped_because)
        return result

    truth = label_series(sub, closure.holdout_column, recipe.holdout_bins)
    result.summary, result.recovery = score_recovery(labels, truth, min_label=MIN_LABEL)
    result.per_cluster = recovery_by_cluster(labels, truth)
    result.inference = predictions(sub, labels, truth, np.asarray(sub["gene_id"]),
                                   min_precision=recipe.min_precision,
                                   min_cluster_labelled=MIN_LABEL)
    if not len(result.inference):
        # An empty answer must say how close it came. "No genes" and "no cluster was 80% pure, the
        # purest was 75%" are different findings, and only the second tells a reader whether the
        # question failed or the threshold did.
        purity = dominant_by_cluster(labels, truth)
        scoreable = purity[purity.n_labelled >= MIN_LABEL]
        result.inference_note = (
            f"no cluster reached the {recipe.min_precision:.2f} purity this recipe asks for; the "
            f"purest scoreable cluster was {scoreable.share_in_cluster.max():.2f} "
            f"({scoreable.loc[scoreable.share_in_cluster.idxmax(), 'label']})"
            if len(scoreable) else
            f"no cluster holds {MIN_LABEL} labelled genes, so none can support a prediction")
        log(result.inference_note)
    if closure.control_column:
        control_truth = label_series(sub, closure.control_column, recipe.control_bins)
        result.control = dominant_by_cluster(labels, control_truth).rename(
            columns={"label": "control_label", "n_labelled": "n_control_labelled"})
        result.control_summary, _per = score_recovery(labels, control_truth, min_label=MIN_LABEL)
        # The agreement goes on the inference's own rows. A control reported in a separate table is a
        # control nobody reads next to the thing it qualifies -- which is the entire point of it.
        if len(result.inference):
            cols = ["cluster", "control_label", "n_control_labelled", "enrichment", "agrees", "why"]
            result.inference = result.inference.merge(
                result.control[cols].rename(columns={"agrees": "control_agrees",
                                                     "why": "control_says"}),
                on="cluster", how="left")
    log(f"{len(result.inference)} gene(s) named"
        + (f", {int(result.inference['control_agrees'].fillna(False).sum())} of them in clusters the "
           f"control corroborates" if "control_agrees" in result.inference else ""))
    return result


# --------------------------------------------------------------------------- keeping them
class RecipeStore:
    """Recipes on disk, so a question asked once can be asked again after the window closes.

    One JSON file per recipe, named from the question rather than a counter: the file list is then
    the question list, and two runs of the same question overwrite rather than accumulate.
    """

    def __init__(self, root: str | None = None):
        self.root = root
        self.recipes: list = []
        if root:
            os.makedirs(root, exist_ok=True)

    @staticmethod
    def filename(recipe: Recipe) -> str:
        """A filesystem-safe name from the question itself."""
        safe = "".join(c if c.isalnum() or c in " -_" else " " for c in recipe.question)
        return "_".join(safe.split())[:120].lower() + ".json"

    def add(self, recipe: Recipe) -> Recipe:
        """Keep a recipe, replacing any earlier one asking the same question."""
        self.recipes = [r for r in self.recipes if r.question != recipe.question] + [recipe]
        if self.root:
            with open(os.path.join(self.root, self.filename(recipe)), "w") as fh:
                json.dump(recipe.to_dict(), fh, indent=2)
        return recipe

    def get(self, question: str) -> Recipe | None:
        """The recipe asking this question, or None."""
        return next((r for r in self.recipes if r.question == question), None)

    def remove(self, question: str) -> bool:
        """Forget a question and delete its file. False if it was not kept here."""
        recipe = self.get(question)
        if recipe is None:
            return False
        self.recipes = [r for r in self.recipes if r.question != question]
        path = os.path.join(self.root, self.filename(recipe)) if self.root else None
        if path and os.path.exists(path):
            os.remove(path)
        return True

    def load_all(self) -> list:
        """Every saved recipe, newest last. A file that is not a recipe is skipped, not fatal."""
        if not self.root or not os.path.isdir(self.root):
            return self.recipes
        for name in sorted(os.listdir(self.root)):
            if not name.endswith(".json"):
                continue
            try:
                with open(os.path.join(self.root, name)) as fh:
                    recipe = Recipe.from_dict(json.load(fh))
            except (ValueError, TypeError, KeyError):
                continue
            if not any(r.question == recipe.question for r in self.recipes):
                self.recipes.append(recipe)
        return self.recipes
