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
from .search import (ABSENCE_LABELS, CLUSTER_GRID, DEFAULT_SEED, TARGETS, all_edge_layers,
                     best_clustering, excluded_detail, excluded_edges, map_quality, score_recovery,
                     tune_umap)

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
    #: How much a label must be CONCENTRATED in a cluster, over its rate across the whole map,
    #: before that cluster may name its unlabelled genes. This is the gate, and precision is not.
    #:
    #: Measured, and the reason is that a fixed purity bar is a different statistical demand for
    #: every question. `compartment` has 24 classes and a typical one covers 3% of labelled genes, so
    #: 0.80 purity there is a 25-fold enrichment; a quantity binned into thirds has 33% per class, so
    #: the same 0.80 is 2.4-fold. Gating on purity therefore asks almost nothing of a rare label and
    #: nearly the impossible of a common one, which is an accident of how many classes the holdout
    #: happens to have rather than a fact about biology. Enrichment means the same thing at any base
    #: rate, and it is already what the CONTROL is judged by -- the primary deserves the same test.
    min_enrichment: float = 2.0
    #: The purity a cluster must reach as well, carried onto every named gene as its confidence.
    #: Kept as a floor rather than the gate because enrichment and reliability are different
    #: questions: a cluster 15% X where X is 3% everywhere is five-fold enriched and still wrong
    #: about five genes in six. Enrichment says the signal is real; precision says how often an
    #: individual prediction will be right, and a reader needs both numbers.
    #:
    #: **0.30 was chosen by measurement, and the measurement was a surprise.** Every shipped question
    #: was run once and 25 candidate gates applied to the same clusterings, scored not by how many
    #: genes they named but by how many of those the INDEPENDENT CONTROL corroborated
    #: (`instructions/done/45_gate_choice.csv`). At enrichment >= 2 the corroboration rate falls
    #: monotonically as this floor rises -- 0.3: 46.0%, 0.4: 41.5%, 0.5: 31.7%, 0.6: 22.3%,
    #: 0.8: 9.9% -- so the purest clusters are precisely the ones the control does NOT back. The
    #: reading: high purity is reached by small tight clusters of rare labels, where a sparse control
    #: has too few genes to agree; the moderately pure, strongly enriched clusters are the large
    #: coherent ones where it can. Raising this floor was selecting against corroborated answers.
    #:
    #: It still does work rather than being switched off: a 24-class holdout at 2x enrichment reaches
    #: only ~6% purity, and naming genes from a cluster that is 6% one label would be wrong 94 times
    #: in a hundred. 0.30 is where reliability stops being free and starts costing corroboration.
    min_precision: float = 0.3
    holdout_bins: int = 0          # >0 bins a QUANTITY into that many classes before scoring
    control_bins: int = 0
    #: How the question is answered. Every method produces a PARTITION -- clusters, or out-of-fold
    #: predicted classes -- and everything after that is identical, so two methods on one question are
    #: compared on the same closure, the same control and the same statistics rather than on their
    #: own separate write-ups.
    method: str = "umap+hdbscan"
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
    #: Graph layers, which are a leak surface of their own. A method that traverses edges rather than
    #: reading columns can recover a label from the layer built out of it, and no column guard can
    #: see that -- `compartment` is 118,712 edges as well as the commonest holdout here.
    layers: tuple = ()                                # the edge layers a graph method may traverse
    excluded_layers: dict = field(default_factory=dict)
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
                 f" {len(self.excluded)} column(s) excluded",
                 f"{len(self.layers)} edge layer(s) may be traversed; "
                 f"{len(self.excluded_layers)} excluded"
                 + (": " + ", ".join(sorted(self.excluded_layers)) if self.excluded_layers else "")]
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
    banned_layers = {L: f"primary: {why}"
                     for L, why in excluded_edges(holdout, recipe.scope).items()}

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
        for L, why in excluded_edges(control, recipe.scope).items():
            banned_layers.setdefault(L, f"control: {why}")
    out.excluded = excluded
    out.excluded_layers = banned_layers
    out.layers = tuple(L for L in all_edge_layers(recipe.organism) if L not in banned_layers)

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


def infer(labels: np.ndarray, truth: pd.Series, gene_index, min_enrichment: float = 2.0,
          min_precision: float = 0.3, min_labelled: int = MIN_LABEL) -> pd.DataFrame:
    """Name the unlabelled genes in clusters where a label is genuinely concentrated.

    Two gates rather than one, because they answer different questions and a single number hides
    whichever it is not. Enrichment says the label is more concentrated here than at large -- that the
    cluster carries signal about it. Precision says how often a prediction drawn from this cluster
    will actually be right, and it rides on every row as the prediction's confidence.

    Built on `dominant_by_cluster`, which is also what scores the control, so "is this label
    concentrated in this cluster" is computed once and means one thing throughout.
    """
    dom = dominant_by_cluster(labels, truth, min_label=min_labelled, enrichment=min_enrichment)
    labels = np.asarray(labels)
    values = truth.astype("object").where(truth.notna(), None).to_numpy()
    unlabelled = np.array([v is None for v in values])
    genes = np.asarray(gene_index)
    out = []
    for _, r in dom.iterrows():
        if not r.agrees or r.share_in_cluster < min_precision:
            continue
        for gene in genes[(labels == r.cluster) & unlabelled]:
            out.append({"gene_id": gene, "predicted": r.label, "cluster": int(r.cluster),
                        "cluster_precision": float(r.share_in_cluster),
                        "enrichment": float(r.enrichment), "base_rate": float(r.base_rate),
                        "n_labelled_in_cluster": int(r.n_labelled)})
    return pd.DataFrame(out)


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
    #: Which measurements a model actually used, when the method is one that can say. Empty for the
    #: clustering path, which cannot.
    coefficients: pd.DataFrame = field(default_factory=pd.DataFrame)
    genes: np.ndarray = field(default_factory=lambda: np.array([], dtype=bool))
    labels: np.ndarray = field(default_factory=lambda: np.array([], dtype=int))
    #: The embedding itself, and the two label columns as the map saw them. Kept because a figure
    #: rebuilt from a seed is a figure nobody checked: the report draws the coordinates that were
    #: actually scored, not a re-run that should agree with them.
    coords: np.ndarray = field(default_factory=lambda: np.zeros((0, 3)))
    truth: pd.Series = field(default_factory=lambda: pd.Series(dtype="object"))
    control_truth: pd.Series = field(default_factory=lambda: pd.Series(dtype="object"))
    #: Other maps the tuner considered, each fully scored. A report that shows only the winner hides
    #: what the choice was between.
    alternatives: list = field(default_factory=list)
    stopped_because: str = ""

    @property
    def ok(self) -> bool:
        return not self.stopped_because

    def tables(self) -> dict:
        """The result as saveable tables, for `results.save_bundle`."""
        return {"recovery": self.recovery, "per_cluster": self.per_cluster,
                "inference": self.inference, "control": self.control}


def build_map(nodes: pd.DataFrame, recipe: Recipe, closure: Closure, settings: dict,
              log=print) -> RecipeResult:
    """Build one map at the given UMAP settings, cluster it, score it, and name genes.

    Split out of `run` so the winning configuration and the alternatives kept for the report go
    through IDENTICAL code. A report whose alternatives were scored by a second, shorter path would
    be comparing the winner against something else.
    """
    spec = EmbeddingSpec(blocks=closure.blocks, random_state=recipe.seed,
                         **{**recipe.spec, **settings})
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
    settings = {**settings, **params}

    genes = np.zeros(len(nodes), dtype=bool)
    genes[np.arange(len(nodes))[kept] if kept is not None else np.arange(len(nodes))] = True
    result = RecipeResult(recipe=recipe, closure=closure, quality=quality, settings=settings,
                          genes=genes, labels=labels, coords=coords)
    if not quality["usable"]:
        # A recipe on a degenerate map returns that and stops. There is no score to report and a
        # number here would be read as one.
        result.stopped_because = f"the map cannot support a claim: {quality['why_not']}"
        log(result.stopped_because)
        return result

    truth = label_series(sub, closure.holdout_column, recipe.holdout_bins)
    result.truth = truth
    result.summary, result.recovery = score_recovery(labels, truth, min_label=MIN_LABEL)
    result.per_cluster = recovery_by_cluster(labels, truth)
    result.inference = infer(labels, truth, np.asarray(sub["gene_id"]),
                             min_enrichment=recipe.min_enrichment,
                             min_precision=recipe.min_precision)
    if not len(result.inference):
        # An empty answer must say how close it came. "No genes" and "no cluster was 80% pure, the
        # purest was 75%" are different findings, and only the second tells a reader whether the
        # question failed or the threshold did.
        purity = dominant_by_cluster(labels, truth)
        scoreable = purity[purity.n_labelled >= MIN_LABEL]
        result.inference_note = (
            f"no cluster reached {recipe.min_enrichment:.1f}x enrichment at "
            f"{recipe.min_precision:.2f} purity; the best scoreable cluster was "
            f"{scoreable.enrichment.max():.1f}x at "
            f"{scoreable.loc[scoreable.enrichment.idxmax(), 'share_in_cluster']:.2f} purity "
            f"({scoreable.loc[scoreable.enrichment.idxmax(), 'label']})"
            if len(scoreable) else
            f"no cluster holds {MIN_LABEL} labelled genes, so none can support a prediction")
        log(result.inference_note)
    if closure.control_column:
        control_truth = label_series(sub, closure.control_column, recipe.control_bins)
        result.control_truth = control_truth
        result.control = dominant_by_cluster(labels, control_truth).rename(
            columns={"label": "control_label", "n_labelled": "n_control_labelled"})
        result.control_summary, _per = score_recovery(labels, control_truth, min_label=MIN_LABEL)
        # The agreement goes on the inference's own rows. A control reported in a separate table is a
        # control nobody reads next to the thing it qualifies -- which is the entire point of it.
        if len(result.inference):
            cols = ["cluster", "control_label", "n_control_labelled", "enrichment", "agrees", "why"]
            result.inference = result.inference.merge(
                result.control[cols].rename(columns={"agrees": "control_agrees",
                                                     "why": "control_says",
                                                     "enrichment": "control_enrichment"}),
                on="cluster", how="left")
    log(f"{len(result.inference)} gene(s) named"
        + (f", {int(result.inference['control_agrees'].fillna(False).sum())} of them in clusters the "
           f"control corroborates" if "control_agrees" in result.inference else ""))
    return result


def build_model(nodes: pd.DataFrame, recipe: Recipe, closure: Closure, log=print) -> RecipeResult:
    """Answer a question with a classifier instead of a clustering, scored out of fold.

    The partition here is the out-of-fold predicted class, so `score_recovery` reports genuine
    cross-validated performance rather than memorisation, and `infer` names the unlabelled genes in
    predicted classes that are concentrated enough to support a claim. Everything downstream --
    including whether the independent control corroborates the answer -- is the code the clustering
    path uses, unchanged.
    """
    from . import methods
    from .embedding import build_matrix
    spec = EmbeddingSpec(blocks=closure.blocks, random_state=recipe.seed, **recipe.spec)
    if recipe.method == "boosted":
        # Raw, with the gaps intact. Every missing-value policy resolves absence before the model
        # sees it, and this method's entire argument is that it should not be resolved: a tree learns
        # a split for the missing branch, so "not measured" becomes evidence.
        X, names = methods.raw_matrix(nodes, closure.columns)
        kept = None
    else:
        X, names, kept = build_matrix(nodes, spec, log=lambda *a: None)
    sub = nodes.loc[kept] if kept is not None else nodes
    truth = label_series(sub, closure.holdout_column, recipe.holdout_bins)

    if recipe.method == "multiplex":
        n_graph = methods.graph_size()
        if n_graph is not None and len(nodes) != n_graph:
            return RecipeResult(recipe=recipe, closure=closure, stopped_because=(
                f"multiplex communities need the full node table the graph was built from: the "
                f"graph has {n_graph:,} nodes and this table has {len(nodes):,}."))
        if not closure.layers:
            return RecipeResult(recipe=recipe, closure=closure, stopped_because=(
                "leakage closure excluded every edge layer, so there is no graph left to cluster"))
        fit = methods.multiplex_communities(closure.layers, len(nodes), seed=recipe.seed)
        sub, truth = nodes, label_series(nodes, closure.holdout_column, recipe.holdout_bins)
    elif recipe.method.startswith("propagation"):
        _, _, layer = recipe.method.partition(":")
        if not layer:
            # A result rather than an exception, like every other refusal here: "this cannot be
            # asked that way" is a finding about the recipe, and a caller that must catch an
            # exception to learn it will eventually report it as a crash.
            return RecipeResult(recipe=recipe, closure=closure, stopped_because=(
                "name the layer to propagate over, as propagation:<layer>. The layers mean different "
                "things and two of them are attention-biased, so there is no sensible merged "
                "default. Available here: " + ", ".join(closure.layers)))
        if layer in closure.excluded_layers:
            # The edge guard, enforced where it matters rather than only reported. Diffusing a label
            # across the layer built out of that label recovers it perfectly and means nothing.
            return RecipeResult(recipe=recipe, closure=closure, stopped_because=(
                f"the {layer!r} layer is excluded for this holdout: "
                f"{closure.excluded_layers[layer]}"))
        if layer not in closure.layers:
            return RecipeResult(recipe=recipe, closure=closure, stopped_because=(
                f"no layer named {layer!r}; this catalogue has "
                + ", ".join(closure.layers)))
        # The graph's node ids are POSITIONS IN THE FULL NODE TABLE, so a walk over any other table
        # maps every edge to the wrong gene -- and silently, since the indices remain in range for a
        # large enough subset. Refused rather than guessed: there is no correct answer from a sample.
        n_graph = methods.graph_size()
        if n_graph is not None and len(nodes) != n_graph:
            return RecipeResult(recipe=recipe, closure=closure, stopped_because=(
                f"propagation needs the full node table the graph was built from: the graph has "
                f"{n_graph:,} nodes and this table has {len(nodes):,}. Edge ids are positions in "
                f"that table, so a subset would walk between the wrong genes."))
        positions = np.arange(len(nodes))[kept] if kept is not None else np.arange(len(nodes))
        fit = methods.propagation(layer, positions, truth, len(nodes), seed=recipe.seed)
    elif recipe.method == "boosted":
        fit = methods.boosted(np.asarray(X), truth, seed=recipe.seed)
    else:
        fit = methods.logistic(np.asarray(X), truth, seed=recipe.seed)
    labels = np.asarray(fit["partition"])
    genes = np.zeros(len(nodes), dtype=bool)
    genes[np.arange(len(nodes))[kept] if kept is not None else np.arange(len(nodes))] = True
    result = RecipeResult(recipe=recipe, closure=closure, settings=fit["settings"],
                          genes=genes, labels=labels, truth=truth,
                          quality=map_quality(labels))
    result.coefficients = (
        methods.importances(fit["model"], list(names), np.asarray(X), truth, recipe.seed)
        if recipe.method == "boosted"
        else methods.coefficients(fit["model"], fit["classes"], list(names)))
    # An unsupervised method has no classes to learn: its partition IS the grouping, exactly as the
    # clustering path's is, so the "nothing to train on" check does not apply to it.
    if not fit.get("unsupervised") and not len(fit["classes"]):
        result.stopped_because = ("too few labelled genes to train and score a classifier on this "
                                  "holdout")
        return result

    result.summary, result.recovery = score_recovery(labels, truth, min_label=MIN_LABEL)
    result.per_cluster = recovery_by_cluster(labels, truth)
    result.inference = infer(labels, truth, np.asarray(sub["gene_id"]),
                             min_enrichment=recipe.min_enrichment,
                             min_precision=recipe.min_precision)
    if len(result.inference) and len(fit["classes"]):
        # The partition's ids are class codes, so the class each one MEANS goes on the row. A table
        # saying "group 3" where the model said "dense granules" would be unreadable. An unsupervised
        # partition has no such mapping: its groups are groups, and `predicted` already names the
        # label they turned out to carry.
        result.inference["predicted_class"] = [fit["classes"][int(c)] if int(c) >= 0 else ""
                                               for c in result.inference.cluster]
    if closure.control_column:
        control_truth = label_series(sub, closure.control_column, recipe.control_bins)
        result.control_truth = control_truth
        result.control = dominant_by_cluster(labels, control_truth).rename(
            columns={"label": "control_label", "n_labelled": "n_control_labelled"})
        result.control_summary, _p = score_recovery(labels, control_truth, min_label=MIN_LABEL)
        if len(result.inference):
            cols = ["cluster", "control_label", "n_control_labelled", "enrichment", "agrees", "why"]
            result.inference = result.inference.merge(
                result.control[cols].rename(columns={"agrees": "control_agrees",
                                                     "why": "control_says",
                                                     "enrichment": "control_enrichment"}),
                on="cluster", how="left")
    log(f"{recipe.method}: {len(fit['classes'])} classes, "
        f"mean F1 {result.summary.get('mean_f1', float('nan')):.3f}, "
        f"{len(result.inference)} gene(s) named")
    return result


def run(nodes: pd.DataFrame, recipe: Recipe, tune: bool = True, sample: int = 1500,
        alternatives: int = 0, log=print) -> RecipeResult:
    """Answer one question: close leakage, tune a map, cluster it, score it, and name genes.

    The order matters and it is the argument of the whole module. Closure runs FIRST and reports
    before anything is built, so a refusal costs no compute and is legible. Tuning runs on what
    survived. The clustering is checked for degeneracy before it is scored, because a recovery number
    computed on a bisection reads exactly like a real one. Only then is the primary holdout scored,
    the unlabelled genes in concentrated clusters named, and the control put on the SAME clusters
    beside them.

    `alternatives` keeps the next best UMAP settings, each built and scored in full, for the report:
    a document that shows only the winner hides what the choice was between. They cost a full build
    each, so the default is none.
    """
    closure = close(nodes, recipe)
    log(closure.report())
    if not closure.ok:
        return RecipeResult(recipe=recipe, closure=closure, stopped_because=closure.refusal)
    if recipe.method != "umap+hdbscan":
        return build_model(nodes, recipe, closure, log=log)

    choices = [{}]
    if tune:
        spec = EmbeddingSpec(blocks=closure.blocks, random_state=recipe.seed, **recipe.spec)
        walk = tune_umap(nodes, spec, sample=sample, seed=recipe.seed, log=lambda *a: None)
        usable = walk[walk["usable"]] if len(walk) else walk
        if not len(usable):
            return RecipeResult(recipe=recipe, closure=closure,
                                stopped_because="no UMAP setting gave a map worth clustering")
        ranked = usable.sort_values(["stage", "score"], ascending=[False, False])
        # Distinct settings only: the walk scores each configuration on a sample and again in full,
        # so the same n_neighbors/min_dist pair appears twice and would otherwise be "the
        # alternative" to itself.
        seen, choices = set(), []
        for _, row in ranked.iterrows():
            key = (int(row["n_neighbors"]), float(row["min_dist"]))
            if key in seen:
                continue
            seen.add(key)
            choices.append({"n_neighbors": key[0], "min_dist": key[1]})
            if len(choices) > alternatives:
                break

    result = build_map(nodes, recipe, closure, choices[0], log=log)
    for other in choices[1:]:
        log(f"building alternative map {other}")
        result.alternatives.append(build_map(nodes, recipe, closure, other, log=lambda *a: None))
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


# --------------------------------------------------------------------------- relevance (43)
#: Genes named from one cluster beyond which more is not better. A claim about forty genes is a
#: result; a claim about four hundred is usually a cluster that swallowed a compartment.
GOOD_REACH = 40

#: Enrichment at which strength saturates. Ten-fold concentration is already a strong claim, and
#: letting a 300-fold rare-label cluster dominate the ranking would rank by rarity rather than by
#: interest.
STRONG_ENRICHMENT = 10.0


def relevance(inference: pd.DataFrame, nodes: pd.DataFrame) -> pd.DataFrame:
    """Rank named genes by how much they are worth following up. A HEURISTIC, and it says so.

    Instruction 43 asked for a biological-relevance score and warned that it must never override the
    raw numbers it summarises, so the four terms are added as COLUMNS and the ranking is one more
    column beside them. A reader who disagrees with the weighting -- and anyone who cares about one
    gene family will -- re-sorts on the term they care about instead of arguing with a number they
    cannot see inside. That is the same contract `interpret.interest` already makes for discovery
    findings, and `novelty` is imported from there rather than reimplemented.

    The four terms:

    * **strength** -- how concentrated the label is in the cluster, saturating at ten-fold, because
      beyond that the ranking would be sorting by how rare the label is;
    * **reach** -- how many genes the cluster names, saturating at forty, because a cluster naming
      four hundred has usually swallowed a compartment rather than found one;
    * **novelty** -- how unstudied those genes are. A gene can be unlocalised and still be one of the
      twenty everybody works on, and predicting its compartment is a smaller contribution than
      predicting one for a gene with no literature at all;
    * **corroborated** -- whether the INDEPENDENT control agrees. This is the term with the most
      claim to be there, and it is the one no other ranking in this project has: it is the difference
      between "the map found structure" and "the map found the structure I asked about".
    """
    from .interpret import novelty
    if inference is None or not len(inference):
        return pd.DataFrame()
    out = inference.copy()
    enrichment = pd.to_numeric(out.get("enrichment", 1.0), errors="coerce").fillna(1.0)
    out["strength"] = np.clip(np.log1p(enrichment) / np.log1p(STRONG_ENRICHMENT), 0.0, 1.0)
    per_cluster = out.groupby("cluster")["gene_id"].transform("size")
    out["reach"] = np.clip(np.log1p(per_cluster) / np.log1p(GOOD_REACH), 0.0, 1.0)
    out["novelty"] = [novelty(nodes, [g]) for g in out["gene_id"]]
    agrees = out.get("control_agrees", pd.Series(False, index=out.index)).fillna(False).astype(bool)
    # A halving rather than a zeroing: an uncorroborated claim is weaker, not void, and several
    # questions here have controls too sparse to corroborate anything at all -- scoring those to zero
    # would rank a question's answers by whether its control happened to be dense.
    out["corroborated"] = agrees
    out["relevance"] = out.strength * out.reach * out.novelty * np.where(agrees, 1.0, 0.5)
    return out.sort_values(["relevance", "strength"], ascending=False).reset_index(drop=True)
