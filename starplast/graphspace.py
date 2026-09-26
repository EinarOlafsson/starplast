"""One integrated neighbour space over the genes, and the nulls that decide whether it means anything.

The thirteen measured edge layers each answer a different question and are deliberately never merged
in this project (`build_graph`). That is right for reading a layer, and wrong for one question a user
keeps asking: *which genes are near this one, and on what evidence?* Answering it needs every
permitted layer and the measurement table in one graph, with a single number per pair and the
per-source breakdown that produced it. This module builds that graph.

It is a CLAIM, not a measurement. Nothing here writes to `ctx.graph`, replaces a measured layer or
adds an edge to one: a :class:`GraphSpace` is a separate object with its own pairs, and every table
it produces says which kind each pair is -- a `kind` of `measured` or `inferred`, and a `measured_in`
naming the layers that record the pair, empty for exactly the inferred ones. An integrated edge is a
hypothesis with a probability attached; a crosslink is an experiment. Anything that loses that
distinction turns a model's output into a dataset, which is the one failure this project is
organised to prevent.

**The evaluation was written before the model**, because the easy number here is worthless. Strategy
16 (`link_prediction`) reports AUROC 0.97 for held-out crosslinks against random non-pairs, and most
of that is the discovery that well-connected genes are well-connected: a random pair of genes is
usually a pair of obscure genes, and *any* feature correlated with degree separates it from an edge
between two hubs. So every number here is reported against three negative sets at once:

    random              non-pairs drawn uniformly from every gene in the table;
    degree_matched      each positive paired with a non-pair of similar evidence degree at BOTH ends;
    configuration       non-pairs drawn from a configuration model -- the same degree sequence,
                        stub-paired at random -- so the negatives are degree-plausible by design.

The gap between the first and the second is reported as its own quantity, `fame_gap`: it is how much
of the performance is fame rather than biology. The shipped model's gap is SMALL -- +0.008 AUROC on
the Toxoplasma table -- and that is the result, not a disappointment: this model is given no degree
and is trained against degree-matched non-pairs, so there was little fame in it to remove. The same
features plus degree, trained the usual way against random non-pairs, have a gap of +0.024 here and
of +0.216 on a preferential-attachment graph (`tests/test_graphspace.py`), which is what the number
looks like when fame IS the model. That diagnostic is fitted on every run, beside the real one, and
never ranks anything. `docs/graphspace.md` carries the full measured table.

Edges are held out by NODE GROUP (`Context.groups`, which is the orthogroup), never by random edge.
A randomly held-out edge can be recovered through a paralog that stayed visible, which scores the
resemblance of a gene to its own copy rather than the inference being tested.

Self-exclusion is the other rule that cannot be relaxed: when the evidence is scored on layer `l`,
`l` is not one of its features. It is enforced in one place -- :meth:`Evidence.features` with
`skip=l`, which drops `l`'s two columns and rebuilds the shared-partner count on the union of the
other layers -- and the integrated weight of a source is the mean of what it was worth in every
fold, counting as zero the fold in which it was the target. A source that can only predict itself
therefore ends at weight zero as arithmetic rather than as good intentions.

The leakage guard is the project's own: `Context.banned` for columns and `Context.banned_layers` for
layers, reached through `Context.measurement_layers` and `Context.features`, so holding out
`compartment` removes the compartment columns from the measurement similarity AND the compartment
layer from the sources, with no second copy of the closure living here.

Typical use::

    from starplast import graphspace, strategies
    ctx = strategies.Context.shipped("Tg")
    sp = graphspace.space(ctx)
    sp.neighbours("TGME49_233460", 10)      # who is near this gene, and on what evidence
    sp.evidence("TGME49_233460", "TGME49_269010")
    sp.gaps(50)                             # what the evidence implies and no layer records
    print(sp.report.summary())              # the three nulls, side by side
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import strategies as S

#: Layers that are never held out as a prediction target. `orthogroup` restates the node groups the
#: split is made from, so predicting it is predicting the split; `domain` and `compartment` are
#: annotations projected onto pairs rather than an experiment that could have missed an edge.
TARGET_SKIP = ("orthogroup", "domain", "compartment")
#: A layer with fewer edges than this cannot be split into a part to learn from and a part to score.
MIN_TARGET_EDGES = 200
#: Fewer held-out edges than this and the fold's AUROC is noise; the fold is reported, not scored.
MIN_TEST_EDGES = 20
#: Neighbours kept per gene per layer when candidate pairs are collected, strongest weight first.
PER_LAYER = 25
#: Nearest neighbours per gene in measurement space, so a pair with no edge anywhere can still be a
#: candidate. This is what lets the space speak about genes the networks have never touched.
KNN = 15
#: The most candidate pairs a space holds. 8,140 genes make 33 million pairs; scoring all of them
#: would cost an hour and 3 GB to say "no support" about 99% of them.
MAX_PAIRS = 400_000
#: Pairs per chunk when features are built. 20,000 x 367 measurement columns is 60 MB of float64,
#: which is the point: the un-chunked version of this loop needed 10 GB on the shipped table.
CHUNK = 20_000
#: Dimensions of the spectral embedding the learned model uses.
EMBED_DIM = 24
#: Degree bins the degree-matched null draws from. Twelve quantile bins put a hub with 200 partners
#: in a different bin from a gene with 3, and still leave enough genes per bin to sample from.
DEGREE_BINS = 12
#: Share of the node GROUPS held out for scoring.
TEST_FRACTION = 0.25
#: The negative sets every number is reported against, hardest last but all three always shown.
NULLS = ("random", "degree_matched", "configuration")
#: The models. `logistic` is the baseline on the pair features and is always fitted and always
#: shipped. `embedding` adds the spectral node vectors and is offered only where it earns its place on
#: the degree-matched null. `fame` is the USUAL RECIPE reproduced deliberately: the same features plus
#: each gene's degree, trained against randomly drawn non-pairs the way a published link predictor is.
#: It exists to be measured and never to rank anything. Its two AUROCs are far apart and the
#: baseline's are close, and that contrast is the point: the fame in a link-prediction number comes
#: from how the model was trained and what it was given, not from the data being dishonest.
MODELS = ("logistic", "embedding", "fame")


def _group(source: str) -> str:
    """The layer or source a feature column belongs to: `xlms` and `xlms_absent` are both `xlms`."""
    return source[:-len("_absent")] if source.endswith("_absent") else source


def _pack(a, b, n: int) -> np.ndarray:
    """Unordered pairs as single integers, for set membership. `a < b` is enforced by the caller."""
    return np.asarray(a, dtype=np.int64) * int(n) + np.asarray(b, dtype=np.int64)


def _ordered(pairs) -> np.ndarray:
    """Pairs with the smaller position first and self-pairs dropped: one row per unordered pair."""
    p = np.asarray(pairs, dtype=int).reshape(-1, 2)
    if not len(p):
        return p.reshape(0, 2)
    lo = np.minimum(p[:, 0], p[:, 1])
    hi = np.maximum(p[:, 0], p[:, 1])
    keep = lo != hi
    return np.column_stack([lo[keep], hi[keep]])


def _unique_pairs(pairs, n: int) -> np.ndarray:
    """`_ordered` with duplicates removed, in a stable sorted order."""
    p = _ordered(pairs)
    if not len(p):
        return p
    keys, first = np.unique(_pack(p[:, 0], p[:, 1], n), return_index=True)
    return p[np.sort(first)]


# --------------------------------------------------------------------------- the evidence
class Evidence:
    """Every permitted source of evidence about a pair of genes, as one matrix with named columns.

    One object per (context, held-out label), because every fold of the evaluation, every candidate
    pair and every question asked of a finished space must read the same features from the same
    closure. Two implementations of "the features of a pair" is how a reported AUROC stops describing
    the model that is shipped.

    The sources are: for each permitted layer, the pair's weight as a percentile among that layer's
    own edge weights, and an indicator that the layer does not record the pair at all; the cosine
    similarity of the two genes' permitted measurements; the Adamic-Adar shared-partner count on the
    union of the permitted layers; and, where the graph carries no such layer and the closure permits
    the column, a shared-orthogroup or shared-Pfam indicator.

    Two choices are worth the words. Weights are replaced by percentiles because a crosslink count of
    5 and a co-expression correlation of 0.95 are not on one scale, and a single weight per source is
    meaningless until they are. Similarity is a cosine on `Context.matrix`, whose columns are already
    rank-scaled to [-0.5, 0.5] with missing values at the median, so it is a rank correlation across
    the measured columns in all but name -- and it is read through `Context.features`, which applies
    the closure, rather than from `ctx.nodes`, which would not.

    Degree is deliberately NOT a source. It is the strongest predictor of a measured edge and the
    worst one: a model given degree learns which genes are famous. The degree-matched null exists to
    measure that contamination, and putting degree in the features would be paying for the null and
    then defeating it.
    """

    def __init__(self, ctx: S.Context, exclude=None, *, chunk: int = CHUNK):
        """Collect the permitted layers and measurement columns for `ctx`, holding out `exclude`."""
        self.ctx = ctx
        self.exclude = exclude
        self.n = int(ctx.n)
        self.chunk = int(chunk)
        self.banned = set(ctx.banned(exclude)) if exclude else set()
        self.banned_layers = set(ctx.banned_layers(exclude)) if exclude else set()
        # measurement_layers already drops the literature layers (co-mention follows attention), the
        # derived ones (unwritten interactions are crosslinks minus co-mentions) and the banned ones.
        self.layers = tuple(ctx.measurement_layers(exclude))
        self.X, self.columns = ctx.features(exclude)
        self.norms = np.linalg.norm(self.X, axis=1) if self.X.shape[1] else np.zeros(self.n)
        self.extras = tuple(self._extra_indicators())
        self._cache: dict = {}

    def _extra_indicators(self) -> list:
        """Shared-orthogroup and shared-domain indicators, where no layer already carries them.

        The shipped graph has an `orthogroup` and a `domain` layer, so on the real tables this adds
        nothing; a user's own table may have the columns and no layers, and then a paralog pair is
        evidence that would otherwise be thrown away. Both are subject to the closure: a label whose
        exclusion removes `orthogroup` removes the indicator with it.
        """
        from .embedding import as_text
        out = []
        ctx = self.ctx
        if "orthogroup" not in self.layers and "orthogroup" in ctx.nodes \
                and "orthogroup" not in self.banned:
            out.append(("same_orthogroup", ctx.groups()))
        if "domain" not in self.layers and "pfam_id" in ctx.nodes and "pfam_id" not in self.banned:
            pfam = as_text(ctx.nodes["pfam_id"]).to_numpy()
            # A gene with no domain must not match every other gene with no domain, so each blank
            # becomes its own value -- the same device `Context.groups` uses for a missing orthogroup.
            own = np.array([f"__{i}" for i in range(self.n)])
            out.append(("same_pfam", np.where(pfam != "", pfam, own)))
        return out

    # ---------------------------------------------------------------- the layers
    def weighted(self, layer: str):
        """One layer as a symmetric sparse matrix of percentile weights in (0, 1]."""
        key = ("pct", layer)
        if key not in self._cache:
            from scipy.stats import rankdata
            import scipy.sparse as sp
            A = self.ctx.adjacency(layer, "w").tocoo()
            data = rankdata(A.data) / max(len(A.data), 1) if len(A.data) else A.data
            W = sp.coo_matrix((data, (A.row, A.col)), shape=(self.n, self.n)).tocsr()
            self._cache[key] = W
        return self._cache[key]

    def union(self, skip=None):
        """The permitted layers as one binary graph, without `skip`: what a shared partner means."""
        key = ("union", skip)
        if key not in self._cache:
            import scipy.sparse as sp
            layers = [l for l in self.layers if l != skip]
            U = sum((self.ctx.adjacency(l, "binary") for l in layers),
                    sp.csr_matrix((self.n, self.n)))
            U = (U > 0).astype(float).tocsr()
            U.setdiag(0)
            U.eliminate_zeros()
            self._cache[key] = U
        return self._cache[key]

    def degree(self, skip=None) -> np.ndarray:
        """Each gene's degree in the union without `skip`: the fame the hard null removes."""
        key = ("deg", skip)
        if key not in self._cache:
            self._cache[key] = np.asarray(self.union(skip).sum(axis=1)).ravel()
        return self._cache[key]

    def _inverse_log_degree(self, skip=None) -> np.ndarray:
        """1/log(degree) per gene, the Adamic-Adar weight: a shared hub is weaker evidence."""
        key = ("aa", skip)
        if key not in self._cache:
            d = self.degree(skip)
            self._cache[key] = np.divide(1.0, np.log(np.maximum(d, 2.0)),
                                         out=np.zeros(self.n), where=d > 0)
        return self._cache[key]

    def records(self, a: int, b: int) -> list:
        """Which permitted layers record this pair. Empty means the pair is inferred, not measured."""
        return [l for l in self.layers if self.ctx.adjacency(l, "binary")[int(a), int(b)] > 0]

    def records_for(self, pairs) -> list:
        """The same for many pairs at once, as one comma-joined string per pair.

        Nine scalar lookups per row turns a 276,000-row table into minutes; nine vectorised lookups
        over the whole column turns it into a second. The single-pair version above stays for the one
        pair a user asks about.
        """
        pairs = np.asarray(pairs, dtype=int).reshape(-1, 2)
        if not len(pairs):
            return []
        hit = {l: np.asarray(self.ctx.adjacency(l, "binary")[pairs[:, 0], pairs[:, 1]]).ravel() > 0
               for l in self.layers}
        return [", ".join(l for l in self.layers if hit[l][i]) for i in range(len(pairs))]

    # ---------------------------------------------------------------- the features
    def names(self, skip=None) -> tuple:
        """The source names, in the order :meth:`features` returns them."""
        out = []
        for l in self.layers:
            if l == skip:
                continue
            out += [l, f"{l}_absent"]
        out += ["measurement_similarity", "shared_partners"]
        return tuple(out + [name for name, _ in self.extras])

    def measurement(self, skip=None) -> tuple:
        """(matrix, row norms, dropped columns): the measurements a similarity may use with `skip` out.

        A layer computed FROM measurements is those measurements' correlation structure, drawn as
        edges: co-translation is the top correlations among the ribosome-footprint columns, and
        co-expression the same among one transcript series. Hiding such a layer's edges while keeping
        its source columns in the cosine hands the model the layer back by a second route. The
        calibration sweep found exactly that -- held-out co-translation "recovered" at AUROC 0.94
        against degree-matched non-pairs -- so the columns a layer was built from
        (`search.LAYER_SOURCES`), and every other column of the experiments they came from, leave
        the similarity along with the layer's edges.
        """
        key = ("measurement", skip)
        if key not in self._cache:
            from .search import layer_measurement_sources
            drop = layer_measurement_sources(skip, self.columns, self.ctx.organism)
            if drop:
                keep = [i for i, c in enumerate(self.columns) if c not in drop]
                X = self.X[:, keep]
                norms = np.linalg.norm(X, axis=1) if X.shape[1] else np.zeros(self.n)
            else:
                X, norms = self.X, self.norms
            self._cache[key] = (X, norms, sorted(drop))
        return self._cache[key]

    def features(self, pairs, skip=None) -> np.ndarray:
        """One row of sources per pair, computed in chunks. `skip` is the self-exclusion rule.

        With `skip=l` the returned matrix has no column for `l`, its shared-partner count is rebuilt
        on the union of the other layers -- a shared partner reached through `l` is `l`'s own
        evidence arriving by a second route -- and the measurement similarity is computed without
        the columns `l` was built from (`measurement`), which is the same evidence arriving by a
        third. All three together are the rule.
        """
        pairs = np.asarray(pairs, dtype=int).reshape(-1, 2)
        names = self.names(skip)
        out = np.zeros((len(pairs), len(names)), dtype=float)
        if not len(pairs):
            return out
        U, inv = self.union(skip), self._inverse_log_degree(skip)
        X, norms, _dropped = self.measurement(skip)
        for start in range(0, len(pairs), self.chunk):
            self.ctx.check()
            sl = slice(start, min(start + self.chunk, len(pairs)))
            a, b = pairs[sl, 0], pairs[sl, 1]
            col = 0
            for l in self.layers:
                if l == skip:
                    continue
                w = np.asarray(self.weighted(l)[a, b]).ravel()
                out[sl, col] = w
                out[sl, col + 1] = (w <= 0).astype(float)
                col += 2
            if X.shape[1]:
                sim = (X[a] * X[b]).sum(axis=1) / (norms[a] * norms[b] + 1e-9)
                out[sl, col] = sim
            col += 1
            out[sl, col] = np.log1p(np.asarray(U[a].multiply(U[b]) @ inv).ravel())
            col += 1
            for _name, values in self.extras:
                out[sl, col] = (values[a] == values[b]).astype(float)
                col += 1
        return out

    def masked(self, pairs, skip: str) -> np.ndarray:
        """Full features with `skip`'s columns set to "this layer does not record the pair".

        Used to score a held-out edge on the INTEGRATED weights, which are defined over every source:
        dropping the columns would change the length of the vector, so instead the pair is presented
        as one its own target layer never saw. The shared-partner count is still rebuilt without
        `skip`, so the substance of the self-exclusion is unchanged; what the masking costs is a
        constant per fold (the weight of one absence indicator), which the calibration absorbs.
        """
        F = self.features(pairs, skip=None)
        names = list(self.names(None))
        if skip in self.layers:
            F[:, names.index(skip)] = 0.0
            F[:, names.index(f"{skip}_absent")] = 1.0
            # And the shared-partner column, rebuilt without the target layer.
            partial = self.features(pairs, skip=skip)
            F[:, names.index("shared_partners")] = partial[
                :, list(self.names(skip)).index("shared_partners")]
        return F

    def edges(self, layer: str) -> np.ndarray:
        """One layer's edges as unordered pairs."""
        import scipy.sparse as sp
        U = sp.triu(self.ctx.adjacency(layer, "binary"), 1).tocoo()
        return np.column_stack([U.row, U.col]).astype(int)

    def all_edges(self) -> np.ndarray:
        """Every pair any permitted layer records: what a negative may never be."""
        if "all_edges" not in self._cache:
            import scipy.sparse as sp
            U = sp.triu(self.union(None), 1).tocoo()
            self._cache["all_edges"] = np.column_stack([U.row, U.col]).astype(int)
        return self._cache["all_edges"]

    def targets(self) -> list:
        """The layers worth holding out: a measurement, not an annotation, with edges enough to
        score. An annotation layer projected onto pairs has no missing edges to find."""
        return [l for l in self.layers if l not in TARGET_SKIP
                and len(self.edges(l)) >= MIN_TARGET_EDGES]


# --------------------------------------------------------------------------- candidate pairs
def candidates(ev: Evidence, *, per_layer: int = PER_LAYER, knn: int = KNN,
               max_pairs: int = MAX_PAIRS) -> tuple:
    """(pairs, note): every gene's neighbours in every permitted layer, plus its nearest genes.

    What the cap costs, said plainly: a pair that is in no layer and in neither gene's `knn` nearest
    neighbours in measurement space is never scored, and so can never appear as a gap. The space is
    therefore not "every pair ranked" but "every pair with a reason to be looked at, ranked" -- on
    the shipped Toxoplasma table 276,384 of the 33 million pairs, which is the difference between
    half a minute and several hours. `per_layer` bounds the contribution of a dense layer (the
    compartment layer alone would otherwise supply 119,000 pairs, all of them the same claim).
    """
    ctx = ev.ctx
    parts, capped = [], []
    for l in ev.layers:
        ctx.check()
        W = ev.weighted(l)
        deg = np.asarray((W > 0).sum(axis=1)).ravel()
        if deg.max(initial=0) <= per_layer:
            parts.append(ev.edges(l))
            continue
        # Only the genes above the cap need the per-row sort, which is the difference between a loop
        # over 1,700 hubs and a loop over 8,140 genes.
        rows = [ev.edges(l)]
        for i in np.flatnonzero(deg > per_layer):
            row = W.getrow(i)
            idx = row.indices[np.argsort(-row.data, kind="stable")[:per_layer]]
            rows.append(np.column_stack([np.full(len(idx), i), idx]))
        keep = np.flatnonzero(deg <= per_layer)
        edges = ev.edges(l)
        inside = np.isin(edges[:, 0], keep) & np.isin(edges[:, 1], keep)
        parts.append(np.vstack([edges[inside]] + rows[1:]))
        capped.append(l)
    if ev.X.shape[1] and knn:
        from sklearn.neighbors import NearestNeighbors
        ctx.say("nearest genes in measurement space")
        k = int(min(knn + 1, ev.n))
        nn = NearestNeighbors(n_neighbors=k, metric="cosine").fit(ev.X)
        idx = nn.kneighbors(ev.X, return_distance=False)
        rows = np.repeat(np.arange(ev.n), idx.shape[1])
        parts.append(np.column_stack([rows, idx.ravel()]))
    pairs = _unique_pairs(np.vstack([p for p in parts if len(p)]) if parts
                          else np.zeros((0, 2), int), ev.n)
    note = (f"{len(capped)} of {len(ev.layers)} layers capped at {per_layer} neighbours per gene "
            f"({', '.join(capped)})") if capped else ""
    if len(pairs) > max_pairs:
        # Truncation keeps the pairs several layers agree on. It is reported rather than silent,
        # because a truncated space answers "no candidates" about pairs it simply never looked at.
        support = np.asarray(ev.union(None)[pairs[:, 0], pairs[:, 1]]).ravel()
        sim = ev.features(pairs)[:, list(ev.names()).index("measurement_similarity")]
        order = np.argsort(-(support + sim), kind="stable")[:max_pairs]
        pairs = pairs[np.sort(order)]
        note = (note + "; " if note else "") + \
            f"truncated to the {max_pairs:,} best-supported of {len(support):,} candidate pairs"
    return pairs, note


# --------------------------------------------------------------------------- the three nulls
def random_negatives(rng, nodes, taken: set, count: int, n: int) -> np.ndarray:
    """Non-pairs drawn uniformly from `nodes`: the easy null, and the one that flatters a model.

    Kept because the gap between it and the degree-matched one is the quantity worth reporting. The
    evaluation passes every gene here, including the ones no layer touches, since that is the null a
    published AUROC is usually against.
    """
    nodes = np.asarray(nodes, dtype=int)
    out, tries = [], 0
    while sum(len(o) for o in out) < count and tries < 40:
        tries += 1
        draw = _ordered(np.column_stack([rng.choice(nodes, size=2 * count),
                                         rng.choice(nodes, size=2 * count)]))
        if not len(draw):
            continue
        keys = _pack(draw[:, 0], draw[:, 1], n)
        fresh = ~np.isin(keys, np.fromiter(taken, dtype=np.int64, count=len(taken))) \
            if taken else np.ones(len(draw), bool)
        draw = draw[fresh]
        if len(draw):
            _, first = np.unique(_pack(draw[:, 0], draw[:, 1], n), return_index=True)
            draw = draw[np.sort(first)]
            out.append(draw)
            taken = taken | set(_pack(draw[:, 0], draw[:, 1], n).tolist())
    got = np.vstack(out) if out else np.zeros((0, 2), int)
    return got[:count]


def degree_matched_negatives(rng, degree: np.ndarray, positives: np.ndarray, taken: set,
                             nodes, n: int, bins: int = DEGREE_BINS) -> np.ndarray:
    """One non-pair per positive, with a gene of similar degree at BOTH ends.

    The point of the whole module. Genes are put into quantile bins of log degree in the evidence
    union; for a positive (a, b) a partner is drawn from a's bin and one from b's bin until the draw
    is not a real edge. What survives this null is what the model knows beyond which genes are
    well connected -- and on real data that is a much smaller number than the usual one.
    """
    positives = np.asarray(positives, dtype=int).reshape(-1, 2)
    nodes = np.asarray(nodes, dtype=int)
    if not len(positives) or not len(nodes):
        return np.zeros((0, 2), int)
    d = np.log1p(degree[nodes])
    edges = np.unique(np.quantile(d, np.linspace(0, 1, int(bins) + 1)))
    which = np.clip(np.searchsorted(edges, d, side="right") - 1, 0, max(len(edges) - 2, 0))
    members = [nodes[which == k] for k in range(max(len(edges) - 1, 1))]
    members = [m if len(m) else nodes for m in members]
    bin_of = np.zeros(n, dtype=int)
    bin_of[nodes] = which
    chosen, blocked = [], set(taken)
    need = np.arange(len(positives))
    for _round in range(30):
        if not len(need):
            break
        a = np.array([rng.choice(members[bin_of[positives[i, 0]]]) for i in need])
        b = np.array([rng.choice(members[bin_of[positives[i, 1]]]) for i in need])
        draw = np.column_stack([np.minimum(a, b), np.maximum(a, b)])
        keys = _pack(draw[:, 0], draw[:, 1], n)
        ok = (draw[:, 0] != draw[:, 1]) & ~np.isin(
            keys, np.fromiter(blocked, dtype=np.int64, count=len(blocked))) if blocked else \
            (draw[:, 0] != draw[:, 1])
        for row in np.flatnonzero(ok):
            key = int(keys[row])
            if key in blocked:
                continue
            blocked.add(key)
            chosen.append(draw[row])
        need = need[~ok]
    return np.array(chosen, dtype=int).reshape(-1, 2)


def configuration_negatives(rng, positives: np.ndarray, taken: set, n: int) -> np.ndarray:
    """Non-pairs from a configuration model on the positives: the same degrees, rewired at random.

    Every endpoint of every positive becomes a stub; the stubs are shuffled and paired. The result
    has (up to the self-loops and duplicates it loses) the positives' own degree sequence and none of
    their topology, so a model that scores these as highly as the real edges is reading the degree
    sequence. Degree-preserving rewiring and stub pairing are the same null; stub pairing is the one
    that is a single shuffle rather than a chain of swaps whose mixing has to be argued for.
    """
    positives = np.asarray(positives, dtype=int).reshape(-1, 2)
    if not len(positives):
        return np.zeros((0, 2), int)
    stubs = np.r_[positives[:, 0], positives[:, 1]]
    rng.shuffle(stubs)
    draw = _unique_pairs(np.column_stack([stubs[0::2], stubs[1::2]]), n)
    if not len(draw):
        return draw
    keys = _pack(draw[:, 0], draw[:, 1], n)
    if taken:
        draw = draw[~np.isin(keys, np.fromiter(taken, dtype=np.int64, count=len(taken)))]
    return draw


# --------------------------------------------------------------------------- metrics
def reliability(prob, positive, bins: int = 5) -> tuple:
    """(mean absolute gap, table): are the predicted probabilities right, decile by decile?

    A ranking can be perfect and the probabilities still wrong, and a probability is what a gap table
    invites a reader to act on. The gap is the mean over bins of |predicted - observed|.
    """
    p = np.asarray(prob, dtype=float)
    y = np.asarray(positive, dtype=bool)
    ok = np.isfinite(p)
    p, y = p[ok], y[ok]
    if not len(p):
        return float("nan"), pd.DataFrame()
    edges = np.unique(np.quantile(p, np.linspace(0, 1, int(bins) + 1)))
    which = np.clip(np.searchsorted(edges, p, side="right") - 1, 0, max(len(edges) - 2, 0))
    rows = []
    for k in range(max(len(edges) - 1, 1)):
        m = which == k
        if not m.any():
            continue
        rows.append({"bin": k, "pairs": int(m.sum()), "predicted": float(p[m].mean()),
                     "observed": float(y[m].mean())})
    table = pd.DataFrame(rows)
    gap = float(np.average(np.abs(table["predicted"] - table["observed"]),
                           weights=table["pairs"])) if len(table) else float("nan")
    return gap, table


def metrics(prob, positive, k: int | None = None) -> dict:
    """AUROC, precision@k, the Brier score and the reliability gap for one scored pair set.

    `k` defaults to the number of positives, capped at 100: precision@k with k the length of the list
    is just the class balance, which is a property of the sampling and not of the model.
    """
    p = np.asarray(prob, dtype=float)
    y = np.asarray(positive, dtype=bool)
    n_pos = int(y.sum())
    if not n_pos or n_pos == len(y):
        return {"auroc": float("nan"), "precision_at_k": float("nan"), "brier": float("nan"),
                "reliability_gap": float("nan"), "positives": n_pos, "negatives": len(y) - n_pos}
    kk = int(min(max(k if k else min(n_pos, 100), 1), len(y)))
    top = np.argsort(-p, kind="stable")[:kk]
    return {"auroc": S.auroc(p, y), "precision_at_k": float(y[top].mean()),
            "brier": float(np.mean((np.clip(p, 0, 1) - y.astype(float)) ** 2)),
            "reliability_gap": reliability(p, y)[0], "positives": n_pos,
            "negatives": int(len(y) - n_pos), "k": kk}


# --------------------------------------------------------------------------- the evaluation
@dataclass
class Fold:
    """One held-out layer: what was hidden, what each model scored, and what each source was worth."""
    layer: str
    sources: tuple
    coefficients: dict
    intercept: float
    n_train: int
    n_test: int
    scores: dict = field(default_factory=dict)
    note: str = ""
    #: (standardised masked features, truth) for the held-out edges and their degree-matched
    #: non-pairs, kept so the integrated model can be CALIBRATED on pairs no fold's own layer could
    #: see. Calibration cannot be done inside the fold, because the integrated weights are not known
    #: until every fold has run.
    calibration: tuple | None = None

    def auroc(self, model: str, null: str) -> float:
        """The AUROC of this fold's held-out edges under one model and one null."""
        return float(self.scores.get((model, null), {}).get("auroc", float("nan")))


@dataclass
class EvaluationReport:
    """Every model against every null, per held-out layer and pooled: the honest comparison table.

    `fame_gap` is the quantity the module exists to report: the baseline's AUROC against random
    non-pairs minus its AUROC against degree-matched ones. It is not an error term. It is the share
    of an impressive-looking link-prediction number that is the observation that hubs are hubs.
    """
    folds: list
    models: tuple = MODELS
    nulls: tuple = NULLS
    settings: dict = field(default_factory=dict)
    seconds: float = 0.0

    @property
    def scored(self) -> list:
        """The folds with enough held-out edges to mean something."""
        return [f for f in self.folds if f.n_test >= MIN_TEST_EDGES]

    def auroc(self, model: str = "logistic", null: str = "degree_matched") -> float:
        """One pooled number: the mean over scored folds, weighted by held-out edges."""
        vals = [(f.auroc(model, null), f.n_test) for f in self.scored]
        vals = [(v, w) for v, w in vals if np.isfinite(v)]
        if not vals:
            return float("nan")
        return float(np.average([v for v, _ in vals], weights=[w for _, w in vals]))

    def metric(self, name: str, model: str = "logistic", null: str = "degree_matched") -> float:
        """Any of the metrics, pooled the same way."""
        vals = [(f.scores.get((model, null), {}).get(name, float("nan")), f.n_test)
                for f in self.scored]
        vals = [(v, w) for v, w in vals if np.isfinite(v)]
        if not vals:
            return float("nan")
        return float(np.average([v for v, _ in vals], weights=[w for _, w in vals]))

    @property
    def fame_gap(self) -> float:
        """Random-null AUROC minus degree-matched AUROC: how much of the skill is fame."""
        return self.auroc("logistic", "random") - self.auroc("logistic", "degree_matched")

    @property
    def learned_gain(self) -> float:
        """What the learned embedding adds on the degree-matched null. Negative means it costs."""
        return self.auroc("embedding", "degree_matched") - self.auroc("logistic", "degree_matched")

    def learned_gain_per_layer(self) -> list:
        """The embedding's gain on the degree-matched null, one number per held-out layer.

        The pooled gain hides the thing worth knowing. On the shipped table the embedding gains
        +0.022 on crosslinks and +0.023 on structure, +0.002 on co-expression, and LOSES 0.005 on
        co-translation: the profile of a model that finds something in the sparse layers and nothing
        in the dense ones. A user whose table is mostly crosslinks should test it themselves.
        """
        out = []
        for f in self.scored:
            a = f.auroc("embedding", "degree_matched") - f.auroc("logistic", "degree_matched")
            if np.isfinite(a):
                out.append(a)
        return out

    @property
    def fame_recipe_gap(self) -> float:
        """The same gap for the diagnostic `fame` model: the usual recipe's own distance from honesty.

        This is the number to quote when someone shows an AUROC of 0.97 for link prediction. It is
        the same data and the same held-out edges; all that changes is that the model was given degree
        and trained against random non-pairs, which is how such figures are made.
        """
        return self.auroc("fame", "random") - self.auroc("fame", "degree_matched")

    @property
    def learned_earns_its_place(self) -> bool:
        """Whether the embedding may be offered for ranking at all.

        Two conditions, both with a reason. The gain on the DEGREE-MATCHED null must exceed 0.01
        AUROC, because a smaller gain is not worth a model that cannot say why it ranked a pair where
        it did -- the attribution a gap table needs is exactly what an embedding does not have. And
        it must exceed twice its own standard error across the held-out layers, because a mean gain
        that is inside its own spread is a coin landing heads. Measured on the shipped table: +0.0075
        pooled, per-layer gains of +0.002, +0.010, -0.005, +0.023, +0.022, a standard error of 0.005
        and so a bar of 0.011 -- it fails both conditions and the baseline ships. A tie is a loss:
        the interpretable model keeps the place by default.
        """
        gains = self.learned_gain_per_layer()
        if len(gains) < 2 or not np.isfinite(self.learned_gain):
            return False
        se = float(np.std(gains, ddof=1) / math.sqrt(len(gains)))
        return bool(self.learned_gain > 0.01 and self.learned_gain > 2 * se)

    def frame(self) -> pd.DataFrame:
        """One row per model and null: the table `docs/graphspace.md` prints."""
        rows = []
        for model in self.models:
            for null in self.nulls:
                if not any((model, null) in f.scores for f in self.scored):
                    continue
                rows.append({"model": model, "null": null,
                             "auroc": self.auroc(model, null),
                             "precision_at_k": self.metric("precision_at_k", model, null),
                             "brier": self.metric("brier", model, null),
                             "reliability_gap": self.metric("reliability_gap", model, null),
                             "held_out_edges": int(sum(f.n_test for f in self.scored))})
        return pd.DataFrame(rows)

    def per_layer(self) -> pd.DataFrame:
        """The same numbers per held-out layer, because one layer can carry a pooled average."""
        rows = []
        for f in self.folds:
            row = {"layer": f.layer, "train_edges": f.n_train, "held_out_edges": f.n_test,
                   "scored": f.n_test >= MIN_TEST_EDGES, "note": f.note}
            for model in self.models:
                for null in self.nulls:
                    if (model, null) in f.scores:
                        row[f"{model}_{null}"] = f.auroc(model, null)
            rows.append(row)
        return pd.DataFrame(rows)

    def weights(self) -> pd.DataFrame:
        """Each source's integrated weight, and the folds it was allowed to speak in."""
        names = sorted({s for f in self.folds for s in f.sources})
        rows = []
        for name in names:
            vals = [f.coefficients.get(name, 0.0) for f in self.scored]
            spoke = [f.layer for f in self.scored if name in f.sources]
            rows.append({"source": name, "weight": float(np.mean(vals)) if vals else float("nan"),
                         "folds_used": len(spoke), "folds_excluded": len(self.scored) - len(spoke)})
        return pd.DataFrame(rows).sort_values("weight", key=abs, ascending=False,
                                              kind="stable").reset_index(drop=True)

    def summary(self) -> str:
        """The three nulls in one line each, and the verdict on the learned model."""
        out = []
        for model in self.models:
            for null in self.nulls:
                a = self.auroc(model, null)
                if np.isfinite(a):
                    out.append(f"{model} vs {null}: AUROC {a:.3f}, "
                               f"precision@k {self.metric('precision_at_k', model, null):.3f}, "
                               f"Brier {self.metric('brier', model, null):.3f}")
        out.append(f"fame gap (random minus degree-matched): baseline {self.fame_gap:+.3f}, "
                   f"the usual recipe {self.fame_recipe_gap:+.3f}")
        gain = self.learned_gain
        verdict = "earns" if self.learned_earns_its_place else "does not earn"
        out.append(f"the learned embedding {verdict} its place: {gain:+.3f} AUROC on the "
                   f"degree-matched null")
        return "; ".join(out) if out else "no fold had enough held-out edges to score"


def group_holdout(ctx: S.Context, fraction: float = TEST_FRACTION, salt: int = 0) -> np.ndarray:
    """A boolean per gene: whether its NODE GROUP is held out. Paralogs are held out together."""
    groups = ctx.groups()
    rng = ctx.rng(900 + int(salt))
    uniq = rng.permutation(np.unique(groups))
    sizes = pd.Series(groups).value_counts()
    want, total, chosen = fraction * ctx.n, 0, set()
    for g in uniq:
        if total >= want:
            break
        chosen.add(g)
        total += int(sizes[g])
    return np.isin(groups, list(chosen))


def _spectral(U, dim: int, seed: int) -> np.ndarray:
    """A truncated-SVD embedding of one graph: the learned model's node vectors.

    A spectral factorisation rather than random walks, because node2vec's dependency (gensim) is not
    installed and this project does not add a dependency to try a model. The adjacency is
    symmetrically degree-normalised first, so a hub does not dominate every singular vector -- which
    is the same fame problem in a different coat.
    """
    import scipy.sparse as sp
    from scipy.sparse.linalg import svds
    d = np.asarray(U.sum(axis=1)).ravel()
    inv = np.divide(1.0, np.sqrt(d), out=np.zeros_like(d), where=d > 0)
    A = sp.diags(inv) @ U @ sp.diags(inv)
    k = int(min(dim, min(A.shape) - 2))
    if k < 2:
        return np.zeros((U.shape[0], 0))
    rng = np.random.default_rng(seed)
    v0 = rng.normal(size=A.shape[0])
    try:
        left, s, _ = svds(A.asfptype(), k=k, v0=v0)
    except Exception:                       # a graph too small or too degenerate to factorise
        return np.zeros((U.shape[0], 0))
    return left * np.sqrt(np.maximum(s, 0.0))


def _hadamard(emb, train_pairs):
    """The pair features of an embedding: the elementwise product of the two genes' vectors.

    Standardised on the training pairs, and that is not a formality. A normalised adjacency's singular
    vectors have entries of order 1/sqrt(n), so their products are of order 1e-3; handed to a logistic
    regression beside standardised features they receive a negligible coefficient and the "learned
    model" silently becomes the baseline, reporting the baseline's AUROC to six decimal places. That
    is what this function existing separately is for. Returns ``None`` when there is no embedding.
    """
    if emb is None or not emb.shape[1]:
        return None
    train_pairs = np.asarray(train_pairs, dtype=int).reshape(-1, 2)
    H = emb[train_pairs[:, 0]] * emb[train_pairs[:, 1]]
    mean, scale = H.mean(axis=0), H.std(axis=0)
    scale[scale < 1e-12] = 1.0
    return lambda p: ((emb[np.asarray(p, dtype=int)[:, 0]] * emb[np.asarray(p, dtype=int)[:, 1]])
                      - mean) / scale


def _degree_features(degree: np.ndarray):
    """Each gene's degree, as the pair features of the model this module exists to disqualify.

    Fitted beside the baseline and reported beside it, never used to rank a pair. A link predictor
    given degree and trained against random non-pairs scores well by answering a different question --
    "are both of these genes well connected?" -- and the difference between its two nulls is the size
    of that substitution. Measuring it is the only way to claim the shipped model, which is given no
    degree and trained against degree-matched non-pairs, does not do the same thing quietly.
    """
    d = np.log1p(np.asarray(degree, dtype=float))

    def features(pairs):
        p = np.asarray(pairs, dtype=int).reshape(-1, 2)
        a, b = d[p[:, 0]], d[p[:, 1]]
        return np.column_stack([a + b, np.minimum(a, b)])
    return features


def _fit(X, y, seed: int = 0):
    """The logistic regression every model here is: balanced classes, mild regularisation."""
    from sklearn.linear_model import LogisticRegression
    return LogisticRegression(C=1.0, max_iter=2000, class_weight="balanced",
                              random_state=int(seed)).fit(X, y)


def evaluate(ctx: S.Context, ev: Evidence | None = None, *, exclude=None, models=MODELS,
             fraction: float = TEST_FRACTION, embed_dim: int = EMBED_DIM,
             scaler=None, salt: int = 0) -> EvaluationReport:
    """Hold out each layer in turn by node group, and score every model against all three nulls.

    One fold per layer worth holding out. In a fold the layer is removed from the features entirely
    (`Evidence.features(skip=layer)`), its edges between held-out groups become the positives, its
    edges between visible groups become the training positives, and the negatives are drawn three
    ways. Edges that straddle the split are used for neither: they have one endpoint the model has
    seen, which is neither a fair test nor a wasted edge worth arguing about.
    """
    t0 = time.monotonic()
    ev = ev or Evidence(ctx, exclude)
    held = group_holdout(ctx, fraction, salt)
    all_keys = set(_pack(ev.all_edges()[:, 0], ev.all_edges()[:, 1], ev.n).tolist())
    scaler = scaler if scaler is not None else _pool_scaler(ev)
    folds = []
    for layer in ev.targets():
        ctx.say(f"graph space: holding out the {layer} layer")
        edges = ev.edges(layer)
        test = edges[held[edges[:, 0]] & held[edges[:, 1]]]
        train = edges[~held[edges[:, 0]] & ~held[edges[:, 1]]]
        nodes = np.flatnonzero(ev.degree(layer) > 0)
        # The layer's POSITION, not hash(layer): Python randomises string hashes per process, so a
        # salt built from one made every fold's negatives different on every run, and two runs of the
        # same evaluation disagreed in the third decimal. The context's promise is repeatability.
        rng = ctx.rng(1300 + 7 * ev.layers.index(layer) + int(salt))
        if len(test) < MIN_TEST_EDGES or len(train) < MIN_TARGET_EDGES // 2 or len(nodes) < 20:
            folds.append(Fold(layer=layer, sources=ev.names(layer), coefficients={},
                              intercept=0.0, n_train=len(train), n_test=len(test),
                              note="too few edges between held-out groups to score"))
            continue
        train_neg = degree_matched_negatives(rng, ev.degree(layer), train, all_keys, nodes, ev.n)
        if len(train_neg) < 10:
            folds.append(Fold(layer=layer, sources=ev.names(layer), coefficients={},
                              intercept=0.0, n_train=len(train), n_test=len(test),
                              note="no degree-matched negatives could be drawn"))
            continue
        negatives = {
            # Uniformly over EVERY gene, not only the ones the layers touch: that is the null a link
            # prediction is usually scored against, and being the flattering one is its whole purpose
            # here. Restricting it to connected genes, which an earlier version of this did, is itself
            # a degree control, and it closed the fame gap to 0.003 before the gap had been measured.
            "random": random_negatives(rng, np.arange(ev.n), all_keys, len(test), ev.n),
            "degree_matched": degree_matched_negatives(rng, ev.degree(layer), test, all_keys,
                                                       nodes, ev.n),
            "configuration": configuration_negatives(rng, test, all_keys, ev.n)}
        names = ev.names(layer)
        keep = [list(ev.names(None)).index(s) for s in names]
        std = lambda F: (F - scaler[0][keep]) / scaler[1][keep]
        pairs_tr = np.vstack([train, train_neg])
        Ftr = std(ev.features(pairs_tr, skip=layer))
        ytr = np.r_[np.ones(len(train), bool), np.zeros(len(train_neg), bool)]
        fitted = {"logistic": _fit(Ftr, ytr, ctx.seed)}
        emb = _spectral(ev.union(layer), embed_dim, ctx.seed) if "embedding" in models else None
        hadamard = _hadamard(emb, pairs_tr)
        if hadamard is not None:
            fitted["embedding"] = _fit(np.column_stack([Ftr, hadamard(pairs_tr)]), ytr, ctx.seed)
        fame = _degree_features(ev.degree(layer)) if "fame" in models else None
        if fame is not None:
            # Trained against RANDOM non-pairs, not degree-matched ones, because that is the recipe
            # being reproduced. Training on the easy null and reporting on it is how a link predictor
            # reaches 0.97; the same model against degree-matched non-pairs is the honest number, and
            # the difference between the two is what a reader of such a figure is not shown.
            easy_neg = random_negatives(rng, np.arange(ev.n), all_keys, len(train), ev.n)
            if len(easy_neg) >= 10:
                fame_pairs = np.vstack([train, easy_neg])
                fitted["fame"] = _fit(
                    np.column_stack([std(ev.features(fame_pairs, skip=layer)),
                                     fame(fame_pairs)]),
                    np.r_[np.ones(len(train), bool), np.zeros(len(easy_neg), bool)], ctx.seed)
            else:
                fame = None
        scores = {}
        for null, neg in negatives.items():
            if len(neg) < MIN_TEST_EDGES // 2:
                continue
            pairs = np.vstack([test, neg])
            y = np.r_[np.ones(len(test), bool), np.zeros(len(neg), bool)]
            F = std(ev.features(pairs, skip=layer))
            for model, m in fitted.items():
                extra = {"logistic": lambda p: np.zeros((len(p), 0)),
                         "embedding": hadamard, "fame": fame}[model]
                design = np.column_stack([F, extra(pairs)])
                scores[(model, null)] = metrics(m.predict_proba(design)[:, 1], y)
        coefs = dict(zip(names, fitted["logistic"].coef_[0]))
        cal_pairs = np.vstack([test, negatives["degree_matched"]]) \
            if len(negatives["degree_matched"]) else test
        cal_y = np.r_[np.ones(len(test), bool),
                      np.zeros(len(cal_pairs) - len(test), bool)]
        cal = ((ev.masked(cal_pairs, layer) - scaler[0]) / scaler[1], cal_y)
        folds.append(Fold(layer=layer, sources=names, coefficients=coefs,
                          intercept=float(fitted["logistic"].intercept_[0]),
                          n_train=len(train), n_test=len(test), scores=scores,
                          calibration=cal))
    return EvaluationReport(folds=folds, models=tuple(models), seconds=time.monotonic() - t0,
                            settings={"fraction": fraction, "embed_dim": embed_dim,
                                      "exclude": exclude, "layers": ev.layers,
                                      "targets": tuple(ev.targets())})


def _pool_scaler(ev: Evidence, size: int = 20000) -> tuple:
    """(mean, scale) per source over a sample of candidate-like pairs, shared by every fold.

    One standardisation for every fold is what makes the folds' coefficients comparable, and
    comparable coefficients are what makes an average of them a weight rather than a number.
    """
    rng = ev.ctx.rng(77)
    edges = ev.all_edges()
    take = edges[rng.choice(len(edges), size=int(min(size, len(edges))), replace=False)] \
        if len(edges) else np.zeros((0, 2), int)
    rand = random_negatives(rng, np.arange(ev.n), set(), int(min(size, ev.n * 2)), ev.n)
    pairs = np.vstack([p for p in (take, rand) if len(p)]) if len(take) or len(rand) else \
        np.zeros((0, 2), int)
    F = ev.features(pairs)
    if not len(F):
        k = len(ev.names(None))
        return np.zeros(k), np.ones(k)
    mean = F.mean(axis=0)
    scale = F.std(axis=0)
    scale[scale < 1e-9] = 1.0
    return mean, scale


# --------------------------------------------------------------------------- the space itself
@dataclass
class GraphSpace:
    """One integrated graph over genes: a probability per pair, and the evidence that produced it.

    Not a layer. `pairs` are candidate pairs, `probability` is a claim about each, and `measured` says
    whether any permitted layer records it -- the distinction every table here keeps, because an
    inferred edge read as a measurement is a fabricated experiment. The measured layers in
    `ctx.graph` are untouched by everything in this class.

    The probability answers one question: *given the permitted evidence, is this pair a real
    relationship of the kind the measured layers record, rather than a degree-matched non-pair?* It is
    calibrated on held-out edges against degree-matched negatives at one negative per positive, so it
    is conditional on that comparison and not an absolute posterior over all 33 million pairs. The
    honest reading is a relative one, and `report.metric("reliability_gap")` says how far even that
    is from perfect.
    """
    ev: Evidence
    pairs: np.ndarray
    standardised: np.ndarray
    probability: np.ndarray
    weights: dict
    intercept: float
    platt: tuple
    report: EvaluationReport
    measured: np.ndarray
    settings: dict = field(default_factory=dict)
    note: str = ""
    seconds: float = 0.0
    cache: dict = field(default_factory=dict)

    # ---------------------------------------------------------------- identity
    @property
    def ctx(self) -> S.Context:
        """The context this space was built over."""
        return self.ev.ctx

    @property
    def sources(self) -> tuple:
        """The source names, in the order `standardised` holds them."""
        return self.ev.names(None)

    @property
    def layers(self) -> tuple:
        """The permitted layers. A banned layer is absent from this tuple, which is the point."""
        return self.ev.layers

    @property
    def columns(self) -> list:
        """The permitted measurement columns behind `measurement_similarity`."""
        return list(self.ev.columns)

    def position(self, gene) -> int:
        """A gene id or a position, as a position. A name this table does not have is an error."""
        if isinstance(gene, (int, np.integer)):
            return int(gene)
        i = self.ctx.index.get(str(gene).upper())
        if i is None:
            raise ValueError(f"{gene!r} is not a gene of this table")
        return int(i)

    @property
    def _row(self) -> dict:
        """packed pair -> row, built once. A dict rather than a search, since a user asks one pair at
        a time and a sorted search over 350,000 pairs is written wrong more often than it is read."""
        if "row" not in self.cache:
            keys = _pack(self.pairs[:, 0], self.pairs[:, 1], self.ev.n)
            self.cache["row"] = {int(k): i for i, k in enumerate(keys)}
        return self.cache["row"]

    def row_of(self, a, b) -> int:
        """The row holding this pair, or -1 when the pair was never a candidate."""
        i, j = self.position(a), self.position(b)
        if i == j:
            return -1
        lo, hi = min(i, j), max(i, j)
        return int(self._row.get(int(_pack(lo, hi, self.ev.n)), -1))

    # ---------------------------------------------------------------- scoring
    def _probability(self, Z) -> np.ndarray:
        """Calibrated probabilities from standardised features."""
        w = np.array([self.weights[s] for s in self.sources], dtype=float)
        A, B = self.platt
        return 1.0 / (1.0 + np.exp(-(A * (np.asarray(Z, dtype=float) @ w + self.intercept) + B)))

    def edge_strength(self, a, b) -> float:
        """The probability for one pair, computed on demand when it was never a candidate.

        A pair outside the candidate set is not refused: the same weights are applied to its own
        features. What the candidate set decides is which pairs the space RANKS, not which it can
        score, so `neighbours` and `gaps` see the candidates and this sees anything.
        """
        i, j = self.position(a), self.position(b)
        if i == j:
            return float("nan")
        row = self.row_of(i, j)
        if row >= 0:
            return float(self.probability[row])
        Z = (self.ev.features(np.array([[min(i, j), max(i, j)]])) - self._scale[0]) / self._scale[1]
        return float(self._probability(Z)[0])

    @property
    def _scale(self) -> tuple:
        """The (mean, scale) the features were standardised with."""
        return self.settings["scaler"]

    def evidence(self, a, b) -> pd.DataFrame:
        """Which sources spoke for this pair and how much each moved the answer.

        `contribution` is on the calibrated logit scale and sums, with `intercept`, to the logit of
        the probability -- so a reader can see that a pair is high because three layers record it, or
        because nothing records it and the measurements are nearly identical. `measured` marks the
        sources that are a layer recording the pair, as against a similarity or a shared partner.
        """
        i, j = self.position(a), self.position(b)
        row = self.row_of(i, j)
        if row >= 0:
            Z = self.standardised[row]
        else:
            Z = ((self.ev.features(np.array([[min(i, j), max(i, j)]])) - self._scale[0])
                 / self._scale[1])[0]
        A = self.platt[0]
        records = set(self.ev.records(min(i, j), max(i, j)))
        raw = (Z * self._scale[1]) + self._scale[0]
        rows = []
        for k, s in enumerate(self.sources):
            rows.append({"source": s, "group": _group(s), "value": float(raw[k]),
                         "standardised": float(Z[k]), "weight": float(self.weights[s]),
                         "contribution": float(A * self.weights[s] * Z[k]),
                         "measured_edge": _group(s) in records})
        table = pd.DataFrame(rows)
        table = table.sort_values("contribution", key=abs, ascending=False, kind="stable")
        return table.reset_index(drop=True)

    def neighbours(self, gene, k: int = 10) -> pd.DataFrame:
        """The genes nearest `gene` in the integrated space, strongest first.

        `measured_in` names the layers that record the pair and is empty for an inferred neighbour,
        and `top_evidence` names the three sources that moved the probability most. Only candidate
        pairs are ranked, so a neighbour that shares no layer and is outside either gene's nearest
        neighbours in measurement space will not appear -- see :func:`candidates`.
        """
        i = self.position(gene)
        mask = (self.pairs[:, 0] == i) | (self.pairs[:, 1] == i)
        rows = np.flatnonzero(mask)
        if not len(rows):
            return pd.DataFrame(columns=["gene_id", "product", "probability", "measured",
                                         "measured_in", "top_evidence"])
        rows = rows[np.argsort(-self.probability[rows], kind="stable")[:int(k)]]
        other = np.where(self.pairs[rows, 0] == i, self.pairs[rows, 1], self.pairs[rows, 0])
        return pd.DataFrame({
            "gene_id": self.ctx.gene_ids[other], "product": self.ctx.product(other),
            "probability": self.probability[rows], "measured": self.measured[rows],
            "measured_in": self.ev.records_for(
                np.column_stack([np.minimum(i, other), np.maximum(i, other)])),
            "top_evidence": [self._top_evidence(r) for r in rows]})

    def _top_evidence(self, row: int, k: int = 3) -> str:
        """The `k` sources that moved this pair's probability most, as text for a table cell.

        Totalled per layer rather than per column, because a layer contributes through two columns
        (its weight and its absence indicator) and on a layer whose edges all weigh the same the two
        are collinear -- the fit splits one effect between them, and printing both as if they were
        separate reasons says the same thing twice.
        """
        A = self.platt[0]
        contrib = A * np.array([self.weights[s] for s in self.sources]) * self.standardised[row]
        totals: dict = {}
        for s, c in zip(self.sources, contrib):
            totals[_group(s)] = totals.get(_group(s), 0.0) + float(c)
        order = sorted(totals, key=lambda g: -abs(totals[g]))[:k]
        return ", ".join(f"{g} {totals[g]:+.2f}" for g in order)

    # ---------------------------------------------------------------- tables
    def to_frame(self, top: int | None = None, measured: bool | None = None) -> pd.DataFrame:
        """Every candidate pair as a table, strongest first. `measured` filters to one kind."""
        keep = np.ones(len(self.pairs), bool) if measured is None else \
            (self.measured == bool(measured))
        rows = np.flatnonzero(keep)
        rows = rows[np.argsort(-self.probability[rows], kind="stable")]
        if top:
            rows = rows[:int(top)]
        a, b = self.pairs[rows, 0], self.pairs[rows, 1]
        return pd.DataFrame({
            "gene_a": self.ctx.gene_ids[a], "product_a": self.ctx.product(a),
            "gene_b": self.ctx.gene_ids[b], "product_b": self.ctx.product(b),
            "probability": self.probability[rows],
            "kind": np.where(self.measured[rows], "measured", "inferred"),
            "measured_in": self.ev.records_for(np.column_stack([a, b])),
            "top_evidence": [self._top_evidence(r) for r in rows]})

    def gaps(self, top: int = 200) -> pd.DataFrame:
        """The pairs the evidence implies and NO permitted layer records: the gaps, ranked.

        A gap's probability rests only on the sources that are not a layer edge -- measurement
        similarity, shared partners, a shared orthogroup or domain -- since every layer weight is
        zero for a pair no layer records. That is the right way round: a gap is an inference from
        everything except a direct observation, and the table says so in `top_evidence`.
        """
        return self.to_frame(top=top, measured=False)

    def summary(self) -> str:
        """One paragraph: what was built, from what evidence, and how it scored where it is hard."""
        n_inferred = int((~self.measured).sum())
        return (f"An integrated neighbour space over {self.ctx.n:,} genes: {len(self.pairs):,} "
                f"candidate pairs from {len(self.layers)} permitted layers "
                f"({', '.join(self.layers)}) and {len(self.columns)} permitted measurements, of "
                f"which {len(self.pairs) - n_inferred:,} are recorded by at least one layer and "
                f"{n_inferred:,} are inferred. Held-out edges score AUROC "
                f"{self.report.auroc('logistic', 'degree_matched'):.3f} against degree-matched "
                f"non-pairs and {self.report.auroc('logistic', 'random'):.3f} against random ones, "
                f"a fame gap of {self.report.fame_gap:+.3f}. "
                + (f"{self.note}. " if self.note else "")
                + "Inferred edges are a claim: no measured layer was changed.")

    def save(self, folder: str) -> list:
        """Write the space, the gaps, the weights and the three-null table as CSVs."""
        import os
        os.makedirs(folder, exist_ok=True)
        paths = []
        for name, table in (("graph_space", self.to_frame()), ("gaps", self.gaps(2000)),
                            ("weights", self.report.weights()), ("nulls", self.report.frame()),
                            ("per_layer", self.report.per_layer())):
            paths.append(f"{folder}/{name}.csv")
            table.to_csv(paths[-1], index=False)
        return paths


def space(ctx: S.Context, *, exclude=None, per_layer: int = PER_LAYER, knn: int = KNN,
          max_pairs: int = MAX_PAIRS, fraction: float = TEST_FRACTION, models=MODELS,
          embed_dim: int = EMBED_DIM, salt: int = 0) -> GraphSpace:
    """Build the integrated neighbour space, evaluating it on the way: ``space(ctx)`` and read it.

    The order matters and is not an implementation detail: the candidate pairs are collected, the
    standardisation is fixed, the evaluation runs (one fold per layer, each blind to its own layer),
    and only then are the folds' coefficients averaged into the weights that score every candidate.
    A source's weight is the mean of what it was worth in each fold, counting as zero the fold where
    it was the target -- so the shipped weights are a held-out quantity, not a fit to everything.

    `exclude` is a held-out label: its closure removes columns (`Context.banned`) and layers
    (`Context.banned_layers`) before anything is computed, which is what makes a space safe to use
    while scoring that label.
    """
    t0 = time.monotonic()
    ev = Evidence(ctx, exclude)
    if not ev.layers:
        raise ValueError("this table has no permitted edge layer to build a neighbour space from")
    ctx.say("collecting candidate pairs")
    pairs, note = candidates(ev, per_layer=per_layer, knn=knn, max_pairs=max_pairs)
    if not len(pairs):
        raise ValueError("no candidate pairs: the permitted layers are empty and the measurement "
                         "table has no columns to find nearest neighbours in")
    scaler = _pool_scaler(ev)
    report = evaluate(ctx, ev, exclude=exclude, models=models, fraction=fraction,
                      embed_dim=embed_dim, scaler=scaler, salt=salt)
    scored = report.scored
    if not scored:
        raise ValueError("no layer had enough edges between held-out node groups to learn a weight "
                         "from; the space cannot be scored on this table")
    names = ev.names(None)
    weights = {s: float(np.mean([f.coefficients.get(s, 0.0) for f in scored])) for s in names}
    intercept = float(np.mean([f.intercept for f in scored]))
    w = np.array([weights[s] for s in names])
    # Platt scaling on the pooled held-out pairs, each scored by the fold that could not see its own
    # layer. Fitting the calibration on training pairs instead is the classic way to ship a model
    # whose probabilities are confident and wrong.
    cal = [f.calibration for f in scored if f.calibration is not None]
    raw = np.concatenate([Z @ w + intercept for Z, _y in cal]) if cal else np.zeros(0)
    y = np.concatenate([yy for _Z, yy in cal]) if cal else np.zeros(0, bool)
    if len(raw) >= 20 and 0 < y.sum() < len(y):
        platt = _fit(raw.reshape(-1, 1), y, ctx.seed)
        A, B = float(platt.coef_[0][0]), float(platt.intercept_[0])
    else:                                   # nothing to calibrate on: report the raw logit as-is
        A, B = 1.0, 0.0
    ctx.say(f"scoring {len(pairs):,} candidate pairs")
    Z = (ev.features(pairs) - scaler[0]) / scaler[1]
    measured = np.asarray(ev.union(None)[pairs[:, 0], pairs[:, 1]]).ravel() > 0
    sp = GraphSpace(ev=ev, pairs=pairs, standardised=Z, probability=np.zeros(len(pairs)),
                    weights=weights, intercept=intercept, platt=(A, B), report=report,
                    measured=measured, note=note, seconds=time.monotonic() - t0,
                    settings={"exclude": exclude, "per_layer": per_layer, "knn": knn,
                              "max_pairs": max_pairs, "fraction": fraction,
                              "embed_dim": embed_dim, "models": tuple(models),
                              "scaler": scaler})
    sp.probability = sp._probability(Z)
    return sp


# --------------------------------------------------------------------------- the self-test
def _scorecard_ranking(score, positive) -> dict:
    from . import scorecard
    return scorecard.ranking(score, positive)


def hidden_edge_test(ctx: S.Context, key: str, *, layer: str | None = None, exclude=None,
                     model: str = "logistic", fraction: float = TEST_FRACTION,
                     embed_dim: int = EMBED_DIM, n_null: int = 10,
                     min_effect: float = 0.05) -> S.TestResult:
    """Hide one layer's edges by node group, ask for them back, and score against the HARD null.

    The pattern-3 shape (hidden pairs, an AUROC, a re-run null) with the two changes this module
    exists to make. The negatives are DEGREE-MATCHED, not random, so the score is not the discovery
    that hubs are hubs. And the null is not permuted identities but `n_null` configuration-model
    draws -- pair sets with the held-out edges' own degree sequence and none of their topology --
    which is exactly what a model reading fame would score as highly as the real thing. The
    easy number is still computed and reported in `numbers` as `auroc_random`, with the gap between
    the two as `fame_gap`, because a reader deserves to see how large that difference is.
    """
    t0 = time.monotonic()
    ev = Evidence(ctx, exclude)
    targets = ev.targets()
    if layer and layer not in ev.layers:
        raise ValueError(f"{layer!r} is not a layer this table permits here; permitted: "
                         f"{', '.join(ev.layers) or 'none'}")
    if layer and layer not in targets:
        raise ValueError(f"the {layer} layer is not one that can be held out: annotation layers "
                         f"({', '.join(TARGET_SKIP)}) and layers with fewer than "
                         f"{MIN_TARGET_EDGES} edges are not targets")
    layer = layer or (targets[0] if targets else None)
    if layer is None:
        return S.judge(key, "AUROC against degree-matched non-pairs", float("nan"), [],
                       min_effect=min_effect, n_hidden=0, hidden="edges",
                       null_kind="configuration model", t0=t0,
                       note="no layer has enough edges to hide and ask back")
    held = group_holdout(ctx, fraction)
    edges = ev.edges(layer)
    test = edges[held[edges[:, 0]] & held[edges[:, 1]]]
    train = edges[~held[edges[:, 0]] & ~held[edges[:, 1]]]
    nodes = np.flatnonzero(ev.degree(layer) > 0)
    all_keys = set(_pack(ev.all_edges()[:, 0], ev.all_edges()[:, 1], ev.n).tolist())
    rng = ctx.rng(1700)
    if len(test) < S.MIN_HIDDEN or len(train) < 20 or len(nodes) < 20:
        return S.judge(key, "AUROC against degree-matched non-pairs", float("nan"), [],
                       min_effect=min_effect, n_hidden=len(test),
                       hidden=f"{layer} edges between held-out orthogroups",
                       null_kind="configuration model", t0=t0,
                       note=f"only {len(test)} {layer} edges fall between held-out orthogroups")
    ctx.say(f"{key}: {len(test):,} {layer} edges hidden, {len(train):,} left to learn from")
    scaler = _pool_scaler(ev)
    names = ev.names(layer)
    keep = [list(ev.names(None)).index(s) for s in names]
    std = lambda F: (F - scaler[0][keep]) / scaler[1][keep]
    train_neg = degree_matched_negatives(rng, ev.degree(layer), train, all_keys, nodes, ev.n)
    if len(train_neg) < 10:
        return S.judge(key, "AUROC against degree-matched non-pairs", float("nan"), [],
                       min_effect=min_effect, n_hidden=len(test), hidden=f"{layer} edges",
                       null_kind="configuration model", t0=t0,
                       note="no degree-matched non-pairs could be drawn to train against")
    pairs_tr = np.vstack([train, train_neg])
    Ftr = std(ev.features(pairs_tr, skip=layer))
    ytr = np.r_[np.ones(len(train), bool), np.zeros(len(train_neg), bool)]
    emb = _spectral(ev.union(layer), embed_dim, ctx.seed) if model == "embedding" else None
    hadamard = _hadamard(emb, pairs_tr) or (lambda p: np.zeros((len(p), 0)))
    fitted = _fit(np.column_stack([Ftr, hadamard(pairs_tr)]), ytr, ctx.seed)
    score = lambda pairs: fitted.predict_proba(
        np.column_stack([std(ev.features(pairs, skip=layer)), hadamard(pairs)]))[:, 1]
    matched = degree_matched_negatives(rng, ev.degree(layer), test, all_keys, nodes, ev.n)
    if len(matched) < S.MIN_HIDDEN:
        return S.judge(key, "AUROC against degree-matched non-pairs", float("nan"), [],
                       min_effect=min_effect, n_hidden=len(test), hidden=f"{layer} edges",
                       null_kind="configuration model", t0=t0,
                       note="no degree-matched non-pairs could be drawn to score against")
    p_test, p_matched = score(test), score(matched)
    y = np.r_[np.ones(len(test), bool), np.zeros(len(matched), bool)]
    observed = S.auroc(np.r_[p_test, p_matched], y)
    easy = random_negatives(rng, nodes, all_keys, len(test), ev.n)
    auroc_random = S.auroc(np.r_[p_test, score(easy)],
                           np.r_[np.ones(len(test), bool), np.zeros(len(easy), bool)]) \
        if len(easy) else float("nan")
    nulls = []
    for i in range(int(n_null)):
        ctx.say(f"{key}: configuration-model null {i + 1} of {n_null}")
        rewired = configuration_negatives(rng, test, all_keys, ev.n)
        if len(rewired) < S.MIN_HIDDEN:
            continue
        nulls.append(S.auroc(np.r_[score(rewired), p_matched],
                             np.r_[np.ones(len(rewired), bool), np.zeros(len(matched), bool)]))
    m = metrics(np.r_[p_test, p_matched], y)
    details = pd.DataFrame([{"layer": layer, "hidden_edges": len(test), "train_edges": len(train),
                             "degree_matched_negatives": len(matched),
                             "sources": len(names), "model": model}])
    return S.judge(key, f"AUROC of hidden {layer} edges against degree-matched non-pairs",
                   observed, nulls, min_effect=min_effect, n_hidden=len(test),
                   hidden=f"the {layer} edges between held-out orthogroups ({len(test):,}), the "
                          f"{layer} layer itself removed from the features",
                   null_kind=f"{n_null} configuration-model rewirings of the hidden edges",
                   t0=t0, details=details,
                   numbers={"auroc_random": auroc_random, "fame_gap": auroc_random - observed,
                            "precision_at_k": m["precision_at_k"], "brier": m["brier"],
                            "reliability_gap": m["reliability_gap"], "model": model,
                            "layer": layer},
                   task="ranking", scorecard=_scorecard_ranking(np.r_[p_test, p_matched], y))
