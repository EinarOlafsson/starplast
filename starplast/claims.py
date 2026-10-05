"""Claims: knowledge generated with a measured certainty, then tested by evidence that is measured to be
independent of what generated it (`instructions/open/63_claims_generate_and_verify.md`).

Every number here comes from held-out genes -- the track record's folds -- never from an assumption.
Two of them are about leakage, because "these datasets do not overlap" is a guess: a signal peptide is
location information though it shares no dataset with LOPIT.

* :func:`shared_mistakes` -- generator vs verifier. When the generator is wrong, how often does the
  verifier make the SAME wrong call, against how often it would if the two were independent given the
  truth? 1.0 is independent; a random forest "confirming" kNN repeats its mistakes 3.5x chance.
* :func:`family_recovery` -- evidence vs label. How much of a label each evidence family recovers on
  its own. A family that alone recovers the label far better than the rest is either the label's own
  experiment restated or a strong biological signal -- worth a person's look either way.

And the two that make a claim:

* :func:`certainty_model` -- a strategy's support mapped to how often calls with that support were
  right when the answer was hidden (isotonic, cross-fitted, with its calibration error).
* :func:`verification` -- what a verifier's agreement and disagreement do to that certainty, measured.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from . import strategies as S
from . import track_record as T

#: Bootstrap resamples for the intervals here; enough for a stable 95% range on a few hundred genes.
BOOT = 300

#: A verifier is an independent test of a generator only when the upper end of the shared-mistake
#: ratio's 95% interval is below this. Measured, Tg compartment: physical partners and shared fold
#: sit at 1.2-1.4 against kNN and the classifier; coexpression diffusion at ~2, a random forest at ~3.5.
INDEPENDENT_BELOW = 1.6

#: A single evidence family that alone recovers this share of what ALL the evidence recovers (mean
#: per-class recall, same folds) is flagged. Tg compartment, 2026-10-04: the strongest family recovers
#: 0.19 against 0.36 for everything together (0.53) -- nothing flagged.
STANDS_OUT = 0.8


def _answered(ledger: pd.DataFrame, strategy: str) -> pd.DataFrame:
    f = ledger[(ledger["mode"] == "together") & (ledger["strategy"].astype(str) == str(strategy))]
    f = f[~f["abstained"].astype(bool)]
    return f.assign(truth=f["truth"].astype(str), prediction=f["prediction"].astype(str))


# --------------------------------------------------------------------------- leakage, measured
def shared_mistakes(ledger: pd.DataFrame, generator: str, verifier: str, boot: int = BOOT,
                    seed: int = 0) -> dict:
    """How much more often the verifier repeats the generator's wrong call than independence predicts.

    `ledger` is one organism and one label. Expected under independence: for each gene the generator
    got wrong (truth t, called p), the verifier's own rate of calling p when the truth is t, taken from
    all its answered genes. The ratio observed/expected is 1 for independent evidence and grows with
    shared information -- whatever its source, named or not.
    """
    g = _answered(ledger, generator).set_index("gene")[["truth", "prediction"]]
    v_all = _answered(ledger, verifier)
    v = v_all.set_index("gene")["prediction"].rename("v")
    j = g.join(v, how="inner")
    wrong = j[j["prediction"] != j["truth"]]
    out = {"generator": generator, "verifier": verifier, "both_answered": int(len(j)),
           "generator_wrong": int(len(wrong))}
    if len(wrong) < T.MIN_FOR_RATE:
        return {**out, "observed": np.nan, "expected": np.nan, "ratio": np.nan,
                "ratio_low": np.nan, "ratio_high": np.nan}
    conf = pd.crosstab(v_all["truth"], v_all["prediction"], normalize="index")
    expect = np.array([conf.at[t, p] if t in conf.index and p in conf.columns else 0.0
                       for t, p in zip(wrong["truth"], wrong["prediction"])])
    same = (wrong["v"].to_numpy() == wrong["prediction"].to_numpy()).astype(float)
    rng = np.random.default_rng(seed)
    ratios = []
    for _ in range(int(boot)):
        i = rng.integers(0, len(same), len(same))
        e = expect[i].mean()
        if e > 0:
            ratios.append(same[i].mean() / e)
    ratio = same.mean() / expect.mean() if expect.mean() > 0 else np.nan
    lo, hi = (np.percentile(ratios, [2.5, 97.5]) if ratios else (np.nan, np.nan))
    return {**out, "observed": float(same.mean()), "expected": float(expect.mean()),
            "ratio": float(ratio), "ratio_low": float(lo), "ratio_high": float(hi)}


def independent(ledger: pd.DataFrame, generator: str, verifier: str) -> bool:
    """True when the verifier's mistakes are measured to be (nearly) independent of the generator's."""
    m = shared_mistakes(ledger, generator, verifier)
    return bool(np.isfinite(m["ratio_high"]) and m["ratio_high"] < INDEPENDENT_BELOW)


def family_recovery(ctx, target: str, k: int = 15, min_columns: int = 1, log=None) -> pd.DataFrame:
    """Each evidence family ALONE, held out by orthogroup fold: how much of the label it recovers.

    A family is a dataset (`datasets.provenance`); columns without one form a family of their own.
    Reported as accuracy and mean per-class recall with what the commonest-class guess scores, so a
    family can be compared with the rest and with chance. Nothing is dropped here: a family that
    recovers a label very well is flagged for a person to read, because strong biology and a restated
    label look the same from the numbers alone.
    """
    from . import datasets
    say = log or (lambda _m: None)
    _X, cols = ctx.features(target)
    fam: dict = {}
    for c in cols:
        d = datasets.provenance(c, ctx.organism)
        fam.setdefault(d.key if d else f"(unregistered) {c}", []).append(c)
    truth = ctx.truth(target)
    kept = truth.where(truth.map(truth.value_counts()) >= S.MIN_CLASS)
    labelled = np.flatnonzero(kept.notna().to_numpy())
    folds = T.folds_by_group(ctx, labelled)
    every = kept.iloc[labelled].astype(str)
    commonest = float(every.value_counts(normalize=True).iloc[0])
    chance = 1.0 / every.nunique()
    rows = []
    for name, members in sorted(fam.items()):
        if len(members) < min_columns:
            continue
        X = ctx.matrix(members)
        pred = pd.Series([None] * len(kept), dtype=object)
        for f in range(int(folds.max()) + 1):
            hidden = labelled[folds == f]
            visible = kept.copy()
            visible.iloc[hidden] = np.nan
            p, _share = S.knn_vote(X, visible, k, query=hidden)
            pred.iloc[hidden] = pd.Series(p).iloc[hidden].to_numpy()
        got = pred.iloc[labelled]
        answered = got.notna()
        right = (got[answered].astype(str).to_numpy() == every[answered.to_numpy()].to_numpy())
        per_class = (pd.Series(right, index=every[answered.to_numpy()].to_numpy())
                     .groupby(level=0).mean().mean()) if answered.any() else np.nan
        rows.append({"family": name, "columns": len(members), "answered": int(answered.sum()),
                     "accuracy": float(right.mean()) if len(right) else np.nan,
                     "per_class": float(per_class), "commonest": commonest, "chance": chance})
        say(f"{name}: {rows[-1]['accuracy']:.2f}")
    out = pd.DataFrame(rows).sort_values("per_class", ascending=False).reset_index(drop=True)
    # The reference: the same kNN on every family together, same folds. A family that ALONE recovers
    # most of what all the evidence recovers is either the label's own experiment restated or one
    # dominant biological signal -- the numbers cannot tell which, so it is flagged for a person.
    X = ctx.matrix(cols)
    pred = pd.Series([None] * len(kept), dtype=object)
    for f in range(int(folds.max()) + 1):
        hidden = labelled[folds == f]
        visible = kept.copy()
        visible.iloc[hidden] = np.nan
        p, _share = S.knn_vote(X, visible, k, query=hidden)
        pred.iloc[hidden] = pd.Series(p).iloc[hidden].to_numpy()
    got = pred.iloc[labelled]
    ok = got.notna().to_numpy()
    right = got[ok].astype(str).to_numpy() == every[ok].to_numpy()
    reference = float(pd.Series(right, index=every[ok].to_numpy()).groupby(level=0).mean().mean())
    out["all_families"] = reference
    out["share_of_all"] = out["per_class"] / reference if reference > 0 else np.nan
    out["stands_out"] = out["share_of_all"] >= STANDS_OUT
    return out


# --------------------------------------------------------------------------- certainty, measured
@dataclass
class CertaintyModel:
    """A strategy's support mapped to the held-out rate of being right, for one organism and label."""
    strategy: str
    target: str
    n: int
    calibration_error: float          # cross-fitted expected calibration error, over 10 bins
    _x: np.ndarray = field(repr=False, default_factory=lambda: np.zeros(0))
    _y: np.ndarray = field(repr=False, default_factory=lambda: np.zeros(0))

    def __call__(self, support) -> np.ndarray:
        s = np.asarray(support, dtype=float)
        if not len(self._x):
            return np.full(s.shape, np.nan)
        return np.interp(np.nan_to_num(s, nan=float(self._x[0])), self._x, self._y)

    @property
    def calibrated(self) -> bool:
        return bool(self.n >= 100 and self.calibration_error <= 0.05)


def _isotonic(x, y):
    from sklearn.isotonic import IsotonicRegression
    m = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0).fit(x, y)
    xs = np.unique(x)
    return xs, m.predict(xs)


def _ece(p, y, bins: int = 10) -> float:
    p, y = np.asarray(p, float), np.asarray(y, float)
    edges = np.quantile(p, np.linspace(0, 1, bins + 1))
    idx = np.clip(np.searchsorted(edges, p, side="right") - 1, 0, bins - 1)
    return float(sum(abs(p[idx == b].mean() - y[idx == b].mean()) * (idx == b).mean()
                     for b in range(bins) if (idx == b).any()))


def certainty_model(ledger: pd.DataFrame, strategy: str) -> CertaintyModel:
    """Fit support -> P(right) on held-out calls, and measure its error out of fold.

    The calibration error is cross-fitted: each fold's calls are scored by a curve fit on the other
    folds, so the error is what a new gene would see, not how well the curve memorised its own points.
    """
    a = _answered(ledger, strategy)
    a = a[np.isfinite(a["support"].astype(float))]
    target = str(a["target"].iloc[0]) if len(a) else ""
    if len(a) < 30:
        return CertaintyModel(strategy, target, int(len(a)), np.nan)
    x = a["support"].astype(float).to_numpy()
    y = (a["prediction"] == a["truth"]).to_numpy().astype(float)
    fold = a["fold"].to_numpy()
    oof = np.full(len(y), np.nan)
    for f in np.unique(fold):
        tr, te = fold != f, fold == f
        if tr.sum() < 20 or not te.any():
            continue
        xs, ys = _isotonic(x[tr], y[tr])
        oof[te] = np.interp(x[te], xs, ys)
    ok = np.isfinite(oof)
    xs, ys = _isotonic(x, y)
    return CertaintyModel(strategy, target, int(len(y)), _ece(oof[ok], y[ok]) if ok.sum() else np.nan,
                          xs, ys)


def verification(ledger: pd.DataFrame, generator: str, verifier: str) -> dict:
    """What the verifier's verdict does to the generator's chance of being right, on held-out genes.

    Three outcomes, each with its rate and Wilson interval: the verifier agrees, disagrees, or says
    nothing. Measured on the same genes, so any shared information between the two is already in
    these numbers -- and `shared_mistakes` says how much of it there is.
    """
    g = _answered(ledger, generator).set_index("gene")[["truth", "prediction"]]
    v = _answered(ledger, verifier).set_index("gene")["prediction"].rename("v")
    j = g.join(v, how="left")
    right = (j["prediction"] == j["truth"])
    verdict = np.where(j["v"].isna(), "silent", np.where(j["v"] == j["prediction"], "agrees", "disagrees"))
    out = {"generator": generator, "verifier": verifier, "base": _rate_of(right)}
    for name in ("agrees", "disagrees", "silent"):
        out[name] = _rate_of(right[verdict == name])
    return out


def _rate_of(right) -> dict:
    """n, right, rate and Wilson interval for a plain right/wrong series."""
    from .calibration import wilson
    right = np.asarray(right, dtype=bool)
    n, k = int(len(right)), int(right.sum())
    low, high = wilson(k, n) if n else (np.nan, np.nan)
    return {"n": n, "right": k, "rate": k / n if n else np.nan, "low": low, "high": high,
            "enough": n >= T.MIN_FOR_RATE}


# --------------------------------------------------------------------------- verified certainty
@dataclass
class VerifiedModel:
    """Generator certainty plus independent verdicts -> P(right), fit and checked on held-out genes."""
    generator: str
    verifiers: tuple
    target: str
    n: int
    calibration_error: float          # cross-fitted, like the certainty model's
    certainty: CertaintyModel
    _model: object = field(repr=False, default=None)

    def design(self, certainty, verdicts: dict) -> np.ndarray:
        """Rows: logit(certainty), then (agrees, disagrees) per verifier; silence is both zero."""
        c = np.clip(np.asarray(certainty, float), 1e-3, 1 - 1e-3)
        cols = [np.log(c / (1 - c))]
        for v in self.verifiers:
            said = np.asarray(verdicts.get(v, ["silent"] * len(c)))
            cols += [(said == "agrees").astype(float), (said == "disagrees").astype(float)]
        return np.column_stack(cols)

    def __call__(self, certainty, verdicts: dict) -> np.ndarray:
        if self._model is None or not len(np.atleast_1d(certainty)):
            return np.asarray(certainty, float)
        return self._model.predict_proba(self.design(certainty, verdicts))[:, 1]

    @property
    def calibrated(self) -> bool:
        return bool(self.n >= 100 and self.calibration_error <= 0.05)


def _verdicts(ledger: pd.DataFrame, generator: str, verifiers, genes) -> dict:
    g = _answered(ledger, generator).set_index("gene")["prediction"]
    out = {}
    for v in verifiers:
        said = _answered(ledger, v).set_index("gene")["prediction"]
        call = said.reindex(genes)
        mine = g.reindex(genes)
        out[v] = np.where(call.isna(), "silent", np.where(call == mine, "agrees", "disagrees"))
    return out


def verified_model(ledger: pd.DataFrame, generator: str, verifiers) -> VerifiedModel:
    """Fit the combination on held-out genes; measure its calibration error out of fold.

    Logistic regression on the generator's (cross-fitted) certainty and each verifier's verdict. It
    never assumes the verifiers are independent of each other or of the generator: their weights are
    whatever the held-out genes say they are worth, shared information included.
    """
    from sklearn.linear_model import LogisticRegression
    verifiers = tuple(verifiers)
    cert = certainty_model(ledger, generator)
    a = _answered(ledger, generator)
    a = a[np.isfinite(a["support"].astype(float))]
    model = VerifiedModel(generator, verifiers, cert.target, int(len(a)), np.nan, cert)
    if len(a) < 100:
        return model
    genes = a["gene"].to_numpy()
    y = (a["prediction"] == a["truth"]).to_numpy().astype(int)
    fold = a["fold"].to_numpy()
    support = a["support"].astype(float).to_numpy()
    verdicts = _verdicts(ledger, generator, verifiers, genes)
    oof = np.full(len(y), np.nan)
    for f in np.unique(fold):
        tr, te = fold != f, fold == f
        xs, ys = _isotonic(support[tr], y[tr])
        c_tr = np.interp(support[tr], xs, ys)
        c_te = np.interp(support[te], xs, ys)
        sub = {v: d[tr] for v, d in verdicts.items()}
        m = LogisticRegression(C=1.0, max_iter=1000).fit(model.design(c_tr, sub), y[tr])
        oof[te] = m.predict_proba(model.design(c_te, {v: d[te] for v, d in verdicts.items()}))[:, 1]
    model.calibration_error = _ece(oof, y)
    model._model = LogisticRegression(C=1.0, max_iter=1000).fit(
        model.design(cert(support), verdicts), y)
    model._oof = oof
    model._y = y
    return model


# --------------------------------------------------------------------------- claims about unknown genes
#: Strategies asked to verify. Any of them qualifies only if measured independent of the generator.
VERIFIER_CANDIDATES = ("physical_partners", "structural_homology", "layer_propagation", "feature_knn",
                       "supervised_classifier", "random_forest", "map_neighbours")


@dataclass
class Recipe:
    """A generator, the verifiers measured to be independent of it, and their fitted combination."""
    organism: str
    target: str
    generator: str
    verifiers: tuple
    independence: pd.DataFrame
    model: VerifiedModel

    @property
    def proven(self) -> bool:
        """Calibrated end to end on held-out genes: its certainties can be taken at face value."""
        return self.model.certainty.calibrated and self.model.calibrated


#: Generators considered. Stacking and graph convolution are left out on purpose: they absorb every
#: kind of evidence, so nothing is left to test them with -- measured, Tg compartment: no verifier's
#: shared-mistake ratio with stacking is below 1.38 even at the low end of its interval.
GENERATOR_CANDIDATES = ("supervised_classifier", "feature_knn", "random_forest", "physical_partners",
                        "structural_homology", "layer_vote", "map_neighbours")

#: A recipe is judged by the held-out calls it gets right at or above this verified certainty.
CONFIDENT = 0.8


def _independence(led, generator) -> pd.DataFrame:
    rows = []
    for v in VERIFIER_CANDIDATES:
        if v == generator:
            continue
        m = shared_mistakes(led, generator, v)
        m["independent"] = bool(np.isfinite(m["ratio_high"]) and m["ratio_high"] < INDEPENDENT_BELOW)
        rows.append(m)
    return pd.DataFrame(rows)


def candidates(organism: str, target: str | None = None,
               ledger: pd.DataFrame | None = None) -> pd.DataFrame:
    """Every generator with the verifiers measured independent of it, and what the pair achieves.

    `confident_right` is the number of held-out genes whose verified certainty reached `CONFIDENT` and
    whose claim was right; `confident_precision` is how often such claims were right. A generator with
    no independent verifier is listed with `testable` False and is never chosen.
    """
    frame = T.shipped(organism) if ledger is None else ledger
    target = target or T.default_target(organism, frame)
    led = frame[frame["target"].astype(str) == str(target)]
    rows = []
    for g in GENERATOR_CANDIDATES:
        cert = certainty_model(led, g)
        if not np.isfinite(cert.calibration_error):
            continue
        table = _independence(led, g)
        chosen = tuple(table[table["independent"]]["verifier"]) if len(table) else ()
        row = {"generator": g, "verifiers": chosen, "testable": bool(chosen),
               "certainty_error": cert.calibration_error, "answered": cert.n}
        if chosen:
            m = verified_model(led, g, chosen)
            if np.isfinite(m.calibration_error):
                sel = m._oof >= CONFIDENT
                row.update({"verified_error": m.calibration_error,
                            "confident": int(sel.sum()),
                            "confident_right": int(m._y[sel].sum()),
                            "confident_precision": float(m._y[sel].mean()) if sel.any() else np.nan,
                            "calibrated": m.calibrated and cert.calibrated})
        rows.append(row)
    out = pd.DataFrame(rows)
    if "confident_right" not in out:
        out["confident_right"] = np.nan
    return out.sort_values("confident_right", ascending=False, na_position="last").reset_index(drop=True)


def recipe(organism: str, target: str | None = None, generator: str | None = None,
           ledger: pd.DataFrame | None = None) -> Recipe:
    """The generator, testable by independent evidence, whose verified claims are most often right.

    Chosen among calibrated, testable generators by `confident_right`; an explicit `generator` is
    honoured, with whatever verifiers are measured independent of it (possibly none).
    """
    frame = T.shipped(organism) if ledger is None else ledger
    target = target or T.default_target(organism, frame)
    led = frame[frame["target"].astype(str) == str(target)]
    if generator is None:
        table = candidates(organism, target, frame)
        ok = table[table["testable"] & table.get("calibrated", False).fillna(False).astype(bool)]
        generator = str(ok["generator"].iloc[0]) if len(ok) else "feature_knn"
    independence = _independence(led, generator)
    chosen = tuple(independence[independence["independent"]]["verifier"]) if len(independence) else ()
    return Recipe(organism, target, generator, chosen, independence,
                  verified_model(led, generator, chosen))


def _answered_rows(ledger, strategy):
    return ledger[(ledger["mode"] == "together") & (ledger["strategy"].astype(str) == str(strategy))]


def generate(ctx, rec: Recipe, log=None) -> pd.DataFrame:
    """Claims about every gene with no label: what the generator says, what each verifier says, and the
    verified certainty -- all labelled genes visible, as in real use."""
    say = log or (lambda _m: None)
    if not rec.proven:
        # Uncalibrated certainty is not shown as certainty anywhere; the binary screen phenotypes
        # once produced thousands of "confident" claims that were the base rate of "no phenotype".
        say(f"{rec.generator} on {rec.target} is not calibrated end to end; no claims")
        return pd.DataFrame(columns=["gene", "gene_id", "claim", "status", "confidence"])
    truth = ctx.truth(rec.target)
    kept = truth.where(truth.map(truth.value_counts()) >= S.MIN_CLASS)
    query = np.flatnonzero(truth.isna().to_numpy())

    def ask(key):
        strategy = S.get(key)
        chosen = T._permitted(ctx, dict(strategy.settings(ctx, target=rec.target)), rec.target)
        say(f"{key}: calling {len(query):,} unlabelled genes")
        pred, support = T._predict(ctx, strategy, kept, query, chosen, rec.target)
        if pred is None:
            return pd.Series([None] * len(query)), pd.Series([np.nan] * len(query))
        pred = pd.Series(pred).reset_index(drop=True).iloc[query].reset_index(drop=True)
        sup = (pd.Series(support).reset_index(drop=True).iloc[query].reset_index(drop=True)
               if support is not None else pd.Series([np.nan] * len(query)))
        return pred, sup

    claim, support = ask(rec.generator)
    said = claim.notna().to_numpy()
    out = pd.DataFrame({"organism": rec.organism, "target": rec.target, "gene": query,
                        "gene_id": np.asarray(ctx.gene_ids)[query],
                        "claim": claim.to_numpy(), "support": support.to_numpy(dtype=float),
                        "generator": rec.generator})[said].reset_index(drop=True)
    if not len(out):                            # every gene already labelled, or the generator silent
        return out.assign(certainty=[], verified_certainty=[], distance=[], in_range=[], status=[],
                          confidence=[])
    out["certainty"] = rec.model.certainty(out["support"].to_numpy())
    verdicts = {}
    for v in rec.verifiers:
        call, _s = ask(v)
        call = call[said].reset_index(drop=True)
        verdicts[v] = np.where(call.isna(), "silent",
                               np.where(call.astype(str) == out["claim"].astype(str), "agrees", "disagrees"))
        out[f"by {v}"] = np.where(call.isna(), None, call.astype(object))
        out[f"{v} verdict"] = verdicts[v]
    out["verified_certainty"] = rec.model(out["certainty"].to_numpy(), verdicts)
    # Where the certainty was measured: distance to the nearest labelled genes in the evidence space,
    # against the same distance among labelled genes. Beyond their 95th percentile no certainty was
    # ever measured -- 42% of Toxoplasma's unlabelled genes on compartment, mostly genes LOPIT could
    # not detect -- so those claims say so instead of carrying a number.
    distance, limit = applicability(ctx, rec.target, out["gene"].to_numpy())
    out["distance"] = distance
    out["in_range"] = distance <= limit
    tested = np.zeros(len(out), dtype=bool)
    for v in rec.verifiers:
        tested |= out[f"{v} verdict"].to_numpy() != "silent"
    out["status"] = np.where(~out["in_range"], "outside tested range",
                             np.where(tested, "tested", "untested"))
    # The one number to read: verified where tested, certainty where not, nothing where unmeasured.
    out["confidence"] = np.where(out["status"] == "tested", out["verified_certainty"],
                                 np.where(out["status"] == "untested", out["certainty"], np.nan))
    # Against the prior: how common the claimed class is among labelled genes. Confidence 0.85 in a
    # class that is 80% of genes says little; in one that is 2% it says a great deal.
    prior = kept.dropna().astype(str).value_counts(normalize=True)
    out["prior"] = out["claim"].astype(str).map(prior).to_numpy(dtype=float)
    out["lift"] = out["confidence"] / out["prior"]
    return out.sort_values(["confidence"], ascending=False, na_position="last").reset_index(drop=True)


#: Neighbours averaged for the applicability distance, and the labelled-gene percentile beyond which a
#: gene is outside the range where certainty was measured.
RANGE_K = 15
RANGE_PERCENTILE = 95


def applicability(ctx, target: str, positions) -> tuple:
    """(distance per gene, limit): mean distance to the `RANGE_K` nearest labelled genes in the
    evidence space, and the `RANGE_PERCENTILE` of that distance among labelled genes themselves."""
    from sklearn.neighbors import NearestNeighbors
    truth = ctx.truth(target)
    X, _cols = ctx.features(target)
    lab = np.flatnonzero(truth.notna().to_numpy())
    nn = NearestNeighbors(n_neighbors=RANGE_K + 1).fit(X[lab])
    own = nn.kneighbors(X[lab])[0][:, 1:].mean(axis=1)           # each labelled gene's self excluded
    limit = float(np.percentile(own, RANGE_PERCENTILE))
    d = nn.kneighbors(X[np.asarray(positions, dtype=int)], n_neighbors=RANGE_K)[0].mean(axis=1)
    return d, limit


# --------------------------------------------------------------------------- the shipped claims
_SHIPPED: dict = {}


def _data(name: str) -> str:
    import os
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", name)


def shipped(organism: str | None = None) -> pd.DataFrame:
    """The built claims (`scripts/build_claims.py`), loaded once. Empty where none were built."""
    import os
    if "all" not in _SHIPPED:
        path = _data("claims.parquet")
        try:
            _SHIPPED["all"] = pd.read_parquet(path) if os.path.exists(path) else pd.DataFrame()
        except Exception:                                   # a half-written or foreign file
            _SHIPPED["all"] = pd.DataFrame()
    frame = _SHIPPED["all"]
    if organism is None or not len(frame):
        return frame
    key = ("organism", str(organism))
    if key not in _SHIPPED:
        _SHIPPED[key] = frame[frame["organism"] == str(organism)].reset_index(drop=True)
    return _SHIPPED[key]


def recipes(organism: str | None = None) -> pd.DataFrame:
    """Which recipe each label uses, whether it is proven, and how many claims it made."""
    import os
    if "recipes" not in _SHIPPED:
        path = _data("claim_recipes.parquet")
        _SHIPPED["recipes"] = pd.read_parquet(path) if os.path.exists(path) else pd.DataFrame()
    frame = _SHIPPED["recipes"]
    return frame[frame["organism"] == organism] if organism and len(frame) else frame


#: How each status reads, as a mark and a phrase.
MARKS = {"tested": "✓", "untested": "○", "outside tested range": "·"}


def verdict_mark(row) -> str:
    """✓ when every independent check that spoke agreed, ✗ when any disagreed, else the status mark."""
    said = [str(row[c]) for c in row.index if str(c).endswith(" verdict") and pd.notna(row[c])]
    if "disagrees" in said:
        return "✗"
    if "agrees" in said:
        return "✓"
    return MARKS.get(str(row.get("status")), "·")


def for_gene(gene_id: str, organism: str) -> pd.DataFrame:
    frame = shipped(organism)
    if not len(frame):
        return frame
    return frame[frame["gene_id"].astype(str) == str(gene_id)]


def gene_html(gene_id: str, organism: str, row: int | None = None) -> str:
    """The gene card's claims, condensed: one line per label, each a click from its reasoning."""
    from html import escape
    from urllib.parse import quote
    mine = for_gene(gene_id, organism)
    if not len(mine):
        return ""
    lines = ""
    for r in mine.itertuples(index=False):
        rr = pd.Series(r._asdict())
        mark = verdict_mark(rr)
        conf = ("not measured" if not np.isfinite(rr["confidence"])
                else f"{rr['confidence']:.0%}")
        link = (f"starplast://claim/{quote(str(gene_id), safe='')}/{quote(str(rr['target']), safe='')}")
        colour = {"✓": "#2e8b57", "✗": "#c0392b"}.get(mark, "#888888")
        lines += (f"<tr><td style='color:{colour}'>{mark}</td>"
                  f"<td>{escape(str(rr['target']).replace('_', ' '))}</td>"
                  f"<td><a href='{link}'>{escape(str(rr['claim']))}</a></td>"
                  f"<td align='right'>{conf}</td></tr>")
    return ("<h4>What Starplast claims</h4>"
            "<p style='color:#888'>Inferred for labels this gene does not have. ✓ an independent check "
            "agreed · ✗ one disagreed · ○ nothing independent reaches it · · unlike any gene the "
            "certainty was measured on. Click a claim for how it was made and tested.</p>"
            f"<table cellspacing='0' cellpadding='3'>{lines}</table>")


def claim_html(gene_id: str, target: str, organism: str) -> str:
    """One claim in full: what was claimed, how certain, how it was tested, and what that is worth."""
    from html import escape
    mine = for_gene(gene_id, organism)
    mine = mine[mine["target"].astype(str) == str(target)]
    if not len(mine):
        return ""
    r = mine.iloc[0]
    rec = recipes(organism)
    rec = rec[rec["target"].astype(str) == str(target)]
    gen = str(r["generator"]).replace("_", " ")
    status = str(r["status"])
    lead = (f"<b>{escape(str(r['claim']))}</b> for {escape(str(target).replace('_', ' '))}, claimed "
            f"by {escape(gen)}.")
    if status == "outside tested range":
        body = ("<p>This gene is unlike every gene the certainty was measured on: it is farther from "
                "the labelled genes than 95% of them are from each other. The claim is shown, but no "
                "certainty is attached to it, because none was ever measured for genes like it.</p>")
    else:
        body = (f"<p>Certainty from the generator alone: <b>{r['certainty']:.0%}</b> -- of held-out "
                f"genes given the same support, that share were right. The class is "
                f"{r['prior']:.0%} of labelled genes, so this is {r['certainty'] / r['prior']:.1f}x "
                f"the prior.</p>")
        checks = ""
        for c in [c for c in mine.columns if str(c).endswith(" verdict")]:
            v = r[c]
            if pd.isna(v):
                continue
            name = str(c)[:-len(" verdict")]
            said = r.get(f"by {name}")
            checks += (f"<tr><td>{escape(name.replace('_', ' '))}</td><td>{escape(str(v))}</td>"
                       f"<td style='color:#888'>{escape(str(said)) if pd.notna(said) else '—'}</td></tr>")
        if checks:
            body += ("<p>Independent checks -- strategies whose mistakes were measured not to repeat "
                     "the generator's:</p><table cellspacing='0' cellpadding='3'><tr><th align='left'>"
                     "check</th><th align='left'>verdict</th><th align='left'>it said</th></tr>"
                     f"{checks}</table>")
        if status == "tested":
            body += (f"<p>Certainty after the checks: <b>{r['verified_certainty']:.0%}</b>, fit and "
                     f"measured on held-out genes.</p>")
        else:
            body += "<p>No independent check reaches this gene, so the claim is untested.</p>"
    if len(rec):
        q = rec.iloc[0]
        if pd.notna(q.get("heldout_confident_precision")):
            body += (f"<p style='color:#888'>This recipe, on held-out genes: of "
                     f"{int(q['heldout_confident']):,} claims at 80% or more after checks, "
                     f"{q['heldout_confident_precision']:.0%} were right. Calibration error "
                     f"{q['certainty_error']:.3f} (certainty), {q['verified_error']:.3f} (after "
                     f"checks).</p>")
    return f"<h4>The claim</h4><p>{lead}</p>{body}"


def discoveries(organism: str, target: str | None = None, status=("tested",),
                min_confidence: float = 0.8, min_lift: float = 2.0) -> pd.DataFrame:
    """The claims worth a person's time: confident, well above the prior, and (by default) tested."""
    frame = shipped(organism)
    if not len(frame):
        return frame
    keep = frame["status"].astype(str).isin(status) & (frame["confidence"] >= min_confidence) \
        & (frame["lift"] >= min_lift)
    if target:
        keep &= frame["target"].astype(str) == str(target)
    return frame[keep].sort_values("confidence", ascending=False).reset_index(drop=True)
