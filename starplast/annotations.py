#!/usr/bin/env python3
"""Annotations: a proposed label for a gene, and the numbers that make it a proposal rather than a claim.

A cluster is mostly one category; the unlabeled genes in it are candidates. Writing those down is a
few lines of code and is the single easiest way for this project to start manufacturing exactly what
it exists to prevent -- because a candidate list looks identical whether it is 90% right or 6%, and
measured on the coarse level-of-detail components of this map it was **1-6%**.

So the store refuses a row that does not carry its own error rate. Every annotation records:

    gene, the label proposed, the configuration and cluster it came from, that cluster's
    composition, the VALIDATED precision and recall for that category, whether the validation
    re-embedded per fold, the date, and the reasoning.

Three rules, and each is a refusal rather than a convention:

* **Never written back into the node table.** A separate file, always. The node table is
  measurement; the moment an inference is stored beside it, every downstream reader has to know
  which columns are which, and one of them eventually will not.
* **Never saved without a validated precision for that category.** Not a warning -- an error. A
  precision from the same clustering the candidate came from is the only honest estimate of how
  often that annotation will be right, and a stored annotation without one is a guess wearing a
  record's clothes.
* **Its own color, used for nothing else.** Measurement, inference and absence already have one
  each in this application. An annotation is a fourth thing, and reading as any of the other three
  is the failure mode -- see `ANNOTATION_COLOR`.

The file is CSV rather than a database: it is small, a person should be able to read it, and the
reasoning field is the part a collaborator will want to argue with.
"""
from __future__ import annotations

import os
from dataclasses import asdict, dataclass, field, fields

import numpy as np
import pandas as pd

from .logging_util import get_logger

_log = get_logger(__name__)

#: The color an annotation is drawn in, and nothing else is. Deliberately outside the categorical
#: palette and outside the greys: a proposal must not read as a measured call, as an inference from
#: a held-out feature, or as absence.
ANNOTATION_COLOR = (0.95, 0.35, 0.85)

#: Below this validated precision, saving is refused outright rather than warned about. Measured on
#: this map's coarse components, annotating from them would be right 1-6% of the time; a store that
#: accepted those rows would be a machine for producing them.
MIN_PRECISION = 0.10


@dataclass
class Annotation:
    """One proposed label, with everything needed to judge it later.

    Every field is here because leaving it out makes the row unfalsifiable. Without the
    configuration nobody can rebuild the map it came from; without the cluster's composition nobody
    can see how mixed it was; without the validation numbers nobody can say how often it is wrong.
    """
    gene_id: str
    proposed: str                    # the label being proposed
    target: str                      # the column that label belongs to
    cluster: int
    blocks: str = ""                 # the configuration, as the search records it
    na_policy: str = ""
    scaling: str = ""
    n_neighbors: int = 0
    min_dist: float = 0.0
    min_cluster_size: int = 0
    seed: int = 0
    cluster_size: int = 0
    cluster_frac_category: float = float("nan")
    cluster_frac_contradicting: float = float("nan")
    enrichment: float = float("nan")
    precision: float = float("nan")          # VALIDATED, for this category
    recall: float = float("nan")
    folds: int = 0
    refit: bool = False
    date: str = ""
    reasoning: str = ""

    def check(self) -> str:
        """Why this annotation cannot be saved, or "" if it can.

        A sentence rather than an exception because the interface has to show it: refusing is the
        feature, and a dialog that says "error" teaches people the store is broken rather than that
        the annotation is unjustified.
        """
        if not self.gene_id or not self.proposed:
            return "an annotation needs a gene and a label"
        if not np.isfinite(self.precision):
            return ("this candidate has no validated precision for its category, so there is no "
                    "estimate of how often it would be right. Run the Validation tab for "
                    f"{self.target or 'that label'} first -- that is what makes it an annotation "
                    "rather than a guess.")
        if self.precision < MIN_PRECISION:
            return (f"validated precision for {self.proposed!r} is {self.precision:.0%}. Annotating "
                    f"from that cluster would be wrong about {1 - self.precision:.0%} of the time, "
                    f"which is not an annotation, it is a coin toss with a record attached.")
        return ""


#: Column order in the file. Fixed, so a file written by one version reads in another and a person
#: opening it in a spreadsheet finds the gene and the label first and the reasoning last.
COLUMNS = [f.name for f in fields(Annotation)]


class AnnotationStore:
    """Annotations on disk, in their own CSV, never in the node table.

    Re-read before every write rather than held in memory: this file is meant to be edited by hand
    and shared, and a store that cached it would silently overwrite whatever a collaborator added
    while the window was open.
    """

    def __init__(self, path: str):
        self.path = path

    def load(self) -> pd.DataFrame:
        """Everything saved so far, with the columns in their fixed order."""
        if not os.path.exists(self.path):
            return pd.DataFrame(columns=COLUMNS)
        df = pd.read_csv(self.path)
        for c in COLUMNS:
            if c not in df.columns:
                df[c] = pd.NA
        return df[COLUMNS + [c for c in df.columns if c not in COLUMNS]]

    def save(self, annotations, replace: bool = False) -> pd.DataFrame:
        """Append annotations, refusing any that cannot justify themselves.

        The refusal is per row and it raises: a partial save that silently dropped the unjustified
        rows would leave the user believing they had saved what they saw.

        A gene already annotated for the same target is REPLACED rather than duplicated -- the
        second opinion is the one someone just formed -- and the replacement is logged, because a
        changed annotation is exactly the kind of thing that should be traceable afterwards.
        """
        rows = [annotations] if isinstance(annotations, Annotation) else list(annotations)
        for a in rows:
            problem = a.check()
            if problem:
                raise ValueError(problem)
        new = pd.DataFrame([asdict(a) for a in rows], columns=COLUMNS)
        old = self.load()
        if len(old) and len(new):
            key = ["gene_id", "target"]
            merged = old.merge(new[key].drop_duplicates(), on=key, how="left", indicator=True)
            superseded = int((merged._merge == "both").sum())
            if superseded:
                _log.info("annotations: replacing %d existing row(s) for the same gene and target",
                          superseded)
                old = old[merged._merge != "both"].copy()
        out = new if replace else pd.concat([old, new], ignore_index=True)
        os.makedirs(os.path.dirname(os.path.abspath(self.path)) or ".", exist_ok=True)
        out.to_csv(self.path, index=False)
        _log.info("annotations: %d saved to %s (%d in the file)", len(new), self.path, len(out))
        return out

    def remove(self, gene_ids, target: str | None = None) -> pd.DataFrame:
        """Delete annotations for these genes. Withdrawing a proposal has to be as easy as making one."""
        df = self.load()
        if df.empty:
            return df
        drop = df.gene_id.isin(list(gene_ids))
        if target is not None:
            drop &= df.target == target
        out = df[~drop]
        out.to_csv(self.path, index=False)
        _log.info("annotations: removed %d row(s)", int(drop.sum()))
        return out

    def mask(self, gene_ids) -> np.ndarray:
        """Which of `gene_ids` carry an annotation, as a boolean array over that order.

        What the map draws its fourth color from. Built here rather than in the window so that
        "is this gene annotated" has one answer.
        """
        have = set(self.load().gene_id.astype(str)) if os.path.exists(self.path) else set()
        return np.array([str(g) in have for g in gene_ids], dtype=bool)


def from_candidates(candidates: pd.DataFrame, target: str, scores, configuration=None,
                    reasoning: str = "", date: str = "") -> list:
    """Turn `validate.candidates` rows into annotations, carrying the validation numbers with them.

    The join between a candidate and its error rate happens HERE, once, rather than at each call
    site: a candidate list and a precision that came from a different category, or a different
    clustering, would be worse than no number at all.
    """
    cfg = dict(configuration or {})
    out = []
    for _, r in candidates.iterrows():
        out.append(Annotation(
            gene_id=str(r.get("gene_id", "")), proposed=str(r.get("proposed", "")), target=target,
            cluster=int(r.get("cluster", -1)),
            blocks=str(cfg.get("blocks", "")), na_policy=str(cfg.get("na_policy", "")),
            scaling=str(cfg.get("scaling", "")),
            n_neighbors=int(float(cfg.get("n_neighbors", 0) or 0)),
            min_dist=float(cfg.get("min_dist", 0.0) or 0.0),
            min_cluster_size=int(float(cfg.get("min_cluster_size", 0) or 0)),
            seed=int(float(cfg.get("seed", 0) or 0)),
            cluster_size=int(r.get("cluster_size", 0)),
            cluster_frac_category=float(r.get("cluster_frac_category", float("nan"))),
            cluster_frac_contradicting=float(r.get("cluster_frac_contradicting", float("nan"))),
            enrichment=float(r.get("enrichment", float("nan"))),
            precision=float(getattr(scores, "precision", float("nan"))),
            recall=float(getattr(scores, "recall", float("nan"))),
            folds=int(getattr(scores, "n_folds", 0) or 0),
            refit=bool(getattr(scores, "refit", False)),
            date=date, reasoning=reasoning))
    return out
