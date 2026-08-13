#!/usr/bin/env python3
"""Clusterings, kept: named, dated, saved with their recipe, and knowing which genes they cover.

A clustering used to be computed, drawn, and lost. The second run replaced the first with no way
back, which makes the one comparison this application exists for -- does this structure hold under
different settings -- impossible to make by looking. Runs are kept here instead.

Three things a kept run has to carry, and each of them was missing:

* **A name.** Auto-generated from the timestamp to the second, so two runs a minute apart are
  distinguishable without anyone typing anything, and renameable afterwards because "hdbscan_60 on
  expression" is what a person will actually look for.
* **Its recipe.** An embedding spec and the clustering parameters. Without them a name is a label on
  nothing: nobody can rebuild the run, and two runs cannot be told apart except by their numbers.
* **Which genes it applies to.** Not how many -- which. A clustering of a subsample has fewer labels
  than the map has genes, and lining them up by position would put cluster 3's colour on whichever
  gene happens to sit at that index. The mask travels with the labels for exactly that reason.

Saved as one `.npz` of arrays plus one `.json` of the recipe, in the same shape as `EmbeddingStore`,
so a run and the embedding it came from are read the same way and neither is a special case.
"""
from __future__ import annotations

import json
import os
import re
from dataclasses import asdict, dataclass, field

import numpy as np
import pandas as pd

#: What a label of -1 means, everywhere: HDBSCAN put this gene in no cluster. Kept as a name because
#: "unclustered" is a finding about the gene, not a gap in the drawing, and the interface says so.
NOISE = -1
NOISE_NAME = "unclustered"


@dataclass
class Run:
    """One clustering, with everything needed to redraw it, name it, and rebuild it."""
    name: str
    labels: np.ndarray                       # one label per gene of `genes`, -1 for unclustered
    genes: np.ndarray                        # boolean mask over the node table
    recipe: dict = field(default_factory=dict)
    created: str = ""

    @property
    def n_clusters(self) -> int:
        lab = np.asarray(self.labels)
        return int(len({int(v) for v in lab if v != NOISE}))

    @property
    def noise_frac(self) -> float:
        lab = np.asarray(self.labels)
        return float((lab == NOISE).mean()) if len(lab) else 0.0

    @property
    def n_genes(self) -> int:
        return int(np.sum(self.genes))

    def describe(self) -> str:
        """One line for a list: the name, how many clusters, and how much it left unclustered."""
        return (f"{self.name}  ·  {self.n_clusters} clusters, "
                f"{self.noise_frac:.0%} unclustered, {self.n_genes:,} genes")

    def values(self, n_nodes: int) -> pd.Series:
        """The run as a categorical column over the WHOLE node table.

        Genes the run does not cover get "", which every part of this application already treats as
        absence -- they are not unclustered, they were not in the map. Conflating the two would put
        6,000 genes into a category that means something quite different.
        """
        out = np.full(n_nodes, "", dtype=object)
        genes = np.asarray(self.genes, dtype=bool)
        lab = np.asarray(self.labels)
        if genes.sum() != len(lab):
            # A run whose mask and labels disagree cannot be placed on any gene, so it is placed on
            # none. Better an empty column than 8,140 confident mislabels.
            return pd.Series(out, dtype=object)
        named = np.where(lab == NOISE, NOISE_NAME, np.char.add("cluster ", lab.astype(str)))
        out[genes] = named
        return pd.Series(out, dtype=object)


def timestamp_name(when=None, prefix: str = "run") -> str:
    """A name from the clock, to the second, so two runs a minute apart are distinguishable."""
    import datetime
    when = when or datetime.datetime.now()
    return f"{prefix}_{when.strftime('%Y%m%d_%H%M%S')}"


class RunStore:
    """Kept clusterings, in memory and on disk, in the order they were made.

    In memory as well as on disk because the panel lists them constantly and the disk copy exists so
    a run survives the window; `load_all` puts the file back into the list on the next start.
    """

    def __init__(self, root: str | None = None):
        self.root = root
        self.runs: list = []
        if root:
            os.makedirs(root, exist_ok=True)

    # ------------------------------------------------------------------ in memory
    def add(self, labels, genes, recipe=None, name: str = "", created: str = "") -> Run:
        """Keep a clustering. Unnamed, it takes a timestamp name."""
        import datetime
        run = Run(name=name or timestamp_name(), labels=np.asarray(labels),
                  genes=np.asarray(genes, dtype=bool), recipe=dict(recipe or {}),
                  created=created or datetime.datetime.now().isoformat(timespec="seconds"))
        self.runs.append(run)
        return run

    def get(self, name: str):
        """The run of that name, or None."""
        for r in self.runs:
            if r.name == name:
                return r
        return None

    def rename(self, old: str, new: str) -> bool:
        """Rename a run, refusing a name already in use.

        Refused rather than silently suffixed: two runs called the same thing in a list you choose
        from is worse than being told to pick another name.
        """
        run = self.get(old)
        if run is None or not new or self.get(new) is not None:
            return False
        if self.root:
            for ext in (".npz", ".json"):
                p = self._path(old, ext)
                if os.path.exists(p):
                    os.remove(p)
        run.name = new
        if self.root:
            self.save(run)
        return True

    def names(self) -> list:
        """Every kept run's name, in the order they were made -- what the colour-by list shows."""
        return [r.name for r in self.runs]

    # ------------------------------------------------------------------ on disk
    def _path(self, name: str, ext: str) -> str:
        safe = re.sub(r"[^A-Za-z0-9._-]", "_", name)
        return os.path.join(self.root, f"{safe}{ext}")

    def save(self, run: Run) -> str:
        """Write a run and its recipe. Without the recipe the name is a label on nothing."""
        if not self.root:
            return ""
        np.savez_compressed(self._path(run.name, ".npz"),
                            labels=np.asarray(run.labels), genes=np.asarray(run.genes))
        json.dump({"name": run.name, "created": run.created, "recipe": run.recipe,
                   "n_clusters": run.n_clusters, "noise_frac": run.noise_frac,
                   "n_genes": run.n_genes},
                  open(self._path(run.name, ".json"), "w"), indent=1)
        return self._path(run.name, ".npz")

    def load_all(self) -> list:
        """Read every saved run back into the list, newest last, skipping any that will not open."""
        if not self.root or not os.path.isdir(self.root):
            return self.runs
        for f in sorted(os.listdir(self.root)):
            if not f.endswith(".json"):
                continue
            try:
                meta = json.load(open(os.path.join(self.root, f)))
                z = np.load(os.path.join(self.root, f[:-5] + ".npz"), allow_pickle=False)
                name = meta.get("name", f[:-5])
                if self.get(name) is None:
                    self.runs.append(Run(name=name, labels=z["labels"], genes=z["genes"],
                                         recipe=meta.get("recipe", {}),
                                         created=meta.get("created", "")))
            except Exception as exc:
                # One unreadable run must not cost the rest of them, and saying which is what makes
                # the file worth going and looking at.
                print(f"starplast: could not read run {f}: {type(exc).__name__}: {exc}")
        return self.runs


def bin_column(values, bins: int = 5, log=None) -> pd.Series:
    """A numeric column as a categorical one, labelled by its ranges.

    Colouring a quantity by a ramp is already offered; binning makes it behave like a category,
    which is what makes it comparable with a clustering -- the comparison this panel exists for.

    Quantile bins, not equal-width: nearly every quantity in this table is heavy-tailed (the fitness
    screens span 64x within themselves), and equal-width bins put 95% of the genes in one of them
    and call the result a colouring. Missing values become "", which is absence everywhere else --
    they are not a low bin.

    **Fewer bins than asked for is an answer, not a failure, and it is reported rather than worked
    around.** `n_publications` is zero for most of this proteome, so its quartile edges are all zero
    and quantile binning yields one bin. Splitting the tie anyway -- by rank, or by equal width --
    would draw four colours over a column that has one level, which is a picture of a distinction
    that does not exist. The honest colouring of a mostly-constant quantity is mostly one colour.
    """
    s = pd.to_numeric(pd.Series(values), errors="coerce")
    out = pd.Series([""] * len(s), dtype=object)
    ok = s.notna().to_numpy()
    if ok.sum() < 2 or bins < 2:
        return out
    try:
        cats = pd.qcut(s[ok], q=int(bins), duplicates="drop", precision=2)
    except ValueError:                      # every value identical: one bin, not an error
        cats = None
    # A constant column does NOT raise: qcut returns a column of NaN categories, and casting those
    # to str gives the literal "nan" -- which would arrive in the interface as a category called
    # "nan", sitting in the legend beside real ones. Both paths end in one honest bin instead.
    if cats is None or cats.isna().all():
        out[ok] = "all one value"
        if log:
            log("every value is the same, so there is one bin")
        return out
    text = cats.astype(str)
    # Any remaining unbinnable value is absence, not a bin named "nan".
    out[ok] = [("" if c == "nan" else c) for c in text]
    got = len(set(out[ok]))
    if log and got < int(bins):
        log(f"{got} bins, not {bins}: the values are too tied for finer quantiles -- "
            f"{(s[ok] == s[ok].mode().iloc[0]).mean():.0%} of them share one value")
    return out
