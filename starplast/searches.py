#!/usr/bin/env python3
"""A search you can put down and pick up: every configuration it tried, and everything each found.

A climb is expensive and its result is not one map -- it is a hundred of them, with their scores,
their metrics and their claims. Keeping only the winner throws away the two things that make a
search evaluable afterwards: the runner-up, and how far behind it came. A single best configuration
with nothing to compare it against cannot be told apart from a lucky one.

## What is written, and what is left to be rebuilt

Written: the manifest (what was asked, of which table, by which version), one row per configuration
with its score and every metric, every finding from every configuration, and the cluster labels
each configuration produced.

Not written: the coordinates. They are the largest thing by an order of magnitude and the one thing
that is exactly reproducible -- a recipe plus a seed rebuilds the same map, which is a promise this
project makes everywhere else and may as well keep here. Labels ARE written, because kmeans and
HDBSCAN are not seeded identically across versions of their libraries and a clustering that came
back different on reload would silently change every claim built on it.

## The fingerprint, which is the point of the manifest

A search is a set of statements about genes identified by their POSITION in a node table. Reloaded
against a different table -- a rebuilt cache, a newer release, a filtered subset -- those positions
address different genes, and every finding in the run quietly becomes a claim about the wrong ones.
So the table is fingerprinted when the search is saved, checked when it is loaded, and a mismatch
is reported rather than papered over. It is not fatal: the configurations and scores remain
readable and comparable. What stops is the mapping from a label to a gene.
"""
from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

#: Columns `optimize.climb` carries for the caller rather than for the table.
ARTEFACTS = ("_labels", "_findings")


def fingerprint(nodes: pd.DataFrame) -> dict:
    """What table this search was about, in enough detail to notice it is a different one.

    The gene identifiers hashed in order, not merely counted: a table with the same number of rows
    in a different order is the commonest way for this to go wrong, and a count would not see it.
    """
    import hashlib
    if nodes is None:
        return {"n_genes": 0, "genes": "", "columns": 0}
    ids = (nodes["gene_id"].astype(str) if "gene_id" in nodes.columns
           else pd.Series(nodes.index.astype(str)))
    digest = hashlib.sha256("\n".join(ids.tolist()).encode()).hexdigest()[:16]
    return {"n_genes": int(len(nodes)), "genes": digest, "columns": int(nodes.shape[1])}


@dataclass
class Search:
    """One climb, complete: what was asked, what was tried, and what each attempt found."""
    manifest: dict = field(default_factory=dict)
    configs: pd.DataFrame = field(default_factory=pd.DataFrame)
    findings: pd.DataFrame = field(default_factory=pd.DataFrame)
    labels: dict = field(default_factory=dict)

    @property
    def name(self) -> str:
        """What this search is called on disk."""
        return str(self.manifest.get("name", "search"))

    @property
    def best(self):
        """The winning row, or None where nothing was evaluated."""
        return None if self.configs.empty else self.configs.iloc[0]

    def describe(self) -> str:
        """One line for a list of saved searches."""
        m = self.manifest
        best = "" if self.best is None else f", best {float(self.best.score):.2f}"
        return (f"{self.name}  ·  {m.get('mode', '?')} / {m.get('layer', '?')}  ·  "
                f"{len(self.configs)} configurations{best}, {len(self.findings)} findings")

    def for_config(self, index: int) -> tuple:
        """(row, its findings, its labels) for one configuration, by position in `configs`.

        By position in the saved table, which is sorted by score -- so 0 is the winner and 1 is the
        runner-up, which is the comparison this exists to make possible.
        """
        if self.configs.empty or not 0 <= index < len(self.configs):
            return None, pd.DataFrame(), None
        row = self.configs.iloc[index]
        key = int(row.get("config", index))
        mine = (self.findings[self.findings.config == key]
                if not self.findings.empty and "config" in self.findings.columns
                else pd.DataFrame())
        return row, mine.reset_index(drop=True), self.labels.get(key)

    def matches(self, nodes: pd.DataFrame) -> bool:
        """Is this the table the search was run against?"""
        saved = self.manifest.get("fingerprint", {})
        return bool(saved) and fingerprint(nodes) == saved

    def report(self, nodes: pd.DataFrame, index: int = 0, top: int = 8) -> str:
        """The reading for one configuration, with a warning where the table has moved on."""
        from . import interpret
        _row, found, _labels = self.for_config(index)
        text = interpret.report(found, nodes, top=top)
        if not self.matches(nodes):
            text = ("> **This search was run against a different node table.** Scores and "
                    "configurations still compare; the genes named below are addressed by position "
                    "and may not be the genes this run was about.\n\n") + text
        return text


def from_climb(result: pd.DataFrame, nodes: pd.DataFrame = None, **manifest) -> Search:
    """Split a climb's return into the part that is a table and the part that is artefacts."""
    import datetime
    from . import __version__
    result = pd.DataFrame() if result is None else result.reset_index(drop=True)
    configs = result.drop(columns=[c for c in ARTEFACTS if c in result.columns], errors="ignore")
    configs = configs.assign(config=range(len(configs)))
    findings, labels = [], {}
    for i, row in result.iterrows():
        found = row.get("_findings") if "_findings" in result.columns else None
        if found is not None and len(found):
            findings.append(pd.DataFrame(found).assign(config=i))
        lab = row.get("_labels") if "_labels" in result.columns else None
        if lab is not None:
            labels[int(i)] = np.asarray(lab)
    return Search(
        manifest={"created": datetime.datetime.now().isoformat(timespec="seconds"),
                  "version": __version__, "fingerprint": fingerprint(nodes), **manifest},
        configs=configs,
        findings=(pd.concat(findings, ignore_index=True) if findings else pd.DataFrame()),
        labels=labels)


class SearchStore:
    """Saved searches on disk, one directory each.

    A directory rather than a single file because the pieces have different natural formats and
    very different sizes -- a manifest that should stay readable in a text editor, two tables, and
    a block of integer arrays -- and because a half-written search should be identifiable as one.
    """

    def __init__(self, root: str):
        self.root = root
        os.makedirs(root, exist_ok=True)

    def path(self, name: str) -> str:
        """Where a search of this name lives, with anything awkward taken out of the name."""
        return os.path.join(self.root, re.sub(r"[^A-Za-z0-9._-]", "_", name))

    def list(self) -> list:
        """Saved searches, newest first, as (name, manifest)."""
        out = []
        for entry in sorted(os.listdir(self.root), reverse=True):
            manifest = os.path.join(self.root, entry, "manifest.json")
            if os.path.exists(manifest):
                try:
                    with open(manifest) as fh:
                        out.append((entry, json.load(fh)))
                except (OSError, ValueError):
                    # A half-written or hand-edited manifest is skipped rather than raised on: one
                    # unreadable directory must not make the whole list unopenable.
                    continue
        return out

    def save(self, search: Search, name: str = None) -> str:
        """Write a search. Returns the directory it went into."""
        from .runs import timestamp_name
        name = name or search.manifest.get("name") or timestamp_name(prefix="search")
        search.manifest["name"] = name
        where = self.path(name)
        os.makedirs(where, exist_ok=True)
        with open(os.path.join(where, "manifest.json"), "w") as fh:
            json.dump(search.manifest, fh, indent=2, default=str)
        # Findings hold lists of gene ids per row, which parquet will not take from an object
        # column without a schema. JSON keeps them as they are and stays readable.
        search.configs.to_csv(os.path.join(where, "configs.csv"), index=False)
        search.findings.to_json(os.path.join(where, "findings.json"), orient="records")
        if search.labels:
            np.savez_compressed(os.path.join(where, "labels.npz"),
                                **{str(k): v for k, v in search.labels.items()})
        return where

    def load(self, name: str) -> Search:
        """Read a search back. Missing pieces come back empty rather than raising."""
        where = self.path(name)
        manifest = {}
        if os.path.exists(os.path.join(where, "manifest.json")):
            with open(os.path.join(where, "manifest.json")) as fh:
                manifest = json.load(fh)
        configs = _read_csv(os.path.join(where, "configs.csv"))
        findings = _read_json(os.path.join(where, "findings.json"))
        labels = {}
        if os.path.exists(os.path.join(where, "labels.npz")):
            with np.load(os.path.join(where, "labels.npz")) as z:
                labels = {int(k): z[k] for k in z.files}
        return Search(manifest=manifest, configs=configs, findings=findings, labels=labels)


def _read_csv(path: str) -> pd.DataFrame:
    return pd.read_csv(path) if os.path.exists(path) else pd.DataFrame()


def _read_json(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        return pd.DataFrame()
    try:
        return pd.read_json(path, orient="records")
    except ValueError:                        # an empty findings file is "[]" or nothing at all
        return pd.DataFrame()
