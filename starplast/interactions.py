#!/usr/bin/env python3
"""Physical binding and structural similarity — and the models that show how the binding happens.

Three edge types come from work already done elsewhere in this tree, all of it already keyed by
`TGME49_` gene id, so none of it needs the identity layer:

* **xlms** — StarPath DSS crosslinking mass-spectrometry. This is the only edge type in starplast that is
  *measured physical proximity*: two residues were covalently joined in a living-cell lysate, so the two
  proteins were within the crosslinker's reach. Weight is the number of crosslinks supporting the pair.
* **ip_ms** — immunoprecipitation of a tagged bait, replicated, with an untagged control. Few, but the
  most direct binding evidence in the map.
* **struct** — Foldseek TM-align over 6,900 Toxoplasma AlphaFold models, TM >= 0.7. Structural similarity
  needs no orthology, so it reaches lineage-specific effectors that homology edges cannot see.

None of these is a co-mention, so none is attention-biased. That matters: `xlms` and `struct` are the
first edge types here that can connect a gene nobody has written about to one everybody has.

**The models.** StarPath ships four Chai-1 predicted complexes per interaction — 12,265 CIF files on
local disk — and `starpath_crosslink_mining` scored 2,411 of them for whether the predicted complex
actually places the crosslinked residues within reach. `crosslink_models()` assembles that into one table
per gene pair: residue positions, paths to the model files, the fraction of crosslinks the model
satisfies, and Chai-1's own confidence. A low `frac_satisfied` is informative rather than fatal -- it says
the model does not explain the measurement, not that the measurement is wrong.
"""
from __future__ import annotations

import json
import os
import re

import numpy as np
import pandas as pd

# Edge types to lift out of the toxonet edge table, and what their weight means.
TOXONET_EDGES = {
    "xlms": "number of DSS crosslinks supporting the pair",
    "ip_ms": "replicated IP-MS of a tagged bait",
    "struct": "Foldseek TM-score (>= 0.7)",
}


def load_toxonet(path: str, node_ids, log=print):
    """Lift the physical-interaction and structural-similarity edges out of the toxonet edge table."""
    if not os.path.exists(path):
        log(f"toxonet edges absent ({path}); binding and structure layers will be empty")
        return {}, pd.DataFrame()

    e = pd.read_parquet(path)
    known = set(node_ids)
    idx = {g: i for i, g in enumerate(node_ids)}
    edges, keep = {}, []
    for label in TOXONET_EDGES:
        s = e[e.edge_type == label]
        if s.empty:
            continue
        s = s[s.src.isin(known) & s.dst.isin(known)]
        # canonical undirected order, and never trust an upstream table to be free of duplicates
        pair = pd.DataFrame({"a": np.minimum(s.src.map(idx), s.dst.map(idx)),
                             "b": np.maximum(s.src.map(idx), s.dst.map(idx)),
                             "w": s.weight.astype(float).to_numpy()})
        pair = pair[pair.a != pair.b].groupby(["a", "b"], as_index=False).w.max()
        if pair.empty:
            continue
        w = pair.w.to_numpy()
        edges[label] = (pair.a.to_numpy(), pair.b.to_numpy(), w, w.copy())
        keep.append(s)
        log(f"{label}: {len(pair):,} edges over {len(set(pair.a) | set(pair.b)):,} genes "
            f"({TOXONET_EDGES[label]})")
    detail = pd.concat(keep, ignore_index=True) if keep else pd.DataFrame()
    return edges, detail


# --------------------------------------------------------------------------- the binding models
def _num(pattern, text, cast=float):
    m = re.search(pattern, str(text))
    return cast(m.group(1)) if m else np.nan


def crosslink_models(base: str, node_ids, log=print) -> pd.DataFrame:
    """One row per crosslinked gene pair: the residues, the model files, and how well they agree.

    StarPath identifies proteins by RH88 accession, whose numbering does NOT correspond to ME49 --
    TGRH88_016370 is TGME49_210408, not TGME49_216370 -- so the mapping comes from the alias column the
    export ships, never from the suffix rule used for GT1 and VEG in identity.py.
    """
    inter_csv = os.path.join(base, "starpath_interactions.csv")
    xl_json = os.path.join(base, "starpath_crosslinks.json")
    sat_csv = os.path.join(base, "starpath_crosslink_mining", "crosslink_satisfaction.csv")
    cif_dir = os.path.join(base, "starpath_dump", "cifs")
    if not (os.path.exists(inter_csv) and os.path.exists(xl_json)):
        log("StarPath export not found; no binding models will be linked")
        return pd.DataFrame()

    inter = pd.read_csv(inter_csv)
    known = set(node_ids)
    rows = {}
    for r in inter.itertuples(index=False):
        a, b = str(r.protA_alias), str(r.protB_alias)
        if a not in known or b not in known or a == b:
            continue
        rows[int(r.edge_id)] = {
            "starpath_id": int(r.edge_id),
            "gene_a": min(a, b), "gene_b": max(a, b),
            "n_crosslinks": int(r.crosslinks_number),
            "identification_score": float(r.identification_score),
            "local_fdr": float(r.local_fdr),
            "loc_a": r.protA_loc, "loc_b": r.protB_loc,
        }

    # residue-level positions: what was actually joined to what
    for rec in json.load(open(xl_json)):
        row = rows.get(int(rec.get("id", -1)))
        if row is None:
            continue
        row["crosslink_positions"] = json.dumps(
            [[c.get("pos_a"), c.get("pos_b")] for c in rec.get("crosslinks", [])])

    # the predicted complexes, and whether they explain the crosslinks
    if os.path.exists(sat_csv):
        for r in pd.read_csv(sat_csv).itertuples(index=False):
            row = rows.get(int(r.id))
            if row is None:
                continue
            row.update(frac_satisfied=float(r.frac_satisfied), min_ca=float(r.min_ca),
                       chai_iptm=float(r.chai_iptm), chai_plddt=float(r.chai_plddt))

    have_cifs = os.path.isdir(cif_dir)
    if have_cifs:
        by_id = {}
        for f in os.listdir(cif_dir):
            head = f.split("_")[0]
            if head.isdigit() and f.endswith(".cif"):
                by_id.setdefault(int(head), []).append(f)
        for sid, files in by_id.items():
            row = rows.get(sid)
            if row is not None:
                row["n_models"] = len(files)
                row["model_dir"] = os.path.relpath(cif_dir, base)
                row["model_files"] = ";".join(sorted(files))

    df = pd.DataFrame(list(rows.values()))
    if df.empty:
        return df
    for c in ("frac_satisfied", "min_ca", "chai_iptm", "chai_plddt", "n_models"):
        if c not in df.columns:
            df[c] = np.nan
    scored = df.frac_satisfied.notna().sum()
    # A model is only worth looking at if Chai-1 is confident about the interface AND the pose puts the
    # crosslinked residues in reach. Most are neither, and saying so is the point: the measurement stands
    # on the crosslink, not on the picture.
    df["model_trustworthy"] = (df.frac_satisfied >= 0.5) & (df.chai_iptm >= 0.5)
    log(f"crosslink models: {len(df):,} gene pairs; {int(df.n_models.notna().sum()):,} with local Chai-1 "
        f"structures; {scored:,} scored for crosslink satisfaction")
    if scored:
        s = df.dropna(subset=["frac_satisfied"])
        log(f"  {(s.frac_satisfied == 0).sum():,} ({(s.frac_satisfied == 0).mean():.0%}) explain no "
            f"crosslink at all; median interface ipTM {s.chai_iptm.median():.2f}; "
            f"{int(df.model_trustworthy.sum()):,} are both confident and crosslink-consistent "
            f"-- treat only those as showing how the binding happens")
    if not have_cifs:
        log(f"  model CIFs not on this machine ({cif_dir}); paths omitted")
    return df.sort_values(["gene_a", "gene_b"]).reset_index(drop=True)


def gene_attributes(edges: dict, models: pd.DataFrame, nodes: pd.DataFrame) -> None:
    """Per-gene counts the evidence panel can show, added in place."""
    n = len(nodes)
    for label in ("xlms", "struct", "ip_ms"):
        col = {"xlms": "n_xlink_partners", "struct": "n_struct_similar",
               "ip_ms": "n_ipms_partners"}[label]
        deg = np.zeros(n, dtype=int)
        if label in edges:
            a, b = edges[label][0], edges[label][1]
            np.add.at(deg, a, 1)
            np.add.at(deg, b, 1)
        nodes[col] = deg
    # best model agreement across a gene's crosslinked pairs: does any predicted complex explain them?
    best = pd.Series(np.nan, index=nodes.gene_id)
    if not models.empty and "frac_satisfied" in models:
        m = models.dropna(subset=["frac_satisfied"])
        if not m.empty:
            s = pd.concat([m[["gene_a", "frac_satisfied"]].rename(columns={"gene_a": "g"}),
                           m[["gene_b", "frac_satisfied"]].rename(columns={"gene_b": "g"})])
            best = s.groupby("g").frac_satisfied.max()
    nodes["best_model_agreement"] = nodes.gene_id.map(best).astype(float)
