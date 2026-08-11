#!/usr/bin/env python3
"""hyperLOPIT localisation: both inference methods, and orthoLOPIT across three species.

hyperLOPIT assigns a compartment to 3,827 of 8,140 *Toxoplasma* genes — 47%. The other 53% are not
"cytosolic by default"; assignment tracks protein abundance, so the unlabelled set is biased toward
low-abundance proteins, which are disproportionately the secreted effectors people care about. Raising
coverage therefore means transferring labels, and transfer must never be confused with measurement.

Three things this module keeps distinct:

* **MAP and MCMC are separate columns.** They are two inferences over the same experiment and they
  disagree for **980 of 3,827 genes (26%)**. Collapsing to one hides that; `lopit_methods_agree` exposes it.
* **orthoLOPIT is a transfer, not a measurement.** For a gene with no native call, the measured labels of
  its orthologs in *P. falciparum* and *C. parvum* are consulted. Compartment vocabularies differ per
  species (Toxo says "micronemes", Plasmodium "Micronemes", Crypto "Microneme"), so everything is mapped
  through `lopit_vocabulary_dictionary.csv` to one of 12 unified categories first. A label transfers only
  when every donor agrees; disagreement is left unlabelled rather than resolved by vote.
* **`compartment` stays measured.** `compartment_best` carries the fallback and `compartment_source`
  records which of the two any given label came from, so nothing downstream can silently mistake an
  inference for an observation.

*C. parvum* deserves a note: `MASTER_parasite_wide_by_orthogroup.csv` carries a `cpar_lopit_native`
column that is **empty for every row**, even though the measured Crypto hyperLOPIT exists on disk with
1,742 assignments. This module joins it from the source file instead of trusting that column.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

# Per-species measured hyperLOPIT: file, gene-id column, label column.
SPECIES = {
    "tgon": ("lopit_toxoplasma_gondii_ME49.csv", "gene_source_id", "MAP_location"),
    "pfal": ("lopit_plasmodium_falciparum_3D7.csv", "gene_source_id", "simplified"),
    "cpar": ("lopit_cryptosporidium_parvum_MEASURED_Guerin2023.csv", "gene_source_id",
             "MAP_location"),
}
# Donors for transfer. Trypanosoma is available in the orthogroup table but excluded: it is not an
# apicomplexan, and its compartment set does not include the apical secretory organelles that matter here.
DONORS = ("pfal", "cpar")

UNMAPPED = "(lineage-specific: unmapped)"


def _vocabulary(ds: str) -> dict:
    """(species, native term) -> unified category."""
    p = os.path.join(ds, "lopit_vocabulary_dictionary.csv")
    if not os.path.exists(p):
        return {}
    v = pd.read_csv(p)
    return {(str(s).strip(), str(t).strip().lower()): str(u).strip()
            for s, t, u in zip(v.species, v.native_term, v.unified_category)}


def _measured(ds: str, species: str, vocab: dict) -> pd.Series:
    """gene id -> unified category, for one species' measured hyperLOPIT."""
    fn, idc, labc = SPECIES[species]
    p = os.path.join(ds, fn)
    if not os.path.exists(p):
        return pd.Series(dtype=object)
    d = pd.read_csv(p, low_memory=False)
    if idc not in d.columns or labc not in d.columns:
        return pd.Series(dtype=object)
    lab = d[labc].astype(str).str.strip()
    lab = lab.where(~lab.str.lower().isin(["nan", "unknown", ""]))
    uni = lab.map(lambda x: vocab.get((species, str(x).strip().lower()))
                  if isinstance(x, str) else None)
    return pd.Series(uni.values, index=d[idc].astype(str)).dropna()


def lopit_labels(ds: str, nodes: pd.DataFrame, log=print) -> pd.DataFrame:
    """Attach measured hyperLOPIT (MAP + MCMC) and orthoLOPIT transfer to `nodes`, in place."""
    src = os.path.join(ds, SPECIES["tgon"][0])
    if os.path.exists(src):
        d = pd.read_csv(src, low_memory=False)
        cols = {"gene_source_id": "gene_id", "MAP_location": "lopit_map",
                "MCMC_location": "lopit_mcmc", "prob_MAP": "lopit_prob_map",
                "prob_MCMC": "lopit_prob_mcmc"}
        d = d[[c for c in cols if c in d.columns]].rename(columns=cols).drop_duplicates("gene_id")
        nodes = nodes.merge(d, on="gene_id", how="left")
        if {"lopit_map", "lopit_mcmc"} <= set(nodes.columns):
            both = nodes.lopit_map.notna() & nodes.lopit_mcmc.notna()
            nodes["lopit_methods_agree"] = np.where(
                both, (nodes.lopit_map == nodes.lopit_mcmc).astype(float), np.nan)
            n_ag = int(np.nansum(nodes.lopit_methods_agree))
            log(f"hyperLOPIT: {int(both.sum()):,} genes assigned by both methods; "
                f"they agree for {n_ag:,} ({n_ag / max(int(both.sum()), 1):.0%}) -- "
                f"the rest is real inference uncertainty, not noise")

    # `compartment` is the MAP call: what the experiment measured for this species.
    nodes["compartment"] = nodes.get("lopit_map")

    vocab = _vocabulary(ds)
    nodes["lopit_unified"] = [
        vocab.get(("tgon", str(x).strip().lower())) if isinstance(x, str) else None
        for x in nodes.compartment]

    # ---- orthoLOPIT: unanimous transfer from measured orthologs in the donor species
    master = os.path.join(ds, "MASTER_parasite_wide_by_orthogroup.csv")
    nodes["ortholopit_label"] = np.nan
    nodes["ortholopit_donors"] = np.nan
    if os.path.exists(master) and vocab:
        m = pd.read_csv(master, low_memory=False)
        m = m[m.tgon_gene_id.notna()]
        donor_maps = {s: _measured(ds, s, vocab) for s in DONORS}
        for s, series in donor_maps.items():
            log(f"orthoLOPIT: {len(series):,} measured labels available from {s}")
        label, who = {}, {}
        for row in m.itertuples(index=False):
            votes = {}
            for s in DONORS:
                gid = getattr(row, f"{s}_gene_id", None)
                if isinstance(gid, str):
                    v = donor_maps[s].get(gid.strip())
                    if isinstance(v, str) and v != UNMAPPED:
                        votes[s] = v
            if not votes:
                continue
            vals = set(votes.values())
            if len(vals) == 1:                       # unanimous, or a single donor
                label[row.tgon_gene_id] = vals.pop()
                who[row.tgon_gene_id] = "+".join(sorted(votes))
        nodes["ortholopit_label"] = nodes.gene_id.map(label)
        nodes["ortholopit_donors"] = nodes.gene_id.map(who)

    # ---- provenance, and the best available label
    measured = nodes.compartment.notna()
    transferred = ~measured & nodes.ortholopit_label.notna()
    nodes["compartment_source"] = np.where(measured, "hyperLOPIT",
                                           np.where(transferred, "orthoLOPIT", "unknown"))
    # unified space, so a transferred coarse category and a measured fine one are never mixed in one column
    nodes["compartment_best"] = nodes.lopit_unified.where(measured)
    nodes["compartment_best"] = nodes.compartment_best.fillna(
        nodes.ortholopit_label).fillna("unassigned")
    # hyperLOPIT assignment tracks abundance: a missing call is unknown, never a 27th compartment
    nodes["compartment"] = nodes.compartment.fillna("unassigned")
    log(f"localisation: {int(measured.sum()):,} measured + {int(transferred.sum()):,} orthoLOPIT "
        f"= {int(measured.sum() + transferred.sum()):,} of {len(nodes):,} genes labelled "
        f"({(measured.sum() + transferred.sum()) / len(nodes):.0%})")
    return nodes
