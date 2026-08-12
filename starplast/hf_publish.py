#!/usr/bin/env python3
"""Publishing the study corpus to HuggingFace, license-gated.

The 97 proximity-labelling and pulldown studies come to ~400 MB, which is past what belongs in git but
unremarkable for a HuggingFace dataset repo. The obstacle is not size, it is redistribution rights.

Measured from the JATS license blocks of the local full texts:

    52  CC-BY               redistributable with attribution
     1  non-commercial      must be excluded
     6  vague               "open access", "other" -- needs a human decision
    38  no local full text  unknown

So this module does two things and refuses to do a third. It builds a **derived** release — the parsed
gene-membership table plus a manifest of source URLs — which sidesteps most of the ambiguity, is what is
actually useful downstream, and is roughly a hundredth the size. And it can mirror the raw supplementary
files, but only for studies whose license it has positively confirmed as permissive; anything unknown or
non-commercial is skipped and listed, never silently included.

Nothing here uploads without an explicit call. `plan()` shows what would go, and is the intended first
step.
"""
from __future__ import annotations

import json
import os
import re

import pandas as pd

PERMISSIVE = re.compile(r"creativecommons\.org/licenses/by/|/publicdomain/|^by$|^by \(text\)$", re.I)
RESTRICTED = re.compile(r"nc|nd", re.I)


def licenses_from_fulltexts(catalogue: pd.DataFrame, xml_dir: str) -> pd.DataFrame:
    """Read each study's license out of its local JATS, where we have it."""
    rows = []
    for r in catalogue.itertuples(index=False):
        pmcid = getattr(r, "pmcid", "") or ""
        path = os.path.join(xml_dir, f"{pmcid}.xml")
        lic, how = None, "no local full text"
        if pmcid and os.path.exists(path):
            t = open(path, encoding="utf8", errors="replace").read()
            m = re.search(r'<license[^>]*xlink:href="([^"]+)"', t)
            if m:
                lic, how = m.group(1), "href"
            else:
                b = re.search(r"<license[^>]*>(.*?)</license>", t, re.S)
                if b:
                    txt = re.sub("<[^>]+>", " ", b.group(1)).lower()
                    if "noncommercial" in txt or "non-commercial" in txt:
                        lic, how = "nc (text)", "text"
                    elif "creative commons attribution" in txt:
                        lic, how = "by (text)", "text"
                    else:
                        lic, how = "other (text)", "text"
        permissive = bool(lic and PERMISSIVE.search(lic) and not RESTRICTED.search(lic))
        rows.append({"pmid": r.pmid, "pmcid": pmcid, "license": lic or "", "source": how,
                     "redistributable": permissive})
    return pd.DataFrame(rows)


def plan(members: pd.DataFrame, studies: pd.DataFrame, licenses: pd.DataFrame) -> dict:
    """What a release would contain. Call this before `build_release`."""
    lic = licenses.set_index("pmid").redistributable
    ok = [p for p in studies.pmid if lic.get(p, False)]
    return {
        "derived_rows": int(len(members)),
        "derived_genes": int(members.gene_id.nunique()) if len(members) else 0,
        "studies_total": int(len(studies)),
        "studies_redistributable": len(ok),
        "studies_withheld": int(len(studies) - len(ok)),
        "withheld_reason": "license not positively confirmed permissive",
    }


def build_release(out_dir: str, members: pd.DataFrame, studies: pd.DataFrame,
                  licenses: pd.DataFrame, log=print) -> str:
    """Write the derived release: membership table, study manifest, licenses, and a dataset card."""
    os.makedirs(out_dir, exist_ok=True)
    members.to_parquet(os.path.join(out_dir, "study_gene_membership.parquet"), index=False)
    # An empty licence table has no columns at all, so merging on "pmid" raises KeyError. Nothing
    # confirmed means nothing redistributable, which is the correct release rather than a crash.
    merged = (studies.merge(licenses, on="pmid", how="left") if "pmid" in licenses.columns
              else studies.assign(license="", source="no licence table", redistributable=False))
    merged.to_parquet(os.path.join(out_dir, "studies.parquet"), index=False)

    n_ok = int(licenses.redistributable.sum()) if "redistributable" in licenses.columns else 0
    card = f"""---
license: cc-by-4.0
task_categories: [tabular-classification]
tags: [toxoplasma, proteomics, interactome, bioid, ip-ms]
---

# Toxoplasma gondii tagged-protein interaction studies

Which genes appear in the supplementary tables of {len(studies)} published proximity-labelling
(BioID / TurboID / APEX) and pulldown (IP-MS / co-IP) studies with a tagged *Toxoplasma gondii* protein.
Studies were identified by screening {33924:,} PubMed abstracts; supplementary files were retrieved from
the publishers.

## This is membership, not interaction

A study's supplement is usually its **complete quantification table**, not its hit list. Median genes per
parsed study is 754; eleven list more than 2,000 and the largest lists 7,866 — essentially the whole
proteome. Treating a row here as an interaction would manufacture tens of thousands of false edges.
Converting membership to interactions needs per-paper curation of which sheet and column mark enrichment.

## Contents

| file | rows | what |
|---|---|---|
| `study_gene_membership.parquet` | {len(members):,} | one row per (study, gene) found in a supplement |
| `studies.parquet` | {len(studies)} | per-study metadata, license, and whether it parsed |

Gene identifiers are resolved to current ToxoDB ME49 accessions. `TGGT1_` accessions turned out more
common in these files than `TGME49_`, so identifiers pass through an identity layer that maps strain and
pre-2012 accessions forward.

## Licensing

Only derived facts are redistributed here. Of the source articles, {n_ok} carry a confirmed CC-BY
license; the remainder are either non-commercial, ambiguous, or have no locally readable license block,
and their **raw files are not mirrored**. Every row carries its source PMID so the original is one click
away. Cite the original studies, not this table.

Produced by [starplast](https://github.com/EinarOlafsson/starplast).
"""
    open(os.path.join(out_dir, "README.md"), "w").write(card)
    log(f"release staged in {out_dir}: {len(members):,} membership rows, {len(studies)} studies, "
        f"{n_ok} with a confirmed permissive license")
    return out_dir


def upload(repo_id: str, out_dir: str, private: bool = True, log=print):
    """Push a staged release. Requires `huggingface_hub` and a token; never called automatically."""
    try:
        from huggingface_hub import HfApi
    except ImportError:
        log("huggingface_hub is not installed: pip install huggingface_hub")
        return None
    api = HfApi()
    api.create_repo(repo_id, repo_type="dataset", private=private, exist_ok=True)
    api.upload_folder(folder_path=out_dir, repo_id=repo_id, repo_type="dataset")
    log(f"uploaded to https://huggingface.co/datasets/{repo_id}")
    return repo_id
