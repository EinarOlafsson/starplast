#!/usr/bin/env python3
"""Fetch the ToxoDB identity tables -> data/toxodb_identity.tsv and data/toxodb_strain_{gt1,veg}.tsv.

The literature does not cite genes by one identifier, so the identity layer needs more than symbols:

* ``gene_name``          -- the symbol (GRA16). ToxoDB has one for ~1,600 of 8,843 ME49 genes.
* ``gene_previous_ids``  -- pre-2012 accessions (TGME49_008830) that older papers still cite.
* GT1 / VEG accessions   -- the strain ids papers use interchangeably with ME49; they map to ME49 by
                            numeric suffix (verified at 99.5% orthogroup agreement, see identity.py).

Run once; the outputs are committed so the build works offline.
"""
import json
import os
import urllib.request

from . import paths

OUT = paths.data_dir()
URL = ("https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByTaxon"
       "/reports/attributesTabular")

# ToxoDB's tabular report ships display names as the header; map them to stable column names.
RENAME = {"Gene ID": "gene_id", "Gene Name or Symbol": "gene_name",
          "Previous ID(s)": "previous_ids", "Ortholog Group": "orthogroup",
          "Product Description": "product"}


def fetch(organism: str, attributes: list) -> str:
    """Retrieve a tabular attribute report from ToxoDB for one organism."""
    body = {"searchConfig": {"parameters": {"organism": json.dumps([organism])}},
            "reportConfig": {"attributes": attributes, "includeHeader": True,
                             "attachmentType": "plain"}}
    req = urllib.request.Request(
        URL, data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "starplast (research)"})
    return urllib.request.urlopen(req, timeout=900).read().decode("utf8", "replace")


def write(txt: str, path: str) -> int:
    """Write a fetched table to disk, creating the directory if needed."""
    lines = txt.splitlines()
    if lines:
        lines[0] = "\t".join(RENAME.get(c.strip(), c.strip()) for c in lines[0].split("\t"))
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"wrote {path}  ({len(lines) - 1} rows)")
    return len(lines) - 1


def main():
    """Fetch the ToxoDB identity and strain accession tables the build joins through."""
    os.makedirs(OUT, exist_ok=True)
    write(fetch("Toxoplasma gondii ME49",
                ["primary_key", "gene_name", "gene_previous_ids", "gene_product"]),
          os.path.join(OUT, "toxodb_identity.tsv"))
    for tag, org in (("gt1", "Toxoplasma gondii GT1"), ("veg", "Toxoplasma gondii VEG")):
        write(fetch(org, ["primary_key", "gene_name"]),
              os.path.join(OUT, f"toxodb_strain_{tag}.tsv"))


if __name__ == "__main__":
    main()
