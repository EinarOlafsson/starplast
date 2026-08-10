#!/usr/bin/env python3
"""Fetch the ToxoDB gene-name (symbol) table -> data/toxodb_gene_names.tsv.

Needed because the OrthoMCL product strings carry a usable symbol for only ~400 genes, while ToxoDB has
symbols for ~2,900. Without this table the literature layer under-counts massively: matching on products
alone found only 234 of 8,140 genes named anywhere in 33,924 abstracts.

Run once; the file is committed so the build works offline.
"""
import json
import os
import urllib.request

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "data", "toxodb_gene_names.tsv")
URL = ("https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByTaxon"
       "/reports/attributesTabular")


def main(organism="Toxoplasma gondii ME49"):
    body = {"searchConfig": {"parameters": {"organism": json.dumps([organism])}},
            "reportConfig": {"attributes": ["primary_key", "gene_name", "gene_product",
                                            "gene_source_id"],
                             "includeHeader": True, "attachmentType": "plain"}}
    req = urllib.request.Request(
        URL, data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "starplast (research)"})
    txt = urllib.request.urlopen(req, timeout=600).read().decode("utf8", "replace")
    with open(OUT, "w") as fh:
        fh.write(txt)
    print(f"wrote {OUT}  ({len(txt.splitlines()) - 1} genes)")


if __name__ == "__main__":
    main()
