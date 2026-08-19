#!/usr/bin/env python3
"""Fetch the VEuPathDB identity tables -> data/toxodb_*.tsv and data/plasmodb_identity.tsv.

The literature does not cite genes by one identifier, so the identity layer needs more than symbols:

* ``gene_name``          -- the symbol (GRA16). ToxoDB has one for ~1,600 of 8,843 ME49 genes.
* ``gene_previous_ids``  -- pre-2012 accessions (TGME49_008830) that older papers still cite.
* GT1 / VEG accessions   -- the strain ids papers use interchangeably with ME49; they map to ME49 by
                            numeric suffix (verified at 99.5% orthogroup agreement, see identity.py).

The Plasmodium arm needs the same thing and for the same reason. It had no identity layer at all
until 2026-08-18, which was survivable only while every Pf source happened to be keyed on current
`PF3D7_` accessions -- the first one that was not (a 2014 ribosome-profiling deposit, keyed on the
pre-2012 chromosome-based ids `PFE0630c` and `PF13_0222`) joined **zero** of its 3,629 genes. That
is the same failure the Toxoplasma layer was built for, met on the second arm.

Run once; the outputs are committed so the build works offline.
"""
import json
import os
import urllib.request

from . import paths

OUT = paths.data_dir()
URL = ("https://toxodb.org/toxo/service/record-types/transcript/searches/GenesByTaxon"
       "/reports/attributesTabular")

#: The same WDK service on the malaria site. One code path, two sites: the report, the parameter
#: shape and the header names are identical, so a second implementation could only differ by being
#: wrong on one arm.
PLASMODB_URL = ("https://plasmodb.org/plasmo/service/record-types/transcript/searches/GenesByTaxon"
                "/reports/attributesTabular")

#: The Plasmodium identity table, beside the Toxoplasma ones and never merged with them.
PLASMODB_IDENTITY = "plasmodb_identity.tsv"

# ToxoDB's tabular report ships display names as the header; map them to stable column names.
RENAME = {"Gene ID": "gene_id", "Gene Name or Symbol": "gene_name",
          "Previous ID(s)": "previous_ids", "Ortholog Group": "orthogroup",
          "Product Description": "product"}


#: VEuPathDB closed this endpoint to anonymous use, discovered 2026-08-19 when a call that had
#: worked hours earlier came back `401 Valid API Key required for this endpoint`. A guest session is
#: not enough -- the service issues one and then answers 403 -- so the key of a registered account is
#: now required, and it goes in this header.
API_KEY_ENV = "VEUPATHDB_API_KEY"
API_KEY_HEADER = "Auth-Key"

#: What to tell someone who has not set one. The URL is the site's own, so it stays right when the
#: page moves; the instruction is the part that is easy to get wrong.
API_KEY_HELP = (
    "VEuPathDB now requires an API key for its tabular reports. Register at toxodb.org or "
    "plasmodb.org, copy the key from Profile > Web Services Access, and set it in the environment:\n"
    f"    export {API_KEY_ENV}=<your key>\n"
    "The identity tables this fetches are COMMITTED, so a build works without it; the key is only "
    "needed to refresh them.")


def api_key() -> str:
    """The configured VEuPathDB key, or an empty string."""
    return os.environ.get(API_KEY_ENV, "").strip()


def fetch(organism: str, attributes: list, url: str = URL) -> str:
    """Retrieve a tabular attribute report from ToxoDB or PlasmoDB for one organizm.

    Raises `PermissionError` with instructions when no API key is configured, rather than letting
    the request come back 401 from three frames down with a message about an endpoint.
    """
    key = api_key()
    if not key:
        raise PermissionError(API_KEY_HELP)
    body = {"searchConfig": {"parameters": {"organism": json.dumps([organism])}},
            "reportConfig": {"attributes": attributes, "includeHeader": True,
                             "attachmentType": "plain"}}
    req = urllib.request.Request(
        url, data=json.dumps(body).encode(),
        headers={"Content-Type": "application/json", "User-Agent": "starplast (research)",
                 API_KEY_HEADER: key})
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
    write(fetch("Plasmodium falciparum 3D7",
                ["primary_key", "gene_name", "gene_previous_ids", "gene_product"],
                url=PLASMODB_URL),
          os.path.join(OUT, PLASMODB_IDENTITY))


if __name__ == "__main__":
    main()
