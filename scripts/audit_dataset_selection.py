"""Audit every registered source and search slot-scoped literature alternatives.

Queries, responses, citation snapshots and code/data hashes are retained with an
executed notebook. Results are proposals: abstract matches and storage coverage
do not admit a new dataset or establish that an incumbent is literature-best.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import platform
import sys
import time
import urllib.parse
import urllib.request

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from starplast import datasets as D, organisms as O, paths, slots as S  # noqa: E402
from starplast.dataset_selection import DatasetCandidate, SelectionPolicy, evaluate  # noqa: E402

SERVICE = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"


class PublicationClient:
    """Cache exact primary-service responses; failed queries remain visible and resumable."""

    def __init__(self, output: Path):
        self.cache = output / "requests"
        self.cache.mkdir(parents=True, exist_ok=True)

    def search(self, query: str, page_size=25, sort="") -> dict:
        """Fetch one core-metadata page, preserving URL, time, body and errors."""
        sorts = {"": "", "CITED desc": "sort_cited:y", "FIRST_PDATE desc": "sort_date:y"}
        if sort not in sorts:
            raise ValueError("Unsupported search ordering")
        params = {"query": query + (" " + sorts[sort] if sort else ""),
                  "format": "json", "resultType": "core", "pageSize": page_size}
        url = SERVICE + "?" + urllib.parse.urlencode(params)
        path = self.cache / (hashlib.sha256(url.encode()).hexdigest() + ".json")
        if path.exists():
            return json.loads(path.read_text())
        result = {"url": url, "retrieved_utc": datetime.now(timezone.utc).isoformat(),
                  "status": "error", "error": "", "response": {}}
        for attempt in range(3):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "starplast-dataset-audit/1"})
                with urllib.request.urlopen(request, timeout=30) as response:
                    result["response"] = json.load(response)
                if "hitCount" not in result["response"]:
                    raise ValueError("Search response has no hit count")
                result.update(status="ok", error="")
                break
            except Exception as exc:
                result["error"] = f"{type(exc).__name__}: {exc}"
                if attempt < 2:
                    time.sleep(attempt + 1)
        path.write_text(json.dumps(result, indent=2) + "\n")
        time.sleep(.25)
        return result


def _results(response):
    return response.get("response", {}).get("resultList", {}).get("result", [])


def _doi(dataset):
    # Only parse recorded identifiers; never invent a DOI from a guessed title.
    for text in (dataset.citation or "", dataset.accession or "", urllib.parse.unquote(dataset.url or "")):
        match = re.search(r"10\.\d{4,9}/[A-Za-z0-9._;/()+:-]+", text)
        if match:
            value = match.group().split("/MediaObjects", 1)[0].rstrip(".;)")
            value = re.sub(r"\.s\d+$", "", value)
            return value
    return ""


def _publication(record):
    types = record.get("pubTypeList", {}).get("pubType", [])
    corrections = record.get("commentCorrectionList", {}).get("commentCorrection", [])
    retracted = any("retracted publication" in str(t).lower() for t in types) or any(
        str(c.get("type", "")).lower() == "retraction in" for c in corrections)
    return {"publication_id": str(record.get("source", "")) + ":" + str(record.get("id", "")),
            "pmid": str(record.get("pmid") or (record.get("id") if record.get("source") == "MED" else "") or ""),
            "doi": record.get("doi", ""), "title": record.get("title", ""),
            "first_publication_date": record.get("firstPublicationDate", ""),
            "citations": record.get("citedByCount"), "retracted": retracted,
            "publication_types": types, "abstract": record.get("abstractText", "")}


def citation_factors(publication: dict, as_of: str) -> dict:
    """Calculate dated bibliometrics alone, without implying dataset suitability."""
    candidate = DatasetCandidate(publication["publication_id"], "bibliometrics_only",
                                 publication["publication_id"], publication["first_publication_date"],
                                 publication["citations"], 0, "Europe PMC", as_of, "accepted",
                                 "Publication metadata calculation only; not biological admission")
    result = evaluate(candidate, as_of)
    return {key: result[key] for key in ("age_years", "citations_per_year", "recency_factor")} | {
        "citation_metadata_status": result["status"], "citation_metadata_gap": result["gap"],
        "citation_provider": "Europe PMC", "citation_snapshot": as_of,
        "comprehensiveness": None, "selection_score": None}


def query_for(slot) -> str:
    """Construct a primary-literature query for a declared organism/context/question.

    Queries intentionally retain different host cell types. Their hits still need
    sample/deposit review; publications naming a stage are not proof of an assay
    in that stage. First pages by citations and date are both retained.
    """
    from generate_slot_table import ASSAY_TERMS

    subject = ""
    if slot.unit == "host_gene":
        host = 'TITLE_ABS:"Anopheles"' if "Anopheles" in slot.name else \
               '(TITLE_ABS:"human" OR TITLE_ABS:"Homo sapiens")' if "human" in slot.name else \
               '(TITLE_ABS:"mouse" OR TITLE_ABS:"Mus musculus")' if "mouse" in slot.name else ""
        question = slot.name.split(" · ", 1)[0]
        assays = {
            "host transcriptome": "transcriptome OR RNA-seq OR expression atlas OR CAGE OR FANTOM5 OR GTEx",
            "host proteome": "proteome OR proteomics OR mass spectrometry",
            "host surface / receptor repertoire": "surfaceome OR cell-surface OR surface proteome",
            "host response to infection": "infection OR host response OR transcriptome",
            "host gene requirement": "CRISPR OR knockout OR genetic screen",
        }
        assay = next((v for k, v in assays.items() if question.startswith(k)), "host OR proteome OR CRISPR")
        tissue = slot.name.split(" · ", 1)[-1].replace("human ", "").replace("mouse ", "").replace("Anopheles ", "")
        synonyms = {"bone-marrow macrophage": ('"bone marrow-derived macrophage"', 'BMDM'),
                    "erythrocyte": ('"red blood cell"', 'RBC'), "hepatocyte": ('hepatocytes',),
                    "fibroblast": ('fibroblasts', 'HFF')}
        tissue_terms = [f'"{tissue}"', *synonyms.get(tissue, ()), '"tissue atlas"', 'FANTOM5', 'GTEx']
        context = 'TITLE_ABS:(' + ' OR '.join(tissue_terms) + ')' if " · " in slot.name else ""
        organism = host or f'TITLE_ABS:"{O.get(slot.organism).species}"'
        if "infection" in question or "gene requirement" in question or "vacuole" in question:
            organism += f' AND TITLE_ABS:"{O.get(slot.organism).species}"'
    else:
        organism = f'TITLE_ABS:"{O.get(slot.organism).species}"'
        assay = ASSAY_TERMS.get(slot.axis, "").strip("()")
        question = slot.name.split(" · ", 1)[0]
        if not assay:
            assay = " OR ".join('"' + word + '"' for word in re.findall(r"[A-Za-z][A-Za-z-]{3,}", question))
        elif question not in {slot.axis, "transcription", "translation", "protein abundance", "fitness"}:
            generic = {"protein", "proteins", "gene", "genes", "measured", "predicted", "derived",
                       "transferred", "sequence", "presence", "counts", "level", "levels", "from", "with"}
            words = [word for word in re.findall(r"[A-Za-z][A-Za-z-]{3,}", question) if word.lower() not in generic]
            if words:
                subject = "TITLE_ABS:(" + " OR ".join('"' + word + '"' for word in words) + ")"
        context = ""
        if slot.context in O.get(slot.organism).contexts:
            context = f'TITLE_ABS:"{slot.context}"'
    # Quoted multiword assay alternatives cannot dissolve into independent words.
    assay = " OR ".join('"' + term.strip().strip('"') + '"' for term in assay.split(" OR ") if term.strip())
    return " AND ".join(part for part in (organism, "TITLE_ABS:(" + assay + ")" if assay else "", subject, context) if part)


def comparison_factors(candidate: DatasetCandidate, as_of: str) -> dict:
    """Expose metadata diagnostics without admitting unreviewed slot outputs.

    A slot can contain several quantities and boolean detection encodings. Its
    regex groups and filled storage values cannot establish biological scope or
    measured assay coverage. Diagnostic scores must never become preferences.
    """
    result = evaluate(candidate, as_of)
    return {**result, "metadata_status": result["status"], "diagnostic_score": result["score"],
            "score": None, "admission": "pending", "status": "pending_scope_review",
            "admission_reason": "Quantity, context, source lineage and assay comprehensiveness require review"}


def audit(output: Path, as_of: str, discover=True) -> dict:
    """Resolve all source identities, compare incumbents and retain literature search gaps."""
    output.mkdir(parents=True, exist_ok=True)
    client = PublicationClient(output)
    catalog = S.all_slots()
    inputs = [ROOT / "starplast" / name for name in
              ("datasets.py", "slots.py", "dataset_selection.py", "data/slots.json")]
    inputs += [Path(__file__), ROOT / "scripts/generate_slot_table.py"]
    tables = {code: pd.read_parquet(O.nodes_path(code)) for code in O.codes()}
    hosts = {code: pd.read_parquet(paths.cache_file(name)) for code, name in O.HOST_TABLES.items()}
    inputs += [Path(O.nodes_path(code)) for code in tables]
    inputs += [Path(paths.cache_file(name)) for name in O.HOST_TABLES.values()]
    hashes = {str(path): hashlib.sha256(path.read_bytes()).hexdigest() for path in inputs}
    pmids = sorted({str(d.pmid) for d in D.REGISTRY if d.pmid} |
                   {str(c[0]) for s in catalog for c in s.candidates if c and str(c[0]).isdigit()})
    publications, id_queries = {}, {}
    batches = [pmids[i:i + 30] for i in range(0, len(pmids), 30)]
    def retrieve(batch):
        query = "(" + " OR ".join("EXT_ID:" + value for value in batch) + ") AND SRC:MED"
        response = client.search(query, page_size=1000)
        return batch, response
    with ThreadPoolExecutor(max_workers=2) as pool:
        for batch, response in pool.map(retrieve, batches):
            for pmid in batch:
                id_queries[pmid] = response
            for record in _results(response):
                publication = _publication(record)
                if publication["pmid"] in batch:
                    publications[publication["pmid"]] = publication
    source_rows = []
    for dataset in D.REGISTRY:
        publication, status, url = {}, "unresolved_publication", ""
        if dataset.pmid:
            publication = publications.get(str(dataset.pmid), {})
            response = id_queries[str(dataset.pmid)]
            status = "resolved_registry_pmid" if publication else "identifier_not_found" if response["status"] == "ok" else "lookup_failed"
            url = response["url"]
        elif _doi(dataset):
            doi = _doi(dataset)
            response = client.search('DOI:"' + doi + '"', page_size=10)
            matches = [r for r in _results(response) if str(r.get("doi", "")).lower() == doi.lower()]
            if len(matches) == 1:
                publication = _publication(matches[0])
                status = "resolved_recorded_doi"
            else:
                status = "doi_ambiguous" if matches else "identifier_not_found" if response["status"] == "ok" else "lookup_failed"
            url = response["url"]
        elif dataset.accession and re.search(r"\bGSE\d+\b", dataset.accession):
            accession = re.search(r"\bGSE\d+\b", dataset.accession).group()
            response = client.search('"' + accession + '"', page_size=10)
            matches = _results(response)
            # A text hit is a candidate association, not a verified originating paper.
            status = "accession_publications_need_review" if matches else "accession_publication_not_found"
            url = response["url"]
        elif dataset.derived_from or "comput" in dataset.kind.lower() or "COMPUTED" in dataset.name:
            status = "derived_or_local_no_own_publication"
        elif dataset.level == "reference":
            status = "reference_resource_publication_unresolved"
        row = {"dataset_id": dataset.key, "organism": dataset.organism, "name": dataset.name,
               "kind": dataset.kind, "columns": json.dumps(dataset.columns), "registry_pmid": dataset.pmid,
               "registry_accession": dataset.accession, "publication_status": status,
               "source_link_status": "requires_source_review" if not dataset.citation or "CONFIRM" in (dataset.citation or "")
                                      or "not verification" in dataset.note else "registry_association_not_independently_revalidated",
               "query_url": url, "citation_provider": "Europe PMC", "citation_snapshot": as_of,
               "publication_id": "", "title": "", "first_publication_date": "", "citations": None,
               "citations_per_year": None, "retracted": False,
               "metadata_status": "unresolved", "metadata_gap": "Publication association unavailable"}
        row.update({k: v for k, v in publication.items() if k not in {"abstract", "publication_types"}})
        if publication:
            candidate = DatasetCandidate(dataset.key, "publication_metadata_only", publication.get("publication_id", ""),
                                         publication.get("first_publication_date", ""), publication.get("citations"),
                                         0, "Europe PMC", as_of, "accepted", "Metadata calculation only, not dataset admission")
            calculated = evaluate(candidate, as_of)
            row["citations_per_year"] = calculated["citations_per_year"]
            row["metadata_status"] = calculated["status"]
            row["metadata_gap"] = calculated["gap"]
        source_rows.append(row)
    source_table = pd.DataFrame(source_rows)
    source_table.to_csv(output / "source_publications.csv", index=False)
    by_source = {r["dataset_id"]: r for r in source_rows}
    comparisons, query_rows, literature = [], [], []
    unique_queries = {s.key: query_for(s) for s in catalog}
    search_results = {}
    if discover:
        queries = sorted(set(unique_queries.values()))
        jobs = [(query, order) for query in queries for order in ("CITED desc", "FIRST_PDATE desc")]
        def search(job):
            query, order = job
            return job, client.search(query, page_size=10, sort=order)
        with ThreadPoolExecutor(max_workers=2) as pool:
            for i, (job, result) in enumerate(pool.map(search, jobs), 1):
                search_results[job] = result
                if i % 40 == 0:
                    print(f"Literature queries: {i}/{len(jobs)}", flush=True)
    for slot in catalog:
        table = tables.get(slot.organism)
        if slot.unit == "host_gene":
            table = S.unit_table(slot, {"host_gene": hosts})
        groups = S._groups(table, slot) if table is not None and slot.unit in {"gene", "host_gene"} else []
        candidates = []
        for pattern, columns in groups:
            source_columns = {}
            for column in columns:
                # Host slot scopes belong to a parasite UI arm; infer no species
                # from an identifier. Match the explicit host table object.
                code = next((code for code, frame in hosts.items() if table is frame), slot.organism)
                source = D.provenance(column, organism=code)
                if source:
                    source_columns.setdefault(source.key, []).append(column)
            for key, declared in source_columns.items():
                metadata = by_source[key]
                coverage = int(table[declared].notna().any(axis=1).sum()) / len(table) if len(table) else None
                candidate = DatasetCandidate(key + ":" + pattern, slot.key, metadata["publication_id"],
                                             metadata["first_publication_date"], metadata["citations"], coverage,
                                             "Europe PMC", as_of, "accepted",
                                             "Provisional comparison of registered slot outputs; assay extent needs source review")
                result = comparison_factors(candidate, as_of)
                candidates.append({**result, "source_id": key, "pattern": pattern,
                                   "coverage_basis": "installed stored-value fraction; assay completeness unverified"})
        scored = sorted([r for r in candidates if r["diagnostic_score"] is not None and not by_source[r["source_id"]]["retracted"]],
                        key=lambda r: r["diagnostic_score"], reverse=True)
        default_patterns = [groups[0][0]] if slot.policy == "one" and groups else []
        provisional = scored[0] if scored else {}
        comparisons.append({"slot": slot.key, "organism": slot.organism, "name": slot.name,
                            "context": slot.context, "unit": slot.unit, "policy": slot.policy,
                            "default_patterns": json.dumps(default_patterns), "candidate_groups": len(groups),
                            "rankable_groups": 0, "provisional_preference": "", "provisional_pattern": "",
                            "change_proposed": False, "diagnostic_groups": len(scored),
                            "diagnostic_source": provisional.get("source_id", ""),
                            "diagnostic_pattern": provisional.get("pattern", ""),
                            "diagnostic_default_change": bool(default_patterns and provisional and provisional["pattern"] not in default_patterns),
                            "decision": "source/assay review required; no automatic promotion",
                            "candidate_factors": json.dumps(candidates, default=str)})
        query = unique_queries[slot.key]
        for order in ("CITED desc", "FIRST_PDATE desc"):
            response = search_results.get((query, order), {})
            found = _results(response)
            query_rows.append({"slot": slot.key, "query": query, "order": order, "query_url": response.get("url", ""),
                               "status": response.get("status", "not_run"), "hit_count": response.get("response", {}).get("hitCount"),
                               "retained": len(found), "truncated": bool(response.get("response", {}).get("hitCount", 0) > len(found)),
                               "error": response.get("error", "")})
            for record in found:
                publication = _publication(record)
                literature.append({"slot": slot.key, "order": order, "query_url": response["url"],
                                   **publication, **citation_factors(publication, as_of),
                                   "admission": "pending_sample_and_deposit_review"})
    pd.DataFrame(comparisons).to_csv(output / "slot_comparisons.csv", index=False)
    pd.DataFrame(query_rows).to_csv(output / "slot_searches.csv", index=False)
    pd.DataFrame(literature).to_json(output / "literature_candidates.json", orient="records", indent=2)
    if len(source_table) != len(D.REGISTRY) or set(source_table.dataset_id) != {d.key for d in D.REGISTRY}:
        raise ValueError("Source audit does not reconcile with the registry")
    if len(comparisons) != len(catalog):
        raise ValueError("Slot comparisons do not reconcile")
    for path, digest in hashes.items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
            raise ValueError("Audit input changed: " + path)
    summary = {"source_count": len(source_table), "publication_statuses": source_table.publication_status.value_counts().to_dict(),
               "slots": len(catalog), "unique_search_queries": len(set(unique_queries.values())),
               "search_jobs": len(search_results), "search_failures": sum(r.get("status") != "ok" for r in search_results.values()),
               "literature_rows": len(literature), "provisional_changes": sum(r["change_proposed"] for r in comparisons),
               "diagnostic_flags": sum(r["diagnostic_default_change"] for r in comparisons),
               "promoted_sources": 0}
    used_urls = {r["url"] for r in id_queries.values()} | {r["query_url"] for r in source_rows if r["query_url"]} | \
                {r["url"] for r in search_results.values()}
    request_hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in client.cache.glob("*.json")
                      if json.loads(p.read_text())["url"] in used_urls}
    manifest = {"as_of": as_of, "retrieved_utc": datetime.now(timezone.utc).isoformat(), "policy": SelectionPolicy().to_dict(),
                "input_sha256": hashes, "summary": summary,
                "software": {"python": platform.python_version(), "pandas": pd.__version__},
                "used_request_sha256": request_hashes,
                "limits": "First ten highest-cited and newest hits; no exhaustive-literature or biological-best assertion; storage coverage proxy",
                "outputs": {name: hashlib.sha256((output / name).read_bytes()).hexdigest() for name in
                            ("source_publications.csv", "slot_comparisons.csv", "slot_searches.csv", "literature_candidates.json")}}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2), flush=True)
    return summary


def main():
    """Execute the dated primary-service audit in a resumable annotated notebook."""
    from notebook_runner import ExecutedNotebook

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--as-of", required=True)
    parser.add_argument("--skip-discovery", action="store_true")
    args = parser.parse_args()
    notebook = ExecutedNotebook("Dataset citation preferences, source identities and slot literature search")
    notebook.md("# Dataset selection audit", "All registry entries and slots remain visible. "
                "Bibliometrics express the user's preference; eligible biology and source validation precede replacement. "
                "Primary-service responses and query URLs are retained. Remote citation coverage is provider-specific.")
    notebook.code("from pathlib import Path", "import sys", f"sys.path.insert(0, {str(ROOT / 'scripts')!r})",
                  "from audit_dataset_selection import audit",
                  f"summary = audit(Path({str(args.out.resolve())!r}), {args.as_of!r}, discover={not args.skip_discovery!r})",
                  "summary")
    notebook.md("## Limits", "A registry PMID is an asserted source link, not independently revalidated assay lineage. "
                "Accession text matches need origin-paper review. Searches retain first pages by citations and newest date; "
                "large result sets remain explicitly truncated. Abstract candidates are not admitted sources. "
                "Incumbent diagnostic scores use stored coverage, not verified assay completeness; all admissions remain pending. "
                "No dataset, calibration or strategy output is replaced by this audit.")
    notebook.write(str(args.out / "audit.ipynb"))


if __name__ == "__main__":
    main()
