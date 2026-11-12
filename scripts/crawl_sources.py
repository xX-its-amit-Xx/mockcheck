#!/usr/bin/env python3
"""Inventory candidate ground-truth sources from public APIs: Zenodo, Europe PMC (literature + preprints), GitHub.
Read-only, polite (1 request/second, small page sizes, descriptive User-Agent), no scraping of sites whose terms forbid it.
Not crawled on purpose: LinkedIn (user agreement forbids automated scraping), Google Patents (no API; use USPTO PatentsView with a key).
Writes results/source_inventory.tsv and results/source_queries.tsv. Re-runnable; every row keeps its URL and date for provenance."""
import csv, json, os, sys, time, urllib.parse, urllib.request

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UA = "mockcheck-research/0.1 (open benchmark project; read-only metadata queries)"

ZENODO = ["\"mock community\" metagenomic", "viral host prediction dataset", "SARS-CoV-2 deep mutational scanning", "pathogen metagenomics benchmark", "wastewater metagenomics sequencing"]
EPMC = ["mock community benchmark metagenomic classifier", "viral host range prediction benchmark", "pathogen detection metagenomic clinical validation", "deep mutational scanning viral escape prediction benchmark"]
GITHUB = ["metagenomic classifier benchmark", "viral host prediction", "protein language model viral fitness benchmark", "mock community kraken2"]


def get(url):
    time.sleep(1.0)
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=40) as r:
        return json.load(r)


def zenodo(q):
    d = get("https://zenodo.org/api/records?" + urllib.parse.urlencode({"q": q, "size": 8, "sort": "bestmatch", "type": "dataset"}))
    rows = [dict(source="zenodo", title=h["metadata"]["title"][:140], date=h["metadata"].get("publication_date", ""),
                 license=(h["metadata"].get("license") or {}).get("id", ""), url=h["links"]["self_html"]) for h in d["hits"]["hits"]]
    return d["hits"]["total"], rows


def epmc(q):
    d = get("https://www.ebi.ac.uk/europepmc/webservices/rest/search?" + urllib.parse.urlencode({"query": q, "format": "json", "pageSize": 8, "sort": "CITED desc"}))
    rows = [dict(source="europepmc", title=r.get("title", "")[:140], date=r.get("firstPublicationDate", ""), license=r.get("license", ""),
                 url="https://europepmc.org/article/%s/%s" % (r.get("source", "MED"), r.get("id", ""))) for r in d["resultList"]["result"]]
    return d["hitCount"], rows


def github(q):
    # unauthenticated search API allows ~10 requests/minute: space calls out
    time.sleep(7.0)
    d = get("https://api.github.com/search/repositories?" + urllib.parse.urlencode({"q": q, "sort": "stars", "per_page": 8}))
    rows = [dict(source="github", title=(r["full_name"] + ": " + (r.get("description") or ""))[:140], date=(r.get("pushed_at") or "")[:10],
                 license=(r.get("license") or {}).get("spdx_id") or "", url=r["html_url"]) for r in d["items"]]
    return d["total_count"], rows


def main():
    inv, qs = [], []
    for fn, queries in ((zenodo, ZENODO), (epmc, EPMC), (github, GITHUB)):
        for q in queries:
            try:
                total, rows = fn(q)
            except Exception as e:  # network or rate limit: record and move on
                qs.append(dict(source=fn.__name__, query=q, total_hits="error", note=str(e)[:80])); continue
            qs.append(dict(source=fn.__name__, query=q, total_hits=total, note=""))
            for r in rows:
                r["query"] = q; inv.append(r)
    for name, rows in (("source_inventory.tsv", inv), ("source_queries.tsv", qs)):
        if rows:
            with open(f"{ROOT}/results/{name}", "w", newline="\n", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    print(f"{len(inv)} sources from {len(qs)} queries")


if __name__ == "__main__":
    main()
