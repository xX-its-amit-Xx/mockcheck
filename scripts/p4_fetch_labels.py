#!/usr/bin/env python3
"""Phase 4 step 1: fetch GenBank records for the label-curation set (public INSDC data, read-only, polite).

Viruses are non-select-agent and the task is classification of existing public sequences only.
For each virus we query region strata (six world regions) and host strata, take a deterministic window of ids,
fetch GenBank flat files, and keep the raw qualifier strings so nothing is lost before curation.

  python scripts/p4_fetch_labels.py [OUTDIR]     # default C:/Temp/mc/p4 ; raw responses cached, resumable
Writes OUTDIR/records.tsv (one row per record, protein translation + raw source qualifiers)."""
import csv, hashlib, json, os, re, sys, time, urllib.parse, urllib.request

OUT = sys.argv[1] if len(sys.argv) > 1 else "C:/Temp/mc/p4"
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
UA = "mockcheck-research/0.1 (open benchmark; read-only metadata and sequence queries)"
PER_QUERY = 120
BATCH = 50

REGIONS = {
    "africa": "Kenya Uganda Nigeria Ghana South_Africa Egypt Ethiopia Tanzania Cameroon Senegal Mali Zambia Zimbabwe Mozambique Madagascar Tunisia Morocco Algeria Burkina_Faso Malawi Rwanda Niger Angola",
    "asia": "China Japan South_Korea India Thailand Vietnam Bangladesh Indonesia Cambodia Taiwan Mongolia Pakistan Iran Israel Malaysia Philippines Nepal Kazakhstan Laos Myanmar",
    "europe": "Germany France Italy Spain United_Kingdom Netherlands Poland Belgium Sweden Denmark Hungary Czech_Republic Austria Switzerland Portugal Greece Finland Norway Ireland Serbia Russia",
    "north_america": "USA Canada Mexico Guatemala Cuba Costa_Rica Panama",
    "south_america": "Brazil Argentina Peru Chile Colombia Venezuela Bolivia Ecuador Uruguay Paraguay",
    "oceania": "Australia New_Zealand Papua_New_Guinea Fiji",
}

# virus: organism term, title/gene term, length filter, product regex for the CDS to keep, host strata queries
VIRUSES = {
    "influenza_a": dict(org='"Influenza A virus"[Organism]', term="hemagglutinin[Title]", slen="1500:2400", prod=r"hemagglutinin", per=250,
                        hosts=["Homo sapiens", "Sus scrofa", "Gallus gallus", "Anas platyrhynchos", "Anatidae", "Equus caballus", "Canis lupus familiaris", "Aves"]),
    "rabies": dict(org='"Rabies lyssavirus"[Organism]', term="nucleoprotein[Title]", slen="1200:1600", prod=r"nucleo",
                   hosts=["Canis lupus familiaris", "Chiroptera", "Vulpes vulpes", "Homo sapiens", "Bos taurus", "Procyon lotor"]),
    "rotavirus_a": dict(org='"Rotavirus A"[Organism]', term="VP7[Title]", slen="900:1100", prod=r"VP7|outer capsid|glycoprotein",
                        hosts=["Homo sapiens", "Bos taurus", "Sus scrofa", "Gallus gallus", "Equus caballus", "Felis catus"]),
}

last = [0.0]


def http(url, data=None):
    wait = 0.45 - (time.time() - last[0])
    if wait > 0:
        time.sleep(wait)
    last[0] = time.time()
    req = urllib.request.Request(url, data=data, headers={"User-Agent": UA})
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            if attempt == 3:
                raise
            time.sleep(2 + 3 * attempt)


def esearch(term, retstart=0, retmax=0):
    q = urllib.parse.urlencode({"db": "nuccore", "retmode": "json", "term": term, "retstart": retstart, "retmax": retmax})
    return json.loads(http(f"{E}/esearch.fcgi?{q}"))["esearchresult"]


def ids_for(term, per=PER_QUERY):
    c = int(esearch(term)["count"])
    if c == 0:
        return c, []
    start = 0 if c <= per else int(hashlib.sha1(term.encode()).hexdigest(), 16) % (c - per)   # deterministic window
    return c, esearch(term, start, per)["idlist"]


def parse_gb(rec):
    acc = re.search(r"^VERSION\s+(\S+)", rec, re.M)
    if not acc:
        return None
    locus = re.search(r"^LOCUS\s+\S+\s+(\d+) bp", rec, re.M)
    dfn = re.search(r"^DEFINITION\s+(.*?)^ACCESSION", rec, re.M | re.S)
    feats = rec.split("\nFEATURES", 1)[1].split("\nORIGIN", 1)[0] if "\nFEATURES" in rec else ""
    items, cur = [], None
    for line in feats.split("\n")[1:]:
        if re.match(r"^ {5}\S", line):
            cur = [line.split()[0], {}]; items.append(cur); lastq = None
        elif cur is not None and line.startswith(" " * 21):
            t = line.strip()
            if t.startswith("/"):
                k, _, v = t[1:].partition("=")
                lastq = k; cur[1][k] = cur[1].get(k, "") + v.strip('"')
            elif lastq:
                cur[1][lastq] += ("" if lastq == "translation" else " ") + t.strip('"')
    src = next((q for k, q in items if k == "source"), {})
    cds = [q for k, q in items if k == "CDS"]
    return dict(accession=acc.group(1), nt_len=int(locus.group(1)) if locus else 0, definition=re.sub(r"\s+", " ", dfn.group(1)).strip() if dfn else "", src=src, cds=cds)


def main():
    os.makedirs(f"{OUT}/raw", exist_ok=True)
    rows = []
    for vname, v in VIRUSES.items():
        qs = [(f"region:{r}", f"{v['org']} AND {v['term']} AND {v['slen']}[SLEN] AND (" + " OR ".join(f"{c.replace('_', ' ')}[Country]" for c in cs.split()) + ")") for r, cs in REGIONS.items()]
        qs += [(f"host:{h}", f"{v['org']} AND {v['term']} AND {v['slen']}[SLEN] AND \"{h}\"[Host]") for h in v["hosts"]]
        ids = {}
        for tag, term in qs:
            cache = f"{OUT}/raw/ids_{vname}_{re.sub(r'[^A-Za-z0-9]+', '_', tag)}.json"
            if os.path.exists(cache):
                c, got = json.load(open(cache))
            else:
                c, got = ids_for(term, v.get("per", PER_QUERY))
                json.dump([c, got], open(cache, "w"))
            for i in got:
                ids.setdefault(i, tag)
            print(f"{vname:12s} {tag:32s} hits={c:7d} taken={len(got)}", flush=True)
        idl = sorted(ids)
        for b in range(0, len(idl), BATCH):
            chunk = idl[b:b + BATCH]
            cache = f"{OUT}/raw/gb_{vname}_{b:05d}.txt"
            if not os.path.exists(cache):
                body = urllib.parse.urlencode({"db": "nuccore", "id": ",".join(chunk), "rettype": "gb", "retmode": "text"}).encode()
                open(cache, "w", encoding="utf-8", newline="\n").write(http(f"{E}/efetch.fcgi", body))
            for rec in open(cache, encoding="utf-8").read().split("\n//\n"):
                p = parse_gb(rec)
                if not p:
                    continue
                cds = [c for c in p["cds"] if re.search(v["prod"], c.get("product", ""), re.I) and c.get("translation")]
                if not cds:
                    continue
                c = max(cds, key=lambda x: len(x["translation"]))
                s = p["src"]
                rows.append(dict(accession=p["accession"], virus=vname, protein=c["translation"], protein_len=len(c["translation"]), nt_len=p["nt_len"],
                                 strain=s.get("strain", ""), isolate=s.get("isolate", ""), host_raw=s.get("host", ""), country_raw=s.get("geo_loc_name", "") or s.get("country", ""),
                                 collection_date=s.get("collection_date", ""), isolation_source=s.get("isolation_source", ""), lab_host=s.get("lab_host", ""),
                                 serotype=s.get("serotype", "") or s.get("serovar", ""), lat_lon=s.get("lat_lon", ""), note=s.get("note", ""),
                                 definition=p["definition"][:200], query_stratum=ids.get(chunk[0], "")))
            print(f"  {vname} fetched {min(b + BATCH, len(idl))}/{len(idl)}", flush=True)
    with open(f"{OUT}/records.tsv", "w", newline="\n", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print("records", len(rows))


if __name__ == "__main__":
    main()
