#!/usr/bin/env python3
"""Phase 4 step 2: turn raw GenBank records into a curated ground-truth set with evidence tiers and holdouts.

Input : OUT/records.tsv from p4_fetch_labels.py
Output: data/groundtruth/p4_groundtruth.tsv.gz (one row per sequence, all labels, tiers, split)
        data/groundtruth/p4_audit_sample.tsv (60 rows with evidence strings, for human review)
        results/p4_curation_summary.tsv (what was kept or dropped, and why)

Label rule (the point of this script): a label enters the GOLD set only when an INDEPENDENT field agrees with the
submitter's field. Example: host from the /host qualifier must agree with the host encoded in the strain name
(A/swine/Iowa/..., RVA/Human-wt/USA/...). Labels with no independent evidence are kept as WEAK (train only).
Labels whose independent evidence CONFLICTS are dropped from the set and counted.
Missing categories are left empty, never guessed."""
import csv, gzip, hashlib, os, re, sys, collections, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from p4_iso3 import ISO3  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = sys.argv[1] if len(sys.argv) > 1 else "C:/Temp/mc/p4"
GT = f"{ROOT}/data/groundtruth"
os.makedirs(GT, exist_ok=True)

# ---------------------------------------------------------------- geography tables
REG = {
    "africa": "kenya uganda nigeria ghana south africa egypt ethiopia tanzania cameroon senegal mali zambia zimbabwe mozambique madagascar tunisia morocco algeria burkina faso malawi rwanda niger angola sudan south sudan libya gabon congo democratic republic of the congo benin togo namibia botswana chad burundi sierra leone liberia guinea mauritius gambia somalia eswatini lesotho cote d'ivoire",
    "asia": "china japan south korea india thailand vietnam viet nam bangladesh indonesia cambodia taiwan mongolia pakistan iran israel malaysia philippines nepal kazakhstan laos myanmar hong kong singapore sri lanka saudi arabia turkey iraq lebanon jordan kuwait oman qatar united arab emirates bhutan uzbekistan afghanistan yemen syria",
    "europe": "germany france italy spain united kingdom netherlands poland belgium sweden denmark hungary czech republic austria switzerland portugal greece finland norway ireland serbia russia romania bulgaria croatia slovenia slovakia ukraine belarus lithuania latvia estonia luxembourg iceland cyprus malta bosnia and herzegovina albania moldova north macedonia montenegro",
    "north_america": "usa canada mexico guatemala cuba costa rica panama honduras el salvador nicaragua haiti dominican republic puerto rico jamaica trinidad and tobago belize",
    "south_america": "brazil argentina peru chile colombia venezuela bolivia ecuador uruguay paraguay french guiana suriname guyana",
    "oceania": "australia new zealand papua new guinea fiji new caledonia samoa tonga vanuatu",
}
ALIAS = {"united states": "usa", "us": "usa", "viet nam": "vietnam", "hong kong": "china", "england": "united kingdom", "scotland": "united kingdom", "wales": "united kingdom",
         "northern ireland": "united kingdom", "uk": "united kingdom", "korea": "south korea", "republic of korea": "south korea", "russian federation": "russia",
         "democratic republic of the congo": "congo drc", "drc": "congo drc", "ivory coast": "cote d'ivoire", "czechia": "czech republic", "taiwan": "taiwan"}
COUNTRY_REGION = {}
for r, names in REG.items():
    # greedy match of multiword names: build from known multiword list
    pass
MULTI = ["south africa", "burkina faso", "south sudan", "sierra leone", "democratic republic of the congo", "south korea", "viet nam", "hong kong", "sri lanka", "saudi arabia",
         "united arab emirates", "united kingdom", "czech republic", "bosnia and herzegovina", "north macedonia", "costa rica", "el salvador", "dominican republic", "puerto rico",
         "trinidad and tobago", "french guiana", "new zealand", "papua new guinea", "new caledonia", "cote d'ivoire"]
for r, names in REG.items():
    rest = names
    for m in sorted(MULTI, key=len, reverse=True):
        if m in rest:
            COUNTRY_REGION[m] = r; rest = rest.replace(m, " ")
    for w in rest.split():
        COUNTRY_REGION.setdefault(w, r)
COUNTRY_REGION["congo"] = "africa"; COUNTRY_REGION["congo drc"] = "africa"; COUNTRY_REGION["guinea"] = "africa"
US_STATES = ("alabama alaska arizona arkansas california colorado connecticut delaware florida hawaii idaho illinois indiana iowa kansas kentucky louisiana maine maryland massachusetts michigan "
             "minnesota mississippi missouri montana nebraska nevada new hampshire new jersey new mexico new york north carolina north dakota ohio oklahoma oregon pennsylvania rhode island "
             "south carolina south dakota tennessee texas utah vermont virginia washington west virginia wisconsin wyoming").split()
US_MULTI = ["new hampshire", "new jersey", "new mexico", "new york", "north carolina", "north dakota", "rhode island", "south carolina", "south dakota", "west virginia"]
PLACE_TO_COUNTRY = {}
rest = " ".join(US_STATES)
for m in US_MULTI:
    PLACE_TO_COUNTRY[m] = "usa"; rest = rest.replace(m, " ")
for w in rest.split():
    PLACE_TO_COUNTRY.setdefault(w, "usa")
for p in "guangdong guangxi hunan hubei zhejiang jiangsu shandong henan hebei beijing shanghai anhui fujian jiangxi sichuan yunnan guizhou shaanxi shanxi liaoning jilin heilongjiang gansu xinjiang tibet chongqing tianjin hainan qinghai ningxia shenzhen wuhan guangzhou nanchang hangzhou nanjing".split():
    PLACE_TO_COUNTRY[p] = "china"
for p in "tokyo osaka hokkaido kyoto".split():
    PLACE_TO_COUNTRY[p] = "japan"
for p in "ontario quebec alberta manitoba saskatchewan british columbia nova scotia".split():
    PLACE_TO_COUNTRY[p] = "canada"
for p in "victoria queensland new south wales tasmania".split():
    PLACE_TO_COUNTRY[p] = "australia"


def canon_country(s):
    s = re.sub(r"[_\-/]+", " ", s.lower()).strip()
    return ALIAS.get(s, s)


# ---------------------------------------------------------------- host tables
HOST_RULES = [  # (class, species label, regex on lowercased host/evidence text)
    ("human", "homo sapiens", r"\bhomo sapiens\b|\bhuman\b|\bhumans\b|\bpatient\b|\bhomo\b"),
    ("swine", "sus scrofa", r"\bsus scrofa\b|\bswine\b|\bpig\b|\bpigs\b|\bporcine\b|\bwild boar\b|\bboar\b"),
    ("avian", "aves", r"\bgallus\b|\bchicken\b|\bduck\b|\banas\b|\bgoose\b|\banser\b|\bmeleagris\b|\bmallard\b|\bteal\b|\bgull\b|\bswan\b|\baves\b|\bavian\b|\bbird\b|\bquail\b|\bpheasant\b|\bpigeon\b|\bcolumba\b|\banatidae\b|\bshorebird\b|\bturnstone\b|\bpelican\b|\bostrich\b|\bpartridge\b|\bguineafowl\b|\bchukar\b|\bwaterfowl\b"),
    ("bovine", "bos", r"\bbos\b|\bcattle\b|\bcow\b|\bbovine\b|\bcalf\b|\bbuffalo\b|\bbubalus\b"),
    ("equine", "equus", r"\bequus\b|\bhorse\b|\bequine\b|\bdonkey\b"),
    ("canine", "canis", r"\bcanis\b|\bdog\b|\bcanine\b"),
    ("feline", "felis", r"\bfelis\b|\bcat\b|\bfeline\b"),
    ("bat", "chiroptera", r"\bchiroptera\b|\bbat\b|\bbats\b|\bmyotis\b|\beptesicus\b|\btadarida\b|\bdesmodus\b|\bpipistrellus\b|\bartibeus\b|\bmolossus\b|\blasiurus\b|\bnyctalus\b|\bvampire\b|\brhinolophus\b|\bpteropus\b|\bmyotis\b"),
    ("other_mammal", "other mammal", r"\bneovison\b|\bmustela\b|\bmink\b|\bodocoileus\b|\bdeer\b|\bcervus\b|\boryctolagus\b|\brabbit\b|\bcamelus\b|\bcamel\b|\bdromedary\b|\bprocyon\b|\braccoon\b|\bvulpes\b|\bfox\b|\bskunk\b|\bmephitis\b|\bmongoose\b|\bjackal\b|\bwolf\b|\bsheep\b|\bovis\b|\bgoat\b|\bcapra\b|\bmouse\b|\brat\b|\bferret\b"),
]
HOST_CLASSES = [c for c, _, _ in HOST_RULES]


GENUS_SYN = {"bovine": "bos", "cattle": "bos", "cow": "bos", "porcine": "sus", "swine": "sus", "pig": "sus", "boar": "sus", "canine": "canis", "dog": "canis", "cat": "felis", "feline": "felis",
             "chicken": "gallus", "goose": "anser", "swan": "cygnus", "duck": "anas", "mallard": "anas", "teal": "anas", "ring-necked duck": "aythya", "horse": "equus", "equine": "equus",
             "human": "homo", "bat": "chiroptera", "turkey": "meleagris", "quail": "coturnix", "pigeon": "columba", "fox": "vulpes", "mink": "neovison", "deer": "odocoileus", "camel": "camelus"}


def host_genus(raw):
    raw = (raw or "").lower().strip()
    return GENUS_SYN.get(raw, raw.split()[0] if raw else "")


def host_class_of(text):
    hits = [c for c, _, rx in HOST_RULES if re.search(rx, text)]
    return hits  # list of classes found


EXCLUDE_TEXT = re.compile(r"vaccine|reassortant|recombinant|\bclone\b|synthetic|reverse genetics|laboratory strain|\bnc[0-9]", re.I)
# passage is flagged, not excluded: it does not change the true host or place of the original sample
PASSAGE = re.compile(r"passage|passaged|\begg\b|\bmdck\b|cell culture|\bvero\b", re.I)
PROT_LEN = dict(influenza_a=(500, 600), rabies=(400, 470), rotavirus_a=(300, 345))   # hev dropped: only 100-199 aa amplicon fragments, no independent evidence


def year_of(s):
    m = re.search(r"(19\d\d|20[0-2]\d)", s or "")
    return int(m.group(1)) if m else None


def evidence(rec):
    """Independent evidence from strain/isolate/definition text: host classes, countries, year."""
    name = " ".join(x for x in (rec["strain"], rec["isolate"]) if x)
    text = " ".join([name, rec["definition"]]).lower()
    t2 = re.sub(r"[_\-/,;()\[\]]+", " ", text)
    hosts, countries, years = set(), set(), set()
    # structured influenza name A/<host>/<place>/<id>/<year>
    m = re.match(r"^\s*a/(.+)$", rec["strain"].strip(), re.I)
    segs = m.group(1).split("/") if m else []
    if rec["virus"] == "influenza_a" and len(segs) >= 3:
        yr = re.sub(r"\(.*$", "", segs[-1])
        years.add((year_of(yr) or (2000 + int(yr) if re.fullmatch(r"[0-4]\d", yr) else 1900 + int(yr) if re.fullmatch(r"[5-9]\d", yr) else 0)))
        place = segs[-3].lower().strip()
        host_tok = segs[0].lower().strip() if len(segs) >= 4 else "human"
        hosts.update(host_class_of(host_tok) if len(segs) >= 4 else ["human"])
        pl = canon_country(place)
        if pl in PLACE_TO_COUNTRY:
            countries.add(PLACE_TO_COUNTRY[pl])
        elif pl in COUNTRY_REGION:
            countries.add(pl)
    elif rec["virus"] == "rotavirus_a" and rec["strain"].startswith("RVA/") and len(rec["strain"].split("/")) >= 5:
        sg = rec["strain"].split("/")
        hosts.update(host_class_of(re.sub(r"[-_]+", " ", sg[1].lower())))
        if sg[2] in ISO3:
            countries.add(ISO3[sg[2]])
        y = year_of(sg[4])
        if y:
            years.add(y)
    else:
        hosts.update(host_class_of(t2))
        for n in sorted(list(COUNTRY_REGION) + list(PLACE_TO_COUNTRY) + list(ALIAS), key=len, reverse=True):
            if re.search(r"\b" + re.escape(n) + r"\b", t2):
                c = canon_country(n)
                countries.add(PLACE_TO_COUNTRY.get(c, c))
        for y in re.findall(r"\b(19\d\d|20[0-2]\d)\b", t2):
            years.add(int(y))
    years.discard(0)
    return hosts, {c for c in countries if c in COUNTRY_REGION or c == "usa"}, years


def tier(field_val, ev):
    """A = independent evidence agrees and does not contradict; B = no independent evidence; X = conflict."""
    if not ev:
        return "B"
    if isinstance(ev, set) and len(ev) == 1 and field_val in ev:
        return "A"
    if field_val in ev and len(ev) > 1:
        return "B"   # ambiguous evidence: do not promote
    return "X"


def kmers(s, k=5):
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def main():
    rows = list(csv.DictReader(open(f"{OUT}/records.tsv", encoding="utf-8"), delimiter="\t"))
    log = collections.Counter()
    kept = []
    for r in rows:
        v = r["virus"]; log[(v, "fetched")] += 1
        if v not in PROT_LEN:
            log[(v, "drop_virus_not_in_scope")] += 1; continue
        lo, hi = PROT_LEN[v]
        p = r["protein"]
        if not (lo <= len(p) <= hi) or p.count("X") > 0.01 * len(p) or "*" in p:
            log[(v, "drop_protein_qc")] += 1; continue
        r["passaged"] = int(bool(PASSAGE.search(" ".join([r["note"], r["isolation_source"], r["strain"], r["isolate"], r["definition"]]))))
        if r["lab_host"] or EXCLUDE_TEXT.search(" ".join([r["note"], r["isolation_source"], r["strain"], r["isolate"], r["definition"]])):
            log[(v, "drop_lab_or_passage_or_construct")] += 1; continue
        hraw = r["host_raw"].lower().split(";")[0].strip()
        hc = host_class_of(hraw)
        r["host_class_field"] = hc[0] if len(hc) == 1 else ""
        r["host_species_field"] = hraw
        craw = r["country_raw"].split(":")[0].strip()
        c = canon_country(craw)
        c = PLACE_TO_COUNTRY.get(c, c)
        r["country_field"] = c if (c in COUNTRY_REGION or c == "usa") else ""
        r["region_field"] = "north_america" if c == "usa" else COUNTRY_REGION.get(c, "")
        r["year_field"] = year_of(r["collection_date"]) or ""
        hosts_ev, countries_ev, years_ev = evidence(r)
        r["host_tier"] = tier(r["host_class_field"], hosts_ev) if r["host_class_field"] else ""
        r["geo_tier"] = tier(r["country_field"], countries_ev) if r["country_field"] else ""
        r["year_tier"] = (("A" if r["year_field"] in years_ev and len(years_ev) == 1 else "X" if years_ev and r["year_field"] not in years_ev else "B") if r["year_field"] else "")
        r["host_evidence"] = ",".join(sorted(hosts_ev)); r["geo_evidence"] = ",".join(sorted(countries_ev)); r["year_evidence"] = ",".join(map(str, sorted(years_ev)))
        for lab, tcol in (("host", "host_tier"), ("geo", "geo_tier"), ("year", "year_tier")):
            if r[tcol] == "X":
                log[(v, f"conflict_{lab}")] += 1
        # a conflicting label is blanked, never kept
        if r["host_tier"] == "X":
            r["host_class_field"] = ""; r["host_tier"] = ""
        if r["geo_tier"] == "X":
            r["country_field"] = ""; r["region_field"] = ""; r["geo_tier"] = ""
        if r["year_tier"] == "X":
            r["year_field"] = ""; r["year_tier"] = ""
        if not (r["host_class_field"] or r["country_field"]):
            log[(v, "drop_no_label")] += 1; continue
        kept.append(r)
    # exact duplicates: same protein and same labels
    seen, dedup = set(), []
    for r in kept:
        key = (r["virus"], r["protein"], r["host_class_field"], r["country_field"], r["year_field"])
        if key in seen:
            log[(r["virus"], "drop_exact_duplicate")] += 1; continue
        seen.add(key); dedup.append(r)
    kept = dedup
    # approximate similarity clusters (5-mer Jaccard >= 0.7, about 96% identity), leader algorithm per virus
    leaders = collections.defaultdict(list)
    for r in sorted(kept, key=lambda x: x["accession"]):
        ks = kmers(r["protein"]); r["_k"] = ks
        for cid, lk in leaders[r["virus"]]:
            if len(ks & lk) / len(ks | lk) >= 0.7:
                r["cluster"] = cid; break
        else:
            cid = f"{r['virus']}_c{len(leaders[r['virus']]):04d}"; leaders[r["virus"]].append((cid, ks)); r["cluster"] = cid
    # splits: blind (sealed, cluster-level), geo-held-out countries, novel host species, time, dev, train
    def h(s, salt):
        return int(hashlib.sha1((salt + s).encode()).hexdigest(), 16) % 1000 / 1000.0
    countries_by_region = collections.defaultdict(set)
    species_by_class = collections.defaultdict(set)
    for r in kept:
        if r["country_field"]:
            countries_by_region[r["region_field"]].add(r["country_field"])
        if r["host_class_field"]:
            species_by_class[r["host_class_field"]].add(host_genus(r["host_species_field"]))
    cnt = collections.Counter(r["country_field"] for r in kept if r["country_field"])
    geo_hold = set()
    for rg, cs in countries_by_region.items():
        big = {c for c, _ in sorted(((c, cnt[c]) for c in cs), key=lambda t: -t[1])[:2]}
        cand = sorted((c for c in cs if c not in big and cnt[c] >= 5 and h(c, "geo") < 0.5), key=lambda c: h(c, "geo"))
        geo_hold.update(cand[:3])
    hcnt = collections.Counter(host_genus(r["host_species_field"]) for r in kept if r["host_class_field"])
    host_hold = set()
    for cl, gs in species_by_class.items():
        big = {x for x, _ in sorted(((x, hcnt[x]) for x in gs), key=lambda t: -t[1])[:2]}
        cand = sorted((x for x in gs if x not in big and hcnt[x] >= 5 and h(x, "host") < 0.6), key=lambda x: h(x, "host"))
        host_hold.update(cand[:2])
    for r in kept:
        r["split"] = "train"
        if h(r["cluster"], "blind") < 0.10:
            r["split"] = "blind"
        elif r["country_field"] in geo_hold:
            r["split"] = "test_geo"
        elif host_genus(r["host_species_field"]) in host_hold:
            r["split"] = "test_host"
        elif isinstance(r["year_field"], int) and r["year_field"] >= 2022:
            r["split"] = "test_time"
        elif h(r["cluster"], "dev") < 0.10:
            r["split"] = "dev"
    # nearest-train similarity (novelty) for every non-train row, per virus
    trains = collections.defaultdict(list)
    for r in kept:
        if r["split"] == "train":
            trains[r["virus"]].append(r["_k"])
    for r in kept:
        if r["split"] == "train":
            r["nn_sim_train"] = ""; continue
        best = 0.0
        for lk in trains[r["virus"]]:
            if len(r["_k"] & lk) / len(r["_k"] | lk) > best:
                best = len(r["_k"] & lk) / len(r["_k"] | lk)
        r["nn_sim_train"] = round(best, 3)
    # derived (synthetic, weak) label: host classes observed among ALL sequences in the same cluster (evidence of compatibility by observation)
    by_cluster = collections.defaultdict(set)
    for r in kept:
        if r["host_class_field"]:
            by_cluster[r["cluster"]].add(r["host_class_field"])
    for r in kept:
        s = by_cluster.get(r["cluster"], set())
        r["derived_cluster_hosts"] = ",".join(sorted(s))
        r["derived_human_in_cluster"] = int("human" in s) if s else ""
    cols = ["accession", "virus", "protein", "protein_len", "host_class_field", "host_species_field", "host_tier", "host_evidence", "country_field", "region_field", "geo_tier", "geo_evidence",
            "year_field", "year_tier", "year_evidence", "serotype", "strain", "isolate", "isolation_source", "split", "cluster", "nn_sim_train", "derived_cluster_hosts", "derived_human_in_cluster", "passaged"]
    hdr = {"host_class_field": "host_class", "host_species_field": "host_species_raw", "country_field": "country", "region_field": "region", "year_field": "year"}
    with gzip.open(f"{GT}/p4_groundtruth.tsv.gz", "wt", newline="\n", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n"); w.writerow([hdr.get(c, c) for c in cols])
        for r in sorted(kept, key=lambda x: (x["virus"], x["accession"])):
            w.writerow([r.get(c, "") for c in cols])
    # audit sample: 60 random gold-tier rows with their evidence, no sequence
    rng = random.Random(7)
    gold = [r for r in kept if r["host_tier"] == "A" or r["geo_tier"] == "A"]
    aud = rng.sample(gold, min(60, len(gold)))
    with open(f"{GT}/p4_audit_sample.tsv", "w", newline="\n", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "virus", "host_class", "host_species_raw", "host_tier", "host_evidence", "country", "geo_tier", "geo_evidence", "year", "strain", "isolate", "reviewer_verdict", "reviewer_note"])
        for r in aud:
            w.writerow([r["accession"], r["virus"], r["host_class_field"], r["host_species_field"], r["host_tier"], r["host_evidence"], r["country_field"], r["geo_tier"], r["geo_evidence"], r["year_field"], r["strain"], r["isolate"], "", ""])
    # summary
    with open(f"{ROOT}/results/p4_curation_summary.tsv", "w", newline="\n") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n"); w.writerow(["virus", "step", "count"])
        for (v, s), n in sorted(log.items()):
            w.writerow([v, s, n])
        for v in sorted({r["virus"] for r in kept}):
            vk = [r for r in kept if r["virus"] == v]
            w.writerow([v, "kept_total", len(vk)])
            for lab, col in (("host", "host_tier"), ("geo", "geo_tier"), ("year", "year_tier")):
                for t in "AB":
                    w.writerow([v, f"{lab}_tier_{t}", sum(r[col] == t for r in vk)])
            for sp in ("train", "dev", "test_geo", "test_host", "test_time", "blind"):
                w.writerow([v, f"split_{sp}", sum(r["split"] == sp for r in vk)])
    print("kept", len(kept), dict(collections.Counter(r["split"] for r in kept)))
    print("geo_hold", sorted(geo_hold)); print("host_hold", sorted(host_hold))


if __name__ == "__main__":
    main()
