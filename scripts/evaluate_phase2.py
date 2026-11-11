#!/usr/bin/env python3
"""Phase 2 scoring: constructed-truth spike-ins in a human background, and real plasma runs.

Inputs (results/raw/phase2): kraken2/*.report.gz (taxa with clade reads > 0), metaphlan/*.profile.txt,
results/spike_manifest.tsv (truth by construction), config/phase2_uganda_runs.tsv.
Outputs: results/spike_detection.tsv, spike_summary.tsv, spike_false_positives.tsv, uganda_summary.tsv,
uganda_top_species.tsv, uganda_tool_agreement.tsv.

Detection uses ABSOLUTE read counts (a species needs >= T read pairs), because in a host-heavy sample
percentages are meaningless: the microbe is a tiny fraction of everything."""
import csv, glob, gzip, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from mockcheck.parsers import norm, genus, load_aliases  # noqa: E402
from mockcheck.report import WATCH_GENERA  # noqa: E402

load_aliases(f"{ROOT}/config/aliases.tsv")
RAW, OUT = f"{ROOT}/results/raw/phase2", f"{ROOT}/results"
THR = [1, 3, 10]
SYN = {"bacillus subtilis": {"bacillus spizizenii"}}   # naming trap, see REPORT.md (reported as relaxed, flagged)
HUMAN = "homo sapiens"


def tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def kraken(path):
    """-> species {name: clade reads}, genera {name: clade reads}, totals dict."""
    sp, ge, tot = {}, {}, dict(root=0, unclassified=0, bacteria=0, viruses=0, eukaryota=0)
    with gzip.open(path, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            reads, rank, name = int(p[1]), p[3], p[5].strip()
            if rank == "U":
                tot["unclassified"] = reads
            elif rank == "R" and name == "root":
                tot["root"] = reads
            elif rank == "D":
                tot[name.lower()] = reads
            elif rank == "S":
                sp[norm(name)] = sp.get(norm(name), 0) + reads
            elif rank == "G":
                ge[norm(name)] = ge.get(norm(name), 0) + reads
    tot["total"] = tot["root"] + tot["unclassified"]
    return sp, ge, tot


def metaphlan(path):
    sp = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or "s__" not in line:
                continue
            p = line.rstrip("\n").split("\t")
            last = p[0].split("|")[-1]
            if last.startswith("s__") and "t__" not in p[0]:
                sp[norm(last[3:])] = float(p[2])
    return sp


def label(path, suffix):
    b = os.path.basename(path)[: -len(suffix)]
    sample, rest = b.split(".", 1)
    return sample, rest


# ---------------------------------------------------------------- spike-ins
manifest = tsv(f"{ROOT}/results/spike_manifest.tsv") if os.path.exists(f"{ROOT}/results/spike_manifest.tsv") else []
kr = {}
for p in glob.glob(f"{RAW}/kraken2/spike_*.report.gz"):
    sample, rest = label(p, ".report.gz")           # spike_P0, 16gb.c0.0
    kr[(sample, "kraken2_" + rest.replace(".", "_", 1).replace(".", "_"))] = kraken(p)
mp = {(label(p, ".profile.txt")[0], "metaphlan4"): metaphlan(p) for p in glob.glob(f"{RAW}/metaphlan/spike_*.profile.txt")}

det = []
for m in manifest:
    if m["organism"] == "(none)":
        continue
    org, sample, exp = norm(m["organism"]), m["sample"], int(m["spiked_pairs"])
    names = {org} | SYN.get(org, set())
    for (s, method), (sp, ge, tot) in kr.items():
        if s != sample:
            continue
        strict, relaxed = sp.get(org, 0), sum(sp.get(n, 0) for n in names)
        det.append(dict(sample=sample, method=method, organism=org, fraction=float(m["fraction"]), expected_pairs=exp,
                        species_reads_strict=strict, species_reads_relaxed=relaxed, genus_reads=ge.get(genus(org), 0),
                        recovery_relaxed=round(relaxed / exp, 3), **{f"detected_ge{t}": int(relaxed >= t) for t in THR}))
    for (s, method), prof in mp.items():
        if s == sample:
            hit = int(any(prof.get(n, 0) > 0 for n in names))
            det.append(dict(sample=sample, method=method, organism=org, fraction=float(m["fraction"]), expected_pairs=exp,
                            species_reads_strict="", species_reads_relaxed="", genus_reads="", recovery_relaxed="",
                            **{f"detected_ge{t}": hit for t in THR}))


def write(name, rows):
    if not rows:
        return
    with open(f"{OUT}/{name}", "w", newline="\n") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)


write("spike_detection.tsv", sorted(det, key=lambda r: (r["method"], -r["fraction"], r["organism"])))

# sensitivity by level, split bacteria vs yeast
YEAST = {"saccharomyces cerevisiae", "cryptococcus neoformans"}
summ = []
for method in sorted({r["method"] for r in det}):
    for frac in sorted({r["fraction"] for r in det}, reverse=True):
        for cls, pick in (("bacteria", lambda o: o not in YEAST), ("yeast", lambda o: o in YEAST)):
            rows = [r for r in det if r["method"] == method and r["fraction"] == frac and pick(r["organism"])]
            if not rows:
                continue
            rec = [r["recovery_relaxed"] for r in rows if r["recovery_relaxed"] != ""]
            summ.append(dict(method=method, fraction=frac, expected_pairs=rows[0]["expected_pairs"], group=cls, n=len(rows),
                             **{f"detected_ge{t}": sum(r[f"detected_ge{t}"] for r in rows) for t in THR},
                             median_recovery=(sorted(rec)[len(rec) // 2] if rec else "")))
write("spike_summary.tsv", summ)

# false positives: species not spiked into that sample (and not human) at >= T reads
truth_by = {}
for m in manifest:
    if m["organism"] != "(none)":
        o = norm(m["organism"]); truth_by.setdefault(m["sample"], set()).update({o} | SYN.get(o, set()))
fps = []
for (sample, method), (sp, ge, tot) in sorted(kr.items()):
    truth = truth_by.get(sample, set())
    for t in THR:
        bad = sorted(((v, k) for k, v in sp.items() if k != HUMAN and k not in truth and v >= t), reverse=True)
        watch = [k for v, k in bad if genus(k) in WATCH_GENERA]
        fps.append(dict(sample=sample, method=method, min_reads=t, false_species=len(bad),
                        false_species_pathogen_genus=len(watch), top3="; ".join(f"{k}:{v}" for v, k in bad[:3]),
                        human_clade_reads=sp.get(HUMAN, 0), total_pairs=tot["total"]))
write("spike_false_positives.tsv", fps)

# ---------------------------------------------------------------- real plasma (no per-sample truth)
contam = {r["genus"] for r in tsv(f"{ROOT}/config/reagent_contaminant_genera.tsv")} if os.path.exists(f"{ROOT}/config/reagent_contaminant_genera.tsv") else set()
runs = {r["run"]: r["sample_alias"] for r in tsv(f"{ROOT}/config/phase2_uganda_runs.tsv")} if os.path.exists(f"{ROOT}/config/phase2_uganda_runs.tsv") else {}
ug, top, agree = [], [], []
mp_ug = {label(p, ".profile.txt")[0]: metaphlan(p) for p in glob.glob(f"{RAW}/metaphlan/ug_*.profile.txt")}
for p in sorted(glob.glob(f"{RAW}/kraken2/ug_*.report.gz")):
    sample, rest = label(p, ".report.gz")
    sp, ge, tot = kraken(p)
    conf = rest.split(".c")[-1]
    run = sample[3:]
    n = tot["total"] or 1
    nonhuman = {k: v for k, v in sp.items() if k != HUMAN}
    for t in (10,):
        called = {k: v for k, v in nonhuman.items() if v >= t}
        cg = sum(v for k, v in called.items() if genus(k) in contam)
        tg = sum(called.values()) or 1
        ug.append(dict(run=run, alias=runs.get(run, ""), confidence=conf, pairs=tot["total"],
                       human_frac=round(sp.get(HUMAN, 0) / n, 4), unclassified_frac=round(tot["unclassified"] / n, 4),
                       bacteria_pairs=tot.get("bacteria", 0), virus_pairs=tot.get("viruses", 0),
                       species_ge10=len(called), reagent_contaminant_genus_share=round(cg / tg, 3),
                       pathogen_genus_species_ge10=sum(genus(k) in WATCH_GENERA for k in called)))
        for v, k in sorted(((v, k) for k, v in called.items()), reverse=True)[:5]:
            top.append(dict(run=run, confidence=conf, species=k, pairs=v, reagent_contaminant_genus=int(genus(k) in contam)))
        if sample in mp_ug and conf == "0.1":
            kset, mset = set(called), {k for k, v in mp_ug[sample].items() if v > 0}
            u = kset | mset
            agree.append(dict(run=run, kraken2_species=len(kset), metaphlan_species=len(mset), shared=len(kset & mset),
                              jaccard=round(len(kset & mset) / len(u), 3) if u else ""))
write("uganda_summary.tsv", ug)
write("uganda_top_species.tsv", top)
write("uganda_tool_agreement.tsv", agree)
print(f"spike rows {len(det)}, summary {len(summ)}, fp rows {len(fps)}, uganda {len(ug)}")
