#!/usr/bin/env python3
"""Parse raw classifier outputs in results/raw and score them against the
expected composition. Writes small tables to results/ and figures to results/fig.

Everything is species level unless a column says genus. Abundances are
percent of the species-level profile (each profile is renormalised to 100)."""
import gzip, glob, os, sys
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW, OUT = f"{ROOT}/results/raw", f"{ROOT}/results"
THRESH = [0.0, 0.01, 0.1, 1.0]
YEASTS = ("saccharomyces cerevisiae", "cryptococcus neoformans")

aliases = pd.read_csv(f"{ROOT}/config/aliases.tsv", sep="\t").set_index("alias")["canonical"].to_dict()


def norm(name):
    n = name.replace("_", " ").strip().lower()
    return aliases.get(n, n)


def genus(name):
    g = name.split()[0]
    return aliases.get(g, g)


def read_kraken(path):
    sp, total, unclass = {}, 0, 0
    with gzip.open(path, "rt") as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            reads, rank, name = int(p[1]), p[3], p[5].strip()
            if rank == "U":
                unclass = reads
            elif rank == "R" and name == "root":
                total = reads
            elif rank == "S":
                sp[norm(name)] = sp.get(norm(name), 0) + reads
    allr = total + unclass
    s = sum(sp.values())
    return ({k: 100 * v / s for k, v in sp.items()} if s else {}), (s / allr if allr else 0)


def read_bracken(path):
    df = pd.read_csv(path, sep="\t")
    sp = {}
    for n, v in zip(df["name"], df["fraction_total_reads"]):
        sp[norm(n)] = sp.get(norm(n), 0) + v
    s = sum(sp.values())
    return ({k: 100 * v / s for k, v in sp.items()} if s else {}), s


def read_metaphlan(path):
    sp = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or "s__" not in line:
                continue
            clade, _, ab = line.split("\t")[:3]
            last = clade.split("|")[-1]
            if last.startswith("s__") and "t__" not in clade:
                sp[norm(last[3:])] = float(ab)
    s = sum(sp.values())
    return ({k: 100 * v / s for k, v in sp.items()} if s else {}), None


profiles, meta = {}, {}
for p in sorted(glob.glob(f"{RAW}/kraken2/*.report.gz")):
    sample, db, conf = os.path.basename(p)[: -len(".report.gz")].split(".", 2)
    m = f"kraken2_{db}_{conf}"
    profiles[(sample, m)], meta[(sample, m)] = read_kraken(p)
for p in sorted(glob.glob(f"{RAW}/bracken/*.bracken")):
    sample, db, conf = os.path.basename(p)[: -len(".bracken")].split(".", 2)
    m = f"bracken_{db}_{conf}"
    profiles[(sample, m)], meta[(sample, m)] = read_bracken(p)
for p in sorted(glob.glob(f"{RAW}/metaphlan/*.profile.txt")):
    sample = os.path.basename(p)[: -len(".profile.txt")]
    profiles[(sample, "metaphlan4")], meta[(sample, "metaphlan4")] = read_metaphlan(p)
methods = sorted({m for _, m in profiles})

exp = pd.read_csv(f"{ROOT}/config/expected.tsv", sep="\t")
exp["organism"] = exp["organism"].map(norm)
samples = pd.read_csv(f"{ROOT}/config/samples.tsv", sep="\t")


def bray(a, b):
    keys = set(a) | set(b)
    num = sum(abs(a.get(k, 0) - b.get(k, 0)) for k in keys)
    den = sum(a.get(k, 0) + b.get(k, 0) for k in keys)
    return num / den if den else float("nan")


def relax(prof):
    """Sensitivity analysis only: merge B. spizizenii (a former B. subtilis subspecies in NCBI)
    into B. subtilis. Reported as *_relaxed columns, never as the headline."""
    r = dict(prof)
    if "bacillus spizizenii" in r:
        r["bacillus subtilis"] = r.get("bacillus subtilis", 0) + r.pop("bacillus spizizenii")
    return r


def prf(pred, truth):
    tp = len(pred & truth)
    fp = len(pred - truth)
    fn = len(truth - pred)
    pr = tp / (tp + fp) if tp + fp else 0.0
    rc = tp / (tp + fn) if tp + fn else 0.0
    return tp, fp, fn, pr, rc, (2 * pr * rc / (pr + rc) if pr + rc else 0.0)


# ---- even mock: quantitative
rows, per = [], []
e = exp[exp["sample"] == "even"].set_index("organism")["expected_percent"].astype(float).to_dict()
# Zymo datasheet: expected_percent is theoretical genomic DNA (matches read share);
# the cell column is the fair truth for tools that report taxon (cell-like) abundance.
ecell = exp[exp["sample"] == "even"].set_index("organism")["expected_cell_percent"].astype(float).to_dict()
CELL_NATIVE = ("metaphlan4",)
bact ={k: v for k, v in e.items() if k not in YEASTS}
sb = sum(bact.values())
bact = {k: 100 * v / sb for k, v in bact.items()}
for m in methods:
    prof = profiles.get(("even", m))
    if prof is None:
        continue
    truth = set(e)
    pb = {k: v for k, v in prof.items() if k in bact}
    spb = sum(pb.values())
    pbn = {k: 100 * v / spb for k, v in pb.items()} if spb else {}
    af = meta[("even", m)]
    for i, t in enumerate(THRESH):
        pred = {k for k, v in prof.items() if v >= t and v > 0}
        tp, fp, fn, pr, rc, f1 = prf(pred, truth)
        rp = relax(prof)
        _, fpr, fnr, prr, rcr, f1r = prf({k for k, v in rp.items() if v >= t and v > 0}, truth)
        r = dict(method=m, threshold_pct=t, tp=tp, fp=fp, fn=fn,
                 precision=round(pr, 3), recall=round(rc, 3), f1=round(f1, 3),
                 fn_relaxed=fnr, recall_relaxed=round(rcr, 3), f1_relaxed=round(f1r, 3))
        if i == 0:
            r["bray_curtis_vs_expected"] = round(bray(prof, e), 3)
            r["bray_curtis_vs_cell_truth"] = round(bray(prof, ecell), 3)
            r["bray_curtis_relaxed_native"] = round(bray(rp, ecell if m.startswith(CELL_NATIVE) else e), 3)
            r["native_basis"] = "cell" if m.startswith(CELL_NATIVE) else "dna"
            r["bray_curtis_native_basis"] = r["bray_curtis_vs_cell_truth"] if r["native_basis"] == "cell" else r["bray_curtis_vs_expected"]
            r["false_positive_mass_pct"] = round(sum(v for k, v in prof.items() if k not in truth), 2)
            r["bray_curtis_8_bacteria_only"] = round(bray(pbn, bact), 3)
            r["species_assigned_frac_of_reads"] = None if af is None else round(af, 3)
        rows.append(r)
    for org, ev in e.items():
        per.append(dict(method=m, organism=org, expected_pct=ev, expected_cell_pct=ecell[org],
                        observed_pct=round(prof.get(org, 0), 3)))
    for v, k in sorted(((v, k) for k, v in prof.items() if k not in truth), reverse=True)[:5]:
        per.append(dict(method=m, organism=f"[not in mock] {k}", expected_pct=0, expected_cell_pct=0, observed_pct=round(v, 3)))
pd.DataFrame(rows).to_csv(f"{OUT}/even_summary.tsv", sep="\t", index=False)
pd.DataFrame(per).to_csv(f"{OUT}/even_per_organism.tsv", sep="\t", index=False)

# ---- negative controls: single isolates
nc = []
for _, s in samples[samples["role"] == "negative_control"].iterrows():
    tr = norm(s["truth"])
    for m in methods:
        prof = profiles.get((s["sample"], m))
        if prof is None:
            continue
        pur = prof.get(tr, 0)
        others = {k: v for k, v in prof.items() if k != tr}
        topk, topv = max(others.items(), key=lambda kv: kv[1]) if others else ("", 0)
        nc.append(dict(sample=s["sample"], method=m, truth=tr, truth_pct=round(pur, 2),
                       other_species_ge_0_1pct=sum(v >= 0.1 for v in others.values()),
                       other_species_ge_1pct=sum(v >= 1 for v in others.values()),
                       other_mass_pct=round(100 - pur, 2),
                       top_other=topk, top_other_pct=round(topv, 2),
                       same_genus_other_pct=round(sum(v for k, v in others.items() if genus(k) == genus(tr)), 2)))
pd.DataFrame(nc).to_csv(f"{OUT}/negative_controls.tsv", sep="\t", index=False)

# ---- log mock: presence only plus the one reliable quantitative anchor
lg = []
ten = set(e)
for m in methods:
    prof = profiles.get(("log", m))
    if prof is None:
        continue
    lg.append(dict(method=m, listeria_pct=round(prof.get("listeria monocytogenes", 0), 2),
                   listeria_expected_pct=89.1,
                   mock_species_detected_ge_0_01pct=len({k for k, v in prof.items() if v >= 0.01} & ten),
                   non_mock_species_ge_0_1pct=len({k for k, v in prof.items() if v >= 0.1} - ten),
                   non_mock_mass_pct=round(sum(v for k, v in prof.items() if k not in ten), 2)))
pd.DataFrame(lg).to_csv(f"{OUT}/log_summary.tsv", sep="\t", index=False)

# ---- timing (even mock only; one classifier call per row)
tp = f"{RAW}/timing.tsv"
if os.path.exists(tp):
    t = pd.read_csv(tp, sep="\t", header=None, names=["label", "wall_s", "max_rss_kb"])
    t = t.drop_duplicates("label", keep="last")  # a failed first attempt may precede the real run
    t = t[t["label"].str.endswith(".even")].copy()
    t["max_rss_gb"] = (t["max_rss_kb"] / 1048576).round(1)
    t[["label", "wall_s", "max_rss_gb"]].to_csv(f"{OUT}/timing_even.tsv", sep="\t", index=False)

# ---- figures
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    os.makedirs(f"{OUT}/fig", exist_ok=True)
    pe = pd.DataFrame(per)
    pe = pe[~pe.organism.str.startswith("[")]
    piv = pe.pivot(index="organism", columns="method", values="observed_pct")
    exp_s = pe.drop_duplicates("organism").set_index("organism")["expected_pct"]
    piv = piv.loc[exp_s.sort_values(ascending=False).index]
    ax = piv.plot.bar(figsize=(11, 5), width=0.85, colormap="tab10")
    ax.scatter(range(len(piv)), [exp_s[o] for o in piv.index], marker="_", s=900, color="k", zorder=5, label="expected")
    ax.set_ylabel("% of species-level profile")
    ax.set_xlabel("")
    ax.set_title("Even mock (D6300): observed vs expected")
    ax.legend(fontsize=7, ncol=2)
    plt.xticks(rotation=40, ha="right")
    plt.tight_layout()
    plt.savefig(f"{OUT}/fig/even_observed_vs_expected.png", dpi=130)
    plt.close()
    ncd = pd.DataFrame(nc)
    if len(ncd):
        g = ncd.groupby("method")["other_mass_pct"].mean().sort_values()
        ax = g.plot.barh(figsize=(7, 4), color="#4c78a8")
        ax.set_xlabel("mean % of profile that is NOT the true organism (10 isolates)")
        plt.tight_layout()
        plt.savefig(f"{OUT}/fig/negative_control_false_mass.png", dpi=130)
        plt.close()
except Exception as ex:
    print("figures skipped:", ex, file=sys.stderr)
print("methods:", methods)
