#!/usr/bin/env python3
"""Phase 4 evaluation suite (CPU only): checks at every stage, from label provenance to host and region predictions.

Stages (each prints tripwires = things that would invalidate or inflate a result):
  S0 label provenance   : evidence tiers per label and split, missingness
  S1 split integrity    : cluster / exact-sequence / country / host-genus / year overlap between train and each split, similarity to train, cell sizes
  S2 shortcut baselines : what virus identity, length, or 3-mers alone achieve (so model gains are measured against them), plus a train-label-permutation control
  S3 model predictions  : per split, per virus, tier-A labels only, calibration (ECE), accuracy vs novelty (similarity to train); clustered bootstrap CIs
Usage: python p4_eval_suite.py [--boot 200] [--out results/p4_eval_suite.json]
Guardrails: classification of existing public sequences only."""
import argparse, collections, csv, glob, gzip, json, os, re
import numpy as np
from sklearn.feature_extraction.text import CountVectorizer, TfidfTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HOSTS = ["avian", "human", "swine", "equine", "canine", "bovine", "feline", "other_mammal", "bat"]
REGIONS = ["africa", "asia", "europe", "north_america", "oceania", "south_america"]
SPLITS = ["dev", "test_geo", "test_host", "test_time", "blind"]
HEADS = {"host": ("host_class", HOSTS, "host"), "region": ("region", REGIONS, "reg")}
MIN_CELL = 30


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, (float, np.floating)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    return o


def load_gt():
    with gzip.open(os.path.join(ROOT, "data/groundtruth/p4_groundtruth.tsv.gz"), "rt") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def genus(raw):
    w = re.sub(r"[^a-z ]", " ", (raw or "").lower()).split()
    return w[0] if w else ""


def acc(y, p):
    return float(np.mean(np.array(y) == np.array(p))) if len(y) else float("nan")


def bacc(y, p):
    y, p = np.array(y), np.array(p)
    cl = [c for c in set(y)]
    return float(np.mean([np.mean(p[y == c] == c) for c in cl])) if cl else float("nan")


def cboot(clusters, y, p, fn, reps, rng):
    """Bootstrap by cluster (sequences in one cluster are not independent)."""
    cl = np.array(clusters); y = np.array(y); p = np.array(p)
    u = sorted(set(cl)); idx = {c: np.where(cl == c)[0] for c in u}
    if len(u) < 3:
        return [float("nan"), float("nan")]
    v = []
    for _ in range(reps):
        ii = np.concatenate([idx[u[j]] for j in rng.integers(0, len(u), len(u))])
        v.append(fn(y[ii], p[ii]))
    return [float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))]


def paired_diff(clusters, ok_a, ok_b, reps, rng):
    cl = np.array(clusters); a = np.array(ok_a, float); b = np.array(ok_b, float)
    u = sorted(set(cl)); idx = {c: np.where(cl == c)[0] for c in u}
    if len(u) < 3:
        return [float("nan"), float("nan")]
    v = []
    for _ in range(reps):
        ii = np.concatenate([idx[u[j]] for j in rng.integers(0, len(u), len(u))])
        v.append(a[ii].mean() - b[ii].mean())
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]


def s0_s1(R):
    out = {"S0_provenance": {}, "S1_integrity": {}, "tripwires": []}
    for s in SPLITS + ["train"]:
        rs = [r for r in R if r["split"] == s]
        n = len(rs)
        out["S0_provenance"][s] = dict(n=n,
            host_tier=dict(collections.Counter(r["host_tier"] or "none" for r in rs)),
            geo_tier=dict(collections.Counter(r["geo_tier"] or "none" for r in rs)))
    tr = [r for r in R if r["split"] == "train"]
    trc, trs = {r["cluster"] for r in tr}, {r["protein"] for r in tr}
    trc_ctry = {r["country"] for r in tr if r["country"]}
    trg = {genus(r["host_species_raw"]) for r in tr if r["host_species_raw"]}
    for s in SPLITS:
        rs = [r for r in R if r["split"] == s]
        sim = [float(r["nn_sim_train"]) for r in rs if r["nn_sim_train"]]
        d = dict(n=len(rs), by_virus=dict(collections.Counter(r["virus"] for r in rs)),
                 clusters_shared_with_train=len({r["cluster"] for r in rs} & trc), n_clusters=len({r["cluster"] for r in rs}),
                 seqs_exact_in_train=sum(r["protein"] in trs for r in rs),
                 nn_sim_quantiles={q: float(np.percentile(sim, q)) for q in (5, 25, 50, 75, 95)} if sim else None,
                 frac_sim_ge_0_9=float(np.mean([x >= 0.9 for x in sim])) if sim else None)
        if s == "test_geo":
            d["countries_shared_with_train"] = sorted({r["country"] for r in rs if r["country"]} & trc_ctry)
        if s == "test_host":
            d["host_genera_shared_with_train"] = sorted({genus(r["host_species_raw"]) for r in rs if r["host_species_raw"]} & trg - {""})
        if s == "test_time":
            ys = [int(r["year"]) for r in rs if r["year"].isdigit()]
            d["min_year"] = min(ys) if ys else None
        out["S1_integrity"][s] = d
        if s == "blind" and (d["clusters_shared_with_train"] or d["seqs_exact_in_train"]):
            out["tripwires"].append(f"LEAK: sealed blind split shares {d['clusters_shared_with_train']} clusters / {d['seqs_exact_in_train']} exact sequences with train")
        if s in ("test_geo", "test_host") and d.get("countries_shared_with_train" if s == "test_geo" else "host_genera_shared_with_train"):
            out["tripwires"].append(f"{s} is not disjoint on its holdout key: {d.get('countries_shared_with_train') or d.get('host_genera_shared_with_train')}")
        if s == "test_time" and d["min_year"] and d["min_year"] < 2022:
            out["tripwires"].append("test_time contains years before 2022")
        if d["seqs_exact_in_train"]:
            out["tripwires"].append(f"{s}: {d['seqs_exact_in_train']} sequences are exact copies of training sequences (inflates any similarity-based method)")
        for v, c in d["by_virus"].items():
            if c < MIN_CELL:
                out["tripwires"].append(f"{s}: only {c} {v} sequences (<{MIN_CELL}); per-virus results there are anecdotes, not evidence")
        missing = [v for v in ("influenza_a", "rabies", "rotavirus_a") if v not in d["by_virus"]]
        if missing:
            out["tripwires"].append(f"{s}: no {', '.join(missing)} sequences at all")
    return out


def feats_lr(rows_tr, rows_te):
    cv = CountVectorizer(analyzer="char", ngram_range=(3, 3), lowercase=False)
    tf = TfidfTransformer()
    Xtr = tf.fit_transform(cv.fit_transform([r["protein"] for r in rows_tr]))
    return Xtr, tf.transform(cv.transform([r["protein"] for r in rows_te]))


def s2(R, reps, rng, restrict=None):
    tr_all = [r for r in R if r["split"] == "train"]
    out, rowp = {}, {}
    for head, (col, classes, _) in HEADS.items():
        tr = [r for r in tr_all if r[col] in classes and (restrict is None or r["accession"] in restrict)]  # same training rows as the models
        te_by = {s: [r for r in R if r["split"] == s and r[col] in classes and (restrict is None or r["accession"] in restrict)] for s in SPLITS}
        allte = [r for s in SPLITS for r in te_by[s]]
        maj = collections.Counter(r[col] for r in tr).most_common(1)[0][0]
        vprior = {v: collections.Counter(r[col] for r in tr if r["virus"] == v).most_common(1)[0][0] for v in {r["virus"] for r in tr}}
        V = sorted({r["virus"] for r in R})
        mat = lambda rs: np.array([[float(r["protein_len"])] + [float(r["virus"] == v) for v in V] for r in rs])
        tree = DecisionTreeClassifier(max_depth=4, random_state=0).fit(mat(tr), [r[col] for r in tr])
        Xtr, Xte = feats_lr(tr, allte)
        lr = LogisticRegression(max_iter=300, C=5.0).fit(Xtr, [r[col] for r in tr])
        perm = np.random.default_rng(0).permutation([r[col] for r in tr])
        lrp = LogisticRegression(max_iter=300, C=5.0).fit(Xtr, perm)
        pred_all = {"majority": [maj] * len(allte), "virus_prior": [vprior.get(r["virus"], maj) for r in allte],
                    "virus+length_tree": list(tree.predict(mat(allte))), "3mer_logreg": list(lr.predict(Xte)),
                    "CONTROL_3mer_logreg_train_labels_permuted": list(lrp.predict(Xte))}
        pos = {id(r): i for i, r in enumerate(allte)}
        out[head] = {}
        rowp[head] = {name: {r["accession"]: pa[pos[id(r)]] for r in allte} for name, pa in pred_all.items()}
        for name, pa in pred_all.items():
            out[head][name] = {}
            for s in SPLITS:
                rs = te_by[s]
                if not rs:
                    continue
                y = [r[col] for r in rs]; p = [pa[pos[id(r)]] for r in rs]; cl = [r["cluster"] for r in rs]
                out[head][name][s] = dict(n=len(rs), acc=acc(y, p), bacc=bacc(y, p), acc_ci95=cboot(cl, y, p, acc, reps, rng),
                                          per_virus={v: dict(n=sum(r["virus"] == v for r in rs),
                                                             acc=acc([a for a, r in zip(y, rs) if r["virus"] == v], [a for a, r in zip(p, rs) if r["virus"] == v]))
                                                     for v in sorted({r["virus"] for r in rs})})
        out[head]["_majority_class"] = maj
        out[head]["_by_similarity_pooled_nondev"] = {}
        for name, pa in pred_all.items():
            bins = {}
            for lo, hi in [(0, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.01)]:
                ii = [i for i, r in enumerate(allte) if r["split"] in ("test_geo", "test_host", "test_time", "blind") and r["nn_sim_train"] and lo <= float(r["nn_sim_train"]) < hi]
                bins[f"{lo}-{min(hi, 1.0)}"] = dict(n=len(ii), acc=float(np.mean([pa[i] == allte[i][col] for i in ii])) if ii else None)
            out[head]["_by_similarity_pooled_nondev"][name] = bins
    return out, rowp


def read_preds(path):
    with gzip.open(path, "rt") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def ece(pm, ok, bins=10):
    pm, ok = np.array(pm, float), np.array(ok, float)
    e = 0.0
    for b in range(bins):
        m = (pm > b / bins) & (pm <= (b + 1) / bins)
        if m.any():
            e += m.mean() * abs(pm[m].mean() - ok[m].mean())
    return float(e)


def s3(R, reps, rng, rowp):
    gt = {r["accession"]: r for r in R}
    conds = sorted({re.sub(r"__seed\d+\.tsv\.gz$", "", os.path.basename(p)) for p in glob.glob(os.path.join(ROOT, "results/p4_preds/*.tsv.gz"))})
    conds = [c for c in conds if not c.startswith("util_single") and c != "majority"]
    out = {}
    for c in conds:
        out[c] = {}
        files = sorted(glob.glob(os.path.join(ROOT, f"results/p4_preds/{c}__seed*.tsv.gz")))
        for head, (col, classes, pre) in HEADS.items():
            out[c][head] = {}
            for s in SPLITS:
                per_seed = []
                for fpath in files:
                    rows = [x for x in read_preds(fpath) if x["split"] == s and x[pre + "_true"] in classes]
                    if not rows:
                        continue
                    y = [x[pre + "_true"] for x in rows]; p = [classes[int(x[pre + "_top"])] for x in rows]
                    ab = int(rows[0]["abstain_idx_" + ("host" if head == "host" else "reg")])
                    ans = [int(x[pre + "_act"]) != ab for x in rows]
                    ok = [a == b for a, b in zip(y, p)]
                    tierA = [gt[x["accession"]][("host_tier" if head == "host" else "geo_tier")] == "A" for x in rows]
                    ti = [classes.index(v) for v in y]
                    util = float(np.mean([0 if int(x[pre + "_act"]) == ab else (1 if int(x[pre + "_act"]) == t else -1) for x, t in zip(rows, ti)]))
                    base = rowp[head]["3mer_logreg"]
                    same = all(x["accession"] in base for x in rows)
                    pdiff = None
                    if len(per_seed) == 0 and same and not c.startswith(("util", "wu")):
                        pdiff = paired_diff([gt[x["accession"]]["cluster"] for x in rows], [a == b for a, b in zip(y, p)],
                                            [base[x["accession"]] == t for x, t in zip(rows, y)], reps, rng)
                    d = dict(n=len(rows), acc=acc(y, p), bacc=bacc(y, p), coverage=float(np.mean(ans)), utility=util, diff_vs_3mer_ci95=pdiff,
                             acc_3mer_same_rows=float(np.mean([base[x["accession"]] == t for x, t in zip(rows, y)])) if same else None,
                             acc_when_answering=float(np.mean([o for o, a in zip(ok, ans) if a])) if any(ans) else float("nan"),
                             acc_tierA_only=float(np.mean([o for o, t in zip(ok, tierA) if t])) if any(tierA) else float("nan"),
                             n_tierA=int(sum(tierA)), ece=ece([float(x[pre + "_pmax"]) for x in rows], ok),
                             ci95=cboot([gt[x["accession"]]["cluster"] for x in rows], y, p, acc, reps if len(per_seed) == 0 else 0, rng) if len(per_seed) == 0 else None,
                             per_virus={v: dict(n=sum(x["virus"] == v for x in rows), acc=acc([a for a, x in zip(y, rows) if x["virus"] == v], [a for a, x in zip(p, rows) if x["virus"] == v]))
                                        for v in sorted({x["virus"] for x in rows})})
                    # novelty curve
                    sims = [(float(x["nn_sim_train"]) if x["nn_sim_train"] else None, o) for x, o in zip(rows, ok)]
                    bins = [(0, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.01)]
                    d["acc_by_similarity_to_train"] = {f"{lo}-{min(hi,1.0)}": dict(n=sum(1 for a, o in sims if a is not None and lo <= a < hi),
                                                       acc=float(np.mean([o for a, o in sims if a is not None and lo <= a < hi])) if any(a is not None and lo <= a < hi for a, o in sims) else None) for lo, hi in bins}
                    per_seed.append(d)
                if per_seed:
                    out[c][head][s] = dict(seeds=len(per_seed), acc_mean=float(np.mean([d["acc"] for d in per_seed])),
                                           acc_range=[min(d["acc"] for d in per_seed), max(d["acc"] for d in per_seed)],
                                           **{k: per_seed[0][k] for k in ("n", "ci95", "per_virus", "n_tierA", "acc_by_similarity_to_train", "diff_vs_3mer_ci95", "acc_3mer_same_rows")},
                                           utility_mean=float(np.mean([d["utility"] for d in per_seed])),
                                           bacc_mean=float(np.mean([d["bacc"] for d in per_seed])), coverage_mean=float(np.mean([d["coverage"] for d in per_seed])),
                                           acc_when_answering_mean=float(np.nanmean([d["acc_when_answering"] for d in per_seed])),
                                           acc_tierA_mean=float(np.nanmean([d["acc_tierA_only"] for d in per_seed])),
                                           ece_mean=float(np.mean([d["ece"] for d in per_seed])))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--boot", type=int, default=200); ap.add_argument("--out", default="results/p4_eval_suite.json")
    a = ap.parse_args()
    R = load_gt(); rng = np.random.default_rng(0)
    res = s0_s1(R)
    print("== S0/S1 tripwires"); [print(" -", t) for t in res["tripwires"]]
    aligned = {x["accession"] for x in read_preds(os.path.join(ROOT, "results/p4_preds/ce_long__seed0.tsv.gz"))}
    res["S2_all_labeled"], _ = s2(R, a.boot, rng)
    res["S2_shortcut_baselines"], rowp = s2(R, a.boot, rng, restrict=aligned)
    res["note_S2"] = "S2_shortcut_baselines uses exactly the rows the model prediction files cover (tier-A, both labels independently evidenced); S2_all_labeled uses every labelled row."
    for head in HEADS:
        print(f"\n== S2 {head} (same rows the models are scored on): accuracy by split (majority = {res['S2_shortcut_baselines'][head]['_majority_class']})")
        for name, d in res["S2_shortcut_baselines"][head].items():
            if name.startswith("_"):
                continue
            print(f"  {name:45s}", " ".join(f"{s}={d[s]['acc']:.2f}(n={d[s]['n']})" for s in SPLITS if s in d))
    res["S3_models"] = s3(R, a.boot, rng, rowp)
    for head in HEADS:
        print(chr(10) + f"== S4 {head}: accuracy vs similarity to train (pooled over test_geo/host/time/blind; n in brackets)")
        for name in ("majority", "virus+length_tree", "3mer_logreg"):
            b = res["S2_shortcut_baselines"][head]["_by_similarity_pooled_nondev"][name]
            print(f"  {name:20s}", " ".join(f"{k}={v['acc']:.2f}[{v['n']}]" if v["acc"] is not None else f"{k}=n/a[0]" for k, v in b.items()))
        for c in ("ce", "ce_long", "hybrid", "knn"):
            tot = {}
            for sp in ("test_geo", "test_host", "test_time", "blind"):
                x = res["S3_models"][c][head].get(sp)
                if not x:
                    continue
                for k, v in x["acc_by_similarity_to_train"].items():
                    if v["n"]:
                        t = tot.setdefault(k, [0, 0.0]); t[0] += v["n"]; t[1] += v["n"] * v["acc"]
            print(f"  {c:20s}", " ".join(f"{k}={t[1]/t[0]:.2f}[{t[0]}]" for k, t in sorted(tot.items())))
    main_c = [c for c in res["S3_models"] if not c.startswith(("wu_single", "util_single"))]
    for head in HEADS:
        print(f"\n== S3 {head}: acc/utility/coverage; (acc minus 3-mer baseline on the same rows, 95% CI by cluster)")
        for c in main_c:
            d = res["S3_models"][c][head]
            def cell(s):
                x = d[s]; ci = x["diff_vs_3mer_ci95"]
                return f"{s}={x['acc_mean']:.2f}/{x['utility_mean']:+.2f}/{x['coverage_mean']:.2f}" + (f"(vs3mer {x['acc_mean']-x['acc_3mer_same_rows']:+.2f}[{ci[0]:+.2f},{ci[1]:+.2f}])" if ci and ci[0] == ci[0] else "")
            print(f"  {c:12s}", " ".join(cell(s) for s in SPLITS if s in d))
    json.dump(jsonable(res), open(os.path.join(ROOT, a.out), "w"), indent=1, allow_nan=False)
    print("\nwrote", a.out)


if __name__ == "__main__":
    main()
