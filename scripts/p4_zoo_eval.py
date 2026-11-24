#!/usr/bin/env python3
"""Model-zoo sandbox harness (CPU): does a pretrained representation add anything beyond a 3-mer baseline on novel clusters?

Input : an embedding file per model (npz with `acc` = accessions, `X` = float32 matrix), made by an extraction job.
        With --emb none it scores only the 3-mer baseline (the reference every model must beat).
Method: the SAME cluster-held-out folds as scripts/p4_robust.py (clusters assigned to 5 folds by rng(1000+seed)); rows = the
        tier-A host+region rows the main experiment uses. Per fold: logistic regression and 1-NN on the embeddings, and the 3-mer
        logistic regression, all trained on the other folds. Out-of-fold accuracy and balanced accuracy for host and region,
        the paired difference vs the 3-mer baseline with a cluster bootstrap (seed 0 folds), mean over seeds.
Use   : embeddings are features only. Classification of existing public sequences; nothing is generated.
  python p4_zoo_eval.py --emb path.npz --name esm2_35m --seeds 0,1,2,3,4 --out results/zoo/esm2_35m.json
  python p4_zoo_eval.py --emb none --name baseline_3mer --out results/zoo/baseline_3mer.json
"""
import argparse, csv, gzip, json, os
import numpy as np
from sklearn.feature_extraction.text import CountVectorizer, TfidfTransformer
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
HOSTS = ["avian", "human", "swine", "equine", "canine", "bovine", "feline", "other_mammal", "bat"]
REGIONS = ["africa", "asia", "europe", "north_america", "oceania", "south_america"]
STRATA = [f"influenza_a|{r}" for r in REGIONS] + [f"rotavirus_a|{r}" for r in ("africa", "asia", "oceania")]
MAXLEN = 600


def load_rows():
    rows = []
    with gzip.open(os.path.join(ROOT, "data/groundtruth/p4_groundtruth.tsv.gz"), "rt") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            st = f"{r['virus']}|{r['region']}"
            if r["host_tier"] == "A" and r["geo_tier"] == "A" and r["host_class"] in HOSTS and st in STRATA:
                rows.append(dict(acc=r["accession"], seq=r["protein"][:MAXLEN], host=r["host_class"], reg=r["region"], cluster=r["cluster"], virus=r["virus"]))
    return rows


def bacc(y, p):
    y, p = np.array(y), np.array(p)
    return float(np.mean([np.mean(p[y == c] == c) for c in set(y)]))


def folds_for(rows, seed, k=5):
    cl = sorted({r["cluster"] for r in rows})
    order = np.random.default_rng(1000 + seed).permutation(len(cl))
    fo = {cl[j]: i % k for i, j in enumerate(order)}
    return np.array([fo[r["cluster"]] for r in rows])


def oof(rows, X3, E, fold, col):
    y = np.array([r[col] for r in rows]); n = len(rows)
    out = {"3mer": np.empty(n, object)}
    if E is not None:
        out["lr"] = np.empty(n, object); out["1nn"] = np.empty(n, object)
    for k in sorted(set(fold)):
        te, tr = np.where(fold == k)[0], np.where(fold != k)[0]
        if len(set(y[tr])) < 2:
            for m in out: out[m][te] = y[tr][0]
            continue
        out["3mer"][te] = LogisticRegression(max_iter=300, C=5.0).fit(X3[tr], y[tr]).predict(X3[te])
        if E is not None:
            sc = StandardScaler().fit(E[tr]); a, b = sc.transform(E[tr]), sc.transform(E[te])
            out["lr"][te] = LogisticRegression(max_iter=300, C=0.1).fit(a, y[tr]).predict(b)
            en = E / (np.linalg.norm(E, axis=1, keepdims=True) + 1e-9)
            out["1nn"][te] = KNeighborsClassifier(1).fit(en[tr], y[tr]).predict(en[te])
    return y, out


def cboot_diff(clusters, ok_a, ok_b, reps=300):
    cl = np.array(clusters); a = np.array(ok_a, float); b = np.array(ok_b, float)
    u = sorted(set(cl)); idx = {c: np.where(cl == c)[0] for c in u}; rng = np.random.default_rng(0); v = []
    for _ in range(reps):
        ii = np.concatenate([idx[u[j]] for j in rng.integers(0, len(u), len(u))]); v.append(a[ii].mean() - b[ii].mean())
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--emb", required=True); ap.add_argument("--name", required=True)
    ap.add_argument("--seeds", default="0,1,2,3,4"); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = load_rows()
    cv = CountVectorizer(analyzer="char", ngram_range=(3, 3), lowercase=False)
    X3 = TfidfTransformer().fit_transform(cv.fit_transform([r["seq"] for r in rows]))
    E = None
    if a.emb != "none":
        z = np.load(a.emb, allow_pickle=True); pos = {x: i for i, x in enumerate(z["acc"])}
        miss = [r["acc"] for r in rows if r["acc"] not in pos]
        if miss:
            raise SystemExit(f"{len(miss)} rows have no embedding (first: {miss[0]})")
        E = z["X"][[pos[r["acc"]] for r in rows]].astype(np.float32)
    res = dict(name=a.name, n_rows=len(rows), n_clusters=len({r["cluster"] for r in rows}), emb_dim=None if E is None else int(E.shape[1]), heads={})
    for head, col in (("host", "host"), ("region", "reg")):
        per = {}
        for seed in map(int, a.seeds.split(",")):
            fold = folds_for(rows, seed); y, o = oof(rows, X3, E, fold, col)
            for m, p in o.items():
                d = per.setdefault(m, dict(acc=[], bacc=[], diff_vs_3mer=[], ci=None))
                d["acc"].append(float(np.mean(p == y))); d["bacc"].append(bacc(y, p))
                if m != "3mer":
                    d["diff_vs_3mer"].append(float(np.mean(p == y) - np.mean(o["3mer"] == y)))
                    if seed == int(a.seeds.split(",")[0]):
                        d["ci"] = cboot_diff([r["cluster"] for r in rows], p == y, o["3mer"] == y)
        res["heads"][head] = {m: dict(acc_mean=float(np.mean(d["acc"])), acc_seeds=d["acc"], bacc_mean=float(np.mean(d["bacc"])),
                                      diff_vs_3mer_mean=(float(np.mean(d["diff_vs_3mer"])) if d["diff_vs_3mer"] else None), diff_ci95_seed0=d["ci"]) for m, d in per.items()}
        print(head, {m: round(v["acc_mean"], 3) for m, v in res["heads"][head].items()}, {m: (round(v["diff_vs_3mer_mean"], 3), v["diff_ci95_seed0"]) for m, v in res["heads"][head].items() if v["diff_vs_3mer_mean"] is not None})
    os.makedirs(os.path.dirname(os.path.join(ROOT, a.out)), exist_ok=True)
    json.dump(res, open(os.path.join(ROOT, a.out), "w"), indent=1, allow_nan=False)
    print("wrote", a.out)


if __name__ == "__main__":
    main()
