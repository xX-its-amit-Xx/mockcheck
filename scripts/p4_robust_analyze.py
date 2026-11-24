#!/usr/bin/env python3
"""Analyse p4_robust outputs (CPU, light). Writes results/p4_robust.json.
 A  cluster-held-out CV: model accuracy/balanced accuracy per seed vs majority and 3-mer LR on the SAME folds; paired cluster bootstrap.
 B  perturbation curves on the original non-train rows (accuracy per kind, per split).
 C  shuffled-training-label control (should sit near chance / majority).
"""
import csv, glob, gzip, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p4_zoo_eval as Z
from sklearn.feature_extraction.text import CountVectorizer, TfidfTransformer

ROOT = Z.ROOT
D = os.path.join(ROOT, "results/p4_robust")
SEEDS = [0, 1, 2, 3, 4]


def read(p):
    return list(csv.DictReader(gzip.open(p, "rt"), delimiter="\t"))


def maj_acc(y_true, y_train_ref):
    vals, cnt = np.unique(y_train_ref, return_counts=True)
    return float(np.mean(np.array(y_true) == vals[cnt.argmax()]))


def main():
    rows = Z.load_rows(); acc = [r["acc"] for r in rows]
    cv = CountVectorizer(analyzer="char", ngram_range=(3, 3), lowercase=False)
    X3 = TfidfTransformer().fit_transform(cv.fit_transform([r["seq"] for r in rows]))
    out = dict(n_rows=len(rows), n_clusters=len({r["cluster"] for r in rows}), cv={}, pert={}, control={})
    for head, col, tcol, pcol in (("host", "host", "host_true", "host_top"), ("region", "reg", "reg_true", "reg_top")):
        per = dict(model=[], bacc=[], m3=[], maj=[], diff=[], ci=None)
        for s in SEEDS:
            R = {r["accession"]: r for r in read(f"{D}/cv__seed{s}.tsv.gz")}
            keep = [i for i, r in enumerate(rows) if r["acc"] in R]
            sub = [rows[i] for i in keep]
            y = np.array([R[r["acc"]][tcol] for r in sub]); p = np.array([R[r["acc"]][pcol] for r in sub])
            fold = Z.folds_for(rows, s)[keep]
            _, o = Z.oof(sub, X3[keep], None, fold, col)
            ok_m, ok_b = p == y, o["3mer"] == y
            per["model"].append(float(ok_m.mean())); per["bacc"].append(Z.bacc(y, p)); per["m3"].append(float(ok_b.mean()))
            per["maj"].append(maj_acc(y, y)); per["diff"].append(float(ok_m.mean() - ok_b.mean()))
            if s == 0:
                per["ci"] = Z.cboot_diff([r["cluster"] for r in sub], ok_m, ok_b)
        out["cv"][head] = dict(model_acc_mean=float(np.mean(per["model"])), model_acc_seeds=per["model"], model_bacc_mean=float(np.mean(per["bacc"])),
                               m3_acc_mean=float(np.mean(per["m3"])), majority_acc_mean=float(np.mean(per["maj"])),
                               diff_vs_3mer_mean=float(np.mean(per["diff"])), diff_ci95_seed0=per["ci"])
        # B
        agg = {}
        for s in SEEDS:
            for r in read(f"{D}/pert__seed{s}.tsv.gz"):
                agg.setdefault((r["fold_or_pert"], r["split"]), []).append(r[tcol] == r[pcol])
        kinds = ["clean", "trunc_50", "trunc_25", "sub_5", "sub_20", "shuffle"]
        out["pert"][head] = {k: {sp: float(np.mean(v)) for (kk, sp), v in agg.items() if kk == k} for k in kinds}
        out["pert"][head]["_all"] = {k: float(np.mean(sum([v for (kk, sp), v in agg.items() if kk == k], []))) for k in kinds}
        # C
        c = []; mj = []
        for s in SEEDS:
            R = read(f"{D}/labelshuf__seed{s}.tsv.gz"); c.append(np.mean([r[tcol] == r[pcol] for r in R]))
            mj.append(maj_acc([r[tcol] for r in R], [r[tcol] for r in R]))
        out["control"][head] = dict(shuffled_acc_mean=float(np.mean(c)), shuffled_acc_seeds=[float(x) for x in c], majority_acc_mean=float(np.mean(mj)))
    json.dump(out, open(os.path.join(ROOT, "results/p4_robust.json"), "w"), indent=1, allow_nan=False)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
