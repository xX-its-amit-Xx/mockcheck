#!/usr/bin/env python3
"""Phase 4 idea test (CPU): a retrieval hybrid built from existing predictions.

Rule: if the nearest training sequence is at least as similar as T1 (5-mer Jaccard), answer with ITS labels (retrieval);
otherwise use the fine-tuned model's answer, abstaining when its top probability is below the model's own threshold.
T1 and the model threshold are tuned on dev only (dev is the 'new lineages' split) and then frozen. No test data is used to choose anything.
Writes results/p4_preds/hybrid__seedN.tsv.gz so p4_analyze.py treats it as one more condition."""
import csv, gzip, os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p4_analyze as A  # noqa: E402

ROOT = A.ROOT
P = f"{ROOT}/results/p4_preds"


def read(path):
    return A.read_preds(path)


def build(ce_rows, knn_rows, t1, thr):
    knn = {r["accession"]: r for r in knn_rows}
    out = []
    for r in ce_rows:
        k = knn[r["accession"]]
        q = dict(r)
        if float(k["host_pmax"]) >= t1:                    # retrieval answers
            for h in ("host", "reg"):
                q[f"{h}_act"] = k[f"{h}_act"]; q[f"{h}_top"] = k[f"{h}_top"]; q[f"{h}_pmax"] = "1.0"
        else:                                              # model answers, with its abstain rule
            for h, K in (("host", len(A.HOSTS)), ("reg", len(A.REGIONS))):
                answered = float(r[f"{h}_pmax"]) >= thr
                q[f"{h}_act"] = r[f"{h}_top"] if answered else str(K); q[f"{h}_top"] = r[f"{h}_top"]; q[f"{h}_pmax"] = "1.0"
        out.append(q)
    return out


def main():
    knn_rows = read(f"{P}/knn__seed0.tsv.gz")
    seeds = sorted(int(f.split("__seed")[1].split(".")[0]) for f in os.listdir(P) if f.startswith("ce__seed"))
    ce = {s: read(f"{P}/ce__seed{s}.tsv.gz") for s in seeds}
    thr = A.tune_threshold(ce[seeds[0]])
    best, bt = -9, 0.9
    for t1 in np.arange(0.5, 0.99, 0.02):
        dev = [r for r in build(ce[seeds[0]], knn_rows, float(t1), thr) if r["split"] == "dev"]
        u, _, _ = A.util_vec(dev)
        if u.mean() > best:
            best, bt = u.mean(), float(t1)
    print(f"tuned on dev: retrieval threshold {bt:.2f}, model abstain threshold {thr:.2f}, dev utility {best:.3f}")
    for s in seeds:
        rows = build(ce[s], knn_rows, bt, A.tune_threshold(ce[s]))     # each seed keeps its own dev-tuned abstain threshold, as ce+abstain does
        with gzip.open(f"{P}/hybrid__seed{s}.tsv.gz", "wt", newline="\n") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
