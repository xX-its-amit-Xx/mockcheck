#!/usr/bin/env python3
"""Phase 4 baselines (CPU): written in the same format as the model predictions so p4_analyze.py treats them alike.
  majority : for each virus, the host class and region most common in the training split (no abstaining)
  knn      : label of the most similar training sequence by shared 5-mers (Jaccard); similarity is its confidence
These show how much of any model's score a trivial rule already gets (the "does it beat a simple baseline" factor)."""
import collections, csv, gzip, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p4_train as T  # noqa: E402  (reuses the exact row filter and label lists)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = f"{ROOT}/results/p4_preds"
os.makedirs(OUT, exist_ok=True)


def kmers(s, k=5):
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def write(path, rows, preds):
    with gzip.open(path, "wt", newline="\n") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["accession", "split", "stratum", "virus", "nn_sim_train", "host_true", "host_act", "host_pmax", "host_top", "reg_true", "reg_act", "reg_pmax", "reg_top", "abstain_idx_host", "abstain_idx_reg"])
        for r, (ha, hp, ra, rp) in zip(rows, preds):
            w.writerow([r["acc"], r["split"], r["stratum"], r["virus"], r["nn"], T.HOSTS[r["y_host"]], ha, f"{hp:.4f}", ha, T.REGIONS[r["y_reg"]], ra, f"{rp:.4f}", ra, len(T.HOSTS), len(T.REGIONS)])


def main():
    rows = T.load_rows(f"{ROOT}/data/groundtruth/p4_groundtruth.tsv.gz")
    train = [r for r in rows if r["split"] == "train"]
    mh = {v: collections.Counter(r["y_host"] for r in train if r["virus"] == v).most_common(1)[0][0] for v in {r["virus"] for r in train}}
    mr = {v: collections.Counter(r["y_reg"] for r in train if r["virus"] == v).most_common(1)[0][0] for v in {r["virus"] for r in train}}
    write(f"{OUT}/majority__seed0.tsv.gz", rows, [(mh[r["virus"]], 1.0, mr[r["virus"]], 1.0) for r in rows])
    tk = [(r, kmers(r["seq"])) for r in train]
    preds = []
    for r in rows:
        ks = kmers(r["seq"]); best, br = -1.0, None
        for t, tks in tk:
            if t["acc"] == r["acc"] or t["virus"] != r["virus"]:
                continue
            j = len(ks & tks) / len(ks | tks)
            if j > best:
                best, br = j, t
        preds.append((br["y_host"], best, br["y_reg"], best))
    write(f"{OUT}/knn__seed0.tsv.gz", rows, preds)
    print("baselines written for", len(rows), "rows")


if __name__ == "__main__":
    main()
