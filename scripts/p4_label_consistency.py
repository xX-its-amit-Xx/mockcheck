#!/usr/bin/env python3
"""Phase 4 label audit (CPU): a sequence-based check that does not rely on the submitter's metadata.

For every host-tier-A sequence, find its most similar other sequence of the same virus (5-mer Jaccard). If the two are
near-identical (J >= T) but carry different host classes, one of two things is true: a real cross-species event, or a label error.
So the disagreement rate among near-identical pairs is an UPPER BOUND on label noise in the neighbourhood (it also counts real
host switches). Reported per virus and per host pair so the biology is visible.
Output: results/p4_label_consistency.tsv"""
import collections, csv, gzip, os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
T = (0.8, 0.9)


def kmers(s, k=5):
    return {s[i:i + k] for i in range(len(s) - k + 1)}


def main():
    rows = [r for r in csv.DictReader(gzip.open(f"{ROOT}/data/groundtruth/p4_groundtruth.tsv.gz", "rt"), delimiter="\t") if r["host_tier"] == "A"]
    out, pairs = [], collections.Counter()
    for v in sorted({r["virus"] for r in rows}):
        vr = [(r, kmers(r["protein"])) for r in rows if r["virus"] == v]
        stat = {t: [0, 0] for t in T}
        for i, (r, ks) in enumerate(vr):
            best, br = 0.0, None
            for j, (q, qs) in enumerate(vr):
                if i == j:
                    continue
                s = len(ks & qs) / len(ks | qs)
                if s > best:
                    best, br = s, q
            for t in T:
                if best >= t:
                    stat[t][0] += 1
                    if br["host_class"] != r["host_class"]:
                        stat[t][1] += 1
                        if t == T[-1]:
                            pairs[(v, " vs ".join(sorted([r["host_class"], br["host_class"]])))] += 1
        for t in T:
            n, d = stat[t]
            out.append(dict(virus=v, similarity_threshold=t, sequences_with_near_neighbour=n, host_disagrees=d, disagreement_rate=round(d / n, 4) if n else ""))
    with open(f"{ROOT}/results/p4_label_consistency.tsv", "w", newline="\n") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(out)
    with open(f"{ROOT}/results/p4_label_consistency_pairs.tsv", "w", newline="\n") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n"); w.writerow(["virus", "host_pair_at_J0.9", "count"])
        for (v, p), n in sorted(pairs.items(), key=lambda t: -t[1]):
            w.writerow([v, p, n])
    print(open(f"{ROOT}/results/p4_label_consistency.tsv").read())
    print(list(pairs.items())[:6])


if __name__ == "__main__":
    main()


def counts():
    """Small table of what is in the set, for the dashboard: virus x split x tier x host class x region."""
    rows = list(csv.DictReader(gzip.open(f"{ROOT}/data/groundtruth/p4_groundtruth.tsv.gz", "rt"), delimiter="\t"))
    c = collections.Counter((r["virus"], r["split"], r["host_class"], r["host_tier"], r["region"], r["geo_tier"]) for r in rows)
    with open(f"{ROOT}/results/p4_gt_counts.tsv", "w", newline="\n") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n"); w.writerow(["virus", "split", "host_class", "host_tier", "region", "geo_tier", "n"])
        for k, n in sorted(c.items()):
            w.writerow(list(k) + [n])


counts()
