#!/usr/bin/env python3
"""Phase 4 analysis (CPU only): turn prediction files into the tables the dashboard shows.

Inputs : results/p4_preds/*.tsv.gz  (from p4_train.py), data/groundtruth/p4_groundtruth.tsv.gz
Outputs: results/p4_eval_by_stratum.tsv   condition x stratum x evalgroup: utility, accuracy when answering, coverage
         results/p4_condition_summary.tsv macro-average utility, worst stratum, with bootstrap 95% CI
         results/p4_transfer_matrix.tsv   train-on-one-stratum rows x evaluate-on-each-stratum columns (held-out utility)
Eval groups: dev (in-distribution), heldout (test_geo + test_host + test_time), blind (sealed; only with --include-blind).
Utility U = mean of +1 right, -1 wrong, 0 abstain, averaged over the host head and the region head.
CE-style conditions (probe, ce) and the nearest-neighbour baseline get a post-hoc abstain rule: abstain when max prob < t, with t chosen on dev."""
import argparse, collections, csv, glob, gzip, os, sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LAM = 1.0
HOSTS = ["avian", "human", "swine", "equine", "canine", "bovine", "feline", "other_mammal", "bat"]
REGIONS = ["africa", "asia", "europe", "north_america", "oceania", "south_america"]
GROUPS = {"dev": {"dev"}, "heldout": {"test_geo", "test_host", "test_time"}, "blind": {"blind"}}


def read_preds(path):
    with gzip.open(path, "rt") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def util_vec(rows, thr_h=None, thr_r=None):
    """Per-sample utility (host+region)/2 and per-head answered/correct flags."""
    out = []
    for r in rows:
        us, ans, cor = [], [], []
        for head, true, K, thr in (("host", HOSTS, len(HOSTS), thr_h), ("reg", REGIONS, len(REGIONS), thr_r)):
            t = true.index(r["host_true"]) if head == "host" else true.index(r["reg_true"])
            act, top, pmax = int(r[f"{head}_act"]), int(r[f"{head}_top"]), float(r[f"{head}_pmax"])
            if thr is not None:                 # CE-style: no abstain action, rule on max prob
                act = top if pmax >= thr else K
            answered = act != K
            correct = answered and act == t
            us.append(1.0 if correct else (-LAM if answered else 0.0)); ans.append(answered); cor.append(correct)
        out.append((sum(us) / 2, ans[0] + ans[1], cor[0] + cor[1]))
    return np.array([o[0] for o in out]), np.array([o[1] for o in out]), np.array([o[2] for o in out])


def tune_threshold(rows):
    dev = [r for r in rows if r["split"] in GROUPS["dev"]]
    best, bt = -9, 0.5
    for t in np.arange(0.30, 0.98, 0.02):
        u, _, _ = util_vec(dev, t, t)
        if u.mean() > best:
            best, bt = u.mean(), float(t)
    return bt


TIERS = [("seen", "in the training set (memorisation reference)"), ("near", "similarity to training >= 0.9"), ("mid", "similarity 0.7 to 0.9"), ("far", "similarity < 0.7")]


def tier_of(r):
    if r["split"] == "train":
        return "seen"
    sim = float(r["nn_sim_train"])
    return "near" if sim >= 0.9 else ("mid" if sim >= 0.7 else "far")


def tiered(data, include_blind, boot, rng):
    """Utility by distance tier, pooled over strata. Tier comes from the sequence's own similarity to training, not from the split name."""
    out = []
    for cond, runs in sorted(data.items()):
        for tname, _ in TIERS:
            per_seed, answered, correct, n = [], [], [], 0
            for seed, rows, thr in runs:
                sub = [r for r in rows if (include_blind or r["split"] != "blind") and tier_of(r) == tname]
                if not sub:
                    continue
                u, an, co = util_vec(sub, thr, thr)
                per_seed.append(u); answered.append(an.mean() / 2); correct.append(co.sum() / max(an.sum(), 1)); n = len(sub)
            if not per_seed or len({len(x) for x in per_seed}) != 1:
                continue
            m = np.mean(per_seed, 0)
            bs = [m[rng.integers(0, len(m), len(m))].mean() for _ in range(boot)]
            out.append(dict(condition=cond, tier=tname, n_sequences=n, seeds=len(per_seed), utility=round(float(m.mean()), 4), ci_lo=round(float(np.percentile(bs, 2.5)), 4),
                            ci_hi=round(float(np.percentile(bs, 97.5)), 4), coverage=round(float(np.mean(answered)), 4), accuracy_when_answering=round(float(np.mean(correct)), 4)))
    if out:
        with open(f"{ROOT}/results/p4_tiered.tsv", "w", newline="\n") as f:
            w = csv.DictWriter(f, fieldnames=list(out[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(out)
    return len(out)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--preds", default=f"{ROOT}/results/p4_preds"); ap.add_argument("--include-blind", action="store_true")
    ap.add_argument("--boot", type=int, default=500); a = ap.parse_args()
    files = sorted(glob.glob(f"{a.preds}/*.tsv.gz"))
    if not files:
        sys.exit("no prediction files")
    data = collections.defaultdict(list)      # cond -> list of (seed, rows, thr)
    for p in files:
        tag = os.path.basename(p)[:-7]; cond, seed = tag.rsplit("__seed", 1)
        rows = read_preds(p)
        thr = tune_threshold(rows) if cond in ("ce", "probe") or cond == "knn" else None
        cond_name = cond + ("+abstain" if thr is not None else "")
        data[cond_name].append((int(seed), rows, thr))
    groups = ["dev", "heldout"] + (["blind"] if a.include_blind else [])
    strata = sorted({r["stratum"] for rows in [v[0][1] for v in data.values()] for r in rows})
    rng = np.random.default_rng(0)
    by = []       # eval_by_stratum rows
    cache = {}
    for cond, runs in sorted(data.items()):
        for g in groups:
            for s in strata:
                us, answ, corr = [], [], []
                for seed, rows, thr in runs:
                    sub = [r for r in rows if r["stratum"] == s and r["split"] in GROUPS[g]]
                    if not sub:
                        continue
                    u, an, co = util_vec(sub, thr, thr)
                    us.append(u.mean()); answ.append(an.mean() / 2); corr.append(co.sum() / max(an.sum(), 1))
                    cache[(cond, g, s, seed)] = u
                if us:
                    by.append(dict(condition=cond, group=g, stratum=s, n_eval=len(sub), seeds=len(us), utility=round(float(np.mean(us)), 4),
                                   accuracy_when_answering=round(float(np.mean(corr)), 4), coverage=round(float(np.mean(answ)), 4)))
    with open(f"{ROOT}/results/p4_eval_by_stratum.tsv", "w", newline="\n") as f:
        w = csv.DictWriter(f, fieldnames=list(by[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(by)
    # condition summary with bootstrap over evaluation samples (resampled within stratum, same draw for all seeds)
    summ = []
    for cond, runs in sorted(data.items()):
        for g in groups:
            per = {}
            for s in strata:
                arrs = [cache[(cond, g, s, sd)] for sd, _, _ in runs if (cond, g, s, sd) in cache]
                if arrs:
                    per[s] = np.mean(arrs, 0)           # per-sample utility averaged over seeds
            if not per:
                continue
            macro = np.mean([v.mean() for v in per.values()]); worst = min(v.mean() for v in per.values())
            bm = []
            for _ in range(a.boot):
                vals = [v[rng.integers(0, len(v), len(v))].mean() for v in per.values()]
                bm.append((np.mean(vals), min(vals)))
            bm = np.array(bm)
            summ.append(dict(condition=cond, group=g, seeds=len(runs), strata=len(per), macro_utility=round(float(macro), 4),
                             macro_ci_lo=round(float(np.percentile(bm[:, 0], 2.5)), 4), macro_ci_hi=round(float(np.percentile(bm[:, 0], 97.5)), 4),
                             worst_stratum_utility=round(float(worst), 4), worst_ci_lo=round(float(np.percentile(bm[:, 1], 2.5)), 4), worst_ci_hi=round(float(np.percentile(bm[:, 1], 97.5)), 4)))
    with open(f"{ROOT}/results/p4_condition_summary.tsv", "w", newline="\n") as f:
        w = csv.DictWriter(f, fieldnames=list(summ[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(summ)
    # transfer matrix on held-out data: rows = trained on stratum s only; columns = evaluated on stratum t
    tm = []
    ref = {r["stratum"]: r["utility"] for r in by if r["condition"] == "util_avg" and r["group"] == "heldout"}
    for r in by:
        if r["condition"].startswith("util_single") and r["group"] == "heldout":
            trained = r["condition"].split("-", 1)[1].replace("_", "|", 1)
            tm.append(dict(trained_on=trained, evaluated_on=r["stratum"], utility=r["utility"], same_stratum=int(trained == r["stratum"]), util_avg_reference=ref.get(r["stratum"], "")))
    if tm:
        with open(f"{ROOT}/results/p4_transfer_matrix.tsv", "w", newline="\n") as f:
            w = csv.DictWriter(f, fieldnames=list(tm[0]), delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(tm)
    nt = tiered(data, a.include_blind, a.boot, rng)
    print(f"tier rows {nt}")
    print(f"conditions {sorted(data)}; eval rows {len(by)}; summary rows {len(summ)}; transfer rows {len(tm)}")


if __name__ == "__main__":
    main()
